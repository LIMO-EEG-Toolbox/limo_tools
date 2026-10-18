"""Extracted/refactored shared-beta GAE; no weighting, inference, or file I/O.

Based on LIMO-EEG-Toolbox/limo_tools, group_outliers commit
2e68021088066be1f8177f9b5f108f3d8b785ff5,
limo_gae_weighting/limo_gae_weighting/gae_weighting.py.
Preserves GCN encoder/linear decoder, SmoothL1, Adam, shared-beta training,
and participant-level folds. Returns predictions instead of deviations.
"""
from __future__ import annotations

import platform
from contextlib import contextmanager

import numpy as np
import torch
from torch import nn
import torch_geometric
from torch_geometric.nn import GCNConv
from torch_geometric.utils import dense_to_sparse


class GraphAE(nn.Module):
    """Original shared model: one disjoint channel graph per beta map."""

    def __init__(self, n_time, edge_index, n_chan, hidden=32, depth=1, dropout=0.3):
        super().__init__()
        dims = [n_time] + [hidden] * depth
        self.convs = nn.ModuleList(GCNConv(dims[i], dims[i + 1]) for i in range(depth))
        self.dec = nn.Sequential(nn.Linear(hidden, hidden), nn.ReLU(), nn.Linear(hidden, n_time))
        self.drop = nn.Dropout(dropout)
        self.register_buffer("ei", edge_index)
        self.n_chan = n_chan

    def forward(self, X):
        batch = X.shape[0]
        h = X.reshape(batch * self.n_chan, -1)
        offsets = (torch.arange(batch, device=X.device) * self.n_chan).repeat_interleave(self.ei.shape[1])
        edges = self.ei.repeat(1, batch) + offsets.unsqueeze(0)
        for i, conv in enumerate(self.convs):
            h = torch.relu(conv(h, edges))
            if i == 0:
                h = self.drop(h)
        return self.dec(h).reshape(batch, self.n_chan, -1)


def _mask_inputs(X, probability):
    mask = torch.rand(X.shape[:2], device=X.device) < probability
    # The source can produce an empty loss. Guarantee one masked node per batch.
    if not mask.any():
        index = torch.randint(mask.numel(), (1,), device=X.device)
        mask.reshape(-1)[index] = True
    hidden = X.clone()
    hidden[mask] = 0
    return hidden, mask


def _train(model, X, *, epochs, masked, mask_probability, batch_size,
           learning_rate, weight_decay):
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    loss_fn = nn.SmoothL1Loss()
    losses, masked_counts = [], []
    step = len(X) if batch_size is None else batch_size
    for _ in range(epochs):
        model.train()
        total, count = 0.0, 0
        for start in range(0, len(X), step):
            target = X[start:start + step]
            if masked:
                inputs, mask = _mask_inputs(target, mask_probability)
                loss = loss_fn(model(inputs)[mask], target[mask])
                n = int(mask.sum().item())
            else:
                loss = loss_fn(model(target), target)
                n = target.shape[0] * target.shape[1]
            if not torch.isfinite(loss):
                raise RuntimeError("Nonfinite training loss; inspect input scale and learning rate.")
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            total += loss.item() * n
            count += n
        losses.append(total / count)
        masked_counts.append(count if masked else 0)
    model.eval()
    return {"loss": losses, "masked_node_counts": masked_counts}


@torch.no_grad()
def _predict_map(model, xs, masked=True, channel_batch_size=16):
    """Return channel x frame; each selected target is zeroed in its own graph."""
    if not masked:
        return model(xs.unsqueeze(0))[0]
    prediction = torch.empty_like(xs)
    for start in range(0, xs.shape[0], channel_batch_size):
        channels = torch.arange(start, min(start + channel_batch_size, xs.shape[0]), device=xs.device)
        rows = torch.arange(len(channels), device=xs.device)
        batch = xs.unsqueeze(0).repeat(len(channels), 1, 1)
        batch[rows, channels] = 0
        prediction[channels] = model(batch)[rows, channels]
    return prediction


def _integer(name, value, minimum=1):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}.")
    return int(value)


def _number(name, value, minimum, maximum, strict_min=False):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(f"Invalid {name}.")
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {name}.") from exc
    if not np.isfinite(value) or value < minimum or value > maximum or (strict_min and value == minimum):
        raise ValueError(f"Invalid {name}.")
    return value


@contextmanager
def _random_state(seed, deterministic, device):
    # Restore caller RNG and determinism settings even after a failed fit.
    devices = list(range(torch.cuda.device_count())) if device.type == "cuda" else []
    previous = torch.are_deterministic_algorithms_enabled()
    warn_only = torch.is_deterministic_algorithms_warn_only_enabled()
    with torch.random.fork_rng(devices=devices):
        torch.random.default_generator.manual_seed(seed)
        if devices:
            torch.cuda.manual_seed_all(seed)
        try:
            torch.use_deterministic_algorithms(deterministic)
            yield
        finally:
            torch.use_deterministic_algorithms(previous, warn_only=warn_only)


def _resolve_device(device):
    """Prefer CUDA automatically; honor explicit CPU/CUDA overrides."""
    if device == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if not isinstance(device, (str, torch.device)):
        raise ValueError("device must be 'auto', 'cpu', 'cuda', or 'cuda:N'.")
    try:
        selected = torch.device(device)
    except (RuntimeError, ValueError) as exc:
        raise ValueError("Invalid device; use 'auto', 'cpu', 'cuda', or 'cuda:N'.") from exc
    if selected.type not in ("cpu", "cuda"):
        raise ValueError("Only CPU and CUDA devices are supported.")
    if selected.type == "cuda":
        if not torch.cuda.is_available():
            raise ValueError("CUDA requested but unavailable; use device='auto' for CPU fallback.")
        if selected.index is not None and selected.index >= torch.cuda.device_count():
            raise ValueError(f"CUDA device index {selected.index} is unavailable.")
    return selected


def limo_gae_predict(betas, adjacency, *, prediction="out-of-fold", folds=4,
                     fold_ids=None, participant_ids=None, masked=True, epochs=60,
                     hidden_dim=32, depth=1, dropout=0.3, batch_size=None,
                     mask_probability=0.3, learning_rate=1e-3, weight_decay=1e-4,
                     seed=0, device="auto", deterministic=True,
                     channel_batch_size=16, return_models=False):
    """Train the existing shared-beta GAE and return aligned predictions.

    betas: finite real [parameter, frame, channel, participant], always 4D.
    adjacency: binary symmetric [channel, channel]. Diagonal is removed;
        GCNConv adds one self-loop per channel during graph convolution.
    fold_ids: optional nonnegative integer label per participant; >=2 unique
        labels in out-of-fold mode. Overrides folds. Default stride assignment
        matches the source. Seed changes model/masks, not default assignment.
    participant_ids: unique strings; supplies stable training order across
        participant permutations when the IDs and fold assignments move too.
    batch_size=None preserves the source's full-batch optimization.
    device='auto' (default) uses CUDA if available, otherwise CPU. Override
        with 'cpu', 'cuda', or 'cuda:N'. Explicit CUDA requests never silently
        fall back. CUDA requires a CUDA-enabled PyTorch installation.
    Returns reconstructed float32 array with the identical 4D shape,
        prediction_info (JSON-compatible), and optional models keyed by fold.
    No group labels, preprocessing, learned scaling, or early stopping.
    """
    B = np.asarray(betas)
    if B.ndim != 4 or 0 in B.shape or B.dtype.kind not in "fiu" or not np.isfinite(B).all():
        raise ValueError("betas must be finite real numeric 4D: parameter x frame x channel x participant.")
    with np.errstate(over="ignore"):
        B = B.astype(np.float32, copy=True)
    if not np.isfinite(B).all():
        raise ValueError("betas exceed float32 range; rescale before training.")
    nb, nt, nc, ns = B.shape
    A = np.asarray(adjacency)
    if A.shape != (nc, nc) or A.dtype.kind not in "bfiu" or not np.isin(A, [0, 1]).all() or not np.array_equal(A, A.T):
        raise ValueError("adjacency must be a binary symmetric channel x channel matrix.")
    A = A.astype(np.float32, copy=True)
    np.fill_diagonal(A, 0)
    for name, value in (("masked", masked), ("deterministic", deterministic), ("return_models", return_models)):
        if not isinstance(value, (bool, np.bool_)):
            raise ValueError(f"{name} must be boolean.")
    epochs = _integer("epochs", epochs)
    hidden_dim = _integer("hidden_dim", hidden_dim)
    depth = _integer("depth", depth)
    channel_batch_size = _integer("channel_batch_size", channel_batch_size)
    seed = _integer("seed", seed, 0)
    if seed >= 2**63 - 1:
        raise ValueError("seed must be less than 2**63-1.")
    if batch_size is not None:
        batch_size = _integer("batch_size", batch_size)
    dropout = _number("dropout", dropout, 0, 1)
    if dropout == 1:
        raise ValueError("dropout must be less than one.")
    mask_probability = _number("mask_probability", mask_probability, 0, 1, True)
    learning_rate = _number("learning_rate", learning_rate, 0, np.inf, True)
    weight_decay = _number("weight_decay", weight_decay, 0, np.inf)
    if prediction not in ("out-of-fold", "in-sample"):
        raise ValueError("prediction must be out-of-fold or in-sample.")
    if participant_ids is None:
        ids = [f"participant-{i:09d}" for i in range(ns)]
    else:
        ids = list(participant_ids)
        if len(ids) != ns or any(not isinstance(x, str) or not x for x in ids) or len(set(ids)) != ns:
            raise ValueError("participant_ids must contain one unique nonempty string per participant.")
    if prediction == "in-sample":
        if fold_ids is not None:
            raise ValueError("fold_ids is not applicable to in-sample prediction.")
        labels = np.zeros(ns, dtype=np.int64)
    elif fold_ids is None:
        folds = _integer("folds", folds, 2)
        if folds > ns:
            raise ValueError("folds cannot exceed participant count.")
        labels = np.arange(ns) % folds
    else:
        labels = np.asarray(fold_ids)
        if labels.shape != (ns,) or labels.dtype.kind not in "iu" or np.any(labels < 0) or len(np.unique(labels)) < 2:
            raise ValueError("fold_ids must be nonnegative integer labels with at least two nonempty folds.")
    selected_device = _resolve_device(device)
    ei = dense_to_sparse(torch.from_numpy(A))[0].to(selected_device)
    reconstructed = np.empty_like(B)
    records, models = [], {}
    for ordinal, label in enumerate(np.unique(labels)):
        test = np.flatnonzero(labels == label)
        train = np.arange(ns) if prediction == "in-sample" else np.flatnonzero(labels != label)
        train = np.array(sorted(train, key=lambda s: ids[s]), dtype=int)
        fold_seed = (seed + ordinal) % (2**63 - 1)
        # np.take avoids NumPy advanced-index axis reordering. All maps from
        # a participant stay together in the fold, then maps are beta-major.
        maps = np.take(B, train, axis=3).transpose(0, 3, 2, 1).reshape(nb * len(train), nc, nt)
        X = torch.tensor(maps, device=selected_device)
        with _random_state(fold_seed, bool(deterministic), selected_device):
            model = GraphAE(nt, ei, nc, hidden_dim, depth, dropout).to(selected_device)
            diagnostics = _train(model, X, epochs=epochs, masked=masked,
                                 mask_probability=mask_probability, batch_size=batch_size,
                                 learning_rate=learning_rate, weight_decay=weight_decay)
            for s in test:
                for b in range(nb):
                    xs = torch.tensor(B[b, :, :, s].T, device=selected_device)
                    pred = _predict_map(model, xs, masked, channel_batch_size)
                    if not torch.isfinite(pred).all():
                        raise RuntimeError("Nonfinite reconstruction.")
                    reconstructed[b, :, :, s] = pred.cpu().numpy().T
        records.append({"fold": int(label), "seed": fold_seed,
                        "train_indices": train.tolist(), "test_indices": test.tolist(),
                        "train_participant_ids": [ids[s] for s in train],
                        "test_participant_ids": [ids[s] for s in test],
                        "training_maps": len(train) * nb, **diagnostics})
        if return_models:
            models[int(label)] = model.cpu()
    info = {"prediction": prediction, "masked": bool(masked), "shared_across_betas": True,
            "layout": "parameter x frame x channel x participant", "input_shape": list(B.shape),
            "index_base": 0, "participant_ids": ids, "fold_ids": labels.tolist(), "folds": records,
            "settings": {"epochs": epochs, "hidden_dim": hidden_dim, "depth": depth, "dropout": dropout,
                         "batch_size": batch_size, "mask_probability": mask_probability,
                         "learning_rate": learning_rate, "weight_decay": weight_decay,
                         "seed": seed, "device_requested": str(device),
                         "device": str(selected_device), "deterministic": bool(deterministic),
                         "channel_batch_size": channel_batch_size, "fold_count": len(records)},
            "self_loops": "input diagonal removed; GCNConv adds self-loops",
            "versions": {"python": platform.python_version(), "numpy": np.__version__,
                         "torch": torch.__version__, "torch_geometric": torch_geometric.__version__}}
    result = {"reconstructed": reconstructed, "prediction_info": info}
    if return_models:
        result["models"] = models
    return result
