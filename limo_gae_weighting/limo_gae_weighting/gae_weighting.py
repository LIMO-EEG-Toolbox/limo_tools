#!/usr/bin/env python3
"""Bounded group-blind quality weights for LIMO's second level.

beta-maps -> masked-channel graph autoencoder (out-of-fold) -> per-(channel, participant)
deviation -> within-channel robust z -> one-sided Tukey bisquare (tau=3, c=4.685, mean 1).

datain [nBeta, nTime, nChan, nSub], adjacency [nChan, nChan] binary symmetric.
Group labels are not used.

    python gae_weighting.py data.npz --beta 0 --out weights.npz
"""
from __future__ import annotations

import os

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse
import numpy as np
import torch
import torch.nn as nn
from torch_geometric.nn import GCNConv
from torch_geometric.utils import dense_to_sparse

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
TAU, BISQUARE_C = 3.0, 4.685
FOLDS, EPOCHS, HIDDEN, MASK_P = 4, 60, 32, 0.3


def build_knn_adjacency(coords, k=6):
    n = coords.shape[0]
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    A = np.zeros((n, n), np.float32)
    for i in range(n):
        A[i, np.argsort(d[i])[:k]] = 1.0
    return np.maximum(A, A.T)


class GraphAE(nn.Module):
    """node = channel, node feature = that channel's beta time series."""

    def __init__(self, n_time, edge_index, n_chan, H=HIDDEN, depth=1, p=0.3):
        super().__init__()
        dims = [n_time] + [H] * depth
        self.convs = nn.ModuleList([GCNConv(dims[i], dims[i + 1]) for i in range(depth)])
        self.dec = nn.Sequential(nn.Linear(H, H), nn.ReLU(), nn.Linear(H, n_time))
        self.drop = nn.Dropout(p)
        self.register_buffer("ei", edge_index)
        self.n_chan = n_chan

    def forward(self, X):                                   # [B, nChan, nTime]
        B = X.shape[0]
        h = X.reshape(B * self.n_chan, -1)
        offs = (torch.arange(B, device=X.device) * self.n_chan).repeat_interleave(self.ei.shape[1])
        eib = self.ei.repeat(1, B) + offs.unsqueeze(0)      # one disjoint graph per map
        for i, conv in enumerate(self.convs):
            h = torch.relu(conv(h, eib))
            if i == 0:
                h = self.drop(h)
        return self.dec(h).reshape(B, self.n_chan, -1)


def _train(model, Xtr, epochs=EPOCHS, mask_p=MASK_P, masked=True):
    """masked=True: zero a random subset of nodes, loss on those nodes only."""
    opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    lossf = nn.SmoothL1Loss()
    for _ in range(epochs):
        model.train()
        if masked:
            msk = torch.rand(Xtr.shape[0], Xtr.shape[1], device=Xtr.device) < mask_p
            xin = Xtr.clone()
            xin[msk] = 0.0
            loss = lossf(model(xin)[msk], Xtr[msk])
        else:
            loss = lossf(model(Xtr), Xtr)
        opt.zero_grad()
        loss.backward()
        opt.step()
    model.eval()
    return model


@torch.no_grad()
def _score_map(model, xs, masked=True):
    """xs [nChan, nTime] -> deviation [nChan]: mask each channel in turn, predict it."""
    n_chan = xs.shape[0]
    if masked:
        b = xs.unsqueeze(0).repeat(n_chan, 1, 1)
        ii = torch.arange(n_chan, device=xs.device)
        b[ii, ii] = 0.0
        return (model(b)[ii, ii] - xs[ii]).abs().mean(1).cpu().numpy()
    return (model(xs.unsqueeze(0))[0] - xs).abs().mean(1).cpu().numpy()


def masked_deviation(datain, adjacency, folds=FOLDS, epochs=EPOCHS, depth=1, seed=0,
                     masked=True, verbose=False):
    """Out-of-fold deviation -> [nBeta, nChan, nSub]."""
    torch.manual_seed(seed)
    nBeta, nTime, nChan, nSub = datain.shape
    ei = dense_to_sparse(torch.tensor(np.asarray(adjacency, np.float32)))[0].to(DEVICE)
    dev = np.zeros((nBeta, nChan, nSub))
    idx = np.arange(nSub)
    for f in range(folds):
        test = idx[f::folds]
        tr = np.array([s for s in idx if s not in test])
        Xtr = torch.tensor(np.concatenate(
            [datain[b][:, :, tr].transpose(2, 1, 0) for b in range(nBeta)], 0),
            dtype=torch.float32, device=DEVICE)
        model = _train(GraphAE(nTime, ei, nChan, depth=depth).to(DEVICE), Xtr, epochs, masked=masked)
        for s in test:
            for b in range(nBeta):
                xs = torch.tensor(datain[b, :, :, s].T, dtype=torch.float32, device=DEVICE)
                dev[b, :, s] = _score_map(model, xs, masked)
        if verbose:
            print(f"  fold {f + 1}/{folds}")
    return dev


def robust_z(x):
    med = np.median(x)
    return (x - med) / (1.4826 * np.median(np.abs(x - med)) + 1e-9)


def dev_to_weight(dev2d, tau=TAU, c=BISQUARE_C, shape="bisquare"):
    """[nChan, nSub] -> W [nChan, nSub]. Flat below tau, zero beyond tau + c."""
    W = np.ones_like(dev2d, dtype=float)
    for ch in range(dev2d.shape[0]):
        e = np.clip(robust_z(dev2d[ch]) - tau, 0, None)
        if shape == "bisquare":
            w = np.where(e < c, (1.0 - (e / c) ** 2) ** 2, 0.0)
        elif shape == "powerlaw":
            w = 1.0 / (1.0 + e ** 2)
        else:
            raise ValueError(shape)
        s = w.sum()
        W[ch] = w * len(w) / s if s > 0 else 1.0
    return W


def compute_weights(datain, adjacency, beta=0, tau=TAU, seeds=(0,), **kw):
    dev = np.mean([masked_deviation(datain, adjacency, seed=s, **kw) for s in seeds], 0)
    return dev_to_weight(dev.mean(0) if beta == "all" else dev[int(beta)], tau)


def weighted_welch_t(Y, groups, W):
    """Y [nTime, nChan, nSub], W [nChan, nSub]. Kish effective sample size in SE and df.
    -> t [nChan, nTime], df, diff (group1 - group0)."""
    def stat(gi):
        Yg, Wg = Y[:, :, gi], W[:, gi]
        den = Wg.sum(1)
        mu = (Yg * Wg[None]).sum(2) / den[None]
        corr = np.maximum(den - (Wg ** 2).sum(1) / den, 1e-9)
        var = (Wg[None] * (Yg - mu[:, :, None]) ** 2).sum(2) / corr[None]
        ess = den ** 2 / (Wg ** 2).sum(1)
        return mu, var / ess[None], ess

    mu0, se0, e0 = stat(groups == 0)
    mu1, se1, e1 = stat(groups == 1)
    t = (mu1 - mu0) / np.sqrt(se0 + se1 + 1e-12)
    df = (se0 + se1) ** 2 / (se0 ** 2 / (e0 - 1)[None] + se1 ** 2 / (e1 - 1)[None] + 1e-12)
    return t.T, df.T, (mu1 - mu0).T


def _tm(x, tr=0.2):
    n = x.shape[-1]; g = int(tr * n); xs = np.sort(x, -1)
    return xs[..., g:n - g].mean(-1)


def _wv(x, tr=0.2):
    n = x.shape[-1]; g = int(tr * n); xs = np.sort(x, -1)
    return np.clip(x, xs[..., g:g + 1], xs[..., n - g - 1:n - g]).var(-1, ddof=1)


def yuen_t(Y, groups, tr=0.2):
    """LIMO 'Trimmed Mean' -> t, diff (group1 - group0)."""
    out = []
    for g in (0, 1):
        x = Y[:, :, groups == g]
        n = x.shape[-1]; h = n - 2 * int(tr * n)
        out.append((_tm(x, tr), _wv(x, tr) * (n - 1) / (h * (h - 1))))
    (m0, d0), (m1, d1) = out
    return ((m1 - m0) / np.sqrt(d0 + d1 + 1e-12)).T, (m1 - m0).T


def main():
    p = argparse.ArgumentParser()
    p.add_argument("npz")
    p.add_argument("--beta", default="0")
    p.add_argument("--tau", type=float, default=TAU)
    p.add_argument("--seeds", type=int, default=1)
    p.add_argument("--out", default="weights.npz")
    a = p.parse_args()
    d = np.load(a.npz, allow_pickle=True)
    datain = d["datain"].astype(np.float32)
    A = d["adjacency"].astype(np.float32) if "adjacency" in d \
        else build_knn_adjacency(d["coords"].astype(float))
    beta = "all" if a.beta == "all" else int(a.beta)
    print(f"datain {datain.shape}, mean degree {A.sum(1).mean():.1f}, {DEVICE}")
    dev = np.mean([masked_deviation(datain, A, seed=s, verbose=True) for s in range(a.seeds)], 0)
    W = dev_to_weight(dev.mean(0) if beta == "all" else dev[beta], a.tau)
    np.savez(a.out, W=W, dev=dev, tau=a.tau, beta=str(beta))
    print(f"flagged {100 * (W < 1).mean():.2f}%, min weight {W.min():.3f} -> {a.out}")


if __name__ == "__main__":
    main()
