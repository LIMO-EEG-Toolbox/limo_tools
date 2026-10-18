"""Small, reproducible CPU smoke example; writes a JSON summary to stdout."""
import json
import numpy as np
import torch

from limo_gae_predict import limo_gae_predict


if __name__ == "__main__":
    torch.set_num_threads(1)
    # Shared smooth time course: synthetic channel x frame x parameter maps.
    rng = np.random.default_rng(17)
    t = np.linspace(0, 2 * np.pi, 12)
    subjects = [np.broadcast_to(np.sin(t)[None, :, None], (4, 12, 2)).copy()
                + rng.normal(0, 0.01, (4, 12, 2)) for _ in range(6)]
    betas = np.stack(subjects, axis=3).transpose(2, 1, 0, 3)
    adjacency = np.ones((4, 4)) - np.eye(4)
    result = limo_gae_predict(betas, adjacency, folds=3, epochs=30, hidden_dim=8,
                              seed=17, device="cpu", dropout=0)
    print(json.dumps({
        "input_shape": list(betas.shape),
        "output_shape": list(result["reconstructed"].shape),
        "finite": bool(np.isfinite(result["reconstructed"]).all()),
        "mean_absolute_error": float(np.abs(betas-result["reconstructed"]).mean()),
        "folds": [{"fold": f["fold"], "train": f["train_indices"], "test": f["test_indices"],
                   "first_loss": f["loss"][0], "last_loss": f["loss"][-1]}
                  for f in result["prediction_info"]["folds"]],
        "versions": result["prediction_info"]["versions"],
    }, indent=2))
