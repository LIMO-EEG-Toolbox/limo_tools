"""CPU tests using small synthetic arrays, not empirical EEG recordings."""
import ast
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np
import torch
from torch import nn
from torch_geometric.nn import GCNConv

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from limo_gae_predict import GraphAE, _mask_inputs, _predict_map, _resolve_device, _train, limo_gae_predict


def fixture(shape=(2, 7, 4, 6)):
    B = np.random.default_rng(19).normal(size=shape).astype(np.float32)
    A = np.ones((shape[2], shape[2]), np.float32) - np.eye(shape[2], dtype=np.float32)
    return B, A


def run(B, A, **kw):
    options = dict(folds=3, epochs=2, hidden_dim=4, device="cpu", seed=11)
    options.update(kw)
    return limo_gae_predict(B, A, **options)


class GAETests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.threads = torch.get_num_threads()
        torch.set_num_threads(1)

    @classmethod
    def tearDownClass(cls):
        torch.set_num_threads(cls.threads)

    def test_wiki_axis_contract(self):
        # Wiki: each participant's beta file is channel x frame x parameter.
        files = [np.fromfunction(lambda c, t, b: 1000*s + 100*b + 10*c + t,
                                 (4, 7, 2), dtype=int) for s in range(6)]
        on_disk = np.stack(files, axis=3)  # channel, frame, parameter, participant
        gae = on_disk.transpose(2, 1, 0, 3)
        self.assertEqual(gae.shape, (2, 7, 4, 6))
        self.assertEqual(gae[1, 5, 3, 4], 4135)
        np.testing.assert_array_equal(gae.transpose(2, 1, 0, 3), on_disk)
        # The same design predicts trial x frame, then Yhat is channel x frame x trial.
        X = np.array([[1, 0], [0, 1], [1, -1]])
        yhat = np.stack([(X @ gae[:, :, c, 0]).T for c in range(4)])
        self.assertEqual(yhat.shape, (4, 7, 3))
        np.testing.assert_array_equal(yhat[3, :, 2], np.full(7, -100))

    def test_shape_folds_metadata_and_no_input_mutation(self):
        B, A = fixture(); before = B.copy(); graph = A.copy()
        result = run(B, A)
        self.assertEqual(result["reconstructed"].shape, B.shape)
        self.assertEqual(result["reconstructed"].dtype, np.float32)
        self.assertTrue(np.isfinite(result["reconstructed"]).all())
        seen = []
        for f in result["prediction_info"]["folds"]:
            self.assertFalse(set(f["train_indices"]) & set(f["test_indices"]))
            self.assertEqual(f["training_maps"], B.shape[0] * len(f["train_indices"]))
            seen.extend(f["test_indices"])
        self.assertEqual(sorted(seen), list(range(B.shape[3])))
        np.testing.assert_array_equal(B, before); np.testing.assert_array_equal(A, graph)
        json.dumps(result["prediction_info"])

    def test_singleton_axes_and_one_subject_in_sample(self):
        B, A = fixture((1, 1, 1, 3))
        self.assertEqual(run(B, A)["reconstructed"].shape, B.shape)
        B, A = fixture((1, 1, 1, 1))
        self.assertEqual(run(B, A, prediction="in-sample")["reconstructed"].shape, B.shape)

    def test_disjoint_graph_batches(self):
        model = GraphAE(3, torch.tensor([[0, 1], [1, 0]]), 2, hidden=4).eval()
        X = torch.arange(12, dtype=torch.float32).reshape(2, 2, 3)
        together = model(X)
        torch.testing.assert_close(together, torch.cat([model(X[i:i+1]) for i in range(2)]))
        X[1] += 1000
        torch.testing.assert_close(together[0], model(X)[0])

    def test_target_mask_and_chunk_order(self):
        class Inspector(nn.Module):
            def forward(self, X):
                # With all original values positive, exactly one zero channel per graph.
                assert torch.all((X == 0).all(2).sum(1) == 1)
                return X + 2
        xs = torch.ones(5, 3)
        torch.testing.assert_close(_predict_map(Inspector(), xs, True, 2), torch.full_like(xs, 2))
        for _ in range(5):
            _, mask = _mask_inputs(xs.unsqueeze(0), 1e-30)
            self.assertGreater(mask.sum().item(), 0)

    def test_held_out_data_cannot_change_model(self):
        B, A = fixture(); labels = [0, 1, 2, 0, 1, 2]
        first = run(B, A, fold_ids=labels, return_models=True)
        changed = B.copy(); changed[:, :, :, [0, 3]] += 50
        second = run(changed, A, fold_ids=labels, return_models=True)
        for name, value in first["models"][0].state_dict().items():
            torch.testing.assert_close(value, second["models"][0].state_dict()[name], rtol=0, atol=0)

    def test_repeatability_and_rng_restoration(self):
        B, A = fixture()
        state = torch.random.get_rng_state().clone()
        settings = torch.are_deterministic_algorithms_enabled()
        first = run(B, A); second = run(B, A)
        np.testing.assert_array_equal(first["reconstructed"], second["reconstructed"])
        torch.testing.assert_close(state, torch.random.get_rng_state(), rtol=0, atol=0)
        self.assertEqual(settings, torch.are_deterministic_algorithms_enabled())

    def test_participant_permutation_with_stable_ids_and_folds(self):
        B, A = fixture(); ids = np.array([f"s{i}" for i in range(6)])
        labels = np.arange(6) % 3; order = [5, 1, 4, 2, 0, 3]
        first = run(B, A, fold_ids=labels, participant_ids=ids)
        second = run(B[:, :, :, order], A, fold_ids=labels[order], participant_ids=ids[order])
        np.testing.assert_allclose(second["reconstructed"], first["reconstructed"][:, :, :, order], rtol=0, atol=0)

    def test_learning_smoke(self):
        B = np.full((2, 5, 4, 6), 0.5, np.float32)
        A = np.ones((4, 4)) - np.eye(4)
        result = run(B, A, epochs=35, learning_rate=0.01, dropout=0)
        for fold in result["prediction_info"]["folds"]:
            self.assertLess(np.mean(fold["loss"][-5:]), np.mean(fold["loss"][:5]))

    def test_matches_source_training_and_deviation(self):
        # Load only the original definitions, avoiding its CLI/environment mutations.
        path = Path(__file__).resolve().parents[2] / "review/gae_source_snapshot/gae_weighting.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        nodes = [node for node in tree.body if isinstance(node, (ast.ClassDef, ast.FunctionDef))
                 and node.name in {"GraphAE", "_train", "_score_map"}]
        env = dict(torch=torch, nn=nn, np=np, GCNConv=GCNConv, HIDDEN=32, EPOCHS=60, MASK_P=0.3)
        exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), env)
        edges = torch.tensor([[0, 1, 1, 2, 2, 3], [1, 0, 2, 1, 3, 2]])
        torch.manual_seed(8); old = env["GraphAE"](7, edges, 4, H=4)
        new = GraphAE(7, edges, 4, hidden=4); new.load_state_dict(old.state_dict())
        X = torch.tensor(fixture()[0].transpose(0, 3, 2, 1).reshape(-1, 4, 7))
        torch.manual_seed(9); env["_train"](old, X, epochs=2)
        torch.manual_seed(9); _train(new, X, epochs=2, masked=True, mask_probability=.3,
                                    batch_size=None, learning_rate=.001, weight_decay=.0001)
        for key, value in old.state_dict().items():
            torch.testing.assert_close(value, new.state_dict()[key], rtol=0, atol=0)
        deviation = (_predict_map(new, X[0], True, 4) - X[0]).abs().mean(1).numpy()
        np.testing.assert_allclose(deviation, env["_score_map"](old, X[0]), rtol=0, atol=0)

    def test_optional_minibatches_and_unmasked(self):
        B, A = fixture()
        result = run(B, A, masked=False, batch_size=3)
        self.assertTrue(np.isfinite(result["reconstructed"]).all())
        self.assertFalse(result["prediction_info"]["masked"])

    def test_rejects_invalid_inputs(self):
        B, A = fixture()
        for kwargs in ({"folds": 7}, {"folds": 1}, {"epochs": 0}, {"mask_probability": 0},
                       {"fold_ids": [0]*6}, {"fold_ids": [0, 1]}, {"participant_ids": ["same"]*6},
                       {"prediction": "unknown"}, {"batch_size": 0}, {"dropout": 1},
                       {"masked": "yes"}, {"seed": -1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                run(B, A, **kwargs)
        for data, graph in [(B[0], A), (B * np.nan, A), (B, A[:2]),
                            (B, A * .5), (B, np.triu(A)), (B.astype(complex), A)]:
            with self.assertRaises(ValueError):
                run(data, graph)

    def test_device_selection_and_overrides(self):
        # Simulate hardware discovery only; do not claim GPU kernel execution.
        with patch("torch.cuda.is_available", return_value=True), patch("torch.cuda.device_count", return_value=2):
            self.assertEqual(str(_resolve_device("auto")), "cuda")
            self.assertEqual(str(_resolve_device("cpu")), "cpu")
            self.assertEqual(str(_resolve_device("cuda:1")), "cuda:1")
            with self.assertRaises(ValueError):
                _resolve_device("cuda:2")
        with patch("torch.cuda.is_available", return_value=False):
            self.assertEqual(str(_resolve_device("auto")), "cpu")
            self.assertEqual(str(_resolve_device("cpu")), "cpu")
            with self.assertRaises(ValueError):
                _resolve_device("cuda")
        for device in ("invalid", "mps", 0, None):
            with self.assertRaises(ValueError):
                _resolve_device(device)

    def test_default_device_cpu_fallback_executes(self):
        B, A = fixture()
        with patch("torch.cuda.is_available", return_value=False):
            # Omit device altogether to exercise the public default.
            result = limo_gae_predict(B, A, folds=3, epochs=1, hidden_dim=4)
        self.assertEqual(result["prediction_info"]["settings"]["device_requested"], "auto")
        self.assertEqual(result["prediction_info"]["settings"]["device"], "cpu")
        self.assertTrue(np.isfinite(result["reconstructed"]).all())


if __name__ == "__main__":
    unittest.main(verbosity=2)
