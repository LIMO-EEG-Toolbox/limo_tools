import numpy as np
import pytest

from limo.limo_glm import _run_glm_handling, _run_glm_bootstrap
from limo.limo_WLS import limo_WLS
from limo.limo_design import flatten_tf, unflatten_tf


def model(x, method="OLS", analysis="Time", level=1, conditions=0, continuous=1):
    return {
        "Analysis": analysis, "Level": level, "Type": "Channels",
        "design": {"X": x, "nb_conditions": conditions, "nb_interactions": 0,
                   "nb_continuous": continuous, "method": method, "status": "to do"},
    }


@pytest.mark.parametrize("method", ["OLS", "WLS", "IRLS"])
@pytest.mark.parametrize("shape", [(2, 12, 60), (1, 12, 60), (2, 1, 60)])
def test_fit_preserves_axes_and_residuals(method, shape):
    rng = np.random.default_rng(42)
    x = np.column_stack((np.linspace(-1, 1, shape[-1]), np.ones(shape[-1])))
    y = rng.normal(size=shape) + 2
    fitted_model, result = _run_glm_handling(model(x, method), {"Yr": y}, variance_estimates="standard")
    assert result["Beta"].shape == (*shape[:-1], 2)
    assert result["Yhat"].shape == shape
    np.testing.assert_allclose(result["Res"], y - result["Yhat"])
    assert np.all(np.isfinite(result["Beta"]))
    expected_weights = shape if method == "IRLS" else (shape[0], shape[-1])
    assert fitted_model["design"]["weights"].shape == expected_weights
    if method == "OLS":
        expected = np.stack([np.linalg.lstsq(x, channel.T, rcond=None)[0].T for channel in y])
        np.testing.assert_allclose(result["Beta"], expected, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("method", ["OLS", "WLS", "IRLS"])
def test_second_level_missing_observations_preserve_subset(method):
    rng = np.random.default_rng(3)
    x = np.column_stack((np.linspace(-1, 1, 65), np.ones(65)))
    y = rng.normal(size=(2, 7, 65)) + 2
    y[0, :, [2, 8, 51]] = np.nan
    fitted_model, result = _run_glm_handling(model(x, method, level=2), {"Yr": y}, variance_estimates="standard")
    keep = np.isfinite(y[0, 0])
    np.testing.assert_allclose(result["Res"][0][:, keep], y[0][:, keep] - result["Yhat"][0][:, keep])
    assert np.all(np.isnan(result["Yhat"][0][:, ~keep]))
    assert np.all(np.isfinite(result["Beta"]))
    if method == "IRLS":
        assert fitted_model["design"]["weights"][0][:, keep].shape == (7, 62)


@pytest.mark.parametrize("level", [1, 2])
def test_irls_bootstrap_retains_frame_trial_weight_axes(level):
    x = np.column_stack((np.linspace(-1, 1, 60), np.ones(60)))
    y = np.random.default_rng(14).normal(size=(2, 3, 60))
    if level == 2:
        y[0, :, [3, 17]] = np.nan
    specification = model(x, "IRLS", level=level)
    specification["design"]["bootstrap"] = 101
    fitted, results = _run_glm_handling(specification, {"Yr": y}, variance_estimates="standard")
    h0 = _run_glm_bootstrap(fitted, results)
    assert h0["H0_Beta"].shape == (2, 3, 2, 101)
    assert np.all(np.isfinite(h0["H0_Beta"]))


def test_wls_other_channels_do_not_change_estimates():
    rng = np.random.default_rng(6)
    x = np.column_stack((np.linspace(-1, 1, 60), np.ones(60)))
    y = rng.normal(size=(2, 12, 60)) + 2
    _, original = _run_glm_handling(model(x, "WLS"), {"Yr": y}, variance_estimates="standard")
    y[1, :, -1] += 30
    _, contaminated = _run_glm_handling(model(x, "WLS"), {"Yr": y}, variance_estimates="standard")
    np.testing.assert_array_equal(original["Beta"][0], contaminated["Beta"][0])


def test_wls_matches_pinned_matlab_reference():
    # Golden betas retained in EEGPrep tests/test_limo_eeglab_tests.py,
    # from MATLAB limo_WLS at bff6d166c5338f05d7ed9b37c1b7615d667492af.
    x = np.column_stack((np.ones(12), np.linspace(-1, 1, 12)))
    y = x @ np.array([[1., 2., 3., 4.], [.5, -1., 1.5, -.25]])
    y[-1] += [5., 4., 6., 3.]
    beta, weights, _ = limo_WLS(x, y)
    np.testing.assert_allclose(beta, [
        [1.00173322768859, 2.00138658215087, 3.00207987322631, 4.00103993661315],
        [.504984652334702, -.996012278132239, 1.50598158280164, -.247009208599179],
    ], rtol=2e-12, atol=2e-12)
    assert weights[-1] == pytest.approx(.04)


@pytest.mark.parametrize("method", ["OLS", "WLS", "IRLS"])
def test_tf_fit_and_effects_use_frequency_fastest_order(method):
    rng = np.random.default_rng(11)
    cat = np.tile([0, 1], 30)
    x = np.column_stack((cat == 0, cat == 1, np.linspace(-1, 1, 60), np.ones(60))).astype(float)
    y = rng.uniform(.5, 3, size=(2, 3, 5, 60))
    fitted_model, result = _run_glm_handling(
        model(x, method, "Time-Frequency", conditions=[2]),
        {"Yr": y}, variance_estimates="standard",
    )
    assert result["Beta"].shape == (2, 3, 5, 4)
    assert result["Condition_effect_1"].shape == (2, 3, 5, 2)
    assert result["Covariate_effect_1"].shape == (2, 3, 5, 2)
    np.testing.assert_allclose(result["Res"], y - result["Yhat"])
    assert np.all(np.isfinite(result["Beta"]))
    if method == "WLS":
        assert fitted_model["design"]["weights"].shape == (2, 3, 60)
        # A TF WLS slice must agree with an independently fitted frequency.
        for freq in range(3):
            _, separate = _run_glm_handling(
                model(x, "WLS", conditions=[2]), {"Yr": y[:, freq]}, variance_estimates="standard",
            )
            np.testing.assert_allclose(result["Beta"][:, freq], separate["Beta"])
            np.testing.assert_allclose(result["Condition_effect_1"][:, freq], separate["Condition_effect_1"])
            np.testing.assert_allclose(result["Covariate_effect_1"][:, freq], separate["Covariate_effect_1"])
    elif method == "IRLS":
        assert fitted_model["design"]["weights"].shape == (2, 15, 60)


def test_tf_reshape_preserves_named_positions():
    y = np.arange(2 * 3 * 5 * 7).reshape(2, 3, 5, 7)
    flat = flatten_tf(y)
    np.testing.assert_array_equal(flat[:, 3 * 2 + 1], y[:, 1, 2])
    np.testing.assert_array_equal(unflatten_tf(flat, y.shape), y)


def test_zero_sum_first_frame_does_not_skip_channel():
    x = np.column_stack((np.linspace(-1, 1, 60), np.ones(60)))
    y = np.random.default_rng(2).normal(size=(2, 7, 60))
    y[:, 0] = np.tile([-1., 1.], 30)
    _, result = _run_glm_handling(model(x), {"Yr": y}, variance_estimates="standard")
    assert np.all(np.isfinite(result["Beta"][:, 1:]))
    np.testing.assert_allclose(result["Beta"][:, 0],
                               np.linalg.lstsq(x, y[:, 0].T, rcond=None)[0].T, atol=1e-12)


def test_tf_wls_multiple_factors_match_separate_frequencies():
    a = np.tile([0, 1], 30)
    b = np.repeat([0, 1], 30)
    x = np.column_stack((a == 0, a == 1, b == 0, b == 1, np.ones(60))).astype(float)
    y = np.random.default_rng(12).uniform(.5, 3, size=(2, 3, 5, 60))
    _, result = _run_glm_handling(model(x, "WLS", "Time-Frequency", conditions=[2, 2], continuous=0),
                                 {"Yr": y}, variance_estimates="standard")
    for freq in range(3):
        _, separate = _run_glm_handling(model(x, "WLS", conditions=[2, 2], continuous=0),
                                       {"Yr": y[:, freq]}, variance_estimates="standard")
        for effect in (1, 2):
            np.testing.assert_allclose(result[f"Condition_effect_{effect}"][:, freq],
                                       separate[f"Condition_effect_{effect}"])
