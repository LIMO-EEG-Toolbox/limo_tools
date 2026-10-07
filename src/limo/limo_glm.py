"""Run the LIMO GLM workflow and save model outputs in HDF5.

This module is the Python counterpart of the MATLAB ``limo_glm*.m`` family:
``limo_glm_handling.m``, ``limo_glm.m``, ``limo_glm_boot.m``, and
``limo_glm_null.m``. The public entry point wraps the file-oriented handling
logic and embeds the core GLM fitting, null-data generation, and bootstrap
under H0 routines in one file.

Primary use:
    Run ``limo_glm(limo_file)`` after ``eeglab_import.py`` and
    ``limo_design.py`` have created ``LIMO.h5`` and ``limo_results.h5``.

Inputs:
    ``limo_file``:
        Path to ``LIMO.h5``.
    ``results_file``:
        Optional path to ``limo_results.h5``. Defaults to ``LIMO.dir`` /
        ``limo_results.h5``.
    ``h0_file``:
        Optional path to the bootstrap-under-H0 output file. Defaults to
        ``LIMO.dir`` / ``limo_H0.h5``.
    ``run_bootstrap``:
        Override whether H0 bootstrap is computed. If omitted, the value in
        ``LIMO.design.bootstrap`` is used.
    ``variance_estimates``:
        Variance estimate mode for the core GLM. ``"standard"`` is the
        documented default. ``"HC4"`` is accepted for the OLS/WLS branch.
    ``tfce_callback``:
        Optional callable used when ``LIMO.design.tfce == 1``. If omitted,
        the module will import the HDF5-aware ``limo_tfce.limo_tfce`` wrapper.

Returned value:
    A tuple ``(updated_limo_path, updated_results_path, h0_path_or_none)`` of
    ``Path`` objects.

Functionality:
    - Loads ``Yr`` and the placeholder outputs created by ``limo_design.py``.
    - Fits the GLM channel-by-channel or component-by-component.
    - Uses ``limo_WLS`` and ``limo_irls`` for robust estimation.
    - Updates ``/Yhat``, ``/Res``, ``/R2``, ``/Beta`` and effect datasets in
      ``limo_results.h5``.
    - Updates ``LIMO.h5`` with model degrees of freedom, trial weights, and
      design completion status.
    - Optionally computes bootstrap distributions under H0 and writes them to
      ``limo_H0.h5``.
        - Optionally dispatches to ``limo_tfce.py`` when the design enables TFCE.

Outputs written to ``limo_results.h5``:
    ``/Yr``:
        The reorganized EEG data used by the model.
    ``/Yhat``:
        The fitted data.
    ``/Res``:
        Residual data.
    ``/R2``:
        Model fit metrics. The last dimension stores ``R2``, ``F``, and ``p``.
    ``/Beta``:
        Beta parameters. For time/frequency analyses this is
        ``channels x frames x parameters``. For time-frequency analyses this
        expands to ``channels x freqs x times x parameters``.
    ``/Condition_effect_*``:
        Factor-wise ``F`` and ``p`` values when categorical effects exist.
    ``/Interaction_effect_*``:
        Interaction ``F`` and ``p`` values when full-factorial interaction
        terms exist.
    ``/Covariate_effect_*``:
        Covariate ``F`` and ``p`` values when continuous predictors exist.

Outputs written to ``limo_H0.h5`` when bootstrapping is enabled:
    ``/H0_R2``, ``/H0_Beta``, ``/boot_table``, and the corresponding H0 effect
    datasets for conditions, interactions, and covariates.

References:
    Christensen, R. (2002). Plane Answers to Complex Questions. 3rd Ed.
    Springer-Verlag.
    Friston, K. et al. (2007). Statistical Parametric Mapping. Academic Press.
    Yandell, B. S. (1997). Practical Data Analysis for Designed Experiments.
    Chapman & Hall.
    Pernet, C. R., et al. LIMO EEG documentation and associated robust GLM
    methods for OLS, WLS, and IRLS.

Notes:
        - TFCE computation is delegated to ``limo_tfce.py`` so the GLM workflow can
            trigger TFCE immediately after writing the observed and H0 statistics.
    - The current implementation follows the ``limo_glm*.m`` statistical path,
      but stores results in HDF5 rather than MATLAB ``.mat`` files.

Command-line usage:
    python limo_glm.py path/to/LIMO.h5
    python limo_glm.py path/to/LIMO.h5 --results-file path/to/limo_results.h5
    python limo_glm.py path/to/LIMO.h5 --run-bootstrap off
"""

from __future__ import annotations

import argparse
import itertools
import json
import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from scipy.stats import f as f_distribution

from .eeglab_import import write_hdf5_structure
from .limo_WLS import limo_WLS
from .limo_irls import limo_irls
from .limo_design import flatten_tf, read_hdf5_structure, unflatten_tf


def limo_glm(
    limo_file: str | Path,
    *,
    results_file: str | Path | None = None,
    h0_file: str | Path | None = None,
    run_bootstrap: bool | None = None,
    variance_estimates: str = "standard",
    tfce_callback: Callable[[Path, Path, Path | None, Mapping[str, Any]], Any] | None = None,
) -> tuple[Path, Path, Path | None]:
    """Run the GLM workflow on ``LIMO.h5`` and update the HDF5 outputs."""

    limo_path = Path(limo_file).expanduser().resolve()
    payload = read_hdf5_structure(limo_path)
    if "LIMO" not in payload or not isinstance(payload["LIMO"], Mapping):
        raise KeyError("The HDF5 file does not contain a top-level 'LIMO' group.")

    limo = dict(payload["LIMO"])
    results_path = (
        Path(results_file).expanduser().resolve()
        if results_file
        else Path(str(limo["data"]["results_file"])).expanduser().resolve()
        if "results_file" in limo.get("data", {})
        else Path(str(limo["dir"])).expanduser().resolve() / "limo_results.h5"
    )
    results_payload = read_hdf5_structure(results_path)

    if str(limo.get("design", {}).get("status", "to do")).lower() == "to do":
        limo, results_payload = _run_glm_handling(limo, results_payload, variance_estimates=variance_estimates)
        write_hdf5_structure(results_path, results_payload)
        write_hdf5_structure(limo_path, {"LIMO": limo})

    do_bootstrap = _resolve_bootstrap_flag(limo, run_bootstrap)
    h0_path: Path | None = None
    if do_bootstrap:
        h0_path = (
            Path(h0_file).expanduser().resolve()
            if h0_file
            else Path(str(limo["dir"])).expanduser().resolve() / "limo_H0.h5"
        )
        h0_payload = _run_glm_bootstrap(limo, results_payload)
        write_hdf5_structure(h0_path, h0_payload)

    _maybe_run_tfce(limo_path, results_path, h0_path, limo, tfce_callback=tfce_callback)

    return limo_path, results_path, h0_path


def limo_glm_handling(
    limo_file: str | Path,
    *,
    results_file: str | Path | None = None,
    h0_file: str | Path | None = None,
    run_bootstrap: bool | None = None,
    variance_estimates: str = "standard",
    tfce_callback: Callable[[Path, Path, Path | None, Mapping[str, Any]], Any] | None = None,
) -> tuple[Path, Path, Path | None]:
    """Compatibility alias for ``limo_glm``."""

    return limo_glm(
        limo_file,
        results_file=results_file,
        h0_file=h0_file,
        run_bootstrap=run_bootstrap,
        variance_estimates=variance_estimates,
        tfce_callback=tfce_callback,
    )


def _maybe_run_tfce(
    limo_path: Path,
    results_path: Path,
    h0_path: Path | None,
    limo: Mapping[str, Any],
    *,
    tfce_callback: Callable[[Path, Path, Path | None, Mapping[str, Any]], Any] | None,
) -> None:
    design = limo.get("design", {}) if isinstance(limo.get("design", {}), Mapping) else {}
    tfce_enabled = int(np.asarray(design.get("tfce", 0)).item()) == 1
    if not tfce_enabled:
        return

    runner = tfce_callback if tfce_callback is not None else _load_tfce_runner()
    if runner is None:
        warnings.warn(
            "LIMO.design.tfce is enabled, but limo_tfce.py could not be imported.",
            RuntimeWarning,
            stacklevel=2,
        )
        return
    runner(limo_path, results_path, h0_path, limo)


def _load_tfce_runner() -> Callable[[Path, Path, Path | None, Mapping[str, Any]], Any] | None:
    try:
        from .limo_tfce import limo_tfce as tfce_runner
    except ImportError:
        return None
    return tfce_runner


def _run_glm_handling(
    limo: Mapping[str, Any],
    results_payload: Mapping[str, Any],
    *,
    variance_estimates: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    updated_limo = _clone_mapping(limo)
    updated_results = _clone_mapping(results_payload)

    yr_native = np.asarray(updated_results["Yr"], dtype=float)
    analysis = str(updated_limo.get("Analysis", ""))
    design = updated_limo.get("design", {}) if isinstance(updated_limo.get("design", {}), Mapping) else {}
    method = str(design.get("method", "OLS"))
    level = int(np.asarray(updated_limo.get("Level", 1)).item())

    if analysis.lower() == "time-frequency":
        n_channels, n_freqs, n_times, n_obs = yr_native.shape
        yr = flatten_tf(yr_native)
    else:
        yr = _ensure_3d(yr_native)
        n_channels, _, n_obs = yr.shape
        n_freqs = n_times = None

    if method.upper() == "IRLS" and n_obs < 50:
        updated_limo["design"]["method"] = "OLS"
        method = "OLS"
        warnings.warn(
            f"With {n_obs} observations detected, IRLS is unlikely to converge; switching to OLS.",
            RuntimeWarning,
            stacklevel=2,
        )

    x_full = np.asarray(design.get("X"), dtype=float)
    nb_conditions = _as_int_list(design.get("nb_conditions"))
    nb_interactions = _as_int_list(design.get("nb_interactions"))
    nb_continuous = int(np.asarray(design.get("nb_continuous", 0)).item())

    working_yhat = np.array(yr, copy=True, dtype=float)
    working_res = np.array(yr, copy=True, dtype=float)
    working_r2 = np.full((yr.shape[0], yr.shape[1], 3), np.nan, dtype=float)
    working_beta = np.full((yr.shape[0], yr.shape[1], x_full.shape[1]), np.nan, dtype=float)

    condition_effect = None
    if _product_or_zero(nb_conditions) != 0:
        condition_effect = np.full((yr.shape[0], yr.shape[1], len(nb_conditions), 2), np.nan, dtype=float)

    interaction_effect = None
    if len(nb_interactions) > 0:
        interaction_effect = np.full((yr.shape[0], yr.shape[1], len(nb_interactions), 2), np.nan, dtype=float)

    covariate_effect = None
    if nb_continuous != 0:
        covariate_effect = np.full((yr.shape[0], yr.shape[1], nb_continuous, 2), np.nan, dtype=float)

    weights = _initialize_weight_container(method=method, analysis=analysis, yr=yr, n_freqs=n_freqs)
    array = _valid_channel_indices(yr, level=level, analysis=analysis)

    model_df = []
    conditions_df = []
    interactions_df = []
    continuous_df = []

    for channel in array:
        if level == 2:
            y_channel = yr[channel]
            index = np.flatnonzero(~np.isnan(y_channel[0]))
            if index.size == 0:
                index = np.arange(y_channel.shape[1])
            y_channel = y_channel[:, index]
            x_channel = x_full[index]
        else:
            y_channel = yr[channel]
            index = np.arange(y_channel.shape[1])
            x_channel = x_full

        model = _fit_limo_glm_model(
            y=y_channel.T,
            x=x_channel,
            nb_conditions=nb_conditions,
            nb_interactions=nb_interactions,
            nb_continuous=nb_continuous,
            method=method,
            analysis=analysis,
            n_freqs=n_freqs,
            n_times=n_times,
            variance_estimates=variance_estimates,
        )

        _collect_model_df(model_df, model.get("df"))
        if "conditions" in model and "df" in model["conditions"]:
            conditions_df.append(np.asarray(model["conditions"]["df"]))
        if "interactions" in model and "df" in model["interactions"]:
            interactions_df.append(np.asarray(model["interactions"]["df"]))
        if "continuous" in model and "df" in model["continuous"]:
            continuous_df.append(np.asarray(model["continuous"]["df"]))

        if analysis.lower() == "time-frequency":
            fitted_flat, beta_flat = _prepare_tf_fitted_outputs(model, x_channel, n_freqs, n_times, method)
            _assign_tf_weights(weights, model.get("W"), channel=channel, index=index, method=method)
        else:
            fitted_flat = x_channel @ np.asarray(model["betas"], dtype=float)
            beta_flat = np.asarray(model["betas"], dtype=float)
            _assign_standard_weights(weights, model.get("W"), channel=channel, index=index, method=method)

        working_yhat[channel][:, index] = fitted_flat.T
        working_res[channel][:, index] = y_channel - fitted_flat.T
        working_r2[channel, :, 0] = np.asarray(model["R2_univariate"], dtype=float)
        working_r2[channel, :, 1] = np.asarray(model["F"], dtype=float)
        working_r2[channel, :, 2] = np.asarray(model["p"], dtype=float)
        working_beta[channel] = beta_flat.T

        if condition_effect is not None and "conditions" in model:
            _assign_effect_array(condition_effect, channel, model["conditions"])
        if interaction_effect is not None and "interactions" in model:
            _assign_effect_array(interaction_effect, channel, model["interactions"])
        if covariate_effect is not None and "continuous" in model:
            _assign_effect_array(covariate_effect, channel, model["continuous"])

    updated_limo = _update_limo_after_glm(
        updated_limo,
        weights=weights,
        model_df=model_df,
        conditions_df=conditions_df,
        interactions_df=interactions_df,
        continuous_df=continuous_df,
    )

    results_update = {
        "Yr": yr_native,
        "Yhat": _restore_analysis_shape(working_yhat, analysis=analysis, n_freqs=n_freqs, n_times=n_times),
        "Res": _restore_analysis_shape(working_res, analysis=analysis, n_freqs=n_freqs, n_times=n_times),
        "R2": _restore_analysis_shape(working_r2, analysis=analysis, n_freqs=n_freqs, n_times=n_times),
        "Beta": _restore_analysis_shape(working_beta, analysis=analysis, n_freqs=n_freqs, n_times=n_times),
    }

    if condition_effect is not None:
        for effect_index in range(condition_effect.shape[2]):
            results_update[f"Condition_effect_{effect_index + 1}"] = _restore_effect_shape(
                condition_effect[:, :, effect_index, :],
                analysis=analysis,
                n_freqs=n_freqs,
                n_times=n_times,
            )
    if interaction_effect is not None:
        for effect_index in range(interaction_effect.shape[2]):
            results_update[f"Interaction_effect_{effect_index + 1}"] = _restore_effect_shape(
                interaction_effect[:, :, effect_index, :],
                analysis=analysis,
                n_freqs=n_freqs,
                n_times=n_times,
            )
    if covariate_effect is not None:
        for effect_index in range(covariate_effect.shape[2]):
            results_update[f"Covariate_effect_{effect_index + 1}"] = _restore_effect_shape(
                covariate_effect[:, :, effect_index, :],
                analysis=analysis,
                n_freqs=n_freqs,
                n_times=n_times,
            )

    updated_results.update(results_update)
    return updated_limo, updated_results


def _fit_limo_glm_model(
    *,
    y: np.ndarray,
    x: np.ndarray,
    nb_conditions: list[int],
    nb_interactions: list[int],
    nb_continuous: int,
    method: str,
    analysis: str,
    n_freqs: int | None,
    n_times: int | None,
    variance_estimates: str,
) -> dict[str, Any]:
    method_upper = method.upper()
    if np.iscomplexobj(y):
        y = np.abs(y) ** 2

    if y.shape[0] != x.shape[0]:
        raise ValueError("The number of events in Y and the design matrix are different")

    if analysis.lower() == "time-frequency" and method_upper == "WLS":
        if n_freqs is None or n_times is None:
            raise ValueError("Time-Frequency WLS requires n_freqs and n_times.")
        return _fit_wls_tf_model(
            y=y,
            x=x,
            nb_conditions=nb_conditions,
            nb_interactions=nb_interactions,
            nb_continuous=nb_continuous,
            n_freqs=n_freqs,
            n_times=n_times,
            variance_estimates=variance_estimates,
        )

    if method_upper == "OLS":
        weights = np.ones((y.shape[0], 1), dtype=float)
        wx = x
        if nb_continuous != 0 and _product_or_zero(nb_conditions) == 0:
            betas = np.linalg.lstsq(wx, y, rcond=None)[0]
        else:
            betas = np.linalg.pinv(wx) @ y
        return _compute_standard_model(
            y=y,
            x=x,
            wx=wx,
            betas=betas,
            weights=weights[:, 0],
            nb_conditions=nb_conditions,
            nb_interactions=nb_interactions,
            nb_continuous=nb_continuous,
            method=method_upper,
            variance_estimates=variance_estimates,
        )

    if method_upper == "WLS":
        betas, weights, _ = limo_WLS(x, y)
        wx = x * weights[:, np.newaxis]
        return _compute_standard_model(
            y=y,
            x=x,
            wx=wx,
            betas=betas,
            weights=weights,
            nb_conditions=nb_conditions,
            nb_interactions=nb_interactions,
            nb_continuous=nb_continuous,
            method=method_upper,
            variance_estimates=variance_estimates,
        )

    if method_upper == "IRLS":
        _, weights = limo_irls(x, y)
        model = _glm_iterate_irls(
            y=y,
            x=x,
            weights=weights,
            nb_conditions=nb_conditions,
            nb_interactions=nb_interactions,
            nb_continuous=nb_continuous,
        )
        model["W"] = weights
        return model

    raise ValueError(f"Unsupported GLM method: {method}")


def _compute_standard_model(
    *,
    y: np.ndarray,
    x: np.ndarray,
    wx: np.ndarray,
    betas: np.ndarray,
    weights: np.ndarray,
    nb_conditions: list[int],
    nb_interactions: list[int],
    nb_continuous: int,
    method: str,
    variance_estimates: str,
) -> dict[str, Any]:
    model: dict[str, Any] = {}
    n_obs = y.shape[0]
    t = (y - np.mean(y, axis=0, keepdims=True)).T @ (y - np.mean(y, axis=0, keepdims=True))
    r = np.eye(n_obs) - wx @ np.linalg.pinv(wx)
    hm = wx @ np.linalg.pinv(wx)
    if variance_estimates.lower() == "hc4":
        h = np.diag(wx @ np.linalg.pinv(wx.T @ wx) @ wx.T)
        d = np.minimum(4, h / np.mean(h))
        residual = r @ y
        hc4 = residual**2 / ((1 - h)[:, np.newaxis] ** d[:, np.newaxis])
        e = np.sum((residual / ((1 - h)[:, np.newaxis] ** d[:, np.newaxis])) ** 2, axis=0)
    else:
        hc4 = None
        e = np.diag(y.T @ r @ y)

    e = np.abs(e)
    df = max(np.linalg.matrix_rank(wx) - 1, 1)
    if method == "OLS":
        dfe = n_obs - np.linalg.matrix_rank(wx)
    else:
        dfe = float(np.trace((np.eye(hm.shape[0]) - hm).T @ (np.eye(hm.shape[0]) - hm)))

    if _product_or_zero(nb_conditions) == 0 and nb_continuous == 0:
        yhat = x @ betas
        h_effect = (yhat - np.mean(yhat, axis=0, keepdims=True)).T @ (y - np.mean(y, axis=0, keepdims=True))
    else:
        c = np.eye(x.shape[1])
        c[:, -1] = 0
        c0 = np.eye(x.shape[1]) - c @ np.linalg.pinv(c)
        x0 = wx @ c0
        r0 = np.eye(n_obs) - x0 @ np.linalg.pinv(x0)
        m = r0 - r
        h_effect = betas.T @ x.T @ m @ x @ betas

    rsquare = np.diag(h_effect) / np.diag(t)
    f_rsquare = (np.diag(h_effect) / df) / (e / dfe)
    p_rsquare = f_distribution.sf(f_rsquare, df, dfe)

    model["W"] = weights
    model["betas"] = betas
    model["betas_se"] = np.array(betas, copy=True)
    for frame in range(y.shape[1]):
        if variance_estimates.lower() == "hc4" and hc4 is not None:
            bread = np.linalg.pinv(wx.T @ wx)
            model["betas_se"][:, frame] = np.diag(bread @ wx.T @ np.diag(hc4[:, frame]) @ wx @ bread)
        else:
            denom = np.sqrt(np.sum((wx - np.mean(wx, axis=0, keepdims=True)) ** 2))
            model["betas_se"][:, frame] = np.sqrt(e[frame] / dfe) / denom
    model["R2_univariate"] = rsquare
    model["F"] = f_rsquare
    model["df"] = np.asarray([df, dfe], dtype=float)
    model["p"] = p_rsquare

    _populate_effects_standard(
        model=model,
        y=y,
        x=x,
        wx=wx,
        betas=betas,
        weights=weights,
        r=r,
        e=e,
        dfe=dfe,
        f_rsquare=f_rsquare,
        p_rsquare=p_rsquare,
        df=df,
        t=t,
        nb_conditions=nb_conditions,
        nb_interactions=nb_interactions,
        nb_continuous=nb_continuous,
    )
    return model


def _populate_effects_standard(
    *,
    model: dict[str, Any],
    y: np.ndarray,
    x: np.ndarray,
    wx: np.ndarray,
    betas: np.ndarray,
    weights: np.ndarray,
    r: np.ndarray,
    e: np.ndarray,
    dfe: float,
    f_rsquare: np.ndarray,
    p_rsquare: np.ndarray,
    df: float,
    t: np.ndarray,
    nb_conditions: list[int],
    nb_interactions: list[int],
    nb_continuous: int,
) -> None:
    nb_factors = len(nb_conditions)
    if nb_factors == 0 and nb_continuous == 0:
        return

    if nb_factors == 1:
        if _product_or_zero(nb_conditions) != 0 and nb_continuous == 0:
            model["conditions"] = {
                "F": f_rsquare,
                "df": np.asarray([df, dfe], dtype=float),
                "p": p_rsquare,
            }
        elif _product_or_zero(nb_conditions) != 0 and nb_continuous != 0:
            c = np.eye(x.shape[1])
            c[:, nb_conditions[0] : x.shape[1]] = 0
            c0 = np.eye(x.shape[1]) - c @ np.linalg.pinv(c)
            x0 = wx @ c0
            r0 = np.eye(y.shape[0]) - x0 @ np.linalg.pinv(x0)
            m = r0 - r
            h = betas.T @ x.T @ m @ x @ betas
            df_conditions = _projection_df(m)
            f_conditions = (np.diag(h) / df_conditions) / (e / dfe)
            p_conditions = f_distribution.sf(f_conditions, df_conditions, dfe)
            model["conditions"] = {
                "F": f_conditions,
                "df": np.asarray([df_conditions, dfe], dtype=float),
                "p": p_conditions,
            }
    elif nb_factors > 1 and len(nb_interactions) == 0:
        model["conditions"] = _compute_multi_factor_main_effects(
            y=y, x=x, wx=wx, betas=betas, r=r, e=e, dfe=dfe, nb_conditions=nb_conditions
        )
    elif nb_factors > 1 and len(nb_interactions) > 0:
        conditions, interactions = _compute_interaction_effects_standard(
            y=y,
            t=t,
            x=x,
            wx=wx,
            betas=betas,
            weights=weights,
            r=r,
            e=e,
            dfe=dfe,
            nb_conditions=nb_conditions,
            nb_interactions=nb_interactions,
            nb_continuous=nb_continuous,
        )
        model["conditions"] = conditions
        model["interactions"] = interactions

    if nb_continuous != 0:
        model["continuous"] = _compute_continuous_effects(
            x=x,
            wx=wx,
            betas=betas,
            r=r,
            e=e,
            dfe=dfe,
            nb_conditions=nb_conditions,
            nb_interactions=nb_interactions,
            nb_continuous=nb_continuous,
        )


def _compute_multi_factor_main_effects(*, y: np.ndarray, x: np.ndarray, wx: np.ndarray, betas: np.ndarray, r: np.ndarray, e: np.ndarray, dfe: float, nb_conditions: list[int]) -> dict[str, Any]:
    df_conditions = np.zeros(len(nb_conditions), dtype=float)
    f_conditions = np.zeros((len(nb_conditions), y.shape[1]), dtype=float)
    p_conditions = np.zeros_like(f_conditions)
    eoi = np.zeros(x.shape[1], dtype=int)
    eoi[: nb_conditions[0]] = np.arange(1, nb_conditions[0] + 1)
    eoni = np.flatnonzero(np.arange(1, x.shape[1] + 1) - eoi)
    for factor in range(len(nb_conditions)):
        c = np.eye(x.shape[1])
        c[:, eoni] = 0
        c0 = np.eye(x.shape[1]) - c @ np.linalg.pinv(c)
        x0 = wx @ c0
        r0 = np.eye(y.shape[0]) - x0 @ np.linalg.pinv(x0)
        m = r0 - r
        h = betas.T @ x.T @ m @ x @ betas
        df_conditions[factor] = _projection_df(m)
        f_conditions[factor] = (np.diag(h) / df_conditions[factor]) / (e / dfe)
        p_conditions[factor] = f_distribution.sf(f_conditions[factor], df_conditions[factor], dfe)
        if factor < len(nb_conditions) - 1:
            update = np.flatnonzero(eoi)[-1]
            eoi = np.zeros(x.shape[1], dtype=int)
            eoi[update + 1 : update + 1 + nb_conditions[factor + 1]] = np.arange(
                update + 2, update + 2 + nb_conditions[factor + 1]
            )
            eoni = np.flatnonzero(np.arange(1, x.shape[1] + 1) - eoi)
    return {"F": f_conditions, "df": np.vstack((df_conditions, np.full_like(df_conditions, dfe))).T, "p": p_conditions}


def _compute_interaction_effects_standard(
    *,
    y: np.ndarray,
    t: np.ndarray,
    x: np.ndarray,
    wx: np.ndarray,
    betas: np.ndarray,
    weights: np.ndarray,
    r: np.ndarray,
    e: np.ndarray,
    dfe: float,
    nb_conditions: list[int],
    nb_interactions: list[int],
    nb_continuous: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    conditions = _compute_multi_factor_main_effects(
        y=y,
        x=np.column_stack((x[:, : sum(nb_conditions)], x[:, sum(nb_conditions) + sum(nb_interactions) : -1], np.ones((x.shape[0], 1)))),
        wx=np.column_stack((wx[:, : sum(nb_conditions)], wx[:, sum(nb_conditions) + sum(nb_interactions) : -1], wx[:, [-1]])),
        betas=np.linalg.pinv(np.column_stack((wx[:, : sum(nb_conditions)], wx[:, sum(nb_conditions) + sum(nb_interactions) : -1], wx[:, [-1]]))) @ y,
        r=np.eye(y.shape[0]) - np.column_stack((wx[:, : sum(nb_conditions)], wx[:, sum(nb_conditions) + sum(nb_interactions) : -1], wx[:, [-1]])) @ np.linalg.pinv(np.column_stack((wx[:, : sum(nb_conditions)], wx[:, sum(nb_conditions) + sum(nb_interactions) : -1], wx[:, [-1]]))),
        e=e,
        dfe=dfe,
        nb_conditions=nb_conditions,
    )

    df_conditions = conditions["df"][:, 0]
    hi = np.full((len(nb_interactions), y.shape[1]), np.nan, dtype=float)
    df_interactions = np.full(len(nb_interactions), np.nan, dtype=float)
    f_interactions = np.full_like(hi, np.nan)
    p_interactions = np.full_like(hi, np.nan)

    covariate_columns = np.arange(sum(nb_conditions) + sum(nb_interactions), x.shape[1] - 1)
    dummy_columns = np.arange(sum(nb_conditions))

    if len(nb_interactions) == 1 and len(nb_conditions) == 2 and nb_continuous == 0:
        hi[0] = np.diag(t) - conditions["F"][0] - conditions["F"][1] - e
        df_interactions[0] = float(np.prod(df_conditions))
        f_interactions[0] = (hi[0] / df_interactions[0]) / (e / dfe)
        p_interactions[0] = f_distribution.sf(f_interactions[0], df_interactions[0], dfe)
    else:
        interaction_map = _interaction_factor_map(len(nb_conditions))
        istart = len(dummy_columns)
        ilowbound = len(dummy_columns)
        main_effects = x[:, dummy_columns]
        cov_and_mean = np.column_stack((x[:, covariate_columns], np.ones((x.shape[0], 1))))
        for interaction_index, width in enumerate(nb_interactions):
            interaction_block = x[:, istart : istart + width]
            current_factors = interaction_map[interaction_index]
            if len(current_factors) == 2:
                x_current = np.column_stack((main_effects, interaction_block, cov_and_mean))
            else:
                isize = sum(nb_interactions[: _first_index_of_factor_size(interaction_map, len(current_factors))])
                ihighbound = len(dummy_columns) + isize
                x_current = np.column_stack((main_effects, x[:, ilowbound:ihighbound], interaction_block, cov_and_mean))
            eoibound = x_current.shape[1] - interaction_block.shape[1] - cov_and_mean.shape[1]
            wx_current = x_current * weights[:, np.newaxis]
            betas_current = np.linalg.pinv(wx_current) @ y
            r_current = np.eye(y.shape[0]) - wx_current @ np.linalg.pinv(wx_current)
            eoi = np.zeros(x_current.shape[1], dtype=int)
            eoi[eoibound : eoibound + width] = np.arange(eoibound + 1, eoibound + 1 + width)
            eoni = np.flatnonzero(np.arange(1, x_current.shape[1] + 1) - eoi)
            c = np.eye(x_current.shape[1])
            c[:, eoni] = 0
            c0 = np.eye(x_current.shape[1]) - c @ np.linalg.pinv(c)
            x0 = wx_current @ c0
            r0 = np.eye(y.shape[0]) - x0 @ np.linalg.pinv(x0)
            m = r0 - r_current
            hi[interaction_index] = np.diag(betas_current.T @ x_current.T @ m @ x_current @ betas_current)
            df_interactions[interaction_index] = float(np.prod(df_conditions[np.asarray(current_factors) - 1]))
            f_interactions[interaction_index] = (hi[interaction_index] / df_interactions[interaction_index]) / (e / dfe)
            p_interactions[interaction_index] = f_distribution.sf(
                f_interactions[interaction_index], df_interactions[interaction_index], dfe
            )
            istart += width
    interactions = {
        "F": f_interactions,
        "df": np.vstack((df_interactions, np.full_like(df_interactions, dfe))).T,
        "p": p_interactions,
    }
    return conditions, interactions


def _compute_continuous_effects(*, x: np.ndarray, wx: np.ndarray, betas: np.ndarray, r: np.ndarray, e: np.ndarray, dfe: float, nb_conditions: list[int], nb_interactions: list[int], nb_continuous: int) -> dict[str, Any]:
    if len(nb_conditions) == 0 and nb_continuous == 1:
        n_frames = betas.shape[1]
        return {
            "F": np.full(n_frames, np.nan, dtype=float),
            "df": np.asarray([1, x.shape[0] - np.linalg.matrix_rank(x)], dtype=float),
            "p": np.full(n_frames, np.nan, dtype=float),
        }

    df_continuous = np.zeros(nb_continuous, dtype=float)
    f_continuous = np.zeros((x.shape[1] if x.ndim == 1 else betas.shape[1], nb_continuous), dtype=float)
    p_continuous = np.zeros_like(f_continuous)
    n_conditions = sum(nb_conditions) + sum(nb_interactions)
    for covariate in range(nb_continuous):
        c = np.zeros((x.shape[1], x.shape[1]))
        c[n_conditions + covariate, n_conditions + covariate] = 1
        c0 = np.eye(x.shape[1]) - c @ np.linalg.pinv(c)
        x0 = wx @ c0
        r0 = np.eye(x.shape[0]) - x0 @ np.linalg.pinv(x0)
        m = r0 - r
        h = betas.T @ x.T @ m @ x @ betas
        df_continuous[covariate] = _projection_df(m)
        f_continuous[:, covariate] = (np.diag(h) / df_continuous[covariate]) / (e / dfe)
        p_continuous[:, covariate] = f_distribution.sf(f_continuous[:, covariate], 1, dfe)
    return {"F": f_continuous, "df": np.asarray([1, dfe], dtype=float), "p": p_continuous}


def _fit_wls_tf_model(
    *,
    y: np.ndarray,
    x: np.ndarray,
    nb_conditions: list[int],
    nb_interactions: list[int],
    nb_continuous: int,
    n_freqs: int,
    n_times: int,
    variance_estimates: str,
) -> dict[str, Any]:
    reshaped = _unpack_tf_trials(y, n_freqs=n_freqs, n_times=n_times)
    betas = np.full((x.shape[1], n_freqs, n_times), np.nan, dtype=float)
    weights = np.full((x.shape[0], n_freqs), np.nan, dtype=float)
    wx_by_freq: list[np.ndarray] = []
    models_by_freq: list[dict[str, Any]] = []
    for freq in range(n_freqs):
        betas_freq, w_freq, _ = limo_WLS(x, reshaped[freq])
        wx_freq = x * w_freq[:, np.newaxis]
        model_freq = _compute_standard_model(
            y=reshaped[freq],
            x=x,
            wx=wx_freq,
            betas=betas_freq,
            weights=w_freq,
            nb_conditions=nb_conditions,
            nb_interactions=nb_interactions,
            nb_continuous=nb_continuous,
            method="WLS",
            variance_estimates=variance_estimates,
        )
        betas[:, freq, :] = betas_freq
        weights[:, freq] = w_freq
        wx_by_freq.append(wx_freq)
        models_by_freq.append(model_freq)

    model = _merge_tf_frequency_models(models_by_freq, betas=betas, weights=weights, n_freqs=n_freqs, n_times=n_times)
    model["WX_by_freq"] = wx_by_freq
    return model


def _merge_tf_frequency_models(models: list[dict[str, Any]], *, betas: np.ndarray, weights: np.ndarray, n_freqs: int, n_times: int) -> dict[str, Any]:
    merged: dict[str, Any] = {
        "W": weights,
        "betas": betas,
        "R2_univariate": _flatten_fortran(np.stack([m["R2_univariate"] for m in models], axis=0)),
        "F": _flatten_fortran(np.stack([m["F"] for m in models], axis=0)),
        "p": _flatten_fortran(np.stack([m["p"] for m in models], axis=0)),
        "df": np.stack([m["df"] for m in models], axis=0),
    }
    if "conditions" in models[0]:
        cond_f = _merge_tf_effect_values(models, "conditions", "F", n_times)
        cond_p = _merge_tf_effect_values(models, "conditions", "p", n_times)
        merged["conditions"] = {
            "F": cond_f,
            "p": cond_p,
            "df": np.stack([m["conditions"]["df"] for m in models], axis=1),
        }
    if "interactions" in models[0]:
        int_f = _merge_tf_effect_values(models, "interactions", "F", n_times)
        int_p = _merge_tf_effect_values(models, "interactions", "p", n_times)
        merged["interactions"] = {
            "F": int_f,
            "p": int_p,
            "df": np.stack([m["interactions"]["df"] for m in models], axis=1),
        }
    if "continuous" in models[0]:
        cont_f = _merge_tf_effect_values(models, "continuous", "F", n_times)
        cont_p = _merge_tf_effect_values(models, "continuous", "p", n_times)
        merged["continuous"] = {
            "F": cont_f,
            "p": cont_p,
            "df": np.stack([m["continuous"]["df"] for m in models], axis=1),
        }
    return merged


def _merge_tf_effect_values(
    models: list[dict[str, Any]], effect: str, statistic: str, n_times: int
) -> np.ndarray:
    """Combine frame-by-effect arrays in frequency-fastest TF order."""
    arrays = []
    for model in models:
        values = np.asarray(model[effect][statistic], dtype=float)
        if values.ndim == 1:
            values = values[:, None]
        elif effect in {"conditions", "interactions"}:
            values = values.T
        if values.ndim != 2 or values.shape[0] != n_times:
            raise ValueError("Expected time-by-effect statistics for each frequency.")
        arrays.append(values)
    stacked = np.stack(arrays, axis=0)
    return stacked.reshape((len(models) * n_times, stacked.shape[-1]), order="F")


def _glm_iterate_irls(
    *,
    y: np.ndarray,
    x: np.ndarray,
    weights: np.ndarray,
    nb_conditions: list[int],
    nb_interactions: list[int],
    nb_continuous: int,
) -> dict[str, Any]:
    nb_factors = len(nb_conditions)
    t = (y - np.mean(y, axis=0, keepdims=True)).T @ (y - np.mean(y, axis=0, keepdims=True))
    betas = np.full((x.shape[1], y.shape[1]), np.nan, dtype=float)
    rsquare = np.full(y.shape[1], np.nan, dtype=float)
    f_rsquare = np.full(y.shape[1], np.nan, dtype=float)
    p_rsquare = np.full(y.shape[1], np.nan, dtype=float)
    dfs = np.full((2, y.shape[1]), np.nan, dtype=float)

    conditions_f = np.full((max(len(nb_conditions), 1), y.shape[1]), np.nan, dtype=float) if nb_factors else None
    conditions_p = np.full_like(conditions_f, np.nan) if conditions_f is not None else None
    interactions_f = np.full((max(len(nb_interactions), 1), y.shape[1]), np.nan, dtype=float) if len(nb_interactions) else None
    interactions_p = np.full_like(interactions_f, np.nan) if interactions_f is not None else None
    continuous_f = np.full((max(nb_continuous, 1), y.shape[1]), np.nan, dtype=float) if nb_continuous else None
    continuous_p = np.full_like(continuous_f, np.nan) if continuous_f is not None else None

    if len(nb_interactions):
        interaction_map = _interaction_factor_map(nb_factors)

    for frame in range(y.shape[1]):
        wx = x * weights[:, frame][:, np.newaxis]
        hm = wx @ np.linalg.pinv(wx)
        r = np.eye(y.shape[0]) - wx @ np.linalg.pinv(wx)
        e = float(y[:, frame].T @ r @ y[:, frame])
        df = _projection_df(hm) - 1
        dfe = float(np.trace((np.eye(hm.shape[0]) - hm).T @ (np.eye(hm.shape[0]) - hm)))
        r_ols = np.eye(y.shape[0]) - x @ np.linalg.pinv(x)
        e_ols = float(y[:, frame].T @ r_ols @ y[:, frame])
        if e < e_ols:
            n = x.shape[0]
            p = np.linalg.matrix_rank(x)
            sigmar = e / (n - p)
            sigmals = e_ols / (n - p)
            mse = (n * sigmar + p**2 * sigmals) / (n + p**2)
            e = mse * dfe

        betas[:, frame] = np.linalg.pinv(wx) @ y[:, frame]
        c = np.eye(x.shape[1])
        c[:, -1] = 0
        c0 = np.eye(x.shape[1]) - c @ np.linalg.pinv(c)
        x0 = wx @ c0
        r0 = np.eye(y.shape[0]) - x0 @ np.linalg.pinv(x0)
        m = r0 - r
        h = float(betas[:, frame].T @ x.T @ m @ x @ betas[:, frame])
        rsquare[frame] = h / t[frame, frame]
        f_rsquare[frame] = (h / df) / (e / dfe)
        p_rsquare[frame] = f_distribution.sf(f_rsquare[frame], df, dfe)
        dfs[:, frame] = np.asarray([df, dfe], dtype=float)

        if nb_factors == 1:
            if _product_or_zero(nb_conditions) != 0 and nb_continuous == 0:
                conditions_f[0, frame] = f_rsquare[frame]
                conditions_p[0, frame] = p_rsquare[frame]
            elif _product_or_zero(nb_conditions) != 0 and nb_continuous != 0:
                c = np.eye(x.shape[1])
                c[:, nb_conditions[0] : x.shape[1]] = 0
                c0 = np.eye(x.shape[1]) - c @ np.linalg.pinv(c)
                x0 = wx @ c0
                r0 = np.eye(y.shape[0]) - x0 @ np.linalg.pinv(x0)
                m = r0 - r
                h_cond = float(betas[:, frame].T @ x.T @ m @ x @ betas[:, frame])
                df_cond = _projection_df(m)
                conditions_f[0, frame] = (h_cond / df_cond) / (e / dfe)
                conditions_p[0, frame] = f_distribution.sf(conditions_f[0, frame], df_cond, dfe)

        elif nb_factors > 1 and len(nb_interactions) == 0:
            eoi = np.zeros(x.shape[1], dtype=int)
            eoi[: nb_conditions[0]] = np.arange(1, nb_conditions[0] + 1)
            eoni = np.flatnonzero(np.arange(1, x.shape[1] + 1) - eoi)
            for factor in range(len(nb_conditions)):
                c = np.eye(x.shape[1])
                c[:, eoni] = 0
                c0 = np.eye(x.shape[1]) - c @ np.linalg.pinv(c)
                x0 = wx @ c0
                r0 = np.eye(y.shape[0]) - x0 @ np.linalg.pinv(x0)
                m = r0 - r
                h_cond = float(betas[:, frame].T @ x.T @ m @ x @ betas[:, frame])
                df_cond = _projection_df(m)
                conditions_f[factor, frame] = (h_cond / df_cond) / (e / dfe)
                conditions_p[factor, frame] = f_distribution.sf(conditions_f[factor, frame], df_cond, dfe)
                if factor < len(nb_conditions) - 1:
                    update = np.flatnonzero(eoi)[-1]
                    eoi = np.zeros(x.shape[1], dtype=int)
                    eoi[update + 1 : update + 1 + nb_conditions[factor + 1]] = np.arange(update + 2, update + 2 + nb_conditions[factor + 1])
                    eoni = np.flatnonzero(np.arange(1, x.shape[1] + 1) - eoi)

        elif nb_factors > 1 and len(nb_interactions) > 0:
            covariate_columns = np.arange(sum(nb_conditions) + sum(nb_interactions), x.shape[1] - 1)
            dummy_columns = np.arange(sum(nb_conditions))
            x_main = np.column_stack((x[:, dummy_columns], x[:, covariate_columns], np.ones((x.shape[0], 1))))
            wx_main = x_main * weights[:, frame][:, np.newaxis]
            betas_main = np.linalg.pinv(wx_main) @ y[:, frame]
            r_main = np.eye(y.shape[0]) - wx_main @ np.linalg.pinv(wx_main)
            eoi = np.zeros(x_main.shape[1], dtype=int)
            eoi[: nb_conditions[0]] = np.arange(1, nb_conditions[0] + 1)
            eoni = np.flatnonzero(np.arange(1, x_main.shape[1] + 1) - eoi)
            df_conditions_frame = np.full(len(nb_conditions), np.nan, dtype=float)
            h_main = np.full(len(nb_conditions), np.nan, dtype=float)
            for factor in range(len(nb_conditions)):
                c = np.eye(x_main.shape[1])
                c[:, eoni] = 0
                c0 = np.eye(x_main.shape[1]) - c @ np.linalg.pinv(c)
                x0 = wx_main @ c0
                r0 = np.eye(y.shape[0]) - x0 @ np.linalg.pinv(x0)
                m = r0 - r_main
                h_main[factor] = float(betas_main.T @ x_main.T @ m @ x_main @ betas_main)
                df_conditions_frame[factor] = _projection_df(m)
                conditions_f[factor, frame] = (h_main[factor] / df_conditions_frame[factor]) / (e / dfe)
                conditions_p[factor, frame] = f_distribution.sf(conditions_f[factor, frame], df_conditions_frame[factor], dfe)
                if factor < len(nb_conditions) - 1:
                    update = np.flatnonzero(eoi)[-1]
                    eoi = np.zeros(x_main.shape[1], dtype=int)
                    eoi[update + 1 : update + 1 + nb_conditions[factor + 1]] = np.arange(update + 2, update + 2 + nb_conditions[factor + 1])
                    eoni = np.flatnonzero(np.arange(1, x_main.shape[1] + 1) - eoi)
            if len(nb_conditions) == 2 and nb_continuous == 0 and len(nb_interactions) == 1:
                hi = t[frame, frame] - np.sum(h_main) - e
                df_i = float(np.prod(df_conditions_frame))
                interactions_f[0, frame] = (hi / df_i) / (e / dfe)
                interactions_p[0, frame] = f_distribution.sf(interactions_f[0, frame], df_i, dfe)
            else:
                main_effects = x[:, dummy_columns]
                cov_and_mean = np.column_stack((x[:, covariate_columns], np.ones((x.shape[0], 1))))
                istart = len(dummy_columns)
                ilowbound = len(dummy_columns)
                for interaction_index, width in enumerate(nb_interactions):
                    interaction_block = x[:, istart : istart + width]
                    factors = interaction_map[interaction_index]
                    if len(factors) == 2:
                        x_current = np.column_stack((main_effects, interaction_block, cov_and_mean))
                    else:
                        isize = sum(nb_interactions[: _first_index_of_factor_size(interaction_map, len(factors))])
                        ihighbound = len(dummy_columns) + isize
                        x_current = np.column_stack((main_effects, x[:, ilowbound:ihighbound], interaction_block, cov_and_mean))
                    wx_current = x_current * weights[:, frame][:, np.newaxis]
                    betas_current = np.linalg.pinv(wx_current) @ y[:, frame]
                    r_current = np.eye(y.shape[0]) - wx_current @ np.linalg.pinv(wx_current)
                    eoibound = x_current.shape[1] - interaction_block.shape[1] - cov_and_mean.shape[1]
                    eoi = np.zeros(x_current.shape[1], dtype=int)
                    eoi[eoibound : eoibound + width] = np.arange(eoibound + 1, eoibound + 1 + width)
                    eoni = np.flatnonzero(np.arange(1, x_current.shape[1] + 1) - eoi)
                    c = np.eye(x_current.shape[1])
                    c[:, eoni] = 0
                    c0 = np.eye(x_current.shape[1]) - c @ np.linalg.pinv(c)
                    x0 = wx_current @ c0
                    r0 = np.eye(y.shape[0]) - x0 @ np.linalg.pinv(x0)
                    m = r0 - r_current
                    hi = float(betas_current.T @ x_current.T @ m @ x_current @ betas_current)
                    df_i = float(np.prod(df_conditions_frame[np.asarray(factors) - 1]))
                    interactions_f[interaction_index, frame] = (hi / df_i) / (e / dfe)
                    interactions_p[interaction_index, frame] = f_distribution.sf(interactions_f[interaction_index, frame], df_i, dfe)
                    istart += width

        if nb_continuous != 0:
            n_conditions = sum(nb_conditions) + sum(nb_interactions)
            for covariate in range(nb_continuous):
                c = np.zeros((x.shape[1], x.shape[1]))
                c[n_conditions + covariate, n_conditions + covariate] = 1
                c0 = np.eye(x.shape[1]) - c @ np.linalg.pinv(c)
                x0 = wx @ c0
                r0 = np.eye(y.shape[0]) - x0 @ np.linalg.pinv(x0)
                m = r0 - r
                h_cont = float(betas[:, frame].T @ x.T @ m @ x @ betas[:, frame])
                df_cont = _projection_df(m)
                continuous_f[covariate, frame] = (h_cont / df_cont) / (e / dfe)
                continuous_p[covariate, frame] = f_distribution.sf(continuous_f[covariate, frame], 1, dfe)

    model = {
        "betas": betas,
        "R2_univariate": rsquare,
        "F": f_rsquare,
        "p": p_rsquare,
        "df": dfs.T,
    }
    if conditions_f is not None:
        model["conditions"] = {"F": conditions_f if conditions_f.shape[0] > 1 else conditions_f[0], "df": np.nan, "p": conditions_p if conditions_p.shape[0] > 1 else conditions_p[0]}
    if interactions_f is not None:
        model["interactions"] = {"F": interactions_f if interactions_f.shape[0] > 1 else interactions_f[0], "df": np.nan, "p": interactions_p if interactions_p.shape[0] > 1 else interactions_p[0]}
    if continuous_f is not None:
        model["continuous"] = {"F": continuous_f.T if continuous_f.shape[0] > 1 else continuous_f[0], "df": np.nan, "p": continuous_p.T if continuous_p.shape[0] > 1 else continuous_p[0]}
    return model


def _run_glm_bootstrap(limo: Mapping[str, Any], results_payload: Mapping[str, Any]) -> dict[str, Any]:
    analysis = str(limo.get("Analysis", ""))
    design = limo.get("design", {}) if isinstance(limo.get("design", {}), Mapping) else {}
    x_full = np.asarray(design.get("X"), dtype=float)
    yr_native = np.asarray(results_payload["Yr"], dtype=float)
    yr = flatten_tf(yr_native) if analysis.lower() == "time-frequency" else _ensure_3d(yr_native)
    n_freqs = yr_native.shape[1] if analysis.lower() == "time-frequency" else None
    n_times = yr_native.shape[2] if analysis.lower() == "time-frequency" else None
    array = _valid_channel_indices(yr, level=int(np.asarray(limo.get("Level", 1)).item()), analysis=analysis)
    nb_conditions = _as_int_list(design.get("nb_conditions"))
    nb_interactions = _as_int_list(design.get("nb_interactions"))
    nb_continuous = int(np.asarray(design.get("nb_continuous", 0)).item())
    method = str(design.get("method", "OLS"))
    nboot = int(np.asarray(design.get("bootstrap", 0)).item())
    if nboot < 800 and nboot != 101:
        nboot = 800

    if int(np.asarray(limo.get("Level", 1)).item()) == 2:
        boot_table = _limo_create_boot_table(yr_native[:, 0, :, :] if analysis.lower() == "time-frequency" else yr_native, nboot)
    else:
        n_obs = yr_native.shape[-1]
        boot_table = np.random.randint(0, n_obs, size=(n_obs, nboot))

    beta_shape = (yr.shape[0], yr.shape[1], x_full.shape[1], nboot)
    r2_shape = (yr.shape[0], yr.shape[1], 3, nboot)
    h0_beta = np.full(beta_shape, np.nan, dtype=float)
    h0_r2 = np.full(r2_shape, np.nan, dtype=float)
    h0_payload: dict[str, Any] = {"H0_Beta": h0_beta, "H0_R2": h0_r2, "boot_table": boot_table}

    if _product_or_zero(nb_conditions) != 0:
        h0_conditions = np.full((yr.shape[0], yr.shape[1], len(nb_conditions), 2, nboot), np.nan, dtype=float)
    else:
        h0_conditions = None
    if len(nb_interactions) > 0:
        h0_interactions = np.full((yr.shape[0], yr.shape[1], len(nb_interactions), 2, nboot), np.nan, dtype=float)
    else:
        h0_interactions = None
    if nb_continuous != 0:
        h0_covariates = np.full((yr.shape[0], yr.shape[1], nb_continuous, 2, nboot), np.nan, dtype=float)
    else:
        h0_covariates = None

    weights = np.asarray(design.get("weights"), dtype=float)
    for channel in array:
        if int(np.asarray(limo.get("Level", 1)).item()) == 2:
            y_channel = yr[channel]
            index = np.flatnonzero(~np.isnan(y_channel[0]))
            if index.size == 0:
                index = np.arange(y_channel.shape[1])
            x_channel = x_full[index]
            y_channel = y_channel[:, index]
            boot_indices = boot_table[channel]
        else:
            index = np.arange(yr.shape[2])
            x_channel = x_full
            y_channel = yr[channel]
            boot_indices = boot_table

        if analysis.lower() == "time-frequency":
            if n_freqs is None or n_times is None:
                raise ValueError("Time-Frequency bootstrap requires native frequency and time sizes.")
            y_channel_tf = _unpack_tf_trials(y_channel.T, n_freqs=n_freqs, n_times=n_times)
            if method.upper() in {"WLS", "OLS"}:
                for freq in range(n_freqs):
                    current_weights = np.ones(index.size) if method.upper() == "OLS" else np.asarray(weights[channel, freq, index], dtype=float)
                    model = _limo_glm_boot(
                        y=y_channel_tf[freq],
                        x=x_channel,
                        weights=current_weights,
                        nb_conditions=nb_conditions,
                        nb_interactions=nb_interactions,
                        nb_continuous=nb_continuous,
                        method=method.upper(),
                        boot_table=boot_indices,
                    )
                    _fill_h0_from_model_freq(
                        model,
                        channel=channel,
                        freq=freq,
                        n_freqs=n_freqs,
                        nboot=nboot,
                        h0_beta=h0_beta,
                        h0_r2=h0_r2,
                        h0_conditions=h0_conditions,
                        h0_interactions=h0_interactions,
                        h0_covariates=h0_covariates,
                    )
            else:
                weights_tf = _reshape_tf_weights_for_boot(weights[channel], n_freqs, n_times, index.size)
                for freq in range(n_freqs):
                    model = _limo_glm_boot(
                        y=y_channel_tf[freq],
                        x=x_channel,
                        weights=weights_tf[freq].T,
                        nb_conditions=nb_conditions,
                        nb_interactions=nb_interactions,
                        nb_continuous=nb_continuous,
                        method=method.upper(),
                        boot_table=boot_indices,
                    )
                    _fill_h0_from_model_freq(
                        model,
                        channel=channel,
                        freq=freq,
                        n_freqs=n_freqs,
                        nboot=nboot,
                        h0_beta=h0_beta,
                        h0_r2=h0_r2,
                        h0_conditions=h0_conditions,
                        h0_interactions=h0_interactions,
                        h0_covariates=h0_covariates,
                    )
        else:
            if method.upper() in {"WLS", "OLS"}:
                current_weights = np.ones(index.size) if method.upper() == "OLS" else np.asarray(weights[channel, index], dtype=float)
            else:
                current_weights = np.asarray(weights[channel][:, index], dtype=float)
            model = _limo_glm_boot(
                y=y_channel.T,
                x=x_channel,
                weights=current_weights,
                nb_conditions=nb_conditions,
                nb_interactions=nb_interactions,
                nb_continuous=nb_continuous,
                method=method.upper(),
                boot_table=boot_indices,
            )
            _fill_h0_from_model(
                model,
                channel=channel,
                nboot=nboot,
                h0_beta=h0_beta,
                h0_r2=h0_r2,
                h0_conditions=h0_conditions,
                h0_interactions=h0_interactions,
                h0_covariates=h0_covariates,
            )

    h0_payload["H0_Beta"] = _restore_analysis_shape(h0_beta, analysis=analysis, n_freqs=n_freqs, n_times=n_times)
    h0_payload["H0_R2"] = _restore_analysis_shape(h0_r2, analysis=analysis, n_freqs=n_freqs, n_times=n_times)
    if h0_conditions is not None:
        for effect_index in range(h0_conditions.shape[2]):
            h0_payload[f"H0_Condition_effect_{effect_index + 1}"] = _restore_effect_shape(h0_conditions[:, :, effect_index, :, :], analysis=analysis, n_freqs=n_freqs, n_times=n_times)
    if h0_interactions is not None:
        for effect_index in range(h0_interactions.shape[2]):
            h0_payload[f"H0_Interaction_effect_{effect_index + 1}"] = _restore_effect_shape(h0_interactions[:, :, effect_index, :, :], analysis=analysis, n_freqs=n_freqs, n_times=n_times)
    if h0_covariates is not None:
        for effect_index in range(h0_covariates.shape[2]):
            h0_payload[f"H0_Covariate_effect_{effect_index + 1}"] = _restore_effect_shape(h0_covariates[:, :, effect_index, :, :], analysis=analysis, n_freqs=n_freqs, n_times=n_times)
    return h0_payload


def _limo_glm_boot(
    *,
    y: np.ndarray,
    x: np.ndarray,
    weights: np.ndarray,
    nb_conditions: list[int],
    nb_interactions: list[int],
    nb_continuous: int,
    method: str,
    boot_table: Any,
) -> dict[str, Any]:
    if isinstance(boot_table, list):
        raise ValueError("Channel-specific boot tables should be resolved before calling _limo_glm_boot.")
    boot_table = np.asarray(boot_table, dtype=int)
    centered_y = _limo_glm_null(y, x, nb_conditions, nb_interactions)
    nboot = boot_table.shape[1]
    model: dict[str, Any] = {"betas": [], "R2_univariate": [], "F": [], "p": []}
    if len(nb_conditions):
        model["conditions"] = {"F": [], "p": []}
    if len(nb_interactions):
        model["interactions"] = {"F": [], "p": []}
    if nb_continuous:
        model["continuous"] = {"F": [], "p": []}

    for boot in range(nboot):
        y_boot = centered_y[boot_table[:, boot]]
        if method == "OLS":
            fitted = _fit_limo_glm_model(
                y=y_boot,
                x=x,
                nb_conditions=nb_conditions,
                nb_interactions=nb_interactions,
                nb_continuous=nb_continuous,
                method="OLS",
                analysis="Time",
                n_freqs=None,
                n_times=None,
                variance_estimates="standard",
            )
        elif method == "WLS":
            w_boot = np.asarray(weights)[boot_table[:, boot]]
            wx = x * w_boot[:, np.newaxis]
            betas = np.linalg.pinv(wx) @ (y_boot * w_boot[:, np.newaxis])
            fitted = _compute_standard_model(
                y=y_boot,
                x=x,
                wx=wx,
                betas=betas,
                weights=w_boot,
                nb_conditions=nb_conditions,
                nb_interactions=nb_interactions,
                nb_continuous=nb_continuous,
                method="WLS",
                variance_estimates="standard",
            )
        else:
            w_boot = np.asarray(weights)[:, boot_table[:, boot]].T
            fitted = _glm_iterate_irls(
                y=y_boot,
                x=x,
                weights=w_boot,
                nb_conditions=nb_conditions,
                nb_interactions=nb_interactions,
                nb_continuous=nb_continuous,
            )

        model["betas"].append(np.asarray(fitted["betas"]).T)
        model["R2_univariate"].append(np.asarray(fitted["R2_univariate"]))
        model["F"].append(np.asarray(fitted["F"]))
        model["p"].append(np.asarray(fitted["p"]))
        if "conditions" in model and "conditions" in fitted:
            model["conditions"]["F"].append(np.asarray(fitted["conditions"]["F"]))
            model["conditions"]["p"].append(np.asarray(fitted["conditions"]["p"]))
        if "interactions" in model and "interactions" in fitted:
            model["interactions"]["F"].append(np.asarray(fitted["interactions"]["F"]))
            model["interactions"]["p"].append(np.asarray(fitted["interactions"]["p"]))
        if "continuous" in model and "continuous" in fitted:
            model["continuous"]["F"].append(np.asarray(fitted["continuous"]["F"]))
            model["continuous"]["p"].append(np.asarray(fitted["continuous"]["p"]))
    return model


def _limo_glm_null(y: np.ndarray, x: np.ndarray, nb_conditions: list[int], nb_interactions: list[int]) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    if _product_or_zero(nb_conditions) == 0:
        return y[np.random.permutation(y.shape[0])]

    null_y = np.full_like(y, np.nan)
    if len(nb_interactions) > 0:
        start_at = sum(nb_conditions) if len(nb_interactions) == 1 else sum(nb_conditions) + sum(nb_interactions[:-1])
        for cell in range(start_at, start_at + nb_interactions[-1]):
            index = np.flatnonzero(x[:, cell])
            null_y[index] = y[index] - np.mean(y[index], axis=0, keepdims=True)
    elif len(nb_conditions) == 1:
        for cell in range(nb_conditions[0]):
            index = np.flatnonzero(x[:, cell])
            null_y[index] = y[index] - np.mean(y[index], axis=0, keepdims=True)
    else:
        tmp_x, interactions = _make_interactions(x[:, : sum(nb_conditions)], nb_conditions)
        start_at = sum(nb_conditions) if len(interactions) == 1 else sum(nb_conditions) + sum(interactions[:-1])
        for cell in range(start_at, start_at + interactions[-1]):
            index = np.flatnonzero(tmp_x[:, cell])
            null_y[index] = y[index] - np.mean(y[index], axis=0, keepdims=True)
    return null_y


def _limo_create_boot_table(data: np.ndarray, nboot: int) -> Any:
    nmin = 3
    data = np.asarray(data, dtype=float)
    chdata = data[0, 0, :] if data.shape[0] == 1 else np.squeeze(data[:, 0, :])
    if data.shape[-1] - 1 <= nmin:
        raise ValueError(f"Not enough subjects in dataset - need at least {nmin} subjects")

    boot_index = np.zeros((data.shape[-1], nboot), dtype=int)
    b = 0
    while b < nboot:
        tmp = np.random.randint(0, data.shape[-1], size=(data.shape[-1],))
        if np.unique(tmp).size >= nmin:
            boot_index[:, b] = tmp
            b += 1

    if data.shape[0] > 1:
        array = np.flatnonzero(np.sum(np.isnan(np.squeeze(data[:, 0, :])), axis=1) < data.shape[-1] - 3)
    else:
        array = np.asarray([0])

    boot_table: list[np.ndarray | None] = [None] * data.shape[0]
    for channel in array:
        tmp = np.squeeze(data[channel])
        bad_subjects = np.flatnonzero(np.isnan(tmp[0]))
        good_subjects = np.flatnonzero(~np.isnan(tmp[0]))
        y = tmp[:, good_subjects]
        if bad_subjects.size:
            boot_index2 = np.zeros((y.shape[1], nboot), dtype=int)
            for column in range(nboot):
                common = np.isin(boot_index[:, column], good_subjects)
                current = boot_index[common, column]
                add = y.shape[1] - current.shape[0]
                if add > 0:
                    new_boot = np.concatenate((current, good_subjects[np.random.randint(0, good_subjects.size, size=add)]))
                else:
                    new_boot = current[: y.shape[1]]
                tmp_boot = new_boot.copy()
                for i, bad in enumerate(bad_subjects):
                    new_boot[tmp_boot > bad] = tmp_boot[tmp_boot > bad] - (i + 1)
                boot_index2[:, column] = new_boot
            boot_table[channel] = boot_index2
        else:
            boot_table[channel] = boot_index
    for index, entry in enumerate(boot_table):
        if entry is None:
            boot_table[index] = boot_index
    return boot_table


def _prepare_tf_fitted_outputs(model: Mapping[str, Any], x: np.ndarray, n_freqs: int | None, n_times: int | None, method: str) -> tuple[np.ndarray, np.ndarray]:
    if n_freqs is None or n_times is None:
        raise ValueError("Time-frequency fitted output requires n_freqs and n_times.")
    method_upper = method.upper()
    if method_upper == "IRLS":
        fitted = np.full((x.shape[0], n_freqs * n_times), np.nan, dtype=float)
        betas = np.asarray(model["betas"], dtype=float)
        weights = np.asarray(model["W"], dtype=float)
        for ft in range(weights.shape[1]):
            wx = x * weights[:, ft][:, np.newaxis]
            fitted[:, ft] = wx @ betas[:, ft]
        return fitted, betas
    if method_upper == "WLS":
        fitted_4d = np.full((1, n_freqs, n_times, x.shape[0]), np.nan, dtype=float)
        betas = np.asarray(model["betas"], dtype=float)
        weights = np.asarray(model["W"], dtype=float)
        for freq in range(n_freqs):
            wx = x * weights[:, freq][:, np.newaxis]
            fitted_4d[0, freq, :, :] = (wx @ betas[:, freq, :]).T
        fitted_flat = flatten_tf(fitted_4d)[0].T
        beta_flat = np.zeros((betas.shape[0], n_freqs * n_times), dtype=float)
        for parameter in range(betas.shape[0]):
            beta_flat[parameter] = np.reshape(betas[parameter], (n_freqs * n_times,), order="F")
        return fitted_flat, beta_flat
    betas = np.asarray(model["betas"], dtype=float)
    return x @ betas, betas


def _assign_standard_weights(weights: np.ndarray, model_weights: Any, *, channel: int, index: np.ndarray, method: str) -> None:
    method_upper = method.upper()
    if method_upper == "WLS":
        weights[channel, index] = np.asarray(model_weights, dtype=float)
    elif method_upper == "IRLS":
        weights[channel][:, index] = np.asarray(model_weights, dtype=float).T


def _assign_tf_weights(weights: np.ndarray, model_weights: Any, *, channel: int, index: np.ndarray, method: str) -> None:
    method_upper = method.upper()
    if method_upper == "WLS":
        weights[channel][:, index] = np.asarray(model_weights, dtype=float).T
    elif method_upper == "IRLS":
        weights[channel][:, index] = np.asarray(model_weights, dtype=float).T


def _assign_effect_array(target: np.ndarray, channel: int, effect_model: Mapping[str, Any]) -> None:
    f_values = np.asarray(effect_model["F"], dtype=float)
    p_values = np.asarray(effect_model["p"], dtype=float)
    if f_values.ndim == 0:
        target[channel, :, 0, 0] = f_values
        target[channel, :, 0, 1] = p_values
    elif f_values.ndim == 1:
        target[channel, :, 0, 0] = f_values
        target[channel, :, 0, 1] = p_values
    else:
        n_frames = target.shape[1]
        n_effects = target.shape[2]
        if f_values.shape == (n_frames, n_effects):
            for effect_index in range(n_effects):
                target[channel, :, effect_index, 0] = f_values[:, effect_index]
                target[channel, :, effect_index, 1] = p_values[:, effect_index]
        elif f_values.shape == (n_effects, n_frames):
            for effect_index in range(n_effects):
                target[channel, :, effect_index, 0] = f_values[effect_index]
                target[channel, :, effect_index, 1] = p_values[effect_index]
        else:
            raise ValueError("Unexpected effect array shape while storing GLM outputs.")


def _update_limo_after_glm(
    limo: Mapping[str, Any],
    *,
    weights: np.ndarray,
    model_df: list[np.ndarray],
    conditions_df: list[np.ndarray],
    interactions_df: list[np.ndarray],
    continuous_df: list[np.ndarray],
) -> dict[str, Any]:
    updated = _clone_mapping(limo)
    updated_design = dict(updated.get("design", {}))
    updated_design["weights"] = weights
    updated_design["status"] = "done"
    if "name" not in updated_design:
        updated_design["name"] = "GLM"
    updated["design"] = updated_design
    updated_model = dict(updated.get("model", {}))
    if model_df:
        updated_model["model_df"] = _stack_or_list(model_df)
    if conditions_df:
        updated_model["conditions_df"] = _stack_or_list(conditions_df)
    if interactions_df:
        updated_model["interactions_df"] = _stack_or_list(interactions_df)
    if continuous_df:
        updated_model["continuous_df"] = _stack_or_list(continuous_df)
    if updated_model:
        updated["model"] = updated_model
    return updated


def _initialize_weight_container(*, method: str, analysis: str, yr: np.ndarray, n_freqs: int | None) -> np.ndarray:
    method_upper = method.upper()
    if method_upper in {"WLS", "OLS"}:
        if analysis.lower() == "time-frequency":
            if n_freqs is None:
                raise ValueError("Time-Frequency weights require n_freqs.")
            return np.ones((yr.shape[0], n_freqs, yr.shape[2]), dtype=float)
        return np.ones((yr.shape[0], yr.shape[2]), dtype=float)
    return np.ones_like(yr, dtype=float)


def _valid_channel_indices(yr: np.ndarray, *, level: int, analysis: str) -> np.ndarray:
    if yr.shape[0] == 1:
        return np.asarray([0])
    if level == 2:
        return np.arange(yr.shape[0])
    return np.flatnonzero(~np.isnan(yr[:, 0, 0]))


def _resolve_bootstrap_flag(limo: Mapping[str, Any], requested: bool | None) -> bool:
    if requested is not None:
        return bool(requested)
    design = limo.get("design", {}) if isinstance(limo.get("design", {}), Mapping) else {}
    return int(np.asarray(design.get("bootstrap", 0)).item()) != 0


def _restore_analysis_shape(data: np.ndarray, *, analysis: str, n_freqs: int | None, n_times: int | None) -> np.ndarray:
    if analysis.lower() != "time-frequency":
        return data
    if n_freqs is None or n_times is None:
        raise ValueError("Time-Frequency restore requires n_freqs and n_times.")
    original_shape = data.shape
    n_channels = original_shape[0]
    trailing = int(np.prod(original_shape[2:])) if data.ndim > 3 else original_shape[2]
    flat = data.reshape((n_channels, original_shape[1], trailing))
    restored = unflatten_tf(flat, (n_channels, n_freqs, n_times, trailing))
    return restored.reshape((n_channels, n_freqs, n_times) + original_shape[2:])


def _restore_effect_shape(data: np.ndarray, *, analysis: str, n_freqs: int | None, n_times: int | None) -> np.ndarray:
    if analysis.lower() != "time-frequency":
        return data
    if n_freqs is None or n_times is None:
        raise ValueError("Time-Frequency effect restore requires n_freqs and n_times.")
    n_channels = data.shape[0]
    trailing = int(np.prod(data.shape[2:]))
    flat = data.reshape((n_channels, data.shape[1], trailing))
    restored = unflatten_tf(flat, (n_channels, n_freqs, n_times, trailing))
    return restored.reshape((n_channels, n_freqs, n_times) + data.shape[2:])


def _unpack_tf_trials(y: np.ndarray, *, n_freqs: int, n_times: int) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    if y.shape[1] != n_freqs * n_times:
        raise ValueError("dimensions disagreement to reshape freq*time")
    reshaped = np.full((n_freqs, y.shape[0], n_times), np.nan, dtype=float)
    for trial in range(y.shape[0]):
        for time_index in range(n_times):
            start = time_index * n_freqs
            stop = start + n_freqs
            reshaped[:, trial, time_index] = y[trial, start:stop]
    return reshaped


def _reshape_tf_weights_for_boot(weights: np.ndarray, n_freqs: int, n_times: int, n_obs: int) -> np.ndarray:
    weights = np.asarray(weights, dtype=float)
    if weights.shape == (n_freqs * n_times, n_obs):
        reshaped = np.full((n_freqs, n_times, n_obs), np.nan, dtype=float)
        for time_index in range(n_times):
            start = time_index * n_freqs
            stop = start + n_freqs
            reshaped[:, time_index, :] = weights[start:stop]
        return reshaped
    if weights.shape == (n_freqs, n_obs):
        return np.repeat(weights[:, np.newaxis, :], n_times, axis=1)
    raise ValueError("Unexpected TF weight shape for bootstrap.")


def _fill_h0_from_model(model: Mapping[str, Any], *, channel: int, nboot: int, h0_beta: np.ndarray, h0_r2: np.ndarray, h0_conditions: np.ndarray | None, h0_interactions: np.ndarray | None, h0_covariates: np.ndarray | None) -> None:
    for boot in range(nboot):
        h0_beta[channel, :, :, boot] = np.asarray(model["betas"][boot], dtype=float)
        h0_r2[channel, :, 0, boot] = np.asarray(model["R2_univariate"][boot], dtype=float)
        h0_r2[channel, :, 1, boot] = np.asarray(model["F"][boot], dtype=float)
        h0_r2[channel, :, 2, boot] = np.asarray(model["p"][boot], dtype=float)
        if h0_conditions is not None and "conditions" in model:
            f_values = np.asarray(model["conditions"]["F"][boot], dtype=float)
            p_values = np.asarray(model["conditions"]["p"][boot], dtype=float)
            if f_values.ndim == 1:
                h0_conditions[channel, :, 0, 0, boot] = f_values
                h0_conditions[channel, :, 0, 1, boot] = p_values
            else:
                for effect in range(f_values.shape[0]):
                    h0_conditions[channel, :, effect, 0, boot] = f_values[effect]
                    h0_conditions[channel, :, effect, 1, boot] = p_values[effect]
        if h0_interactions is not None and "interactions" in model:
            f_values = np.asarray(model["interactions"]["F"][boot], dtype=float)
            p_values = np.asarray(model["interactions"]["p"][boot], dtype=float)
            if f_values.ndim == 1:
                h0_interactions[channel, :, 0, 0, boot] = f_values
                h0_interactions[channel, :, 0, 1, boot] = p_values
            else:
                for effect in range(f_values.shape[0]):
                    h0_interactions[channel, :, effect, 0, boot] = f_values[effect]
                    h0_interactions[channel, :, effect, 1, boot] = p_values[effect]
        if h0_covariates is not None and "continuous" in model:
            f_values = np.asarray(model["continuous"]["F"][boot], dtype=float)
            p_values = np.asarray(model["continuous"]["p"][boot], dtype=float)
            if f_values.ndim == 1:
                h0_covariates[channel, :, 0, 0, boot] = f_values
                h0_covariates[channel, :, 0, 1, boot] = p_values
            else:
                for effect in range(f_values.shape[-1]):
                    h0_covariates[channel, :, effect, 0, boot] = f_values[:, effect]
                    h0_covariates[channel, :, effect, 1, boot] = p_values[:, effect]


def _fill_h0_from_model_freq(model: Mapping[str, Any], *, channel: int, freq: int, n_freqs: int, nboot: int, h0_beta: np.ndarray, h0_r2: np.ndarray, h0_conditions: np.ndarray | None, h0_interactions: np.ndarray | None, h0_covariates: np.ndarray | None) -> None:
    for boot in range(nboot):
        beta_values = np.asarray(model["betas"][boot], dtype=float)
        r2_values = np.asarray(model["R2_univariate"][boot], dtype=float)
        f_values_model = np.asarray(model["F"][boot], dtype=float)
        p_values_model = np.asarray(model["p"][boot], dtype=float)
        flat_indices = freq + np.arange(beta_values.shape[0]) * n_freqs
        h0_beta[channel, flat_indices, :, boot] = beta_values
        h0_r2[channel, flat_indices, 0, boot] = r2_values
        h0_r2[channel, flat_indices, 1, boot] = f_values_model
        h0_r2[channel, flat_indices, 2, boot] = p_values_model
        if h0_conditions is not None and "conditions" in model:
            f_values = np.asarray(model["conditions"]["F"][boot], dtype=float)
            p_values = np.asarray(model["conditions"]["p"][boot], dtype=float)
            if f_values.ndim == 1:
                h0_conditions[channel, flat_indices, 0, 0, boot] = f_values
                h0_conditions[channel, flat_indices, 0, 1, boot] = p_values
            else:
                for effect in range(f_values.shape[0]):
                    h0_conditions[channel, flat_indices, effect, 0, boot] = f_values[effect]
                    h0_conditions[channel, flat_indices, effect, 1, boot] = p_values[effect]
        if h0_interactions is not None and "interactions" in model:
            f_values = np.asarray(model["interactions"]["F"][boot], dtype=float)
            p_values = np.asarray(model["interactions"]["p"][boot], dtype=float)
            if f_values.ndim == 1:
                h0_interactions[channel, flat_indices, 0, 0, boot] = f_values
                h0_interactions[channel, flat_indices, 0, 1, boot] = p_values
            else:
                for effect in range(f_values.shape[0]):
                    h0_interactions[channel, flat_indices, effect, 0, boot] = f_values[effect]
                    h0_interactions[channel, flat_indices, effect, 1, boot] = p_values[effect]
        if h0_covariates is not None and "continuous" in model:
            f_values = np.asarray(model["continuous"]["F"][boot], dtype=float)
            p_values = np.asarray(model["continuous"]["p"][boot], dtype=float)
            if f_values.ndim == 1:
                h0_covariates[channel, flat_indices, 0, 0, boot] = f_values
                h0_covariates[channel, flat_indices, 0, 1, boot] = p_values
            else:
                for effect in range(f_values.shape[-1]):
                    h0_covariates[channel, flat_indices, effect, 0, boot] = f_values[:, effect]
                    h0_covariates[channel, flat_indices, effect, 1, boot] = p_values[:, effect]


def _interaction_factor_map(nb_factors: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    for size in range(2, nb_factors + 1):
        out.extend(itertools.combinations(range(1, nb_factors + 1), size))
    return out


def _first_index_of_factor_size(interaction_map: list[tuple[int, ...]], factor_size: int) -> int:
    for index, interaction in enumerate(interaction_map):
        if len(interaction) == factor_size:
            return index
    return 0


def _make_interactions(x: np.ndarray, nb_conditions: list[int]) -> tuple[np.ndarray, list[int]]:
    if not nb_conditions:
        return x, []
    factors = []
    index = 0
    for levels in nb_conditions:
        factors.append(x[:, index : index + levels])
        index += levels
    tmp_x = x.copy()
    interactions: list[int] = []
    for interaction_size in range(2, len(nb_conditions) + 1):
        for combination in itertools.combinations(range(len(nb_conditions)), interaction_size):
            current = factors[combination[0]]
            for next_index in combination[1:]:
                blocks = []
                for column in range(current.shape[1]):
                    blocks.append(current[:, [column]] * factors[next_index])
                current = np.concatenate(blocks, axis=1)
                current = current[:, np.sum(current, axis=0) != 0]
            interactions.append(current.shape[1])
            tmp_x = np.concatenate((tmp_x, current), axis=1)
    return tmp_x, interactions


def _projection_df(m: np.ndarray) -> float:
    numer = np.trace(m.T @ m) ** 2
    denom = np.trace((m.T @ m) @ (m.T @ m))
    if denom == 0:
        return np.nan
    return float(numer / denom)


def _flatten_fortran(array: np.ndarray) -> np.ndarray:
    return np.reshape(array, (-1,), order="F")


def _flatten_effect_tensor(array: np.ndarray) -> np.ndarray:
    if array.ndim == 2:
        return array
    leading = array.shape[0]
    out = np.zeros((leading, array.shape[1] * array.shape[2]), dtype=float)
    for index in range(leading):
        out[index] = np.reshape(array[index], (-1,), order="F")
    return out


def _stack_or_list(items: list[np.ndarray]) -> Any:
    try:
        return np.stack(items)
    except ValueError:
        return items


def _collect_model_df(store: list[np.ndarray], value: Any) -> None:
    if value is None:
        return
    store.append(np.asarray(value))


def _product_or_zero(values: list[int]) -> int:
    if not values:
        return 0
    out = 1
    for value in values:
        out *= int(value)
    return out


def _ensure_3d(array: np.ndarray) -> np.ndarray:
    array = np.asarray(array, dtype=float)
    if array.ndim == 2:
        return array[:, :, np.newaxis]
    if array.ndim == 3:
        return array
    raise ValueError("Expected a 2D or 3D array.")


def _as_int_list(value: Any) -> list[int]:
    if value is None:
        return []
    array = np.asarray(value)
    if array.size == 0:
        return []
    if array.shape == ():
        scalar = int(array.item())
        return [] if scalar == 0 else [scalar]
    return [int(item) for item in array.tolist()]


def _clone_mapping(mapping: Mapping[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in mapping.items():
        if isinstance(value, Mapping):
            out[key] = _clone_mapping(value)
        elif isinstance(value, list):
            out[key] = [item for item in value]
        elif isinstance(value, np.ndarray):
            out[key] = np.array(value, copy=True)
        else:
            out[key] = value
    return out


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the LIMO GLM workflow on HDF5 outputs")
    parser.add_argument("limo_file", type=Path, help="Path to LIMO.h5")
    parser.add_argument("--results-file", type=Path, default=None, help="Optional explicit path to limo_results.h5")
    parser.add_argument("--h0-file", type=Path, default=None, help="Optional explicit path to limo_H0.h5")
    parser.add_argument(
        "--run-bootstrap",
        choices=["on", "off"],
        default=None,
        help="Override the bootstrap flag stored in LIMO.design.bootstrap",
    )
    parser.add_argument(
        "--variance-estimates",
        choices=["standard", "HC4"],
        default="standard",
        help="Variance estimate mode for the OLS/WLS GLM branches",
    )
    return parser


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()
    limo_path, results_path, h0_path = limo_glm(
        args.limo_file,
        results_file=args.results_file,
        h0_file=args.h0_file,
        run_bootstrap=None if args.run_bootstrap is None else args.run_bootstrap == "on",
        variance_estimates=args.variance_estimates,
    )
    print(json.dumps({"LIMO": str(limo_path), "results": str(results_path), "H0": None if h0_path is None else str(h0_path)}, indent=2))


if __name__ == "__main__":
    main()
