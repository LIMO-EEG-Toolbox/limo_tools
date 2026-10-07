"""Evaluate LIMO contrasts in the HDF5 workflow.

This module is the Python counterpart of MATLAB ``limo_contrast_execute.m``
for the HDF5 workflow introduced by the other Python ports in this
repository. It also embeds the contrast-normalization and validity checks
from ``limo_contrast_checking.m``.

Scope:
    - First-level contrasts and non-repeated second-level GLM contrasts for
      mass-univariate analyses.
    - OLS, WLS, and IRLS branches, including bootstrap contrasts under H0
      when ``limo_H0.h5`` is available.
    - TFCE dispatch through ``limo_tfce.py`` after contrast datasets are
      written.

Not yet ported:
    - Repeated-measures contrast execution from the MATLAB branches that
      depend on ``limo_rep_anova`` and ``limo_robust_rep_anova``.
    - The MATLAB multivariate Hotelling branch.
    - Generalized Welch ANOVA contrast delegation to ``limo_random_robust``.

Primary use:
    Call ``limo_contrast(limo_file, handles)`` where ``handles`` provides the
    contrast vector or matrix ``C`` and the test selector ``F`` (`0` for a
    T contrast, non-zero for an F contrast).

Inputs:
    ``limo_file``:
        Path to ``LIMO.h5``.
    ``handles``:
        Mapping or object with fields ``C`` and ``F``.
    ``results_file``:
        Optional path to ``limo_results.h5``. Defaults to the path recorded in
        ``LIMO.h5`` or to ``LIMO.dir / limo_results.h5``.
    ``h0_file``:
        Optional path to ``limo_H0.h5``.
    ``compute_bootstrap``:
        Override whether contrast bootstrap datasets are computed. If omitted,
        bootstrap contrasts are computed when ``LIMO.design.bootstrap`` is
        non-zero and ``limo_H0.h5`` is available.
    ``contrast``:
        Optional explicit contrast matrix. When omitted, ``handles.C`` is used.
    ``checkfile``:
        Compatibility flag passed to TFCE handling. Defaults to ``"no"`` to
        match the non-interactive Python workflow.

Outputs written to ``limo_results.h5``:
    ``/con_<n>``:
        T-contrast output with last dimension
        ``contrast_value, standard_error, dfe, t_value, p_value``.
    ``/ess_<n>``:
        F-contrast output with last dimension
        ``C*Beta rows, mean_square_error, df, F_value, p_value``.

Outputs written to ``limo_H0.h5``:
    ``/H0_con_<n>`` or ``/H0_ess_<n>`` with bootstrap T/F and p values.

Updates written to ``LIMO.h5``:
    ``LIMO.contrast`` is extended with the evaluated contrast definition and
    test type.

References:
    LIMO Team MATLAB functions ``limo_contrast_execute.m``,
    ``limo_contrast.m``, and ``limo_contrast_checking.m``.
    Christensen, R. (2002). Plane Answers to Complex Questions.
    Friston, K. et al. (2007). Statistical Parametric Mapping.

Command-line usage:
    python limo_contrast.py path/to/LIMO.h5 --contrast "1 0 -1 0"
    python limo_contrast.py path/to/LIMO.h5 --contrast "1 0 -1 0" --test T
    python limo_contrast.py path/to/LIMO.h5 --contrast-file path/to/contrast.txt --test F
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from scipy.stats import f as f_distribution
from scipy.stats import t as t_distribution

from .eeglab_import import write_hdf5_structure
from .limo_design import read_hdf5_structure


def limo_contrast(
    limo_file: str | Path,
    handles: Mapping[str, Any] | Any,
    *,
    results_file: str | Path | None = None,
    h0_file: str | Path | None = None,
    compute_bootstrap: bool | None = None,
    contrast: Any = None,
    checkfile: str = "no",
) -> tuple[Path, Path, Path | None, str]:
    """Execute a contrast analysis and write the outputs to the HDF5 files."""

    limo_path = Path(limo_file).expanduser().resolve()
    payload = read_hdf5_structure(limo_path)
    if "LIMO" not in payload or not isinstance(payload["LIMO"], Mapping):
        raise KeyError("The HDF5 file does not contain a top-level 'LIMO' group.")

    limo = _clone_mapping(payload["LIMO"])
    results_path = _resolve_results_path(limo, results_file)
    results_payload = _clone_mapping(read_hdf5_structure(results_path))

    design = limo.get("design", {}) if isinstance(limo.get("design", {}), Mapping) else {}
    design_name = str(design.get("name", ""))
    if "Repeated" in design_name:
        raise NotImplementedError(
            "Repeated-measures contrast execution is not implemented yet in limo_contrast.py. "
            "It still depends on the unported repeated-measures ANOVA functions."
        )
    if str(design.get("type_of_analysis", "Mass-univariate")) != "Mass-univariate":
        raise NotImplementedError(
            "The multivariate Hotelling contrast branch is not implemented yet in limo_contrast.py."
        )
    if "Generalized Welch's method" in str(design.get("method", "")):
        raise NotImplementedError(
            "Generalized Welch ANOVA contrast delegation is not implemented in limo_contrast.py."
        )

    requested_contrast = np.asarray(_extract_handle_field(handles, "C") if contrast is None else contrast, dtype=float)
    requested_contrast = np.atleast_2d(requested_contrast)
    test_flag = int(_extract_handle_field(handles, "F", default=0))
    x = np.asarray(design.get("X"), dtype=float)
    corrected = np.asarray(limo_contrast_checking(limo_path.parent, x, requested_contrast), dtype=float)
    if not limo_contrast_checking(corrected, x):
        raise ValueError("Invalid contrast as input.")

    if corrected.shape[0] > 1 and test_flag == 0:
        warnings.warn("The requested contrast has multiple rows; switching to an F contrast.", RuntimeWarning, stacklevel=2)
        test_flag = 1

    contrast_entries = list(limo.get("contrast", []))
    contrast_entries.append({"C": corrected, "V": "T" if test_flag == 0 else "F"})
    contrast_index = len(contrast_entries)
    limo["contrast"] = contrast_entries

    yr_native = np.asarray(results_payload["Yr"], dtype=float)
    beta_native = np.asarray(results_payload["Beta"], dtype=float)
    res_native = np.asarray(results_payload["Res"], dtype=float)
    analysis = str(limo.get("Analysis", ""))
    method = str(design.get("method", "OLS"))
    weights = np.asarray(design.get("weights"), dtype=float) if "weights" in design else None
    dfe = _extract_dfe(limo, yr_native, analysis)

    if analysis.lower() == "time-frequency" and method.upper() == "WLS":
        output = _compute_tf_wls_contrast(
            yr=yr_native,
            beta=beta_native,
            res=res_native,
            x=x,
            contrast=corrected,
            test_flag=test_flag,
            weights=weights,
            dfe=dfe,
        )
    else:
        yr_working, tf_shape = _flatten_native_data(yr_native, analysis)
        beta_working, _ = _flatten_native_data(beta_native, analysis)
        res_working, _ = _flatten_native_data(res_native, analysis)
        dfe_working = _flatten_dfe(dfe, yr_working.shape[:2], tf_shape)
        output_flat = _compute_standard_contrast(
            beta=beta_working,
            res=res_working,
            x=x,
            contrast=corrected,
            test_flag=test_flag,
            method=method,
            weights=weights,
            dfe=dfe_working,
        )
        output = _restore_native_data(output_flat, tf_shape)

    result_name = f"con_{contrast_index}" if test_flag == 0 else f"ess_{contrast_index}"
    results_payload[result_name] = output
    write_hdf5_structure(results_path, results_payload)

    h0_path = _resolve_h0_path(limo, h0_file)
    do_bootstrap = _resolve_compute_bootstrap(limo, h0_path, compute_bootstrap)
    if do_bootstrap:
        h0_payload = _clone_mapping(read_hdf5_structure(h0_path)) if h0_path.exists() else {}
        boot_output = _compute_bootstrap_contrast(
            limo=limo,
            results_payload=results_payload,
            h0_payload=h0_payload,
            contrast=corrected,
            test_flag=test_flag,
        )
        h0_payload[f"H0_{result_name}"] = boot_output
        write_hdf5_structure(h0_path, h0_payload)
    else:
        h0_path = h0_path if h0_path.exists() else None

    write_hdf5_structure(limo_path, {"LIMO": limo})

    if int(np.asarray(design.get("tfce", 0)).item()) == 1:
        try:
            from .limo_tfce import limo_tfce as limo_tfce_handling
        except ImportError:
            warnings.warn("TFCE is enabled, but limo_tfce.py could not be imported.", RuntimeWarning, stacklevel=2)
        else:
            limo_tfce_handling(limo_path, results_path, h0_path, limo, stat_name=result_name, checkfile=checkfile)

    return limo_path, results_path, h0_path, result_name


def limo_contrast_execute(
    limo_file: str | Path,
    handles: Mapping[str, Any] | Any,
    **kwargs: Any,
) -> tuple[Path, Path, Path | None, str]:
    """Compatibility alias matching the MATLAB entry-point name."""

    return limo_contrast(limo_file, handles, **kwargs)


def limo_contrast_checking(*args: Any) -> Any:
    """Port of the main normalization and validity checks from MATLAB."""

    if len(args) == 0:
        raise ValueError("limo_contrast_checking expects 1, 2, or 3 arguments")

    diagonalize = False
    if len(args) in {1, 3}:
        limo_path = Path(args[0]).expanduser().resolve()
        if not limo_path.exists():
            warnings.warn(f"{limo_path} does not exist, updating to local dir", RuntimeWarning, stacklevel=2)
            limo_path = Path.cwd()

        if len(args) == 1:
            payload = read_hdf5_structure(limo_path / "LIMO.h5")
            limo = payload["LIMO"]
            x = np.asarray(limo["design"]["X"], dtype=float)
            contrast_entries = list(limo.get("contrast", []))
            if not contrast_entries:
                raise ValueError("No contrast found in LIMO.h5")
            c = np.asarray(contrast_entries[-1]["C"], dtype=float)
        else:
            x = np.asarray(args[1], dtype=float)
            c = np.asarray(args[2]["C"], dtype=float) if isinstance(args[2], Mapping) and "C" in args[2] else np.asarray(args[2], dtype=float)

        c = np.atleast_2d(c)
        rows, cols = c.shape
        if cols != x.shape[1]:
            if rows == x.shape[1]:
                c = c.T
                rows, cols = c.shape
            if rows == cols:
                diagonalize = True
                c = np.diag(c)[np.newaxis, :]
                rows, cols = c.shape
            if cols < x.shape[1]:
                tmp = np.zeros((rows, x.shape[1]), dtype=float)
                tmp[:, :cols] = c
                c = tmp
                cols = c.shape[1]
            if cols != x.shape[1]:
                raise ValueError("c must have the same number of columns as X")
            if diagonalize:
                c = np.diag(c.squeeze())
        if len(args) == 1:
            payload = read_hdf5_structure(limo_path / "LIMO.h5")
            limo = _clone_mapping(payload["LIMO"])
            contrast_entries = list(limo.get("contrast", []))
            contrast_entries[-1]["C"] = c
            limo["contrast"] = contrast_entries
            write_hdf5_structure(limo_path / "LIMO.h5", {"LIMO": limo})
        return c

    if len(args) == 2:
        contrast = np.atleast_2d(np.asarray(args[0], dtype=float))
        x = np.asarray(args[1], dtype=float)
        if contrast.shape[1] != x.shape[1]:
            raise ValueError("the length of the contrast must be equal to the number of regressors in the design")

        for row in contrast:
            n = np.sum(x[:, row != 0], axis=0)
            if np.sum(row) == np.count_nonzero(row):
                valid = True
            elif row[row != 0].size and np.all(_matlab_int16(row[row != 0] - n / np.sum(n)) == 0):
                valid = True
            else:
                projection = x @ np.linalg.pinv(x)
                lam = x @ row[:, np.newaxis]
                check = _matlab_int16(projection @ lam) == _matlab_int16(lam)
                valid = bool(np.sum(check) == x.shape[0])
            if not valid:
                return 0
            if row[-1] != 0:
                raise ValueError("the contrast requested include the constant term, which is not possible")
        return 1

    raise ValueError("the number of arguments must be 2 or 3")


def _matlab_int16(values: Any) -> np.ndarray:
    """Match MATLAB integer conversion: round ties away from zero, saturate."""
    values = np.asarray(values, dtype=float)
    rounded = np.copysign(np.floor(np.abs(values) + 0.5), values)
    return np.clip(rounded, -32768, 32767).astype(np.int16)


def _compute_standard_contrast(
    *,
    beta: np.ndarray,
    res: np.ndarray,
    x: np.ndarray,
    contrast: np.ndarray,
    test_flag: int,
    method: str,
    weights: np.ndarray | None,
    dfe: np.ndarray,
) -> np.ndarray:
    n_channels, n_frames, _ = beta.shape
    if test_flag == 0:
        out = np.full((n_channels, n_frames, 5), np.nan, dtype=float)
    else:
        out = np.full((n_channels, n_frames, contrast.shape[0] + 4), np.nan, dtype=float)

    array = np.flatnonzero(~np.isnan(res[:, 0, 0])) if res.ndim == 3 else np.arange(n_channels)
    contrast_scale_unweighted = float((contrast @ np.linalg.pinv(x.T @ x) @ contrast.T).squeeze()) if test_flag == 0 else None

    for channel in array:
        residual = np.asarray(res[channel], dtype=float)
        beta_channel = np.asarray(beta[channel], dtype=float)
        valid_obs = ~np.isnan(residual[0])
        if not np.any(valid_obs):
            continue

        if method.upper() in {"OLS", "WLS"}:
            wx = x if method.upper() == "OLS" else x * np.asarray(weights[channel], dtype=float)[:, np.newaxis]
            pinv_term = np.linalg.pinv(wx.T @ wx)
            contrast_scale = float((contrast @ pinv_term @ contrast.T).squeeze()) if test_flag == 0 else None

            sse = np.sum(residual * residual, axis=1) / np.asarray(dfe[channel], dtype=float)
            if test_flag == 0:
                estimate = (contrast @ beta_channel.T).squeeze()
                se = np.sqrt(sse * contrast_scale)
                t_values = estimate / se
                p_values = (t_distribution.sf(np.abs(t_values), np.asarray(dfe[channel], dtype=float))) * 2
                out[channel, :, 0] = estimate
                out[channel, :, 1] = se
                out[channel, :, 2] = np.asarray(dfe[channel], dtype=float)
                out[channel, :, 3] = t_values
                out[channel, :, 4] = p_values
            else:
                df = _f_contrast_df(contrast)
                c_matrix = _as_f_contrast_matrix(contrast)
                c0 = np.eye(c_matrix.shape[1]) - c_matrix @ np.linalg.pinv(c_matrix)
                x0 = x @ c0
                r = np.eye(x.shape[0]) - wx @ np.linalg.pinv(wx)
                r0 = np.eye(x.shape[0]) - x0 @ np.linalg.pinv(x0)
                m = r0 - r
                contrast_values = contrast @ beta_channel.T
                h = np.einsum("fp,pq,fq->f", beta_channel, x.T @ m @ x, beta_channel)
                f_values = (h / df) / sse
                p_values = f_distribution.sf(f_values, df, np.asarray(dfe[channel], dtype=float))
                out[channel, :, : contrast.shape[0]] = contrast_values.T
                out[channel, :, contrast.shape[0]] = sse
                out[channel, :, contrast.shape[0] + 1] = df
                out[channel, :, contrast.shape[0] + 2] = f_values
                out[channel, :, contrast.shape[0] + 3] = p_values
        else:
            for frame in range(n_frames):
                residual_frame = residual[frame]
                valid = ~np.isnan(residual_frame)
                if not np.any(valid):
                    continue
                wx = x[valid] * np.asarray(weights[channel, frame, valid], dtype=float)[:, np.newaxis]
                beta_frame = beta_channel[frame]
                dfe_frame = float(np.asarray(dfe[channel, frame]).item())
                sse = float(np.sum(residual_frame[valid] ** 2) / dfe_frame)
                if test_flag == 0:
                    estimate = float((contrast @ beta_frame[:, np.newaxis]).squeeze())
                    se = float(np.sqrt(sse * (contrast @ np.linalg.pinv(wx.T @ wx) @ contrast.T).squeeze()))
                    t_value = estimate / se
                    p_value = float((t_distribution.sf(abs(t_value), dfe_frame)) * 2)
                    out[channel, frame] = np.asarray([estimate, se, dfe_frame, t_value, p_value], dtype=float)
                else:
                    df = _f_contrast_df(contrast)
                    c_matrix = _as_f_contrast_matrix(contrast)
                    c0 = np.eye(c_matrix.shape[1]) - c_matrix @ np.linalg.pinv(c_matrix)
                    x0 = x[valid] @ c0
                    r = np.eye(np.sum(valid)) - wx @ np.linalg.pinv(wx)
                    r0 = np.eye(np.sum(valid)) - x0 @ np.linalg.pinv(x0)
                    m = r0 - r
                    h = float(beta_frame @ x[valid].T @ m @ x[valid] @ beta_frame)
                    f_value = (h / df) / sse
                    p_value = float(f_distribution.sf(f_value, df, dfe_frame))
                    out[channel, frame, : contrast.shape[0]] = (contrast @ beta_frame[:, np.newaxis]).T
                    out[channel, frame, contrast.shape[0] : contrast.shape[0] + 4] = np.asarray([sse, df, f_value, p_value])
    return out


def _compute_tf_wls_contrast(
    *,
    yr: np.ndarray,
    beta: np.ndarray,
    res: np.ndarray,
    x: np.ndarray,
    contrast: np.ndarray,
    test_flag: int,
    weights: np.ndarray | None,
    dfe: np.ndarray,
) -> np.ndarray:
    n_channels, n_freqs, n_times, _ = yr.shape
    if test_flag == 0:
        out = np.full((n_channels, n_freqs, n_times, 5), np.nan, dtype=float)
    else:
        out = np.full((n_channels, n_freqs, n_times, contrast.shape[0] + 4), np.nan, dtype=float)
    for channel in np.flatnonzero(~np.isnan(yr[:, 0, 0, 0])):
        for freq in range(n_freqs):
            wx = x * np.asarray(weights[channel, freq], dtype=float)[:, np.newaxis]
            pinv_term = np.linalg.pinv(wx.T @ wx)
            df_freq = np.asarray(dfe[channel, freq], dtype=float)
            for time in range(n_times):
                sse = float(np.sum(res[channel, freq, time] ** 2) / df_freq)
                beta_frame = np.asarray(beta[channel, freq, time], dtype=float)
                if test_flag == 0:
                    estimate = float((contrast @ beta_frame[:, np.newaxis]).squeeze())
                    se = float(np.sqrt(sse * (contrast @ pinv_term @ contrast.T).squeeze()))
                    t_value = estimate / se
                    p_value = float((t_distribution.sf(abs(t_value), df_freq)) * 2)
                    out[channel, freq, time] = np.asarray([estimate, se, df_freq, t_value, p_value])
                else:
                    df = _f_contrast_df(contrast)
                    c_matrix = _as_f_contrast_matrix(contrast)
                    c0 = np.eye(c_matrix.shape[1]) - c_matrix @ np.linalg.pinv(c_matrix)
                    x0 = x @ c0
                    r = np.eye(x.shape[0]) - wx @ np.linalg.pinv(wx)
                    r0 = np.eye(x.shape[0]) - x0 @ np.linalg.pinv(x0)
                    m = r0 - r
                    h = float(beta_frame @ x.T @ m @ x @ beta_frame)
                    f_value = (h / df) / sse
                    p_value = float(f_distribution.sf(f_value, df, df_freq))
                    out[channel, freq, time, : contrast.shape[0]] = (contrast @ beta_frame[:, np.newaxis]).T
                    out[channel, freq, time, contrast.shape[0] : contrast.shape[0] + 4] = np.asarray([sse, df, f_value, p_value])
    return out


def _compute_bootstrap_contrast(
    *,
    limo: Mapping[str, Any],
    results_payload: Mapping[str, Any],
    h0_payload: Mapping[str, Any],
    contrast: np.ndarray,
    test_flag: int,
) -> np.ndarray:
    design = limo.get("design", {}) if isinstance(limo.get("design", {}), Mapping) else {}
    method = str(design.get("method", "OLS"))
    analysis = str(limo.get("Analysis", ""))
    x = np.asarray(design.get("X"), dtype=float)
    yr_native = np.asarray(results_payload["Yr"], dtype=float)
    h0_beta_native = np.asarray(h0_payload["H0_Beta"], dtype=float)
    boot_table = h0_payload.get("boot_table")
    if boot_table is None:
        raise KeyError("limo_H0.h5 does not contain /boot_table, required for contrast bootstraps.")

    if analysis.lower() == "time-frequency" and method.upper() == "WLS":
        return _compute_bootstrap_tf_wls_contrast(limo, yr_native, h0_beta_native, contrast, test_flag, x, boot_table)

    yr_working, tf_shape = _flatten_native_data(yr_native, analysis)
    h0_beta_working, _ = _flatten_native_data(h0_beta_native, analysis)
    centered_data = np.full_like(yr_working, np.nan, dtype=float)
    for channel in range(yr_working.shape[0]):
        centered_data[channel] = _limo_glm_null(
            np.asarray(yr_working[channel], dtype=float).T,
            x,
            _as_int_list(design.get("nb_conditions")),
            _as_int_list(design.get("nb_interactions")),
        ).T

    weights = np.asarray(design.get("weights"), dtype=float) if "weights" in design else None
    dfe = _flatten_dfe(_extract_dfe(limo, yr_native, analysis), yr_working.shape[:2], tf_shape)
    nboot = h0_beta_working.shape[-1]
    out = np.full((yr_working.shape[0], yr_working.shape[1], 2, nboot), np.nan, dtype=float)

    for channel in np.flatnonzero(~np.isnan(yr_working[:, 0, 0])):
        channel_boot_table = _select_channel_boot_table(boot_table, channel)
        for boot in range(nboot):
            resampling_index = np.asarray(channel_boot_table[:, boot], dtype=int)
            y_boot = centered_data[channel, :, resampling_index].T
            if method.upper() in {"OLS", "WLS"}:
                trials_to_keep = ~np.isnan(y_boot[:, 0])
                y_boot = y_boot[trials_to_keep]
                x_boot = x[trials_to_keep]
                if method.upper() == "OLS":
                    weight_vector = np.ones(x_boot.shape[0], dtype=float)
                else:
                    weight_vector = np.asarray(weights[channel], dtype=float)[trials_to_keep]
                if np.any(~trials_to_keep) and int(np.asarray(design.get("nb_continuous", 0)).item()) != 0 and int(np.asarray(design.get("zscore", 0)).item()) == 1:
                    x_boot = _rezscore_covariates(x_boot, design)
                wx = x_boot * weight_vector[:, np.newaxis]
                r = np.eye(y_boot.shape[0]) - wx @ np.linalg.pinv(wx)
                residual = r @ y_boot
                beta_boot = np.asarray(h0_beta_working[channel, :, :, boot], dtype=float)
                if test_flag == 0:
                    sse = np.sum(residual * residual, axis=0) / np.asarray(dfe[channel], dtype=float)
                    t_values = (contrast @ beta_boot.T).squeeze() / np.sqrt(sse * (contrast @ np.linalg.pinv(wx.T @ wx) @ contrast.T).squeeze())
                    p_values = t_distribution.sf(t_values, np.asarray(dfe[channel], dtype=float))
                    out[channel, :, 0, boot] = t_values
                    out[channel, :, 1, boot] = p_values
                else:
                    e = np.sum(residual * residual, axis=0)
                    c_matrix = _as_f_contrast_matrix(contrast)
                    c0 = np.eye(c_matrix.shape[1]) - c_matrix @ np.linalg.pinv(c_matrix)
                    x0 = wx @ c0
                    r0 = np.eye(y_boot.shape[0]) - x0 @ np.linalg.pinv(x0)
                    m = r0 - r
                    h = np.einsum("fp,pq,fq->f", beta_boot, x_boot.T @ m @ x_boot, beta_boot)
                    df = _f_contrast_df(contrast)
                    f_values = (h / df) / (e / np.asarray(dfe[channel], dtype=float))
                    p_values = f_distribution.sf(f_values, df, np.asarray(dfe[channel], dtype=float))
                    out[channel, :, 0, boot] = f_values
                    out[channel, :, 1, boot] = p_values
            else:
                beta_boot = np.asarray(h0_beta_working[channel, :, :, boot], dtype=float)
                for frame in range(yr_working.shape[1]):
                    valid = ~np.isnan(y_boot[:, frame])
                    if not np.any(valid):
                        continue
                    x_boot = x[valid]
                    weights_frame = np.asarray(weights[channel, frame], dtype=float)[valid]
                    y_frame = y_boot[valid, frame]
                    if np.any(~valid) and int(np.asarray(design.get("nb_continuous", 0)).item()) != 0 and int(np.asarray(design.get("zscore", 0)).item()) == 1:
                        x_boot = _rezscore_covariates(x_boot, design)
                    wx = x_boot * weights_frame[:, np.newaxis]
                    hm = wx @ np.linalg.pinv(wx)
                    r = np.eye(y_frame.shape[0]) - hm
                    dfe_frame = float(np.asarray(dfe[channel, frame]).item())
                    if test_flag == 0:
                        sse = float((r @ y_frame).T @ (r @ y_frame) / dfe_frame)
                        t_value = float((contrast @ beta_boot[frame][:, np.newaxis]).squeeze() / np.sqrt(sse * (contrast @ np.linalg.pinv(x_boot.T @ x_boot) @ contrast.T).squeeze()))
                        p_value = float(t_distribution.sf(t_value, dfe_frame))
                        out[channel, frame, 0, boot] = t_value
                        out[channel, frame, 1, boot] = p_value
                    else:
                        e = float(y_frame.T @ r @ y_frame)
                        c_matrix = _as_f_contrast_matrix(contrast)
                        c0 = np.eye(c_matrix.shape[1]) - c_matrix @ np.linalg.pinv(c_matrix)
                        x0 = wx @ c0
                        r0 = np.eye(y_frame.shape[0]) - x0 @ np.linalg.pinv(x0)
                        m = r0 - r
                        h = float(beta_boot[frame] @ x_boot.T @ m @ x_boot @ beta_boot[frame])
                        df = _f_contrast_df(contrast)
                        f_value = (h / df) / (e / dfe_frame)
                        p_value = float(f_distribution.sf(f_value, df, dfe_frame))
                        out[channel, frame, 0, boot] = f_value
                        out[channel, frame, 1, boot] = p_value

    return _restore_native_data(out, tf_shape)


def _compute_bootstrap_tf_wls_contrast(
    limo: Mapping[str, Any],
    yr_native: np.ndarray,
    h0_beta_native: np.ndarray,
    contrast: np.ndarray,
    test_flag: int,
    x: np.ndarray,
    boot_table: Any,
) -> np.ndarray:
    design = limo.get("design", {}) if isinstance(limo.get("design", {}), Mapping) else {}
    weights = np.asarray(design.get("weights"), dtype=float)
    dfe = _extract_dfe(limo, yr_native, "Time-Frequency")
    centered_data = np.full_like(yr_native, np.nan, dtype=float)
    for channel in range(yr_native.shape[0]):
        flattened = _flatten_tf_native(yr_native[channel : channel + 1])[0]
        centered = _limo_glm_null(flattened.T, x, _as_int_list(design.get("nb_conditions")), _as_int_list(design.get("nb_interactions"))).T
        centered_data[channel] = _restore_tf_native(centered[np.newaxis, ...], (1, yr_native.shape[1], yr_native.shape[2], yr_native.shape[3]))[0]

    nboot = h0_beta_native.shape[-1]
    out = np.full((yr_native.shape[0], yr_native.shape[1], yr_native.shape[2], 2, nboot), np.nan, dtype=float)
    for channel in np.flatnonzero(~np.isnan(yr_native[:, 0, 0, 0])):
        channel_boot_table = _select_channel_boot_table(boot_table, channel)
        for boot in range(nboot):
            resampling_index = np.asarray(channel_boot_table[:, boot], dtype=int)
            for freq in range(yr_native.shape[1]):
                y_boot = centered_data[channel, freq, :, resampling_index].T
                valid = ~np.isnan(y_boot[:, 0])
                if not np.any(valid):
                    continue
                x_boot = x[valid]
                wx = x_boot * np.asarray(weights[channel, freq], dtype=float)[valid][:, np.newaxis]
                residual = (np.eye(y_boot[valid].shape[0]) - wx @ np.linalg.pinv(wx)) @ y_boot[valid]
                beta_boot = np.asarray(h0_beta_native[channel, freq, :, :, boot], dtype=float)
                sse = np.sum(residual * residual, axis=0) / np.asarray(dfe[channel, freq], dtype=float)
                if test_flag == 0:
                    t_values = (contrast @ beta_boot.transpose(1, 0)).squeeze() / np.sqrt(sse * (contrast @ np.linalg.pinv(wx.T @ wx) @ contrast.T).squeeze())
                    p_values = t_distribution.sf(t_values, np.asarray(dfe[channel, freq], dtype=float))
                    out[channel, freq, :, 0, boot] = t_values
                    out[channel, freq, :, 1, boot] = p_values
                else:
                    c_matrix = _as_f_contrast_matrix(contrast)
                    c0 = np.eye(c_matrix.shape[1]) - c_matrix @ np.linalg.pinv(c_matrix)
                    x0 = x_boot @ c0
                    r = np.eye(y_boot[valid].shape[0]) - wx @ np.linalg.pinv(wx)
                    r0 = np.eye(y_boot[valid].shape[0]) - x0 @ np.linalg.pinv(x0)
                    m = r0 - r
                    h = np.einsum("tp,pq,tq->t", beta_boot, x_boot.T @ m @ x_boot, beta_boot)
                    df = _f_contrast_df(contrast)
                    f_values = (h / df) / sse
                    p_values = f_distribution.sf(f_values, df, np.asarray(dfe[channel, freq], dtype=float))
                    out[channel, freq, :, 0, boot] = f_values
                    out[channel, freq, :, 1, boot] = p_values
    return out


def _extract_dfe(limo: Mapping[str, Any], yr_native: np.ndarray, analysis: str) -> np.ndarray:
    model = limo.get("model", {}) if isinstance(limo.get("model", {}), Mapping) else {}
    if "model_df" in model:
        df = np.asarray(model["model_df"], dtype=float)
        if df.shape and df.shape[-1] >= 2:
            dfe = df[..., 1]
        else:
            dfe = df
    else:
        x = np.asarray(limo.get("design", {}).get("X"), dtype=float)
        dfe = np.asarray(yr_native.shape[-1] - np.linalg.matrix_rank(x), dtype=float)

    if analysis.lower() == "time-frequency":
        n_channels, n_freqs, n_times, _ = yr_native.shape
        if np.ndim(dfe) == 0:
            return np.full((n_channels, n_freqs), float(dfe), dtype=float)
        dfe = np.asarray(dfe, dtype=float)
        if dfe.shape == (n_channels,):
            return np.repeat(dfe[:, np.newaxis], n_freqs, axis=1)
        if dfe.shape == (n_channels, n_freqs):
            return dfe
        if dfe.shape == (n_channels, n_freqs, n_times):
            return dfe
    else:
        n_channels = yr_native.shape[0]
        if np.ndim(dfe) == 0:
            return np.full(n_channels, float(dfe), dtype=float)
        dfe = np.asarray(dfe, dtype=float)
        if dfe.shape == (n_channels,):
            return dfe
        if dfe.shape == (n_channels, yr_native.shape[1]):
            return dfe
    return np.asarray(dfe, dtype=float)


def _flatten_dfe(dfe: np.ndarray, working_shape: tuple[int, int], tf_shape: tuple[int, int, int, int] | None) -> np.ndarray:
    dfe = np.asarray(dfe, dtype=float)
    if tf_shape is None:
        if dfe.ndim == 0:
            return np.full(working_shape[0], float(dfe), dtype=float)
        return dfe
    n_channels, n_freqs, n_times, _ = tf_shape
    if dfe.ndim == 0:
        return np.full((n_channels, n_freqs * n_times), float(dfe), dtype=float)
    if dfe.shape == (n_channels,):
        return np.repeat(dfe[:, np.newaxis], n_freqs * n_times, axis=1)
    if dfe.shape == (n_channels, n_freqs):
        return np.repeat(dfe, n_times, axis=1)
    if dfe.shape == (n_channels, n_freqs, n_times):
        return _flatten_tf_native(dfe[..., np.newaxis])[..., 0]
    return dfe


def _resolve_results_path(limo: Mapping[str, Any], results_file: str | Path | None) -> Path:
    if results_file is not None:
        return Path(results_file).expanduser().resolve()
    data = limo.get("data", {}) if isinstance(limo.get("data", {}), Mapping) else {}
    if "results_file" in data:
        return Path(str(data["results_file"])).expanduser().resolve()
    return Path(str(limo["dir"])).expanduser().resolve() / "limo_results.h5"


def _resolve_h0_path(limo: Mapping[str, Any], h0_file: str | Path | None) -> Path:
    if h0_file is not None:
        return Path(h0_file).expanduser().resolve()
    return Path(str(limo["dir"])).expanduser().resolve() / "limo_H0.h5"


def _resolve_compute_bootstrap(limo: Mapping[str, Any], h0_path: Path, compute_bootstrap: bool | None) -> bool:
    if compute_bootstrap is not None:
        return bool(compute_bootstrap)
    design = limo.get("design", {}) if isinstance(limo.get("design", {}), Mapping) else {}
    return int(np.asarray(design.get("bootstrap", 0)).item()) != 0 and h0_path.exists()


def _extract_handle_field(handles: Mapping[str, Any] | Any, name: str, default: Any = None) -> Any:
    if isinstance(handles, Mapping):
        return handles.get(name, default)
    return getattr(handles, name, default)


def _as_f_contrast_matrix(contrast: np.ndarray) -> np.ndarray:
    if contrast.shape[0] == 1:
        return np.diag(contrast.squeeze())
    return np.asarray(contrast, dtype=float)


def _f_contrast_df(contrast: np.ndarray) -> int:
    c_matrix = _as_f_contrast_matrix(contrast)
    df = int(np.linalg.matrix_rank(c_matrix) - 1)
    return 1 if df <= 0 else df


def _rezscore_covariates(x: np.ndarray, design: Mapping[str, Any]) -> np.ndarray:
    x = np.array(x, copy=True, dtype=float)
    nb_conditions = int(np.sum(_as_int_list(design.get("nb_conditions"))))
    nb_interactions = int(np.sum(_as_int_list(design.get("nb_interactions"))))
    start = nb_conditions + nb_interactions
    stop = x.shape[1] - 1
    if start < stop:
        subset = x[:, start:stop]
        means = np.mean(subset, axis=0, keepdims=True)
        stds = np.std(subset, axis=0, ddof=1, keepdims=True)
        stds[stds == 0] = 1
        x[:, start:stop] = (subset - means) / stds
    return x


def _select_channel_boot_table(boot_table: Any, channel: int) -> np.ndarray:
    if isinstance(boot_table, list):
        return np.asarray(boot_table[channel], dtype=int)
    if isinstance(boot_table, np.ndarray) and boot_table.dtype == object:
        return np.asarray(boot_table[channel], dtype=int)
    return np.asarray(boot_table, dtype=int)


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


def _make_interactions(x: np.ndarray, nb_conditions: list[int]) -> tuple[np.ndarray, list[int]]:
    factors = []
    start = 0
    for levels in nb_conditions:
        factors.append(x[:, start : start + levels])
        start += levels
    tmp_x = np.array(x, copy=True)
    interactions: list[int] = []
    if not factors:
        return tmp_x, interactions
    import itertools

    for size in range(2, len(factors) + 1):
        for combo in itertools.combinations(range(len(factors)), size):
            current = factors[combo[0]]
            for next_index in combo[1:]:
                blocks = [current[:, [column]] * factors[next_index] for column in range(current.shape[1])]
                current = np.concatenate(blocks, axis=1)
                current = current[:, np.sum(current, axis=0) != 0]
            interactions.append(current.shape[1])
            tmp_x = np.concatenate((tmp_x, current), axis=1)
    return tmp_x, interactions


def _flatten_native_data(data: np.ndarray, analysis: str) -> tuple[np.ndarray, tuple[int, int, int, int] | None]:
    array = np.asarray(data, dtype=float)
    if analysis.lower() != "time-frequency":
        return array, None
    return _flatten_tf_native(array), array.shape[:4]


def _restore_native_data(data: np.ndarray, tf_shape: tuple[int, int, int, int] | None) -> np.ndarray:
    if tf_shape is None:
        return data
    return _restore_tf_native(data, tf_shape)


def _flatten_tf_native(data: np.ndarray) -> np.ndarray:
    array = np.asarray(data, dtype=float)
    if array.ndim < 4:
        raise ValueError("Time-Frequency arrays must have at least 4 dimensions.")
    rest = array.shape[3:]
    transposed = np.transpose(array, (0, 2, 1) + tuple(range(3, array.ndim)))
    return transposed.reshape((array.shape[0], array.shape[1] * array.shape[2]) + rest)


def _restore_tf_native(data: np.ndarray, tf_shape: tuple[int, int, int, int]) -> np.ndarray:
    array = np.asarray(data, dtype=float)
    n_channels, n_freqs, n_times, _ = tf_shape
    rest = array.shape[2:]
    reshaped = array.reshape((n_channels, n_times, n_freqs) + rest)
    return np.transpose(reshaped, (0, 2, 1) + tuple(range(3, reshaped.ndim)))


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


def _product_or_zero(values: list[int]) -> int:
    if not values:
        return 0
    out = 1
    for value in values:
        out *= int(value)
    return out


def _clone_mapping(mapping: Mapping[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in mapping.items():
        if isinstance(value, Mapping):
            out[key] = _clone_mapping(value)
        elif isinstance(value, list):
            out[key] = [(_clone_mapping(item) if isinstance(item, Mapping) else np.array(item, copy=True) if isinstance(item, np.ndarray) else item) for item in value]
        elif isinstance(value, np.ndarray):
            out[key] = np.array(value, copy=True)
        else:
            out[key] = value
    return out


def _load_contrast_text(path: Path) -> np.ndarray:
    return np.atleast_2d(np.loadtxt(path, dtype=float))


def _parse_inline_contrast(raw: str) -> np.ndarray:
    rows = []
    for row in raw.split(";"):
        row = row.strip()
        if not row:
            continue
        rows.append([float(value) for value in row.replace(",", " ").split()])
    if not rows:
        raise ValueError("No contrast values were provided.")
    return np.asarray(rows, dtype=float)


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Execute LIMO contrasts on HDF5 GLM outputs")
    parser.add_argument("limo_file", type=Path, help="Path to LIMO.h5")
    parser.add_argument("--results-file", type=Path, default=None, help="Optional explicit path to limo_results.h5")
    parser.add_argument("--h0-file", type=Path, default=None, help="Optional explicit path to limo_H0.h5")
    parser.add_argument("--contrast", type=str, default=None, help="Inline contrast values, rows separated by ';'")
    parser.add_argument("--contrast-file", type=Path, default=None, help="Optional text file with contrast values")
    parser.add_argument("--test", choices=["T", "F"], default="T", help="Contrast test family")
    parser.add_argument("--compute-bootstrap", choices=["on", "off"], default=None, help="Override contrast bootstrap computation")
    return parser


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()
    if args.contrast is None and args.contrast_file is None:
        raise ValueError("Provide either --contrast or --contrast-file")
    contrast = _parse_inline_contrast(args.contrast) if args.contrast is not None else _load_contrast_text(args.contrast_file)
    handles = {"C": contrast, "F": 0 if args.test == "T" else 1}
    limo_path, results_path, h0_path, result_name = limo_contrast(
        args.limo_file,
        handles,
        results_file=args.results_file,
        h0_file=args.h0_file,
        compute_bootstrap=None if args.compute_bootstrap is None else args.compute_bootstrap == "on",
    )
    print(
        json.dumps(
            {
                "LIMO": str(limo_path),
                "results": str(results_path),
                "H0": None if h0_path is None else str(h0_path),
                "dataset": result_name,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
