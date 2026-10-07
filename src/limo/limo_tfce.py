"""Apply Threshold-Free Cluster Enhancement to LIMO HDF5 results.

This module combines the roles of MATLAB ``limo_tfce_handling.m`` and
``limo_tfce.m`` for the HDF5-based Python workflow introduced in this
repository.

Primary use:
    Run ``limo_tfce(limo_file, results_file, h0_file)`` after
    ``limo_glm.py`` has written ``LIMO.h5``, ``limo_results.h5``, and
    optionally ``limo_H0.h5``.

Public entry points:
    ``limo_tfce(...)``
        File-oriented wrapper modeled after ``limo_tfce_handling.m``. It
        selects the correct statistic from the results files, computes TFCE
        scores for the observed data, optionally computes TFCE scores for the
        null bootstrap distribution, and saves the outputs under ``tfce`` and
        ``H0`` folders.
    ``limo_tfce_transform(...)``
        Direct array-oriented port of ``limo_tfce.m``. It accepts the MATLAB
        ``type`` convention, a statistic map or stack of maps, an optional
        neighbouring-channel matrix, and the TFCE parameters.

Inputs for ``limo_tfce``:
    ``limo_file``:
        Path to ``LIMO.h5``.
    ``results_file``:
        Optional path to ``limo_results.h5``. Defaults to the results file
        recorded in ``LIMO.h5`` or to ``LIMO.dir / limo_results.h5``.
    ``h0_file``:
        Optional path to ``limo_H0.h5``.
    ``limo``:
        Optional preloaded ``LIMO`` structure. When omitted, the function
        reads it from ``LIMO.h5``.
    ``stat_name``:
        Optional single statistic name to process. When omitted, all TFCE-
        eligible GLM outputs in ``limo_results.h5`` are processed.
    ``checkfile``:
        Compatibility flag from MATLAB. ``"yes"`` is the documented default.
        Because this Python port is non-interactive, existing TFCE files are
        overwritten with a warning instead of prompting.
    ``return_thresholded_maps``:
        If ``True``, return the intermediate thresholded maps that correspond
        to the ``thresholded_maps`` output of MATLAB ``limo_tfce_handling``.
    ``E``, ``H``, ``dh``:
        TFCE parameters. Defaults match MATLAB: ``E=0.5``, ``H=2``,
        ``dh=0.1``.

Returned value from ``limo_tfce``:
    A dictionary keyed by statistic name. Each value contains:
        ``tfce_score_path``:
            Path to the observed TFCE HDF5 file.
        ``h0_tfce_score_path``:
            Path to the H0 TFCE HDF5 file, or ``None`` when no H0 file is
            available.
        ``tfce_score``:
            Observed TFCE score array.
        ``thresholded_maps``:
            Only included when ``return_thresholded_maps=True``.

Outputs on disk:
    ``tfce/<stat_name>_tfce.h5``:
        HDF5 file with dataset ``tfce_score``.
    ``H0/<stat_name>_tfce_H0.h5``:
        HDF5 file with dataset ``tfce_H0_score`` when an H0 file exists.

References:
    Smith, S. M., & Nichols, T. E. (2009). Threshold-free cluster enhancement:
    Addressing problems of smoothing, threshold dependence and localisation in
    cluster inference. NeuroImage, 44(1), 83-98.
    Pernet, C., Latinus, M., Nichols, T. E., & Rousselet, G. A. (2015).
    Cluster-based computational methods for mass univariate analyses of event-
    related brain potentials/fields: a simulation study. Journal of
    Neuroscience Methods, 250, 85-93.
    Pernet, C., & Rousselet, G. (2014). Type 1 error rate using TFCE for ERP.
    figshare. http://dx.doi.org/10.6084/m9.figshare.1008325

Notes:
    - The direct numerical TFCE transform is implemented in
      ``limo_tfce_transform`` and follows the MATLAB branching by data type,
      sign handling, and threshold integration.
    - The HDF5 wrapper stores one TFCE output file per statistic rather than
      one MATLAB ``.mat`` file per statistic.

Command-line usage:
    python limo_tfce.py path/to/LIMO.h5
    python limo_tfce.py path/to/LIMO.h5 --stat-name R2
    python limo_tfce.py path/to/LIMO.h5 --results-file path/to/limo_results.h5 --h0-file path/to/limo_H0.h5
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from scipy import ndimage

from .eeglab_import import write_hdf5_structure
from .limo_design import read_hdf5_structure


def limo_tfce(
    limo_file: str | Path,
    results_file: str | Path | None = None,
    h0_file: str | Path | None = None,
    limo: Mapping[str, Any] | None = None,
    *,
    stat_name: str | None = None,
    checkfile: str = "yes",
    return_thresholded_maps: bool = False,
    E: float = 0.5,
    H: float = 2.0,
    dh: float = 0.1,
) -> dict[str, dict[str, Any]]:
    """HDF5-aware wrapper equivalent to ``limo_tfce_handling.m``."""

    limo_path = Path(limo_file).expanduser().resolve()
    limo_payload = read_hdf5_structure(limo_path)
    if limo is None:
        if "LIMO" not in limo_payload or not isinstance(limo_payload["LIMO"], Mapping):
            raise KeyError("The HDF5 file does not contain a top-level 'LIMO' group.")
        limo_data = _clone_mapping(limo_payload["LIMO"])
    else:
        limo_data = _clone_mapping(limo)

    results_path = _resolve_results_path(limo_data, results_file)
    results_payload = read_hdf5_structure(results_path)

    h0_path = None if h0_file is None else Path(h0_file).expanduser().resolve()
    h0_payload = read_hdf5_structure(h0_path) if h0_path is not None and h0_path.exists() else {}

    design = limo_data.get("design", {}) if isinstance(limo_data.get("design", {}), Mapping) else {}
    design["tfce"] = 1
    limo_data["design"] = design
    write_hdf5_structure(limo_path, {"LIMO": limo_data})

    tfce_dir = Path(str(limo_data["dir"])).expanduser().resolve() / "tfce"
    tfce_dir.mkdir(parents=True, exist_ok=True)
    h0_dir = Path(str(limo_data["dir"])).expanduser().resolve() / "H0"
    h0_dir.mkdir(parents=True, exist_ok=True)

    stat_names = [stat_name] if stat_name is not None else _discover_tfce_statistics(results_payload)
    outputs: dict[str, dict[str, Any]] = {}
    for current_name in stat_names:
        if current_name not in results_payload:
            raise KeyError(f"Statistic '{current_name}' is not present in {results_path}.")

        observed_source = np.asarray(results_payload[current_name], dtype=float)
        analysis = str(limo_data.get("Analysis", ""))
        neighbouring_matrix = _resolve_neighbouring_matrix(limo_data, observed_source, analysis)

        kind, observed_map = _select_observed_statistic(current_name, observed_source, analysis)
        tfce_score, thresholded_maps = limo_tfce_transform(
            kind,
            observed_map,
            neighbouring_matrix,
            1,
            E,
            H,
            dh,
        )

        tfce_score_path = tfce_dir / f"{current_name}_tfce.h5"
        _prepare_output_path(tfce_score_path, checkfile)
        write_hdf5_structure(tfce_score_path, {"tfce_score": tfce_score})

        entry: dict[str, Any] = {
            "tfce_score_path": tfce_score_path,
            "h0_tfce_score_path": None,
            "tfce_score": tfce_score,
        }

        h0_name = f"H0_{current_name}"
        if h0_name in h0_payload:
            h0_source = np.asarray(h0_payload[h0_name], dtype=float)
            h0_map = _select_h0_statistic(current_name, h0_source, analysis)
            h0_tfce_score, h0_thresholded_maps = limo_tfce_transform(
                kind,
                h0_map,
                neighbouring_matrix,
                0,
                E,
                H,
                dh,
            )
            h0_tfce_score_path = h0_dir / f"{current_name}_tfce_H0.h5"
            _prepare_output_path(h0_tfce_score_path, checkfile)
            write_hdf5_structure(h0_tfce_score_path, {"tfce_H0_score": h0_tfce_score})
            entry["h0_tfce_score_path"] = h0_tfce_score_path
            if return_thresholded_maps:
                entry["thresholded_maps"] = [thresholded_maps, h0_thresholded_maps]
        elif return_thresholded_maps:
            entry["thresholded_maps"] = thresholded_maps

        outputs[current_name] = entry

    return outputs


def limo_tfce_transform(
    data_type: int,
    data: np.ndarray,
    channeighbstructmat: np.ndarray | None,
    updatebar: int = 1,
    E: float = 0.5,
    H: float = 2.0,
    dh: float = 0.1,
) -> tuple[np.ndarray, Any]:
    """Direct array-oriented port of MATLAB ``limo_tfce.m``."""

    if data_type not in {1, 2, 3}:
        raise ValueError("type must be 1, 2, or 3")

    array = np.asarray(data, dtype=float)
    thresholded_maps: Any = []

    if data_type == 1:
        if array.ndim == 1:
            tfce_score, thresholded_maps = _tfce_single_map(
                array,
                _label_1d,
                positive_nonnegative=True,
                E=E,
                H=H,
                dh=dh,
            )
            return tfce_score, thresholded_maps
        if array.ndim != 2:
            raise ValueError("type 1 expects a vector or a 2D bootstrap stack")
        return _tfce_bootstrap_stack(
            array,
            _label_1d,
            positive_nonnegative=False,
            E=E,
            H=H,
            dh=dh,
        )

    if data_type == 2:
        labeler = _build_type2_labeler(channeighbstructmat)
        if array.ndim == 2:
            return _tfce_single_map(array, labeler, positive_nonnegative=False, E=E, H=H, dh=dh)
        if array.ndim != 3:
            raise ValueError("type 2 expects a 2D map or a 3D bootstrap stack")
        return _tfce_bootstrap_stack(
            array,
            labeler,
            positive_nonnegative=False,
            E=E,
            H=H,
            dh=dh,
        )

    labeler = _build_type3_labeler(channeighbstructmat)
    if array.ndim == 3:
        return _tfce_single_map(array, labeler, positive_nonnegative=False, E=E, H=H, dh=dh)
    if array.ndim != 4:
        raise ValueError("type 3 expects a 3D map or a 4D bootstrap stack")
    return _tfce_bootstrap_stack(
        array,
        labeler,
        positive_nonnegative=False,
        E=E,
        H=H,
        dh=dh,
    )


def _resolve_results_path(limo: Mapping[str, Any], results_file: str | Path | None) -> Path:
    if results_file is not None:
        return Path(results_file).expanduser().resolve()
    data = limo.get("data", {}) if isinstance(limo.get("data", {}), Mapping) else {}
    if "results_file" in data:
        return Path(str(data["results_file"])).expanduser().resolve()
    return Path(str(limo["dir"])).expanduser().resolve() / "limo_results.h5"


def _discover_tfce_statistics(results_payload: Mapping[str, Any]) -> list[str]:
    preferred = []
    for key in results_payload:
        if key == "R2":
            preferred.append(key)
        elif key.startswith("Condition_effect_"):
            preferred.append(key)
        elif key.startswith("Interaction_effect_"):
            preferred.append(key)
        elif key.startswith("Covariate_effect_"):
            preferred.append(key)
        elif "semi_partial" in key:
            preferred.append(key)
        elif key.startswith("con"):
            preferred.append(key)
    return preferred


def _resolve_neighbouring_matrix(limo: Mapping[str, Any], observed_source: np.ndarray, analysis: str) -> np.ndarray | None:
    if observed_source.shape[0] == 1 and analysis.lower() == "time-frequency":
        return None
    data = limo.get("data", {}) if isinstance(limo.get("data", {}), Mapping) else {}
    matrix = data.get("neighbouring_matrix")
    if matrix is None:
        if observed_source.shape[0] == 1:
            return None
        raise ValueError("TFCE requires LIMO.data.neighbouring_matrix for multi-channel data.")
    array = np.asarray(matrix)
    return None if array.size == 0 else array.astype(bool)


def _prepare_output_path(path: Path, checkfile: str) -> None:
    if path.exists() and str(checkfile).lower() == "yes":
        warnings.warn(f"Overwriting existing TFCE file: {path}", RuntimeWarning, stacklevel=2)


def _select_observed_statistic(stat_name: str, data: np.ndarray, analysis: str) -> tuple[int, np.ndarray]:
    kind = _tfce_type_for_data(data, analysis)
    if stat_name == "R2" or "semi_partial" in stat_name:
        return kind, np.asarray(data[..., 1], dtype=float)
    if _is_t_statistic(stat_name, data):
        return kind, np.asarray(data[..., -2], dtype=float)
    return kind, np.asarray(data[..., 0], dtype=float)


def _select_h0_statistic(stat_name: str, data: np.ndarray, analysis: str) -> np.ndarray:
    _ = analysis
    if stat_name == "R2" or "semi_partial" in stat_name:
        return np.asarray(data[..., 1, :], dtype=float)
    if _is_t_statistic(stat_name, data):
        return np.asarray(data[..., -2, :], dtype=float)
    return np.asarray(data[..., 0, :], dtype=float)


def _is_t_statistic(stat_name: str, data: np.ndarray) -> bool:
    return stat_name.startswith("con") or data.shape[-1] >= 5


def _tfce_type_for_data(data: np.ndarray, analysis: str) -> int:
    if data.shape[0] == 1:
        return 2 if analysis.lower() == "time-frequency" else 1
    return 3 if analysis.lower() == "time-frequency" else 2


def _tfce_bootstrap_stack(
    data: np.ndarray,
    labeler: Any,
    *,
    positive_nonnegative: bool,
    E: float,
    H: float,
    dh: float,
) -> tuple[np.ndarray, list[np.ndarray]]:
    score = np.full(data.shape, np.nan, dtype=float)
    thresholded_maps: list[np.ndarray] = []
    for boot in range(data.shape[-1]):
        boot_score, boot_maps = _tfce_single_map(
            np.asarray(data[..., boot], dtype=float),
            labeler,
            positive_nonnegative=positive_nonnegative,
            E=E,
            H=H,
            dh=dh,
        )
        score[..., boot] = boot_score
        thresholded_maps.append(boot_maps)
    return score, thresholded_maps


def _tfce_single_map(
    data: np.ndarray,
    labeler: Any,
    *,
    positive_nonnegative: bool,
    E: float,
    H: float,
    dh: float,
) -> tuple[np.ndarray, np.ndarray]:
    data = np.asarray(data, dtype=float)
    data_range = float(np.nanmax(data) - np.nanmin(data)) if data.size else 0.0
    if data_range == 0 or np.all(np.isnan(data)):
        return np.zeros_like(data, dtype=float), np.zeros(data.shape + (0,), dtype=float)

    increment = _compute_increment(data_range, dh)
    minimum = float(np.nanmin(data))
    use_positive_only = minimum >= 0 if positive_nonnegative else minimum > 0
    if use_positive_only:
        tfce = _integrate_branch(data, increment, labeler, E, H)
        return np.nansum(tfce, axis=-1), _trim_thresholded_maps(tfce)

    pos_data = np.where(data > 0, data, 0.0)
    neg_data = np.abs(np.where(data < 0, data, 0.0))

    pos_tfce = _integrate_branch_signed(pos_data, increment, labeler, E, H)
    neg_tfce = _integrate_branch_signed(neg_data, increment, labeler, E, H, negative=True)
    score = np.nansum(pos_tfce, axis=-1) + np.nansum(neg_tfce, axis=-1)
    combined = _combine_signed_thresholds(pos_tfce, neg_tfce)
    return score, _trim_thresholded_maps(combined)


def _compute_increment(data_range: float, dh: float) -> float:
    if data_range > 1:
        precision = round(data_range / dh)
        if precision > 200:
            return data_range / 200
        if precision == 0:
            return data_range
        return data_range / precision
    increment = data_range * dh
    return increment if increment > 0 else data_range or dh


def _integrate_branch(data: np.ndarray, increment: float, labeler: Any, E: float, H: float) -> np.ndarray:
    heights = _matlab_range(float(np.nanmin(data)), increment, float(np.nanmax(data)))
    tfce = np.full(data.shape + (len(heights),), np.nan, dtype=float)
    for index, height in enumerate(heights):
        clustered_map, num = labeler(np.asarray(data > height, dtype=bool))
        extent_map = _integrate_clusters(clustered_map, num)
        tfce[..., index] = (extent_map**E) * (height**H) * increment
    return tfce


def _integrate_branch_signed(
    data: np.ndarray,
    increment: float,
    labeler: Any,
    E: float,
    H: float,
    *,
    negative: bool = False,
) -> np.ndarray:
    if np.nanmax(data) == np.nanmin(data):
        return np.zeros(data.shape + (0,), dtype=float)
    n_levels = len(_matlab_range(float(np.nanmin(data)), increment, float(np.nanmax(data))))
    if negative:
        n_levels = max(n_levels - 1, 0)
    if n_levels == 0:
        return np.zeros(data.shape + (0,), dtype=float)
    branch_increment = (float(np.nanmax(data)) - float(np.nanmin(data))) / n_levels
    heights = _matlab_range(float(np.nanmin(data)), branch_increment, float(np.nanmax(data)))
    if negative and heights.size > n_levels:
        heights = heights[:n_levels]
    tfce = np.full(data.shape + (len(heights),), np.nan, dtype=float)
    for index, height in enumerate(heights):
        clustered_map, num = labeler(np.asarray(data > height, dtype=bool))
        extent_map = _integrate_clusters(clustered_map, num)
        tfce[..., index] = (extent_map**E) * (height**H) * increment
    return tfce


def _combine_signed_thresholds(pos_tfce: np.ndarray, neg_tfce: np.ndarray) -> np.ndarray:
    out = np.full(pos_tfce.shape[:-1] + (neg_tfce.shape[-1] + pos_tfce.shape[-1],), np.nan, dtype=float)
    if neg_tfce.shape[-1] > 0:
        out[..., : neg_tfce.shape[-1]] = neg_tfce[..., ::-1]
    if pos_tfce.shape[-1] > 0:
        out[..., neg_tfce.shape[-1] :] = pos_tfce
    return out


def _trim_thresholded_maps(thresholded_maps: np.ndarray) -> np.ndarray:
    if thresholded_maps.size == 0:
        return thresholded_maps
    reduce_axes = tuple(range(thresholded_maps.ndim - 1))
    keep = np.nansum(thresholded_maps, axis=reduce_axes) != 0
    return thresholded_maps[..., keep]


def _label_1d(mask: np.ndarray) -> tuple[np.ndarray, int]:
    mask = np.asarray(mask, dtype=bool).ravel()
    labels = np.zeros(mask.shape[0], dtype=int)
    cluster_id = 0
    active = False
    for index, value in enumerate(mask):
        if value and not active:
            cluster_id += 1
            active = True
        elif not value:
            active = False
        if value:
            labels[index] = cluster_id
    return labels.reshape(mask.shape), cluster_id


def _build_type2_labeler(channeighbstructmat: np.ndarray | None):
    if channeighbstructmat is None or np.asarray(channeighbstructmat).size == 0:
        structure = ndimage.generate_binary_structure(2, 1)
        return lambda mask: ndimage.label(np.asarray(mask, dtype=bool), structure=structure)

    neighbour = np.asarray(channeighbstructmat, dtype=bool)

    def labeler(mask: np.ndarray) -> tuple[np.ndarray, int]:
        cluster, num = _limo_findcluster(np.asarray(mask, dtype=bool)[:, :, np.newaxis], neighbour, 2)
        return cluster[:, :, 0], num

    return labeler


def _build_type3_labeler(channeighbstructmat: np.ndarray | None):
    if channeighbstructmat is None or np.asarray(channeighbstructmat).size == 0:
        structure = ndimage.generate_binary_structure(3, 1)
        return lambda mask: ndimage.label(np.asarray(mask, dtype=bool), structure=structure)

    neighbour = np.asarray(channeighbstructmat, dtype=bool)

    def labeler(mask: np.ndarray) -> tuple[np.ndarray, int]:
        return _limo_findcluster(np.asarray(mask, dtype=bool), neighbour, 2)

    return labeler


def _limo_findcluster(onoff: np.ndarray, spatdimneighbstructmat: np.ndarray, minnbchan: int = 2) -> tuple[np.ndarray, int]:
    onoff = np.asarray(onoff, dtype=bool).copy()
    spatdimlength, nfreq, ntime = onoff.shape
    neighbour = np.asarray(spatdimneighbstructmat, dtype=bool)
    if neighbour.shape != (spatdimlength, spatdimlength):
        raise ValueError("invalid dimension of spatdimneighbstructmat")

    if minnbchan > 0:
        selectmat = np.asarray(neighbour | neighbour.T, dtype=float)
        nremoved = 1
        while nremoved > 0:
            nsigneighb = (selectmat @ onoff.reshape(spatdimlength, nfreq * ntime).astype(float)).reshape(onoff.shape)
            remove = (onoff.astype(float) * nsigneighb) < minnbchan
            nremoved = int(np.count_nonzero(remove & onoff))
            onoff[remove] = False

    labelmat = np.zeros(onoff.shape, dtype=int)
    total = 0
    structure = ndimage.generate_binary_structure(2, 1)
    for spatdimlev in range(spatdimlength):
        labels, num = ndimage.label(onoff[spatdimlev], structure=structure)
        labels[labels != 0] += total
        labelmat[spatdimlev] = labels
        total += num
    if total == 0:
        return np.zeros(onoff.shape, dtype=int), 0

    flat = labelmat.reshape(spatdimlength, nfreq * ntime)
    replaceby = np.arange(1, total + 1, dtype=int)

    for spatdimlev in range(spatdimlength):
        neighbours = np.flatnonzero(neighbour[spatdimlev])
        for nbindx in neighbours:
            overlap = np.flatnonzero((flat[spatdimlev] != 0) & (flat[nbindx] != 0))
            for index in overlap:
                a = flat[spatdimlev, index]
                b = flat[nbindx, index]
                rep_a = replaceby[a - 1]
                rep_b = replaceby[b - 1]
                if rep_a == rep_b:
                    continue
                if rep_a < rep_b:
                    replaceby[replaceby == rep_b] = rep_a
                else:
                    replaceby[replaceby == rep_a] = rep_b

    cluster = np.zeros(flat.shape, dtype=int)
    num = 0
    for label in np.unique(replaceby):
        num += 1
        here = np.flatnonzero(replaceby == label) + 1
        cluster[np.isin(flat, here)] = num
    return cluster.reshape(onoff.shape), num


def _integrate_clusters(clustered_map: np.ndarray, num: int) -> np.ndarray:
    if num == 0:
        return np.zeros_like(clustered_map, dtype=float)
    flat = np.asarray(clustered_map, dtype=int).ravel()
    counts = np.bincount(flat, minlength=num + 1)
    extent = counts[flat]
    extent[flat == 0] = 0
    return extent.reshape(clustered_map.shape).astype(float)


def _matlab_range(start: float, step: float, stop: float) -> np.ndarray:
    if step <= 0 or not np.isfinite(step):
        return np.asarray([start], dtype=float)
    count = int(np.floor((stop - start) / step + 1e-12)) + 1
    values = start + step * np.arange(max(count, 1), dtype=float)
    if values.size == 0:
        return np.asarray([start], dtype=float)
    if values[-1] < stop - (step * 1e-9):
        values = np.append(values, stop)
    return values


def _clone_mapping(mapping: Mapping[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in mapping.items():
        if isinstance(value, Mapping):
            out[key] = _clone_mapping(value)
        elif isinstance(value, list):
            out[key] = list(value)
        elif isinstance(value, np.ndarray):
            out[key] = np.array(value, copy=True)
        else:
            out[key] = value
    return out


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Apply TFCE to LIMO HDF5 statistic outputs")
    parser.add_argument("limo_file", type=Path, help="Path to LIMO.h5")
    parser.add_argument("--results-file", type=Path, default=None, help="Optional path to limo_results.h5")
    parser.add_argument("--h0-file", type=Path, default=None, help="Optional path to limo_H0.h5")
    parser.add_argument("--stat-name", type=str, default=None, help="Optional single statistic to process")
    parser.add_argument("--checkfile", choices=["yes", "no"], default="yes", help="Compatibility overwrite flag")
    parser.add_argument("--return-thresholded-maps", action="store_true", help="Include thresholded maps in the JSON summary")
    parser.add_argument("--E", type=float, default=0.5, help="TFCE extent exponent")
    parser.add_argument("--H", type=float, default=2.0, help="TFCE height exponent")
    parser.add_argument("--dh", type=float, default=0.1, help="TFCE integration step")
    return parser


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()
    outputs = limo_tfce(
        args.limo_file,
        args.results_file,
        args.h0_file,
        stat_name=args.stat_name,
        checkfile=args.checkfile,
        return_thresholded_maps=args.return_thresholded_maps,
        E=args.E,
        H=args.H,
        dh=args.dh,
    )
    printable = {
        key: {
            inner_key: (str(value) if isinstance(value, Path) else value if inner_key == "h0_tfce_score_path" and value is None else str(value) if isinstance(value, Path) else value)
            for inner_key, value in item.items()
            if inner_key != "thresholded_maps"
        }
        for key, item in outputs.items()
    }
    print(json.dumps(printable, indent=2))


if __name__ == "__main__":
    main()