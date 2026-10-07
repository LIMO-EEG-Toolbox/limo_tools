"""Create the LIMO design matrix and placeholder result arrays in HDF5.

Use ``limo_design(limo_file)`` after ``eeglab_import.py`` has exported
``LIMO.h5`` from an EEGLAB ``.set`` file.

The function reads ``LIMO.h5``, reloads the source EEGLAB dataset, rebuilds the
analysis-specific observation matrix ``Y``, creates the design matrix ``X``,
updates the ``LIMO`` structure, and writes a separate ``limo_results.h5`` file.

Inputs:
    ``limo_file``: path to the exported ``LIMO.h5`` file.
    ``results_path``: optional explicit path for ``limo_results.h5``.

Returned value:
    A tuple ``(updated_limo_path, results_path)`` of ``Path`` objects.

Outputs on disk:
    ``LIMO.h5``: updated with ``LIMO.design.X``, the factor counts, design
    name, analysis status, and time-frequency size metadata when relevant.
    ``limo_results.h5``: contains root datasets ``/Yr``, ``/Yhat``, ``/Res``,
    ``/R2``, and ``/Beta``.

Result datasets in ``limo_results.h5``:
    ``/Yr``: the EEG data from the ``.set`` file, reorganized to fit ``X``;
    when categorical regressors are present, observations are grouped by
    condition.
    ``/Yhat``: the predicted data, with the same shape as ``/Yr``.
    ``/Res``: the residual, non-modeled data, with the same shape as ``/Yr``.
    ``/R2``: the model fit array.
    ``/Beta``: the beta values, with shape
    ``channels x frames x number_of_parameters`` for time/frequency analyses,
    or the corresponding time-frequency extension for TF data.

Command-line usage:
    python limo_design.py path/to/LIMO.h5
    python limo_design.py path/to/LIMO.h5 --results-output path/to/limo_results.h5
"""

from __future__ import annotations

import argparse
import itertools
import json
import re
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

try:
    import h5py
except ImportError:  # pragma: no cover
    h5py = None  # type: ignore[assignment]

from .eeglab_import import write_hdf5_structure
from .read_setfile import load_matlab_file, read_set_file


SIGNAL_FIELD_PATTERN = re.compile(r"^(chan|comp)(\d+)$", re.IGNORECASE)


def limo_design(
    limo_file: str | Path,
    *,
    results_path: str | Path | None = None,
) -> tuple[Path, Path]:
    """Update ``LIMO.h5`` with the design matrix and export ``limo_results.h5``."""

    limo_path = Path(limo_file).expanduser().resolve()
    payload = read_hdf5_structure(limo_path)
    if "LIMO" not in payload or not isinstance(payload["LIMO"], Mapping):
        raise KeyError("The HDF5 file does not contain a top-level 'LIMO' group.")

    limo = dict(payload["LIMO"])
    eeg = _load_eeg_for_limo(limo)
    y_native = _load_analysis_data(limo, eeg)
    analysis = str(limo.get("Analysis", ""))

    if analysis == "Time-Frequency":
        if y_native.ndim != 4:
            raise ValueError("Time-Frequency data must be 4D: electrodes x freqs x times x observations.")
        size4d = tuple(int(value) for value in y_native.shape)
        y_working = flatten_tf(y_native)
    else:
        y_working = ensure_3d(y_native)
        size4d = None

    design = _build_design_from_observations(y_working, limo)
    yr_working = design["Yr"]
    yr_working = _apply_expected_chanlocs(yr_working, limo)

    if analysis == "Time-Frequency":
        if size4d is None:
            raise ValueError("Missing 4D size metadata for Time-Frequency analysis.")
        size4d = (yr_working.shape[0], size4d[1], size4d[2], yr_working.shape[2])
        size3d = (yr_working.shape[0], yr_working.shape[1], yr_working.shape[2])
        yr_to_store = unflatten_tf(yr_working, size4d)
    else:
        size3d = tuple(int(value) for value in yr_working.shape)
        yr_to_store = yr_working

    results_output = (
        Path(results_path).expanduser().resolve()
        if results_path
        else Path(str(limo["dir"])).expanduser().resolve() / "limo_results.h5"
    )
    results_output.parent.mkdir(parents=True, exist_ok=True)

    results_payload = _build_results_payload(
        yr=yr_to_store,
        x=design["X"],
        analysis=analysis,
        type_of_analysis=str(limo.get("design", {}).get("type_of_analysis", "Mass-univariate")),
    )
    write_hdf5_structure(results_output, results_payload)

    updated_limo = _update_limo_structure(
        limo=limo,
        design=design,
        results_output=results_output,
        size3d=size3d,
        size4d=size4d,
    )
    write_hdf5_structure(limo_path, {"LIMO": updated_limo})
    return limo_path, results_output


def read_hdf5_structure(file_path: str | Path) -> dict[str, Any]:
    """Read a generic HDF5 structure written by ``write_hdf5_structure``."""

    if h5py is None:  # pragma: no cover
        raise ImportError("h5py is required to read LIMO HDF5 files. Install h5py first.")

    path = Path(file_path).expanduser().resolve()
    with h5py.File(path, "r") as handle:
        return {key: _read_hdf5_item(handle[key]) for key in handle.keys()}


def flatten_tf(data_4d: np.ndarray) -> np.ndarray:
    """Match ``limo_tf_4d_reshape`` for 4D -> 3D conversions."""

    data_4d = np.asarray(data_4d)
    if data_4d.ndim != 4:
        raise ValueError("flatten_tf expects a 4D array.")

    n_electrodes, n_freqs, n_times, n_obs = data_4d.shape
    flattened = np.empty((n_electrodes, n_freqs * n_times, n_obs), dtype=data_4d.dtype)
    for obs in range(n_obs):
        for time_index in range(n_times):
            start = time_index * n_freqs
            stop = start + n_freqs
            flattened[:, start:stop, obs] = data_4d[:, :, time_index, obs]
    return flattened


def unflatten_tf(data_3d: np.ndarray, size4d: Iterable[int]) -> np.ndarray:
    """Match ``limo_tf_4d_reshape`` for 3D -> 4D conversions."""

    data_3d = ensure_3d(data_3d)
    n_electrodes, n_freqs, n_times, n_obs = tuple(int(value) for value in size4d)
    if data_3d.shape != (n_electrodes, n_freqs * n_times, n_obs):
        raise ValueError("3D and 4D shapes disagree for frequency-time reshaping.")

    restored = np.empty((n_electrodes, n_freqs, n_times, n_obs), dtype=data_3d.dtype)
    for obs in range(n_obs):
        for time_index in range(n_times):
            start = time_index * n_freqs
            stop = start + n_freqs
            restored[:, :, time_index, obs] = data_3d[:, start:stop, obs]
    return restored


def _load_eeg_for_limo(limo: Mapping[str, Any]) -> dict[str, Any]:
    set_path = Path(str(limo["data"]["data_dir"])) / str(limo["data"]["data"])
    return read_set_file(set_path, load_data=False)


def _load_analysis_data(limo: Mapping[str, Any], eeg: dict[str, Any]) -> np.ndarray:
    analysis = str(limo.get("Analysis", ""))
    data_type = str(limo.get("Type", ""))
    data = limo["data"]

    if data_type == "Components" and not _is_empty_value(data.get("cluster")) and int(np.asarray(data.get("cluster")).item()) != 0:
        raise NotImplementedError(
            "Component-cluster reordering from EEGLAB STUDY metadata is not implemented in limo_design.py yet."
        )

    trim1 = int(np.asarray(data.get("trim1", 1)).item()) - 1
    trim2 = int(np.asarray(data.get("trim2", 0)).item())

    if analysis == "Time":
        if data_type == "Channels":
            signal = _load_signal_source_for_analysis(eeg, ["daterp"], fallback_keys=["data"])
        else:
            signal = _load_signal_source_for_analysis(eeg, ["icaerp"], allow_ica_fallback=True)
        signal = ensure_3d(signal)
        return signal[:, trim1:trim2, :]

    if analysis == "Frequency":
        if data_type == "Channels":
            signal = _load_signal_source_for_analysis(eeg, ["datspec"], fallback_keys=["specdata"])
        else:
            signal = _load_signal_source_for_analysis(eeg, ["icaspec"], fallback_keys=["specicaact"])
        signal = ensure_3d(signal)
        return signal[:, trim1:trim2, :]

    if analysis == "Time-Frequency":
        trim_lowf = int(np.asarray(data.get("trim_lowf", 1)).item()) - 1
        trim_highf = int(np.asarray(data.get("trim_highf", 0)).item())

        if data_type == "Channels":
            signal = _load_signal_source_for_analysis(eeg, ["dattimef", "datersp"])
        else:
            signal = _load_signal_source_for_analysis(eeg, ["icatimef", "icaersp"])

        signal = ensure_4d(signal)
        return np.abs(signal[:, trim_lowf:trim_highf, trim1:trim2, :]) ** 2

    raise ValueError("Unsupported analysis type in LIMO: expected Time, Frequency, or Time-Frequency.")


def _load_signal_source_for_analysis(
    eeg: dict[str, Any],
    datafile_keys: list[str],
    *,
    fallback_keys: list[str] | None = None,
    allow_ica_fallback: bool = False,
) -> np.ndarray:
    set_path = Path(str(eeg["set_path"]))
    etc = eeg.get("etc", {}) if isinstance(eeg.get("etc", {}), Mapping) else {}
    datafiles = etc.get("datafiles", {}) if isinstance(etc.get("datafiles", {}), Mapping) else {}

    for key in datafile_keys:
        source = datafiles.get(key)
        if _is_empty_value(source):
            continue
        materialized = _materialize_signal_source(source, set_path.parent)
        if materialized is not None:
            return np.asarray(materialized)

    if fallback_keys:
        for key in fallback_keys:
            if key in eeg and not _is_empty_value(eeg[key]):
                return np.asarray(eeg[key])
            if key == "data":
                eeg_with_data = read_set_file(set_path, load_data=True)
                if key in eeg_with_data and not _is_empty_value(eeg_with_data[key]):
                    return np.asarray(eeg_with_data[key])

    if allow_ica_fallback:
        eeg_with_data = read_set_file(set_path, load_data=True)
        return _compute_ica_activations(eeg_with_data)

    raise FileNotFoundError(
        f"Could not locate EEGLAB measure data for keys {', '.join(datafile_keys)}."
    )


def _materialize_signal_source(source: Any, base_dir: Path) -> Any:
    if source is None:
        return None

    if isinstance(source, np.ndarray):
        return source

    if isinstance(source, Mapping):
        return _stack_signal_payload(source, base_dir)

    if isinstance(source, (list, tuple)):
        pieces = []
        for item in source:
            materialized = _materialize_signal_source(item, base_dir)
            if materialized is not None:
                pieces.append(np.asarray(materialized))
        if not pieces:
            return None
        if len(pieces) == 1:
            return pieces[0]
        return np.concatenate(pieces, axis=pieces[0].ndim - 1)

    if isinstance(source, (str, Path)):
        path = Path(source)
        if not path.is_absolute():
            path = (base_dir / path).resolve()
        if not path.exists():
            raise FileNotFoundError(f"Referenced EEGLAB data file does not exist: {path}")
        payload = load_matlab_file(path)
        return _stack_signal_payload(payload, path.parent)

    return np.asarray(source)


def _stack_signal_payload(payload: Mapping[str, Any], base_dir: Path) -> Any:
    keys = []
    for key in payload.keys():
        match = SIGNAL_FIELD_PATTERN.match(str(key))
        if match and "base" not in str(key).lower() and "boot" not in str(key).lower() and "label" not in str(key).lower():
            keys.append((int(match.group(2)), str(key)))

    if keys:
        arrays = [np.asarray(_materialize_signal_source(payload[key], base_dir)) for _, key in sorted(keys)]
        first = arrays[0]
        return np.stack(arrays, axis=0) if first.ndim in {1, 2, 3} else np.asarray(arrays)

    if "data" in payload:
        return np.asarray(payload["data"])

    numeric_values = [np.asarray(value) for value in payload.values() if isinstance(value, np.ndarray)]
    if len(numeric_values) == 1:
        return numeric_values[0]

    if len(payload) == 1:
        return _materialize_signal_source(next(iter(payload.values())), base_dir)

    raise ValueError("Could not infer the signal matrix from the EEGLAB payload.")


def _compute_ica_activations(eeg: Mapping[str, Any]) -> np.ndarray:
    data = np.asarray(eeg.get("data"))
    weights = np.asarray(eeg.get("icaweights"))
    sphere = np.asarray(eeg.get("icasphere"))

    if data.size == 0 or weights.size == 0 or sphere.size == 0:
        raise ValueError("ICA fallback requires EEG.data, EEG.icaweights, and EEG.icasphere.")

    unmixing = weights @ sphere
    data = np.asarray(data)
    if data.ndim == 2:
        return unmixing @ data
    if data.ndim == 3:
        n_channels, n_frames, n_trials = data.shape
        flattened = data.reshape(n_channels, n_frames * n_trials, order="F")
        activations = unmixing @ flattened
        return activations.reshape(activations.shape[0], n_frames, n_trials, order="F")
    raise ValueError("EEG.data must be 2D or 3D to compute ICA activations.")


def _build_design_from_observations(y: np.ndarray, limo: Mapping[str, Any]) -> dict[str, Any]:
    y = ensure_3d(y)
    design_info = limo.get("design", {}) if isinstance(limo.get("design", {}), Mapping) else {}
    cat = _normalize_regressor_matrix(limo.get("data", {}).get("Cat"), y.shape[2], mode="cat")
    cont = _normalize_regressor_matrix(limo.get("data", {}).get("Cont"), y.shape[2], mode="cont")

    full_factorial = int(np.asarray(design_info.get("fullfactorial", 0)).item())
    zscore = int(np.asarray(design_info.get("zscore", 0)).item())

    if cat.size and np.any(np.isnan(cat)):
        keep = ~np.any(np.isnan(cat), axis=1)
        cat = cat[keep]
        if cont.size:
            cont = cont[keep]
        y = y[:, :, keep]

    if cont.size:
        keep = ~np.any(np.isnan(cont), axis=1)
        if not np.all(keep):
            cont = cont[keep]
            if cat.size:
                cat = cat[keep]
            y = y[:, :, keep]

    _check_dimensions(y, cat, cont)

    if full_factorial == 1 and cont.size:
        raise ValueError(
            "LIMO does not compute full factorial ANCOVAs with covariates. "
            "Use a factorial ANCOVA without interaction terms."
        )

    if cont.size:
        if cont.shape[1] + 1 >= y.shape[2]:
            raise ValueError("There are too many regressors for the number of observations.")
        if cont.shape[1] > 1:
            xtx = cont.T @ cont
            if np.linalg.matrix_rank(xtx) < xtx.shape[0]:
                raise ValueError("The regression matrix is singular; at least one predictor depends on others.")
            if not np.isfinite(np.linalg.cond(xtx, p=1)):
                raise ValueError("The regression matrix is close to singular.")

    if cont.size and zscore == 1:
        cont = _zscore_columns(cont)

    nb_conditions: list[int] = []
    nb_interactions: list[int] = []
    nb_continuous = int(cont.shape[1]) if cont.size else 0
    yr = y.copy()

    if cat.size:
        sort_index = np.lexsort(tuple(cat[:, column] for column in range(cat.shape[1] - 1, -1, -1)))
        cat = cat[sort_index]
        yr = yr[:, :, sort_index]
        if cont.size:
            cont = cont[sort_index]

        indices_conditions: list[np.ndarray] = []
        for column in range(cat.shape[1]):
            unique_values = np.unique(cat[:, column])
            unique_values = unique_values[~np.isnan(unique_values)]
            nb_conditions.append(int(unique_values.size))
            for value in unique_values:
                indices_conditions.append(np.flatnonzero(cat[:, column] == value))

        main_design = np.zeros((yr.shape[2], len(indices_conditions)), dtype=float)
        for column_index, indices in enumerate(indices_conditions):
            main_design[indices, column_index] = 1.0

        design_matrix = main_design
        if full_factorial == 1 and cat.shape[1] > 1:
            design_matrix, nb_interactions = _make_interactions(main_design, nb_conditions)

        pieces = [design_matrix]
        if cont.size:
            pieces.append(cont)
        pieces.append(np.ones((yr.shape[2], 1), dtype=float))
        x = np.column_stack(pieces)

        if full_factorial == 1 and nb_interactions:
            x, yr, full_factorial = _balance_full_factorial_design(
                x=x,
                yr=yr,
                nb_conditions=nb_conditions,
                nb_interactions=nb_interactions,
                nb_continuous=nb_continuous,
                main_design=main_design,
            )
            if full_factorial == 0:
                nb_interactions = []

    elif cont.size:
        x = np.column_stack((cont, np.ones((y.shape[2], 1), dtype=float)))
        yr = y
    else:
        x = np.ones((y.shape[2], 1), dtype=float)
        yr = y

    return {
        "X": x,
        "Yr": yr,
        "nb_conditions": 0 if not nb_conditions else np.asarray(nb_conditions, dtype=int),
        "nb_interactions": 0 if not nb_interactions else np.asarray(nb_interactions, dtype=int),
        "nb_continuous": nb_continuous,
        "fullfactorial": full_factorial,
    }


def _balance_full_factorial_design(
    *,
    x: np.ndarray,
    yr: np.ndarray,
    nb_conditions: list[int],
    nb_interactions: list[int],
    nb_continuous: int,
    main_design: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
    interaction_offset = int(sum(nb_conditions) + sum(nb_interactions[:-1]))
    interaction_width = int(nb_interactions[-1])
    higher_interaction = x[:, interaction_offset : interaction_offset + interaction_width]

    if higher_interaction.shape[1] != int(np.prod(nb_conditions)):
        fallback = np.column_stack((main_design, np.ones((yr.shape[2], 1), dtype=float)))
        return fallback, yr, 0

    cell_counts = higher_interaction.sum(axis=0).astype(int)
    if np.unique(cell_counts).size <= 1:
        return x, yr, 1

    sample_to_n = int(cell_counts.min())
    if sample_to_n <= 1:
        fallback = np.column_stack((main_design, np.ones((yr.shape[2], 1), dtype=float)))
        return fallback, yr, 0

    rng = np.random.default_rng(0)
    keep_indices = []
    for column in range(higher_interaction.shape[1]):
        indices = np.flatnonzero(higher_interaction[:, column])
        if indices.size > sample_to_n:
            indices = np.sort(rng.choice(indices, size=sample_to_n, replace=False))
        keep_indices.append(indices)

    sampled_index = np.sort(np.concatenate(keep_indices))
    return x[sampled_index], yr[:, :, sampled_index], 1


def _make_interactions(x: np.ndarray, nb_conditions: list[int]) -> tuple[np.ndarray, list[int]]:
    if not nb_conditions:
        return x, []

    factor_slices = []
    index = 0
    for levels in nb_conditions:
        factor_slices.append(x[:, index : index + levels])
        index += levels

    interactions: list[int] = []
    tmp_x = x.copy()
    n_factors = len(nb_conditions)

    for interaction_size in range(2, n_factors + 1):
        for combination in itertools.combinations(range(n_factors), interaction_size):
            current = factor_slices[combination[0]]
            for factor_index in combination[1:]:
                next_factor = factor_slices[factor_index]
                blocks = []
                for column in range(current.shape[1]):
                    blocks.append(current[:, [column]] * next_factor)
                current = np.concatenate(blocks, axis=1)
                current = current[:, current.sum(axis=0) != 0]
            interactions.append(int(current.shape[1]))
            tmp_x = np.concatenate((tmp_x, current), axis=1)

    return tmp_x, interactions


def _apply_expected_chanlocs(yr: np.ndarray, limo: Mapping[str, Any]) -> np.ndarray:
    data = limo.get("data", {}) if isinstance(limo.get("data", {}), Mapping) else {}
    expected = data.get("expected_chanlocs")
    chanlocs = data.get("chanlocs")

    if _is_empty_value(expected) or yr.shape[0] <= 1:
        return yr

    current_labels = _extract_chanloc_labels(chanlocs)
    expected_labels = _extract_chanloc_labels(expected)
    if not current_labels or not expected_labels:
        return yr

    matched = np.full((len(expected_labels), yr.shape[1], yr.shape[2]), np.nan, dtype=float)
    label_to_index = {label.lower(): index for index, label in enumerate(current_labels)}
    for expected_index, label in enumerate(expected_labels):
        current_index = label_to_index.get(label.lower())
        if current_index is not None:
            matched[expected_index] = yr[current_index]
    return matched


def _extract_chanloc_labels(chanlocs: Any) -> list[str]:
    if isinstance(chanlocs, Mapping):
        chanlocs = [chanlocs]
    if not isinstance(chanlocs, list):
        return []
    labels = []
    for entry in chanlocs:
        if isinstance(entry, Mapping) and "labels" in entry:
            labels.append(str(entry["labels"]))
    return labels


def _build_results_payload(
    *,
    yr: np.ndarray,
    x: np.ndarray,
    analysis: str,
    type_of_analysis: str,
) -> dict[str, Any]:
    spatial_shape = yr.shape[:-1]
    if analysis == "Time-Frequency":
        yhat = np.zeros_like(yr, dtype=np.float32)
        res = np.zeros_like(yr, dtype=np.float32)
        beta = np.zeros(spatial_shape + (x.shape[1],), dtype=np.float32)
        r2 = np.zeros(spatial_shape + (3,), dtype=np.float32)
    else:
        yhat = np.full_like(yr, np.nan, dtype=np.float32)
        res = np.full_like(yr, np.nan, dtype=np.float32)
        beta = np.full(spatial_shape + (x.shape[1],), np.nan, dtype=np.float32)
        r2 = np.full(spatial_shape + (3,), np.nan, dtype=np.float32)

    if type_of_analysis != "Mass-univariate":
        r2[:] = np.nan

    return {
        "Yr": np.asarray(yr),
        "Yhat": yhat,
        "Res": res,
        "R2": r2,
        "Beta": beta,
    }


def _update_limo_structure(
    *,
    limo: Mapping[str, Any],
    design: Mapping[str, Any],
    results_output: Path,
    size3d: tuple[int, ...],
    size4d: tuple[int, ...] | None,
) -> dict[str, Any]:
    updated = dict(limo)
    updated_design = dict(updated.get("design", {}))
    updated_data = dict(updated.get("data", {}))

    updated_design["X"] = np.asarray(design["X"], dtype=float)
    updated_design["nb_conditions"] = design["nb_conditions"]
    updated_design["nb_interactions"] = design["nb_interactions"]
    updated_design["nb_continuous"] = int(design["nb_continuous"])
    has_interactions = len(_as_int_list(design["nb_interactions"])) > 0
    updated_design["fullfactorial"] = int(design["fullfactorial"] and has_interactions)
    updated_design["name"] = _design_name(
        nb_conditions=design["nb_conditions"],
        nb_continuous=int(design["nb_continuous"]),
    )
    updated_design["status"] = "to do"

    updated_data["results_file"] = str(results_output)
    if size4d is not None:
        updated_data["size4D"] = np.asarray(size4d, dtype=int)
        updated_data["size3D"] = np.asarray(size3d, dtype=int)

    updated["design"] = updated_design
    updated["data"] = updated_data
    return updated


def _design_name(*, nb_conditions: Any, nb_continuous: int) -> str:
    condition_counts = _as_int_list(nb_conditions)

    if condition_counts and nb_continuous == 0:
        if len(condition_counts) == 1:
            if condition_counts[0] == 2:
                return f"Categorical: T-test i.e. {condition_counts[0]} conditions"
            return f"Categorical: 1 way ANOVA with {condition_counts[0]} conditions"
        return f"Categorical: N way ANOVA with {len(condition_counts)} factors"

    if not condition_counts and nb_continuous > 0:
        if nb_continuous == 1:
            return "Continuous: Simple Regression"
        return f"Continuous: Multiple Regression with {nb_continuous} continuous variables"

    if condition_counts and nb_continuous > 0:
        if len(condition_counts) == 1:
            return (
                f"AnCOVA with {condition_counts[0]} conditions and {nb_continuous} continuous variable(s)"
            )
        return f"AnCOVA with {len(condition_counts)} factors and {nb_continuous} continuous variable(s)"

    return "Mean"


def _check_dimensions(y: np.ndarray, cat: np.ndarray, cont: np.ndarray) -> None:
    n_obs = y.shape[2]
    if cat.size and cat.shape[0] != n_obs:
        raise ValueError("The number of categorical rows does not match the number of observations.")
    if cont.size and cont.shape[0] != n_obs:
        raise ValueError("The number of continuous rows does not match the number of observations.")


def _normalize_regressor_matrix(value: Any, n_obs: int, *, mode: str) -> np.ndarray:
    if value is None:
        return np.empty((0, 0), dtype=float)

    array = np.asarray(value)
    if array.size == 0:
        return np.empty((0, 0), dtype=float)

    if array.shape == ():
        scalar = array.item()
        if scalar == 0:
            return np.empty((0, 0), dtype=float)
        array = np.asarray([[scalar]], dtype=float)
    elif array.ndim == 1:
        array = array.reshape(-1, 1)
    elif array.ndim > 2:
        raise ValueError(f"{mode} regressors must be 1D or 2D.")

    array = np.asarray(array, dtype=float)
    if array.shape[1] == n_obs and array.shape[0] != n_obs:
        array = array.T
    return array


def _zscore_columns(cont: np.ndarray) -> np.ndarray:
    centered = cont.astype(float).copy()
    for column in range(centered.shape[1]):
        std = centered[:, column].std(ddof=0)
        if std == 0:
            raise ValueError("Continuous regressors contain a constant column; z-scoring would divide by zero.")
        centered[:, column] = (centered[:, column] - centered[:, column].mean()) / std
    return centered


def ensure_3d(array: Any) -> np.ndarray:
    array = np.asarray(array)
    if array.ndim == 2:
        return array[:, :, np.newaxis]
    if array.ndim == 3:
        return array
    raise ValueError("Expected a 2D or 3D array.")


def ensure_4d(array: Any) -> np.ndarray:
    array = np.asarray(array)
    if array.ndim == 3:
        return array[:, :, :, np.newaxis]
    if array.ndim == 4:
        return array
    raise ValueError("Expected a 3D or 4D array.")


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


def _is_empty_value(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, (list, tuple, dict, str, bytes)):
        return len(value) == 0
    if isinstance(value, np.ndarray):
        return value.size == 0
    return False


def _read_hdf5_item(node: Any) -> Any:
    if h5py is None:  # pragma: no cover
        raise ImportError("h5py is required to read LIMO HDF5 files.")

    if isinstance(node, h5py.Dataset):
        data = node[()]
        return _decode_hdf5_dataset(data)

    if node.attrs.get("__none__"):
        return None

    kind = _decode_hdf5_dataset(node.attrs.get("__kind__"))
    items = {key: _read_hdf5_item(node[key]) for key in node.keys()}
    if kind == "list":
        return [items[key] for key in sorted(items)]
    return items


def _decode_hdf5_dataset(data: Any) -> Any:
    if isinstance(data, bytes):
        return data.decode("utf-8", errors="replace")
    if isinstance(data, np.ndarray):
        if data.shape == ():
            return _decode_hdf5_dataset(data.item())
        if data.dtype.kind in {"S", "O", "U"}:
            return np.asarray(data).astype(str)
        return np.asarray(data)
    if isinstance(data, np.generic):
        return data.item()
    return data


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Update LIMO.h5 and export limo_results.h5")
    parser.add_argument("limo_file", type=Path, help="Path to the LIMO.h5 file")
    parser.add_argument(
        "--results-output",
        type=Path,
        default=None,
        help="Optional explicit output path for limo_results.h5",
    )
    return parser


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()
    limo_path, results_path = limo_design(args.limo_file, results_path=args.results_output)
    print(json.dumps({"LIMO": str(limo_path), "results": str(results_path)}, indent=2))


if __name__ == "__main__":
    main()
