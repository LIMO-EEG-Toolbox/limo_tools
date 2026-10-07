"""Utilities for exporting LIMO metadata from an EEGLAB .set file.

Use ``export_limo_h5(set_file, cat=None, cont=None, defaults=...)`` and
save the result as ``LIMO.h5``.

The ``defaults`` argument must provide the LIMO options needed to build the
output structure, either as a Python mapping or as a path to a JSON file. The
optional ``cat`` and ``cont`` arguments can be numeric values, arrays, or paths
to ``.txt`` or ``.mat`` regressor files.

The returned variable is a ``Path`` pointing to the exported ``LIMO.h5`` file,
for example ``output = export_limo_h5(...)``. The file contains a top-level
``LIMO`` group with the nested ``data`` and ``design`` fields assembled from
the EEGLAB dataset and the provided defaults.

Command-line usage:
    python eeglab_import.py path/to/file.set --defaults path/to/defaults.json
    python eeglab_import.py path/to/file.set --defaults path/to/defaults.json --cat cat.txt --cont cont.mat
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

try:
    import h5py
except ImportError:  # pragma: no cover
    h5py = None  # type: ignore[assignment]

from .read_setfile import load_matlab_file, read_set_file


DEFAULT_LIMO_OPTIONS = {
    "type": "Channels",
    "method": "WLS",
    "type_of_analysis": "Mass-univariate",
    "fullfactorial": 0,
    "zscore": 0,
    "bootstrap": 0,
    "tfce": 0,
}


def export_limo_h5(
    set_file: str | Path,
    cat: Any = None,
    cont: Any = None,
    defaults: Mapping[str, Any] | str | Path | None = None,
    output_path: str | Path | None = None,
) -> Path:

    options = _coerce_defaults(defaults)
    set_path = Path(set_file).expanduser().resolve()
    eeg = read_set_file(set_path, load_data=False)
    limo = build_limo_structure(eeg=eeg, set_path=set_path, cat=cat, cont=cont, defaults=options)

    output = Path(output_path).expanduser().resolve() if output_path else Path(options["name"]).expanduser().resolve() / "LIMO.h5"
    output.parent.mkdir(parents=True, exist_ok=True)
    write_hdf5_structure(output, {"LIMO": limo})
    return output


def import_eeglab_to_limo_h5(
    set_file: str | Path,
    cat: Any = None,
    cont: Any = None,
    defaults: Mapping[str, Any] | str | Path | None = None,
    output_path: str | Path | None = None,
) -> Path:
    """Backward-compatible alias for export_limo_h5()."""

    return export_limo_h5(
        set_file=set_file,
        cat=cat,
        cont=cont,
        defaults=defaults,
        output_path=output_path,
    )


def build_limo_structure(
    *,
    eeg: Mapping[str, Any],
    set_path: Path,
    cat: Any,
    cont: Any,
    defaults: Mapping[str, Any],
) -> dict[str, Any]:
    options = dict(DEFAULT_LIMO_OPTIONS)
    options.update(defaults)

    if not options.get("analysis"):
        raise KeyError("defaults['analysis'] is required.")
    if not options.get("name"):
        raise KeyError("defaults['name'] is required.")

    limo = {
        "dir": str(Path(options["name"]).expanduser().resolve()),
        "Analysis": options["analysis"],
        "Type": options["type"],
        "Level": 1,
        "data": {
            "data": set_path.name,
            "data_dir": str(set_path.parent),
            "sampling_rate": _scalar(eeg.get("srate")),
            "Cat": _load_regressor(cat, mode="cat"),
            "Cont": _load_regressor(cont, mode="cont"),
        },
        "design": {
            "zscore": options["zscore"],
            "method": options["method"],
            "type_of_analysis": options["type_of_analysis"],
            "fullfactorial": options["fullfactorial"],
            "bootstrap": options["bootstrap"],
            "tfce": options["tfce"],
            "status": "to do",
        },
    }

    if "labels" in options:
        limo["design"]["labels"] = options["labels"]

    if "icaclustering" in options:
        limo["data"]["cluster"] = options["icaclustering"]

    if "chanlocs" in options:
        limo["data"]["chanlocs"] = options["chanlocs"]
    else:
        limo["data"]["chanlocs"] = eeg.get("chanlocs", [])

    if "neighbouring_matrix" in options:
        limo["data"]["neighbouring_matrix"] = options["neighbouring_matrix"]

    if "studyinfo" in options:
        limo["data"]["studyinfo"] = options["studyinfo"]

    analysis = str(options["analysis"])
    etc = eeg.get("etc", {}) if isinstance(eeg.get("etc", {}), Mapping) else {}

    if analysis == "Time":
        timevect = _pick_vector(etc.get("timeerp"), eeg.get("times"), field_name="time vector")
        _populate_time_window(limo["data"], timevect, options.get("start"), options.get("end"))
        limo["data"]["timevect"] = limo["data"].pop("_selected_vector")

    elif analysis == "Frequency":
        freqvect = _pick_vector(etc.get("freqspec"), eeg.get("freqs"), field_name="frequency vector")
        _populate_frequency_window(limo["data"], freqvect, options.get("lowf"), options.get("highf"))
        limo["data"]["freqlist"] = limo["data"].pop("_selected_vector")

    elif analysis == "Time-Frequency":
        timevect = _pick_vector(etc.get("timeersp"), eeg.get("times"), field_name="time vector")
        freqvect = _pick_vector(etc.get("freqersp"), eeg.get("freqs"), field_name="frequency vector")
        _populate_time_window(limo["data"], timevect, options.get("start"), options.get("end"))
        limo["data"]["tf_times"] = limo["data"].pop("_selected_vector")
        _populate_tf_frequency_window(limo["data"], freqvect, options.get("lowf"), options.get("highf"))
        limo["data"]["tf_freqs"] = limo["data"].pop("_selected_vector")

    else:
        raise ValueError(
            "Unsupported analysis type. Expected 'Time', 'Frequency', or 'Time-Frequency'."
        )

    return limo


def write_hdf5_structure(output_path: str | Path, payload: Mapping[str, Any]) -> None:
    if h5py is None:  # pragma: no cover
        raise ImportError("h5py is required to write LIMO HDF5 files. Install h5py first.")

    path = Path(output_path).expanduser().resolve()
    with h5py.File(path, "w") as handle:
        handle.attrs["creator"] = "eeglab_import.py"
        handle.attrs["format"] = "LIMO-HDF5"
        for key, value in payload.items():
            _write_hdf5_item(handle, key, value)


def _coerce_defaults(defaults: Mapping[str, Any] | str | Path | None) -> dict[str, Any]:
    if defaults is None:
        raise ValueError("defaults must be provided as a mapping or a JSON file path.")
    if isinstance(defaults, Mapping):
        return dict(defaults)

    path = Path(defaults).expanduser().resolve()
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _load_regressor(value: Any, *, mode: str) -> np.ndarray:
    if value is None:
        return np.asarray([])

    if isinstance(value, np.ndarray):
        return np.asarray(value)

    if isinstance(value, (list, tuple)):
        return np.asarray(value)

    if isinstance(value, (int, float, np.integer, np.floating, bool)):
        return np.asarray(value)

    if isinstance(value, (str, Path)):
        path = Path(value).expanduser().resolve()
        suffix = path.suffix.lower()
        if suffix == ".txt":
            return np.loadtxt(path, ndmin=1)
        if suffix == ".mat":
            payload = load_matlab_file(path)
            if not payload:
                return np.asarray([])
            if mode == "cont" and path.stem in payload:
                return np.asarray(payload[path.stem])
            first_key = next(iter(payload))
            return np.asarray(payload[first_key])
        raise ValueError("Regressor files must be .txt or .mat files.")

    return np.asarray(value)


def _pick_vector(*candidates: Any, field_name: str) -> np.ndarray:
    for candidate in candidates:
        if candidate is None:
            continue
        array = np.asarray(candidate).squeeze()
        if array.size == 0:
            continue
        return np.asarray(array, dtype=float)
    raise KeyError(f"Could not find the required {field_name} in the EEGLAB dataset.")


def _populate_time_window(data: dict[str, Any], vector: np.ndarray, start: Any, end: Any) -> None:
    start_index = _nearest_index(vector, start, lower_bound=vector[0])
    end_index = _nearest_index(vector, end, upper_bound=vector[-1])

    data["start"] = float(vector[start_index])
    data["trim1"] = int(start_index + 1)
    data["end"] = float(vector[end_index])
    data["trim2"] = int(end_index + 1)
    data["_selected_vector"] = vector[start_index : end_index + 1]


def _populate_frequency_window(data: dict[str, Any], vector: np.ndarray, lowf: Any, highf: Any) -> None:
    start_index = _nearest_index(vector, lowf, lower_bound=vector[0])
    end_index = _nearest_index(vector, highf, upper_bound=vector[-1])

    data["start"] = float(vector[start_index])
    data["trim1"] = int(start_index + 1)
    data["end"] = float(vector[end_index])
    data["trim2"] = int(end_index + 1)
    data["_selected_vector"] = vector[start_index : end_index + 1]


def _populate_tf_frequency_window(data: dict[str, Any], vector: np.ndarray, lowf: Any, highf: Any) -> None:
    start_index = _nearest_index(vector, lowf, lower_bound=vector[0])
    end_index = _nearest_index(vector, highf, upper_bound=vector[-1])

    data["lowf"] = float(vector[start_index])
    data["trim_lowf"] = int(start_index + 1)
    data["highf"] = float(vector[end_index])
    data["trim_highf"] = int(end_index + 1)
    data["_selected_vector"] = vector[start_index : end_index + 1]


def _nearest_index(
    vector: np.ndarray,
    requested: Any,
    *,
    lower_bound: float | None = None,
    upper_bound: float | None = None,
) -> int:
    if requested is None or requested == []:
        if lower_bound is not None:
            return 0
        if upper_bound is not None:
            return int(vector.size - 1)

    value = float(np.asarray(requested).item())

    if lower_bound is not None and value < lower_bound:
        return 0
    if upper_bound is not None and value > upper_bound:
        return int(vector.size - 1)

    return int(np.argmin(np.abs(vector - value)))


def _scalar(value: Any) -> Any:
    if isinstance(value, np.ndarray) and value.shape == ():
        return value.item()
    if isinstance(value, np.generic):
        return value.item()
    return value


def _write_hdf5_item(parent: Any, name: str, value: Any) -> None:
    key = str(name)

    if value is None:
        group = parent.create_group(key)
        group.attrs["__none__"] = True
        return

    if isinstance(value, Path):
        _write_string_dataset(parent, key, str(value))
        return

    if isinstance(value, Mapping):
        group = parent.create_group(key)
        group.attrs["__kind__"] = "dict"
        for child_key, child_value in value.items():
            _write_hdf5_item(group, child_key, child_value)
        return

    if isinstance(value, np.ndarray):
        _write_array_dataset(parent, key, value)
        return

    if isinstance(value, (list, tuple)):
        if all(isinstance(item, Mapping) for item in value):
            group = parent.create_group(key)
            group.attrs["__kind__"] = "list"
            for index, item in enumerate(value):
                _write_hdf5_item(group, f"{index:04d}", item)
            return

        if all(not isinstance(item, (Mapping, list, tuple)) for item in value):
            _write_array_dataset(parent, key, np.asarray(value))
            return

        group = parent.create_group(key)
        group.attrs["__kind__"] = "list"
        for index, item in enumerate(value):
            _write_hdf5_item(group, f"{index:04d}", item)
        return

    if isinstance(value, str):
        _write_string_dataset(parent, key, value)
        return

    if isinstance(value, (bool, int, float, np.integer, np.floating)):
        parent.create_dataset(key, data=value)
        return

    _write_string_dataset(parent, key, json.dumps(_json_ready(value)))


def _write_array_dataset(parent: Any, name: str, array: np.ndarray) -> None:
    if h5py is None:  # pragma: no cover
        raise ImportError("h5py is required to write LIMO HDF5 files.")

    array = np.asarray(array)
    if array.dtype.kind in {"U", "S"}:
        dtype = h5py.string_dtype(encoding="utf-8")
        parent.create_dataset(name, data=array.astype(dtype), dtype=dtype)
        return

    if array.dtype == object:
        _write_string_dataset(parent, name, json.dumps(_json_ready(array.tolist())))
        return

    parent.create_dataset(name, data=array)


def _write_string_dataset(parent: Any, name: str, value: str) -> None:
    if h5py is None:  # pragma: no cover
        raise ImportError("h5py is required to write LIMO HDF5 files.")
    dtype = h5py.string_dtype(encoding="utf-8")
    parent.create_dataset(name, data=value, dtype=dtype)


def _json_ready(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Mapping):
        return {key: _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Export a LIMO.h5 file from an EEGLAB .set")
    parser.add_argument("set_file", type=Path, help="Path to the EEGLAB .set file")
    parser.add_argument(
        "--defaults",
        required=True,
        type=Path,
        help="JSON file containing the defaults structure expected by LIMO",
    )
    parser.add_argument("--cat", type=Path, default=None, help="Optional categorical regressor file")
    parser.add_argument("--cont", type=Path, default=None, help="Optional continuous regressor file")
    parser.add_argument("--output", type=Path, default=None, help="Optional explicit LIMO.h5 path")
    return parser


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()
    output = export_limo_h5(
        set_file=args.set_file,
        cat=args.cat,
        cont=args.cont,
        defaults=args.defaults,
        output_path=args.output,
    )
    print(output)


if __name__ == "__main__":
    main()