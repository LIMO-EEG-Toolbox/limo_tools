"""Utilities for reading EEGLAB .set files in Python.

Use ``read_set_file(set_file)`` when you want the EEGLAB dataset as a Python
dictionary built from the documented MATLAB ``EEG`` structure. This path handles
both embedded data and sidecar ``.fdt`` float32 files.

The returned variable is a dictionary named by the caller, for example
``eeg = read_set_file(...)``. It contains the EEGLAB fields such as ``srate``,
``times``, ``chanlocs``, ``event``, ``etc``, and ``data`` when ``load_data`` is
left at its default value of ``True``. The helper also adds ``set_path`` and
fills in ``filename`` and ``filepath`` when those fields are missing.

Use ``read_set_file_with_mne(set_file)`` when you want an MNE ``Raw`` object
and the dataset is compatible with ``mne.io.read_raw_eeglab``.

Command-line usage:
    python read_setfile.py path/to/file.set
    python read_setfile.py path/to/file.set --no-data
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat

try:
    from scipy.io.matlab import mat_struct
except ImportError:  # pragma: no cover
    mat_struct = None  # type: ignore[assignment]

try:
    import h5py
except ImportError:  # pragma: no cover
    h5py = None  # type: ignore[assignment]


def read_set_file(
    set_file: str | Path,
    *,
    load_data: bool = True,
) -> dict[str, Any]:
    """Read an EEGLAB .set file using its documented open MATLAB layout.

    EEGLAB documents .set files as MATLAB files that either contain an ``EEG``
    struct or the fields of that struct directly. Signal samples may be embedded
    in the .set file or stored in a sidecar float32 file such as ``.fdt``.
    """

    set_path = Path(set_file).expanduser().resolve()
    eeg = load_matlab_file(set_path)
    eeg = _extract_eeg_root(eeg)

    # loadmat squeezes singleton axes. Restore embedded data from the EEGLAB
    # dimension metadata, including for metadata-only reads used by design.
    if "data" in eeg and not isinstance(eeg["data"], str):
        embedded = np.asarray(eeg["data"])
        if embedded.size:
            nbchan = int(np.asarray(eeg["nbchan"]).item())
            pnts = int(np.asarray(eeg["pnts"]).item())
            trials = int(np.asarray(eeg.get("trials", 1)).item())
            shape = (nbchan, pnts, trials) if trials > 1 else (nbchan, pnts)
            if embedded.size != int(np.prod(shape)):
                raise ValueError("Embedded EEG data size disagrees with nbchan, pnts and trials.")
            eeg["data"] = embedded.reshape(shape)

    if load_data:
        eeg["data"] = _resolve_data_array(eeg, set_path)

    eeg.setdefault("filename", set_path.name)
    eeg.setdefault("filepath", str(set_path.parent))
    eeg["set_path"] = str(set_path)
    return eeg


def read_set_file_with_mne(
    set_file: str | Path,
    *,
    preload: bool = False,
    uint16_codec: str | None = None,
    montage_units: str = "auto",
) -> Any:
    """Read a continuous EEGLAB .set file through MNE if it is installed.

    This mirrors the documented MNE entry point ``mne.io.read_raw_eeglab``.
    MNE expects any referenced sidecar ``.fdt`` file to live next to the .set.
    """

    try:
        import mne
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "MNE is required for read_set_file_with_mne(). Install mne first."
        ) from exc

    return mne.io.read_raw_eeglab(
        str(Path(set_file).expanduser().resolve()),
        preload=preload,
        uint16_codec=uint16_codec,
        montage_units=montage_units,
    )


def load_matlab_file(file_path: str | Path) -> dict[str, Any]:
    """Load a MATLAB file saved in classic MAT or v7.3 HDF5 form."""

    path = Path(file_path).expanduser().resolve()
    try:
        payload = loadmat(path, squeeze_me=True, struct_as_record=False)
        return {
            key: _coerce_mat_value(value)
            for key, value in payload.items()
            if not key.startswith("__")
        }
    except NotImplementedError as exc:
        if h5py is None:
            raise ImportError(
                "h5py is required to read MATLAB v7.3 files. Install h5py first."
            ) from exc
        return _load_hdf5_matlab_file(path)


def _extract_eeg_root(payload: dict[str, Any]) -> dict[str, Any]:
    if "EEG" in payload:
        eeg = payload["EEG"]
        if isinstance(eeg, dict):
            return eeg
        raise TypeError("The MATLAB file contains an EEG variable with an unsupported type.")

    required_fields = {"srate", "nbchan", "pnts"}
    if required_fields.issubset(payload):
        return payload

    raise KeyError("Could not find an EEGLAB EEG structure in the .set file.")


def _coerce_mat_value(value: Any) -> Any:
    if mat_struct is not None and isinstance(value, mat_struct):
        return {
            field_name: _coerce_mat_value(getattr(value, field_name))
            for field_name in value._fieldnames
        }

    if isinstance(value, np.ndarray):
        if value.shape == ():
            return _coerce_mat_value(value.item())

        if value.dtype.names:
            items = []
            for entry in value.flat:
                items.append(
                    {
                        field_name: _coerce_mat_value(entry[field_name])
                        for field_name in value.dtype.names
                    }
                )
            return items[0] if len(items) == 1 else items

        if value.dtype == object:
            items = [_coerce_mat_value(item) for item in value.flat]
            return items[0] if len(items) == 1 else items

        if value.dtype.kind in {"U", "S"}:
            if value.ndim == 1:
                return "".join(str(item) for item in value.tolist())
            return value.astype(str).tolist()

        return np.asarray(value)

    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")

    if isinstance(value, np.generic):
        return value.item()

    return value


def _load_hdf5_matlab_file(path: Path) -> dict[str, Any]:
    if h5py is None:  # pragma: no cover
        raise ImportError("h5py is required to read MATLAB v7.3 files.")

    with h5py.File(path, "r") as handle:
        return {
            key: _read_hdf5_node(handle, handle[key])
            for key in handle.keys()
            if not key.startswith("#")
        }


def _read_hdf5_node(handle: Any, node: Any) -> Any:
    if h5py is None:  # pragma: no cover
        raise ImportError("h5py is required to read MATLAB v7.3 files.")

    if isinstance(node, h5py.Dataset):
        refs_dtype = h5py.check_dtype(ref=node.dtype)
        if refs_dtype is not None:
            data = node[()]
            refs = [_read_reference(handle, ref) for ref in np.ravel(data)]
            return refs[0] if len(refs) == 1 else refs

        matlab_class = _decode_attr(node.attrs.get("MATLAB_class"))
        data = node[()]

        if matlab_class == "char":
            return _decode_matlab_chars(data)

        if isinstance(data, np.ndarray) and data.shape == ():
            data = data.item()

        if isinstance(data, np.ndarray):
            return np.asarray(data)

        if isinstance(data, bytes):
            return data.decode("utf-8", errors="replace")

        return data

    matlab_class = _decode_attr(node.attrs.get("MATLAB_class"))
    fields = {name: _read_hdf5_node(handle, child) for name, child in node.items()}

    if matlab_class == "struct":
        return _collapse_struct_fields(fields)

    return fields


def _read_reference(handle: Any, ref: Any) -> Any:
    if not ref:
        return None
    return _read_hdf5_node(handle, handle[ref])


def _collapse_struct_fields(fields: dict[str, Any]) -> Any:
    list_lengths = [len(value) for value in fields.values() if isinstance(value, list)]
    if not list_lengths:
        return fields

    target_length = max(list_lengths)
    if not all(
        (not isinstance(value, list) and target_length == 1)
        or (isinstance(value, list) and len(value) == target_length)
        for value in fields.values()
    ):
        return fields

    entries = []
    for index in range(target_length):
        entry = {}
        for field_name, field_value in fields.items():
            entry[field_name] = field_value[index] if isinstance(field_value, list) else field_value
        entries.append(entry)
    return entries[0] if len(entries) == 1 else entries


def _decode_attr(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.ndarray):
        if value.shape == ():
            return _decode_attr(value.item())
        return "".join(_decode_attr(item) or "" for item in value.flat)
    return str(value)


def _decode_matlab_chars(data: Any) -> str:
    array = np.asarray(data)
    if array.size == 0:
        return ""
    if array.dtype.kind in {"U", "S"}:
        return "".join(array.astype(str).ravel(order="F")).rstrip("\x00")
    if array.dtype.kind in {"u", "i"}:
        return "".join(chr(int(value)) for value in array.ravel(order="F") if int(value) != 0)
    return str(array)


def _resolve_data_array(eeg: dict[str, Any], set_path: Path) -> np.ndarray | None:
    data_field = eeg.get("data")

    if data_field is None:
        return None

    if isinstance(data_field, str):
        sidecar = (set_path.parent / data_field).resolve()
        if not sidecar.exists():
            raise FileNotFoundError(f"Referenced sidecar data file does not exist: {sidecar}")
        return _load_fdt_data(
            sidecar,
            nbchan=int(np.asarray(eeg["nbchan"]).item()),
            pnts=int(np.asarray(eeg["pnts"]).item()),
            trials=int(np.asarray(eeg.get("trials", 1)).item()),
        )

    if isinstance(data_field, list):
        return np.asarray(data_field)

    if isinstance(data_field, np.ndarray):
        return np.asarray(data_field)

    return np.asarray(data_field)


def _load_fdt_data(sidecar: Path, *, nbchan: int, pnts: int, trials: int) -> np.ndarray:
    expected_values = nbchan * pnts * trials
    data = np.fromfile(sidecar, dtype="<f4")

    if data.size != expected_values:
        raise ValueError(
            f"Unexpected data size in {sidecar}. Expected {expected_values} float32 values, "
            f"found {data.size}."
        )

    # EEGLAB floatwrite stores the channel-by-time-by-trial matrix in
    # MATLAB column-major order, with channels varying fastest.
    data = data.reshape((nbchan, pnts, trials), order="F")
    return data[:, :, 0] if trials == 1 else data


def _json_ready(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {key: _json_ready(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Inspect an EEGLAB .set file.")
    parser.add_argument("set_file", type=Path, help="Path to the .set file")
    parser.add_argument(
        "--no-data",
        action="store_true",
        help="Skip loading the embedded or sidecar signal array",
    )
    return parser


def main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()
    eeg = read_set_file(args.set_file, load_data=not args.no_data)
    print(json.dumps(_json_ready(eeg), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
