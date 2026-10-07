from pathlib import Path

import numpy as np
import pytest
from scipy.io import savemat
from scipy.stats import t as t_distribution

from limo.eeglab_import import export_limo_h5
from limo.limo_design import limo_design, read_hdf5_structure
from limo.limo_glm import limo_glm
from limo.limo_contrast import limo_contrast, limo_contrast_checking, _matlab_int16
from limo.read_setfile import read_set_file


def write_set(path, signal, sidecar=False):
    nchan, nframe, ntrial = signal.shape
    data = signal
    if sidecar:
        data = path.with_suffix(".fdt").name
        signal.astype("<f4").ravel(order="F").tofile(path.with_suffix(".fdt"))
    savemat(path, {"EEG": {
        "data": data, "srate": 100., "times": np.arange(nframe) * 10.,
        "xmin": 0., "xmax": (nframe - 1) / 100.,
        "trials": ntrial, "nbchan": nchan, "pnts": nframe,
        "etc": {"timeerp": np.arange(nframe) * 10.},
        "chanlocs": np.array([{"labels": f"C{i}"} for i in range(nchan)], dtype=object),
    }})


@pytest.mark.parametrize("method", ["OLS", "WLS", "IRLS"])
def test_set_to_design_glm_and_t_contrast(tmp_path, method):
    rng = np.random.default_rng(17)
    categorical = np.tile([1, 2], 30)
    y = rng.normal(size=(2, 12, 60))
    y[:, :, categorical == 2] += 2
    set_path = tmp_path / "sub-01.set"
    write_set(set_path, y)
    path = export_limo_h5(set_path, categorical, None,
                          {"name": str(tmp_path / method), "analysis": "Time", "method": method})
    path, results = limo_design(path)
    metadata = read_hdf5_structure(path)["LIMO"]
    stored = read_hdf5_structure(results)
    order = np.argsort(categorical, kind="stable")
    np.testing.assert_array_equal(stored["Yr"], y[..., order])
    np.testing.assert_array_equal(metadata["design"]["X"][:, :2],
                                  np.column_stack((categorical[order] == 1, categorical[order] == 2)))
    limo_glm(path)
    before_contrast = read_hdf5_structure(results)
    _, _, _, result_name = limo_contrast(path, {"C": [1., -1., 0.], "F": 0})
    output = read_hdf5_structure(results)[result_name]
    assert output.shape == (2, 12, 5)
    np.testing.assert_allclose(output[..., 0], before_contrast["Beta"][..., 0] - before_contrast["Beta"][..., 1])
    assert np.all(np.isfinite(output))
    if method == "OLS":
        x = metadata["design"]["X"]
        df = 60 - np.linalg.matrix_rank(x)
        c = np.array([1., -1., 0.])
        variance = np.sum(before_contrast["Res"]**2, axis=-1) / df
        expected_se = np.sqrt(variance * (c @ np.linalg.pinv(x.T @ x) @ c))
        np.testing.assert_allclose(output[..., 1], expected_se)
        np.testing.assert_allclose(output[..., 4], 2 * t_distribution.sf(np.abs(output[..., 3]), df))
    # The fitted file can be reopened without recomputing the completed model.
    limo_glm(path)
    np.testing.assert_array_equal(read_hdf5_structure(results)["Beta"], before_contrast["Beta"])


@pytest.mark.parametrize("shape", [(2, 4, 5), (1, 4, 5), (2, 1, 5)])
@pytest.mark.parametrize("sidecar", [False, True])
def test_set_reader_preserves_channel_time_trial_axes(tmp_path, shape, sidecar):
    y = np.arange(np.prod(shape), dtype=float).reshape(shape)
    path = tmp_path / "signal.set"
    write_set(path, y, sidecar)
    loaded = read_set_file(path)["data"]
    assert loaded.shape == shape
    np.testing.assert_array_equal(loaded, y)


def test_contrast_conversion_uses_matlab_rounding_and_saturation():
    np.testing.assert_array_equal(_matlab_int16([-.5, .5, 1.999999999, -1.999999999, 99999, -99999]),
                                  [-1, 1, 2, -2, 32767, -32768])
    x = np.column_stack((np.tile([1., 0.], 30), np.tile([0., 1.], 30), np.ones(60)))
    assert limo_contrast_checking([1., -1., 0.], x) == 1
