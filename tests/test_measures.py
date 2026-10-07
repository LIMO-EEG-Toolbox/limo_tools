import numpy as np
import pytest
from scipy.io import savemat

from limo.eeglab_import import export_limo_h5
from limo.limo_design import limo_design, read_hdf5_structure
from limo.limo_glm import limo_glm
from limo.limo_contrast import limo_contrast


@pytest.mark.parametrize("method", ["OLS", "WLS", "IRLS"])
@pytest.mark.parametrize("analysis", ["Frequency", "Time-Frequency"])
def test_precomputed_measure_file_workflow(tmp_path, method, analysis):
    rng = np.random.default_rng(23)
    n = 60
    cat = np.tile([1, 2], 30)
    times = np.arange(5) * 10.
    freqs = np.array([5., 10., 15.])
    if analysis == "Frequency":
        signal = rng.uniform(.5, 3., size=(2, 3, n))
        key = "datspec"
        etc = {"freqspec": freqs}
        expected = signal
    else:
        signal = rng.normal(size=(2, 3, 5, n)) + 1j * rng.normal(size=(2, 3, 5, n))
        key = "dattimef"
        etc = {"freqersp": freqs, "timeersp": times}
        expected = np.abs(signal)**2
    measure = tmp_path / f"sub-01.{key}"
    savemat(measure, {"chan1": signal[0], "chan2": signal[1]})
    etc["datafiles"] = {key: measure.name}
    set_file = tmp_path / "sub-01.set"
    savemat(set_file, {"EEG": {
        "data": rng.normal(size=(2, 5, n)), "nbchan": 2, "pnts": 5,
        "trials": n, "srate": 100., "times": times, "etc": etc,
        "chanlocs": np.array([{"labels": "C1"}, {"labels": "C2"}], dtype=object),
    }})
    path = export_limo_h5(set_file, cat, None,
                          {"name": str(tmp_path / method), "analysis": analysis, "method": method})
    path, results = limo_design(path)
    order = np.argsort(cat, kind="stable")
    np.testing.assert_allclose(read_hdf5_structure(results)["Yr"], expected[..., order])
    limo_glm(path)
    fitted = read_hdf5_structure(results)
    assert fitted["Beta"].shape == (*expected.shape[:-1], 3)
    np.testing.assert_allclose(fitted["Res"], fitted["Yr"] - fitted["Yhat"])
    assert np.all(np.isfinite(fitted["Beta"]))
    _, _, _, name = limo_contrast(path, {"C": [1., -1., 0.], "F": 0})
    output = read_hdf5_structure(results)[name]
    assert output.shape == (*expected.shape[:-1], 5)
    assert np.all(np.isfinite(output))
    np.testing.assert_allclose(output[..., 0], fitted["Beta"][..., 0] - fitted["Beta"][..., 1])
