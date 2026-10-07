"""Cluster-based multiple-comparison correction for LIMO EEG results.

This module is the Python counterpart of the MATLAB clustering path centered
on ``limo_clustering.m``. It embeds the dependent routines used directly by
that entry point:

- ``limo_cluster_test.m`` for spatial-temporal cluster-mass thresholding
- ``limo_ecluster_make.m`` for 1D temporal bootstrap thresholds
- ``limo_ecluster_test.m`` for 1D temporal observed-data testing
- ``limo_findcluster.m`` for channel-aware connected-component labeling

It also adds a Python neighbour-matrix path inspired by MNE-Python. When
channel locations are available and a neighbouring matrix is missing, the
module can:

- call ``mne.channels.find_ch_adjacency`` when MNE is installed, or
- fall back to a Delaunay-triangulation adjacency built from channel
  positions, following the same broad idea described in the MNE adjacency
  documentation for EEG sensors.

Credits and provenance:
    The statistical procedures implemented here are ports of the LIMO EEG
    MATLAB functions listed above, from the LIMO Team.
    The optional adjacency-building behavior is informed by the public
    MNE-Python APIs and documentation for:
        ``mne.channels.find_ch_adjacency``
        ``mne.stats.combine_adjacency``
        ``mne.stats.spatio_temporal_cluster_test``
    This file does not copy MNE-Python source code verbatim. Instead, it uses
    the public MNE API directly when available and otherwise provides a local
    fallback implementation.

Primary use:
    Call ``limo_clustering(M, P, bootM, bootP, LIMO, MCC, p)`` with the same
    data layout as MATLAB:

    - ``M``: observed F or squared-t statistic map
    - ``P``: observed p-value map
    - ``bootM``: bootstrap statistic maps under H0
    - ``bootP``: bootstrap p-value maps under H0
    - ``LIMO``: structure or mapping containing at least ``data.chanlocs`` and
      optionally ``data.neighbouring_matrix``
    - ``MCC``: ``2`` for spatial-temporal clustering, ``3`` for temporal
      clustering
    - ``p``: uncorrected cluster-forming threshold

Inputs:
    ``M``, ``P``:
        2D matrices for electrode x frame analyses. For single-channel data,
        the first dimension must still be present with length 1.
    ``bootM``, ``bootP``:
        3D matrices with bootstrap samples in the last dimension.
    ``LIMO``:
        Mapping or object-like structure with a ``data`` field. The function
        uses ``data.neighbouring_matrix`` when present. If it is missing and
        channel locations are available, a neighbouring matrix can be built
        automatically.
    ``MCC``:
        ``2`` for spatial-temporal clustering and ``3`` for temporal-only
        clustering. As in MATLAB, one-channel data force temporal clustering.
    ``p``:
        Cluster-forming threshold.
    ``fig``:
        When ``1``, plot the bootstrap cluster-mass distribution if no
        significant cluster survives correction. Requires ``matplotlib``.
    ``build_neighbours``:
        If ``True``, try to build a neighbouring matrix when the LIMO
        structure does not already contain one.
    ``prefer_mne``:
        If ``True``, try MNE first for channel adjacency before falling back
        to the local Delaunay path.
    ``neighbourdist``:
        Optional distance threshold for the pure distance-based fallback.
        When omitted, the Delaunay fallback is used.

Returned value:
    A tuple ``(mask, cluster_pval, max_th)`` matching MATLAB.

Outputs:
    ``mask``:
        Labelled cluster mask with the same first two dimensions as ``M``.
    ``cluster_pval``:
        Corrected p-value matrix, matching the shape of ``mask``.
    ``max_th``:
        Cluster-mass threshold controlling the family-wise error rate.

References:
    Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of
    EEG- and MEG-data. Journal of Neuroscience Methods, 164(1), 177-190.
    Rousselet, G. A., Pernet, C. R., and LIMO Team clustering functions in the
    LIMO EEG toolbox.
    MNE-Python documentation for sensor adjacency and spatio-temporal cluster
    tests: https://mne.tools/stable/

Command-line usage:
    This file primarily provides a library function. No standalone CLI is
    exposed because the clustering inputs are already in-memory statistic
    arrays in the MATLAB API.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
from scipy import ndimage
from scipy.spatial import Delaunay, distance


def limo_clustering(
    M: np.ndarray,
    P: np.ndarray,
    bootM: np.ndarray,
    bootP: np.ndarray,
    LIMO: Mapping[str, Any] | Any,
    MCC: int,
    p: float,
    fig: int | None = None,
    *,
    build_neighbours: bool = True,
    prefer_mne: bool = True,
    neighbourdist: float | None = None,
) -> tuple[np.ndarray, np.ndarray | list[Any], float | np.ndarray | None]:
    """Port of MATLAB ``limo_clustering.m``."""

    observed = _ensure_2d(M)
    observed_p = _ensure_2d(P)
    boot_stat = _ensure_3d(bootM)
    boot_p = _ensure_3d(bootP)

    if observed.shape != observed_p.shape:
        raise ValueError("M and P must have the same shape.")
    if boot_stat.shape != boot_p.shape:
        raise ValueError("bootM and bootP must have the same shape.")
    if observed.shape != boot_stat.shape[:2]:
        raise ValueError("Observed and bootstrap arrays must agree on the first two dimensions.")

    if observed.shape[0] == 1:
        MCC = 3

    nboot = boot_stat.shape[2]
    cluster_pval: np.ndarray | list[Any] = np.asarray([])
    mask = np.asarray([])
    max_th: float | np.ndarray | None = None
    boot_maxclustersum: np.ndarray | None = None
    maxval = 0.0

    if MCC == 2 and boot_stat.shape[0] > 1:
        minnbchan = 2
        channeighbstructmat = _resolve_neighbouring_matrix(
            LIMO,
            n_channels=observed.shape[0],
            build_neighbours=build_neighbours,
            prefer_mne=prefer_mne,
            neighbourdist=neighbourdist,
        )
        boot_maxclustersum = np.zeros(nboot, dtype=float)
        for boot in range(nboot):
            labels, n_clusters = _limo_findcluster(boot_p[:, :, boot] <= p, channeighbstructmat, 2)
            if n_clusters != 0:
                tmp = np.zeros(n_clusters, dtype=float)
                current_boot = boot_stat[:, :, boot]
                for cluster in range(1, n_clusters + 1):
                    tmp[cluster - 1] = np.sum(current_boot[labels == cluster])
                boot_maxclustersum[boot] = np.max(tmp)
            else:
                boot_maxclustersum[boot] = 0.0

        mask, cluster_pval, maxval, max_th = _limo_cluster_test(
            observed,
            observed_p,
            boot_maxclustersum,
            channeighbstructmat,
            minnbchan,
            p,
        )

    if (MCC == 2 and boot_stat.shape[0] == 1) or MCC == 3:
        th, boot_values = _limo_ecluster_make(np.squeeze(boot_stat), np.squeeze(boot_p), p)
        max_th = float(np.max(np.atleast_1d(th["elec"])))
        sigcluster, cluster_pval, maxval_channel = _limo_ecluster_test(
            np.squeeze(observed),
            np.squeeze(observed_p),
            th,
            p,
            boot_values,
        )
        mask = np.asarray(sigcluster["elec_mask"], dtype=float)
        boot_maxclustersum = np.max(np.atleast_2d(boot_values), axis=0)
        maxval = float(np.max(np.atleast_1d(maxval_channel)))

    if np.size(mask) and np.sum(mask) == 0 and fig == 1 and boot_maxclustersum is not None and max_th is not None:
        _plot_bootstrap_distribution(boot_maxclustersum, max_th, maxval)

    return np.asarray(mask), np.asarray(cluster_pval), max_th


def _limo_cluster_test(
    ori_f: np.ndarray,
    ori_p: np.ndarray,
    boot_maxclustersum: np.ndarray,
    channeighbstructmat: np.ndarray,
    minnbchan: int = 2,
    alphav: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    labels, n_clusters = _limo_findcluster(ori_p <= alphav, channeighbstructmat, minnbchan)
    nboot = boot_maxclustersum.shape[0]
    sort_clustermax = np.sort(np.asarray(boot_maxclustersum, dtype=float))
    n_nan = int(np.isnan(sort_clustermax).sum())
    if n_nan == nboot:
        raise ValueError("no cluster max value found - check boot_maxclustersum parameter")
    if n_nan:
        clean = sort_clustermax[~np.isnan(sort_clustermax)]
        sort_clustermax = np.concatenate((np.full(n_nan, np.nan), clean))
    threshold_index = int(np.clip(round((1 - alphav) * nboot) - 1, 0, nboot - 1))
    maxclustersum_th = float(sort_clustermax[threshold_index])

    mask = np.zeros_like(ori_f, dtype=float)
    cluster_label = 1
    cluster_masses: list[float] = []
    if n_clusters != 0:
        for cluster in range(n_clusters, 0, -1):
            mass = float(np.sum(ori_f[labels == cluster]))
            cluster_masses.append(mass)
            if mass >= maxclustersum_th:
                mask[labels == cluster] = cluster_label
                cluster_label += 1

    pval = np.ones_like(mask, dtype=float)
    mask2 = mask.astype(bool)
    if np.any(mask2):
        filtered_labels = labels * mask2
        cluster_ids = [int(value) for value in np.unique(filtered_labels) if value != 0]
        for cluster_id in cluster_ids:
            cluster_mass = float(np.sum(ori_f[filtered_labels == cluster_id]))
            p_cluster = 1 - np.sum(cluster_mass >= sort_clustermax) / nboot
            if p_cluster == 0:
                p_cluster = 1 / nboot
            tmp = np.ones_like(mask, dtype=float)
            tmp[filtered_labels == cluster_id] = p_cluster
            pval *= tmp
    pval[pval == 1] = np.nan
    if np.any(pval[~np.isnan(pval)] > alphav):
        raise ValueError(
            "some corrected p-values are above the set alpha value, which should not happen"
        )
    maxval = float(max(cluster_masses)) if cluster_masses else 0.0
    return mask, pval, maxval, maxclustersum_th


def _limo_ecluster_make(
    bootf: np.ndarray,
    bootp: np.ndarray,
    alphav: float = 0.05,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    bootf = np.asarray(bootf, dtype=float)
    bootp = np.asarray(bootp, dtype=float)
    if bootf.ndim == 3:
        n_electrodes = bootf.shape[0]
        b = bootf.shape[2]
        u = int(np.clip(round((1 - alphav) * b) - 1, 0, b - 1))
        boot_values = np.zeros((n_electrodes, b), dtype=float)
        for electrode in range(n_electrodes):
            for boot in range(b):
                labels, n_clusters = _label_temporal_clusters(np.squeeze(bootp[electrode, :, boot]) <= alphav)
                if n_clusters != 0:
                    tmp = np.zeros(n_clusters, dtype=float)
                    for cluster in range(1, n_clusters + 1):
                        tmp[cluster - 1] = np.sum(np.squeeze(bootf[electrode, :, boot])[labels == cluster])
                    boot_values[electrode, boot] = np.max(tmp)
                else:
                    boot_values[electrode, boot] = 0.0
        sort_sc = np.sort(boot_values, axis=1)
        th = {
            "elec": sort_sc[:, u],
            "max": np.sort(np.max(boot_values, axis=0))[u],
        }
        return th, boot_values

    if bootf.ndim == 2:
        b = bootf.shape[1]
        u = int(np.clip(round((1 - alphav) * b) - 1, 0, b - 1))
        boot_values = np.zeros(b, dtype=float)
        for boot in range(b):
            labels, n_clusters = _label_temporal_clusters(np.squeeze(bootp[:, boot]) <= alphav)
            if n_clusters != 0:
                tmp = np.zeros(n_clusters, dtype=float)
                for cluster in range(1, n_clusters + 1):
                    tmp[cluster - 1] = np.sum(np.squeeze(bootf[:, boot])[labels == cluster])
                boot_values[boot] = np.max(tmp)
            else:
                boot_values[boot] = 0.0
        sort_sc = np.sort(boot_values)
        th = {"elec": np.asarray(sort_sc[u], dtype=float)}
        return th, boot_values

    raise ValueError("bootf should be 2D or 3D")


def _limo_ecluster_test(
    orif: np.ndarray,
    orip: np.ndarray,
    th: Mapping[str, Any],
    alpha_value: float = 0.05,
    boot_maxclustersum: np.ndarray | None = None,
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    orif = np.asarray(orif, dtype=float)
    orip = np.asarray(orip, dtype=float)
    if orif.ndim == 1:
        orif = orif[np.newaxis, :]
        orip = orip[np.newaxis, :]

    n_electrodes, n_frames = orif.shape
    sigcluster: dict[str, np.ndarray] = {}
    pval = np.asarray([])
    maxval = np.asarray([])

    if "max" in th:
        sigcluster["max_mask"] = np.zeros((n_electrodes, n_frames), dtype=float)
        for electrode in range(n_electrodes):
            labels, n_clusters = _label_temporal_clusters(orip[electrode] <= alpha_value)
            current_max = np.zeros(n_clusters, dtype=float)
            cluster_label = 1
            for cluster in range(1, n_clusters + 1):
                current_max[cluster - 1] = np.sum(np.abs(orif[electrode, labels == cluster]))
                if current_max[cluster - 1] >= float(np.asarray(th["max"]).item()):
                    sigcluster["max_mask"][electrode, labels == cluster] = cluster_label
                    cluster_label += 1
            maxval = np.asarray([np.max(current_max) if current_max.size else 0.0])

    if "elec" in th:
        sigcluster["elec_mask"] = np.zeros((n_electrodes, n_frames), dtype=float)
        maxval_channel = np.zeros(n_electrodes, dtype=float)
        pval = np.full((n_electrodes, n_frames), np.nan, dtype=float)
        cluster_label = 1

        if boot_maxclustersum is None:
            boot_vector = np.asarray([], dtype=float)
        else:
            boot_vector = np.asarray(boot_maxclustersum, dtype=float)
            if boot_vector.ndim > 1:
                boot_vector = np.max(boot_vector, axis=0)

        thresholds = np.asarray(th["elec"], dtype=float)
        if thresholds.ndim == 0:
            thresholds = np.repeat(thresholds.item(), n_electrodes)

        for channel in range(n_electrodes):
            labels, n_clusters = _label_temporal_clusters(orip[channel] <= alpha_value)
            current_max = np.zeros(n_clusters, dtype=float)
            for cluster in range(1, n_clusters + 1):
                current_max[cluster - 1] = np.sum(np.abs(orif[channel, labels == cluster]))
                if current_max[cluster - 1] >= thresholds[channel]:
                    sigcluster["elec_mask"][channel, labels == cluster] = cluster_label
                    cluster_label += 1
                    if boot_vector.size:
                        p_cluster = 1 - np.sum(current_max[cluster - 1] >= boot_vector) / boot_vector.size
                        if p_cluster == 0:
                            p_cluster = 1 / boot_vector.size
                        pval[channel, labels == cluster] = p_cluster
            maxval_channel[channel] = np.max(current_max) if current_max.size else 0.0
        maxval = maxval_channel

    return sigcluster, pval, maxval


def _limo_findcluster(
    onoff: np.ndarray,
    spatdimneighbstructmat: np.ndarray,
    minnbchan: int = 2,
) -> tuple[np.ndarray, int]:
    onoff = np.asarray(onoff, dtype=bool)
    if onoff.ndim == 2:
        onoff = onoff[:, :, np.newaxis]
        squeeze_last = True
    elif onoff.ndim == 3:
        squeeze_last = False
    else:
        raise ValueError("onoff must be 2D or 3D")

    spatdimlength, nfreq, ntime = onoff.shape
    neighbour = np.asarray(spatdimneighbstructmat, dtype=bool)
    if neighbour.shape != (spatdimlength, spatdimlength):
        raise ValueError("invalid dimension of spatdimneighbstructmat")

    working = onoff.copy()
    if minnbchan > 0:
        selectmat = np.asarray(neighbour | neighbour.T, dtype=float)
        nremoved = 1
        while nremoved > 0:
            nsigneighb = (selectmat @ working.reshape(spatdimlength, nfreq * ntime).astype(float)).reshape(working.shape)
            remove = (working.astype(float) * nsigneighb) < minnbchan
            nremoved = int(np.count_nonzero(remove & working))
            working[remove] = False

    labelmat = np.zeros(working.shape, dtype=int)
    total = 0
    structure = ndimage.generate_binary_structure(2, 1)
    for spatdimlev in range(spatdimlength):
        labels, num = ndimage.label(working[spatdimlev], structure=structure)
        labels[labels != 0] += total
        labelmat[spatdimlev] = labels
        total += int(num)

    if total == 0:
        output = np.zeros_like(labelmat, dtype=int)
        return (output[:, :, 0], 0) if squeeze_last else (output, 0)

    flat = labelmat.reshape(spatdimlength, nfreq * ntime)
    replaceby = np.arange(1, total + 1, dtype=int)
    for spatdimlev in range(spatdimlength):
        neighbours = np.flatnonzero(neighbour[spatdimlev])
        for neighbour_index in neighbours:
            overlap = np.flatnonzero((flat[spatdimlev] != 0) & (flat[neighbour_index] != 0))
            for idx in overlap:
                a = flat[spatdimlev, idx]
                b = flat[neighbour_index, idx]
                rep_a = replaceby[a - 1]
                rep_b = replaceby[b - 1]
                if rep_a == rep_b:
                    continue
                if rep_a < rep_b:
                    replaceby[replaceby == rep_b] = rep_a
                else:
                    replaceby[replaceby == rep_a] = rep_b

    cluster = np.zeros(flat.shape, dtype=int)
    num_clusters = 0
    for label in np.unique(replaceby):
        num_clusters += 1
        label_ids = np.flatnonzero(replaceby == label) + 1
        cluster[np.isin(flat, label_ids)] = num_clusters
    cluster = cluster.reshape(spatdimlength, nfreq, ntime)
    return (cluster[:, :, 0], num_clusters) if squeeze_last else (cluster, num_clusters)


def _resolve_neighbouring_matrix(
    LIMO: Mapping[str, Any] | Any,
    *,
    n_channels: int,
    build_neighbours: bool,
    prefer_mne: bool,
    neighbourdist: float | None,
) -> np.ndarray:
    data = _extract_data_mapping(LIMO)
    for key in ("neighbouring_matrix", "channeighbstructmat"):
        if key in data:
            matrix = np.asarray(data[key], dtype=bool)
            if matrix.shape == (n_channels, n_channels):
                return matrix

    if not build_neighbours:
        raise KeyError("No neighbouring matrix found in LIMO.data and automatic building is disabled.")

    chanlocs = data.get("chanlocs") or data.get("expected_chanlocs")
    if chanlocs is None:
        raise KeyError("No neighbouring matrix found and no channel locations are available to build one.")

    labels, positions = _extract_chanloc_positions(chanlocs)
    if len(labels) != n_channels:
        raise ValueError("The number of channel locations does not match the first dimension of the data.")

    if prefer_mne:
        adjacency = _build_adjacency_with_mne(labels, positions)
        if adjacency is not None:
            return adjacency

    if neighbourdist is not None:
        return _build_distance_adjacency(positions, neighbourdist)
    return _build_delaunay_adjacency(positions)


def _build_adjacency_with_mne(labels: list[str], positions: np.ndarray) -> np.ndarray | None:
    try:
        import mne
    except ImportError:
        return None

    try:
        ch_pos = {label: positions[index] for index, label in enumerate(labels)}
        montage = mne.channels.make_dig_montage(ch_pos=ch_pos, coord_frame="head")
        info = mne.create_info(ch_names=labels, sfreq=1.0, ch_types="eeg")
        info.set_montage(montage, on_missing="ignore")
        adjacency, adjacency_labels = mne.channels.find_ch_adjacency(info, ch_type="eeg")
    except Exception:
        return None

    if hasattr(adjacency, "toarray"):
        adjacency = adjacency.toarray()
    adjacency = np.asarray(adjacency, dtype=bool)
    return _reorder_adjacency_matrix(adjacency, list(adjacency_labels), labels)


def _build_distance_adjacency(positions: np.ndarray, neighbourdist: float) -> np.ndarray:
    distances = distance.squareform(distance.pdist(positions))
    adjacency = (distances <= neighbourdist) & (distances > 0)
    return adjacency.astype(bool)


def _build_delaunay_adjacency(positions: np.ndarray) -> np.ndarray:
    projected = _project_positions_to_2d(positions)
    triangulation = Delaunay(projected)
    adjacency = np.zeros((projected.shape[0], projected.shape[0]), dtype=bool)
    for simplex in triangulation.simplices:
        for start in simplex:
            for stop in simplex:
                if start != stop:
                    adjacency[start, stop] = True
    return adjacency


def _reorder_adjacency_matrix(adjacency: np.ndarray, source_labels: list[str], target_labels: list[str]) -> np.ndarray:
    lookup = {label.lower(): index for index, label in enumerate(source_labels)}
    indices = []
    for label in target_labels:
        key = label.lower()
        if key not in lookup:
            raise KeyError(f"Could not find adjacency entry for channel '{label}'.")
        indices.append(lookup[key])
    indices_array = np.asarray(indices, dtype=int)
    return adjacency[np.ix_(indices_array, indices_array)]


def _extract_chanloc_positions(chanlocs: Any) -> tuple[list[str], np.ndarray]:
    labels: list[str] = []
    positions: list[np.ndarray] = []
    for chanloc in _iter_chanlocs(chanlocs):
        label = str(_get_field(chanloc, "labels", default="")).strip()
        if not label:
            raise ValueError("Each channel location must have a label.")
        xyz = []
        for key in ("X", "Y", "Z"):
            value = _get_field(chanloc, key, default=np.nan)
            xyz.append(float(value) if value not in (None, "") else np.nan)
        xyz_array = np.asarray(xyz, dtype=float)
        if np.isnan(xyz_array).any():
            theta = _safe_float(_get_field(chanloc, "theta", default=np.nan))
            radius = _safe_float(_get_field(chanloc, "radius", default=np.nan))
            if np.isnan(theta) or np.isnan(radius):
                raise ValueError("Channel locations must provide either XYZ or theta/radius coordinates.")
            angle = np.deg2rad(theta)
            xyz_array = np.asarray([radius * np.cos(angle), radius * np.sin(angle), 0.0], dtype=float)
        labels.append(label)
        positions.append(xyz_array)
    return labels, np.vstack(positions)


def _project_positions_to_2d(positions: np.ndarray) -> np.ndarray:
    positions = np.asarray(positions, dtype=float)
    centered = positions - np.mean(positions, axis=0, keepdims=True)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    basis = vt[:2].T
    return centered @ basis


def _plot_bootstrap_distribution(boot_maxclustersum: np.ndarray, max_th: float, maxval: float) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    mass = np.sort(np.asarray(boot_maxclustersum, dtype=float))
    fig, ax = plt.subplots(num="Correction by clustering: results under H0")
    ax.plot(mass, linewidth=3)
    ax.grid(True)
    threshold_index = int(np.argmin(np.abs(mass - max_th)))
    ax.plot(threshold_index, max_th, "r*")
    ax.text(threshold_index, max_th, f"bootstrap threshold {max_th} ->", fontsize=10, ha="right")
    observed_index = int(np.argmin(np.abs(mass - maxval)))
    ax.plot(observed_index, maxval, "r*")
    ax.text(observed_index, maxval, f"biggest observed cluster mass: {maxval} ->", fontsize=10, ha="right")
    ax.set_title("Cluster-mass Maxima under H0 - no significant results", fontsize=12)
    ax.set_xlabel("sorted bootstrap iterations", fontsize=12)
    ax.set_ylabel("Freq.", fontsize=12)
    ax.set_axisbelow(False)
    plt.show()


def _label_temporal_clusters(mask: np.ndarray) -> tuple[np.ndarray, int]:
    labels, num = ndimage.label(np.asarray(mask, dtype=bool), structure=np.asarray([1, 1, 1], dtype=int))
    return labels, int(num)


def _extract_data_mapping(LIMO: Mapping[str, Any] | Any) -> Mapping[str, Any]:
    if isinstance(LIMO, Mapping):
        data = LIMO.get("data", {})
        if isinstance(data, Mapping):
            return data
    if hasattr(LIMO, "data"):
        data = getattr(LIMO, "data")
        if isinstance(data, Mapping):
            return data
    raise TypeError("LIMO must provide a mapping-like 'data' field.")


def _iter_chanlocs(chanlocs: Any) -> list[Any]:
    if isinstance(chanlocs, np.ndarray):
        return [item for item in chanlocs.flat]
    if isinstance(chanlocs, list | tuple):
        return list(chanlocs)
    return [chanlocs]


def _get_field(obj: Any, name: str, default: Any = None) -> Any:
    if isinstance(obj, Mapping):
        return obj.get(name, default)
    return getattr(obj, name, default)


def _safe_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _ensure_2d(array: np.ndarray) -> np.ndarray:
    array = np.asarray(array, dtype=float)
    if array.ndim == 1:
        return array[np.newaxis, :]
    if array.ndim != 2:
        raise ValueError("Expected a 1D or 2D array.")
    return array


def _ensure_3d(array: np.ndarray) -> np.ndarray:
    array = np.asarray(array, dtype=float)
    if array.ndim == 2:
        return array[np.newaxis, :, :]
    if array.ndim != 3:
        raise ValueError("Expected a 2D or 3D array.")
    return array