"""Python port of ``limo_WLS.m``.

Use ``limo_WLS(X, Y)`` to compute the LIMO weighted least-squares solution
for a design matrix ``X`` and trial-by-frame data matrix ``Y``.

Inputs:
    ``X``: design matrix with shape ``trials x parameters``.
    ``Y``: EEG data with shape ``trials x frames``.

Outputs:
    ``b``: beta coefficients with shape ``parameters x frames``.
    ``W``: one weight per trial.
    ``rf``: reduction factor returned by the principal-components projection.

This port follows the MATLAB algorithm in ``limo_WLS.m`` and includes a local
implementation.

About ``limo_pcout``:
    ``limo_pcout`` is the principal-components projection method used by LIMO
    to identify multivariate trial outliers and convert them into robust trial
    weights. It computes location and scatter weights in principal-component
    space, then combines them into the final distance-based weights used by
    ``limo_WLS``.

Reference:
    Filzmoser, P., Maronna, R., and Werner, M. (2008). Outlier
    identification in high dimensions. Computational Statistics and Data
    Analysis, 52, 1694-1711.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.stats import chi2


def limo_WLS(X: Any, Y: Any) -> tuple[np.ndarray, np.ndarray, int]:
    """Replicate the MATLAB ``limo_WLS`` weighted least-squares routine."""

    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)

    if X.ndim != 2 or Y.ndim != 2:
        raise ValueError("X and Y must both be 2D arrays.")

    rows, cols = X.shape
    if rows <= cols:
        raise ValueError("WLS cannot be computed, there is not enough trials for this design")

    median_profile = np.median(np.abs(Y - np.median(Y, axis=0, keepdims=True)), axis=0)
    if Y[:, median_profile > 1e-6].size == 0:
        raise ValueError(
            "WLS cannot be computed, for at least 1 condition, all trials have the same values"
        )

    hat = np.diag(X @ np.linalg.pinv(X.T @ X) @ X.T)
    hat = np.minimum(hat, 1.0)

    adjfactor = 1.0 / np.sqrt(1.0 - hat)
    adjfactor[np.isinf(adjfactor)] = 1.0

    b = np.linalg.pinv(X) @ Y

    tune = 4.685
    res = Y - X @ b
    resadj = res * adjfactor[:, np.newaxis]

    re = np.median(np.abs(resadj), axis=0) / 0.6745
    re[re < 1e-5] = 1e-5
    r = resadj / (tune * re)[np.newaxis, :]

    W, _, rf, _, _ = limo_pcout(r, figure="off")
    WY = Y * W[:, np.newaxis]
    WX = X * W[:, np.newaxis]
    b = np.linalg.pinv(WX) @ WY
    return b, W, rf


def limo_pcout(
    x: Any,
    *,
    downsample: str = "on",
    figure: str = "off",
    weightsas: str = "Kernel",
    xaxis: Any | None = None,
) -> tuple[np.ndarray, np.ndarray, int, np.ndarray, np.ndarray]:
    """Python port of the non-plotting core of ``limo_pcout``.

    The plotting-related arguments are accepted for interface compatibility,
    but no figure is produced.
    """

    del figure, weightsas, xaxis

    x = np.asarray(x, dtype=float)
    if x.ndim != 2:
        raise ValueError("x must be a 2D array")

    x = x[:, _mad(x, axis=0) > 1e-6]
    if x.size == 0:
        raise ValueError("WLS cannot be computed, for at least 1 frame, all trials have the same values")

    n, p = x.shape
    if n < p and downsample.lower() == "on":
        f = p / n
        if int(np.ceil(f)) == 2:
            x = x[:, ::2]
            n, p = x.shape

    madx = _mad(x, axis=0) * 1.4826
    madx[madx == 0] = np.finfo(float).eps
    x2 = (x - np.median(x, axis=0, keepdims=True)) / madx[np.newaxis, :]
    x3 = x2 - np.mean(x2, axis=0, keepdims=True)

    singular_values = np.linalg.svd(x3, compute_uv=False)
    a = singular_values**2 / (n - 1)
    cumulative = np.cumsum(a) / np.sum(a)
    above_threshold = np.flatnonzero(cumulative > 0.99)
    p1 = int(above_threshold[0] + 1) if above_threshold.size else len(a)
    rf = int(p - p1)

    _, _, vh = np.linalg.svd(x3, full_matrices=False)
    xpc = x2 @ vh[:p1, :].T
    madxpc = _mad(xpc, axis=0) * 1.4826
    madxpc[madxpc == 0] = np.finfo(float).eps
    xpcsc = (xpc - np.median(xpc, axis=0, keepdims=True)) / madxpc[np.newaxis, :]

    wp = np.abs(np.mean(xpcsc**4, axis=0) - 3.0)
    wp_sum = np.sum(wp)
    if wp_sum == 0:
        wp_scaled = np.ones_like(wp) / max(wp.size, 1)
    else:
        wp_scaled = wp / wp_sum
    xpcwsc = xpcsc @ np.diag(wp_scaled)
    xpcnorm = np.sqrt(np.sum(xpcwsc**2, axis=1))

    chi_half = np.sqrt(chi2.ppf(0.5, p1))
    xdist1 = (xpcnorm * chi_half) / np.median(xpcnorm)
    # MATLAB's default quantile uses midpoint plotting positions (Hazen),
    # unlike NumPy's default inclusive interpolation.
    M1 = np.quantile(xdist1, 1.0 / 3.0, method="hazen")
    const1 = np.median(xdist1) + 2.5 * _mad(xdist1, axis=0) * 1.4826
    w1 = _translated_biweight(xdist1, M1, const1)

    xpcnorm = np.sqrt(np.sum(xpcsc**2, axis=1))
    xdist2 = (xpcnorm * chi_half) / np.median(xpcnorm)
    M2 = np.sqrt(chi2.ppf(0.25, p1))
    const2 = np.sqrt(chi2.ppf(0.99, p1))
    w2 = _translated_biweight(xdist2, M2, const2)

    s = 0.25
    dist = (w1 + s) * (w2 + s) / ((1 + s) ** 2)
    out = np.round(dist + s)
    return dist, out, rf, w1, w2


def _translated_biweight(xdist: np.ndarray, m: float, const: float) -> np.ndarray:
    if np.isclose(const, m):
        weights = np.ones_like(xdist)
        weights[xdist > const] = 0.0
        return weights

    weights = (1 - ((xdist - m) / (const - m)) ** 2) ** 2
    weights = np.asarray(weights, dtype=float)
    weights[xdist < m] = 1.0
    weights[xdist > const] = 0.0
    return weights


def _mad(x: np.ndarray, axis: int) -> np.ndarray:
    median = np.median(x, axis=axis, keepdims=True)
    return np.median(np.abs(x - median), axis=axis)
