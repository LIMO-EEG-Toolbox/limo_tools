"""Python port of ``limo_IRLS.m``.

Use ``limo_irls(X, Y)`` to compute the LIMO iteratively reweighted least-
squares solution for a design matrix ``X`` and trial-by-frame data matrix
``Y``.

Inputs:
    ``X``: design matrix with shape ``trials x parameters``.
    ``Y``: EEG data with shape ``trials x frames``.
    ``tune``: bisquare tuning constant. Lower values downweight large
    residuals more strongly. Default is ``4.685``.
    ``figure``: ``"on"`` or ``"off"``. When ``"on"``, the routine attempts
    to plot residual convergence across iterations.

Outputs:
    ``b``: beta coefficients with shape ``parameters x frames``.
    ``w``: weights with shape ``trials x frames``.

References:
    GroB, J. (2003). Linear Regression, pp. 191-215.
    Street, J. O., Carroll, R. J., and Ruppert, D. (1988). A note on
    computing robust regression estimates via iteratively reweighted least
    squares. The American Statistician.
    Wager, T. D., Keller, M. C., Lacey, S. C., and Jonides, J. (2005).
    Increased sensitivity in neuroimaging analyses using robust regression.
    NeuroImage, 26(1), 99-113.

This file follows the MATLAB ``limo_IRLS.m`` algorithm directly and also
provides a ``limo_IRLS(...)`` compatibility alias.
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np


def limo_irls(
    X: Any,
    Y: Any,
    *options: Any,
    tune: float = 4.685,
    figure: str = "off",
    iterlim: int = 100,
) -> tuple[np.ndarray, np.ndarray]:
    """Replicate the MATLAB ``limo_IRLS`` routine."""

    if options:
        tune, figure = _parse_option_pairs(options, tune=tune, figure=figure)

    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)

    if X.ndim != 2 or Y.ndim != 2:
        raise ValueError("X and Y must both be 2D arrays.")

    rows, cols = X.shape
    if rows <= cols:
        raise ValueError("IRLS cannot be computed, there is not enough trials for this design")

    b = np.linalg.pinv(X) @ Y

    hat = np.diag(X @ np.linalg.pinv(X.T @ X) @ X.T)
    adjfactor = 1.0 / np.sqrt(1.0 - hat)
    adjfactor[np.isinf(adjfactor)] = 1.0

    numiter = 0
    old_res = 1.0
    new_res = 10.0
    convergence_history = np.full(iterlim, np.nan, dtype=float)

    while np.max(np.abs(old_res - new_res)) > 1e-4:
        numiter += 1
        old_res = new_res

        if numiter > iterlim:
            warnings.warn(
                "limo_irls could not converge after "
                f"{iterlim} iterations. iteration limit can be adjusted, "
                "with no guarantee it improves convergence.",
                RuntimeWarning,
                stacklevel=2,
            )
            break

        res = Y - X @ b
        resadj = res * adjfactor[:, np.newaxis]

        re = np.median(np.abs(resadj), axis=0) / 0.6745
        re[re < 1e-5] = 1e-5
        r = resadj / (tune * re)[np.newaxis, :]

        w = (np.abs(r) < 1.0) * (1.0 - r**2) ** 2
        w = np.sqrt(w)
        yw = Y * w

        for frame in range(Y.shape[1]):
            xw = X * w[:, frame][:, np.newaxis]
            b[:, frame] = np.linalg.pinv(xw) @ yw[:, frame]

        new_res = float(np.sum(res.ravel() ** 2))
        convergence_history[numiter - 1] = abs(old_res - new_res)

        if figure.lower() == "on":
            _plot_convergence(
                iteration=numiter,
                convergence=convergence_history,
                iterlim=iterlim,
            )

    return b, w


def limo_IRLS(
    X: Any,
    Y: Any,
    *options: Any,
    tune: float = 4.685,
    figure: str = "off",
    iterlim: int = 100,
) -> tuple[np.ndarray, np.ndarray]:
    """Compatibility alias for ``limo_irls``."""

    return limo_irls(X, Y, *options, tune=tune, figure=figure, iterlim=iterlim)


def _parse_option_pairs(
    options: tuple[Any, ...],
    *,
    tune: float,
    figure: str,
) -> tuple[float, str]:
    if len(options) % 2 != 0:
        raise ValueError("Optional arguments must be provided as keyword/value pairs.")

    parsed_tune = tune
    parsed_figure = figure

    for index in range(0, len(options), 2):
        key = str(options[index])
        value = options[index + 1]
        if "fig" in key.lower():
            parsed_figure = str(value)
        elif key.lower() == "tune":
            parsed_tune = float(value)

    return parsed_tune, parsed_figure


def _plot_convergence(*, iteration: int, convergence: np.ndarray, iterlim: int) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        warnings.warn(
            "matplotlib is not installed, so the IRLS convergence figure cannot be shown.",
            RuntimeWarning,
            stacklevel=2,
        )
        return

    current = convergence[:iteration]
    plt.figure("Residual Mean Squares")
    plt.plot(iteration, current[-1], "ro", linewidth=3)
    ymax = np.nanmax(current)
    ymin = np.nanmin(current)
    plt.axis([1, iterlim + 0.5, -0.1, ymax + 0.1 * ymin if np.isfinite(ymax) and np.isfinite(ymin) else 1])
    if iteration > 2:
        plt.title(f"iteration {iteration} convergence {current[-1]}")
    plt.grid(True)
    plt.draw()
    if iteration == iterlim:
        plt.axis([1, iterlim + 0.5, -0.1 * ymax, ymax + 0.1 * ymin])
