"""Spatial balance diagnostics and a spatially unaware baseline design."""

from __future__ import annotations

import numpy as np

__all__ = ["spatial_balance", "systematic_pps"]


def systematic_pps(pi, rng=None, return_mask=False):
    """Randomised systematic pps sampling -- the spatially *unaware* baseline.

    Permutes the units at random, walks the cumulative sum of pi and takes the
    units straddling ``u, u+1, ..., u+n-1`` for ``u ~ Uniform(0, 1)``.  Gives
    exactly n units and respects the inclusion probabilities, but ignores the
    coordinates entirely -- which is exactly what LPM is meant to improve on.
    """
    if rng is None:
        rng = np.random.default_rng()
    p = np.asarray(pi, dtype=float).ravel()
    N = p.size
    n = int(round(p.sum()))
    order = rng.permutation(N)
    cum = np.cumsum(p[order])
    marks = rng.random() + np.arange(n, dtype=float)
    picks = order[np.searchsorted(cum, marks, side="right")]
    mask = np.zeros(N, dtype=bool)
    mask[picks] = True
    if return_mask:
        return mask
    return np.flatnonzero(mask)


# --------------------------------------------------------------------------
# spatial balance
# --------------------------------------------------------------------------


def spatial_balance(pi, coords, selected, *, chunk=4096):
    """Voronoi spatial balance B (Stevens & Olsen); lower is better.

    Every unit of the frame is assigned to its nearest SELECTED unit.  For each
    selected unit i, ``v_i`` is the sum of the inclusion probabilities of the
    units assigned to it.  A perfectly balanced sample has every ``v_i == 1``
    (each selected unit "represents" exactly one unit's worth of probability),
    so::

        B = mean_i (v_i - 1)^2

    Parameters
    ----------
    pi : array_like, shape (N,)
    coords : array_like, shape (N, d) or (N,)
    selected : array_like
        Indices of the selected units, or a boolean mask of length N.
    chunk : int
        Rows of the frame processed per block, to bound memory at O(chunk * n).

    Returns
    -------
    float
        B, or ``nan`` if nothing is selected.
    """
    p = np.asarray(pi, dtype=float).ravel()
    coords = np.asarray(coords, dtype=float)
    if coords.ndim == 1:
        coords = coords.reshape(-1, 1)
    sel = np.asarray(selected)
    if sel.dtype == bool:
        sel = np.flatnonzero(sel)
    sel = sel.astype(np.intp).ravel()
    if sel.size == 0:
        return float("nan")

    S = coords[sel]
    v = np.zeros(sel.size, dtype=float)
    for a in range(0, coords.shape[0], chunk):
        block = coords[a:a + chunk]
        # squared Euclidean distance, block x selected
        d2 = (
            np.einsum("ij,ij->i", block, block)[:, None]
            - 2.0 * block @ S.T
            + np.einsum("ij,ij->i", S, S)[None, :]
        )
        owner = np.argmin(d2, axis=1)
        v += np.bincount(owner, weights=p[a:a + chunk], minlength=sel.size)
    return float(np.mean((v - 1.0) ** 2))
