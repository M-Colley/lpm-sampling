"""Inclusion probabilities proportional to a size measure."""

from __future__ import annotations

import numpy as np

__all__ = ["pi_from_size"]

EPS = 1e-12


def pi_from_size(sizes, n, *, eps=1e-12):
    """Inclusion probabilities proportional to size, capped at 1.

    Starts from ``pi = n * sizes / sum(sizes)`` and then applies the standard
    iterative capping: any unit whose pi would exceed 1 becomes a *certainty*
    unit (pi = 1) and its excess is redistributed proportionally among the
    remaining units, repeated until no further unit exceeds 1.

    Parameters
    ----------
    sizes : array_like, shape (N,)
        Non-negative size measure.  Units with size 0 get pi = 0.
    n : int
        Target sample size.  Must satisfy ``0 <= n <= (# units with size > 0)``.

    Returns
    -------
    pi : ndarray, shape (N,), float
        Inclusion probabilities, ``sum(pi) == n`` (up to float error).
    certainty : ndarray, int
        Indices of the units that were capped to 1.

    Raises
    ------
    ValueError
        If sizes are negative, or n exceeds the number of positive-size units.
    """
    sizes = np.asarray(sizes, dtype=float).ravel()
    if np.any(sizes < 0) or not np.all(np.isfinite(sizes)):
        raise ValueError("sizes must be finite and non-negative")
    n = int(n)
    N = sizes.size
    if n < 0:
        raise ValueError("n must be non-negative")

    pi = np.zeros(N, dtype=float)
    if n == 0:
        return pi, np.empty(0, dtype=int)

    positive = sizes > 0
    n_pos = int(positive.sum())
    if n > n_pos:
        raise ValueError(
            f"n={n} exceeds the number of units with positive size ({n_pos})"
        )

    certainty = np.zeros(N, dtype=bool)
    while True:
        free = positive & ~certainty
        n_free = n - int(certainty.sum())
        if n_free <= 0:
            # All of the sample is taken by certainty units.
            pi[free] = 0.0
            break
        total = sizes[free].sum()
        pi_free = n_free * sizes[free] / total
        over = pi_free > 1.0 + eps
        if not np.any(over):
            pi[free] = np.minimum(pi_free, 1.0)
            break
        idx_free = np.flatnonzero(free)
        certainty[idx_free[over]] = True

    pi[certainty] = 1.0
    pi[~positive] = 0.0
    pi = np.clip(pi, 0.0, 1.0)
    return pi, np.flatnonzero(certainty)


# --------------------------------------------------------------------------
# the pivotal pair update
# --------------------------------------------------------------------------
