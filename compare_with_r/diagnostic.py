"""Deterministic half of the comparison: do both sides compute the same functions?

Tests 1-3 in compare.py are distributional, because the two samplers cannot
produce identical samples (different RNGs, different tie-breaking).  The two
checks here have no such limitation.  Both cover functions that are
deterministic given their inputs, so the same inputs must give the same outputs
on both sides -- any difference beyond floating point is a real disagreement.

(a) spatial balance -- BalancedSampling::sb vs lpm_sampling.spatial_balance.
    This has to pass before test 3 means anything: test 3 compares the two
    samplers *using spatial balance as the yardstick*, which is only a fair
    comparison if both sides agree on what the yardstick measures.

(b) pps probabilities -- BalancedSampling::getPips vs lpm_sampling.pi_from_size.
    The certainty-unit iteration (cap a unit at 1, redistribute its excess,
    repeat) is the fiddliest arithmetic in the package, and it feeds every
    design weight, so it is checked against the reference directly.

Samples are exchanged as 0/1 masks rather than index lists, so the 0-based /
1-based convention gap between Python and R cannot silently shift a sample.
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from lpm_sampling import lpm, pi_from_size, spatial_balance  # noqa: E402

FRAME_DIR = os.path.join(HERE, "frames")
RESULTS = os.path.join(HERE, "results")

# Gamma shapes from mild to extreme skew, so the capping iteration runs zero,
# one and many rounds across the cases.
PPS_SHAPES = (2.0, 0.5, 0.2, 0.05)
PPS_SIZES = (10, 50, 100)
PPS_UNITS = 400


def build_balance_cases(pi, coords, cases=200, seed=11):
    """Fixed samples spanning well-spread to deliberately clumped.

    Includes non-LPM samples on purpose: a purely random sample and a
    corner-clustered one push the Voronoi counts far from 1, which is where
    two differing definitions of balance would diverge most visibly.
    """
    rng = np.random.default_rng(seed)
    n = int(round(pi.sum()))
    N = len(pi)
    masks = []

    for index in range(cases):
        kind = index % 4

        if kind == 0:
            selected = lpm(pi, coords, method="lpm2", rng=rng)
        elif kind == 1:
            selected = lpm(pi, coords, method="lpm1", rng=rng)
        elif kind == 2:
            selected = rng.choice(N, size=n, replace=False)
        else:
            # Worst case: n units packed into one corner of the frame.
            corner = np.argsort(coords[:, 0] + coords[:, 1])[:3 * n]
            selected = rng.choice(corner, size=n, replace=False)

        mask = np.zeros(N, dtype=int)
        mask[selected] = 1
        masks.append(mask)

    return np.array(masks)


def build_pps_cases(seed=3):
    """Size vectors and target sample sizes covering a range of skew."""
    rng = np.random.default_rng(seed)
    sizes, targets = [], []

    for shape in PPS_SHAPES:
        for n in PPS_SIZES:
            # The offset keeps every size strictly positive: getPips documents
            # its input as positive numbers, and a zero would be a different
            # question (units that cannot be sampled) than the one asked here.
            sizes.append(rng.gamma(shape, 1000.0, size=PPS_UNITS) + 1e-3)
            targets.append(n)

    return np.column_stack(sizes), np.array(targets)


def main():
    frame = np.genfromtxt(os.path.join(FRAME_DIR, "frame.csv"), delimiter=",", names=True)
    coords = np.stack([frame["x"], frame["y"]], axis=1)
    pi = np.asarray(frame["pi"], dtype=float)
    os.makedirs(RESULTS, exist_ok=True)

    # (a) spatial balance
    masks = build_balance_cases(pi, coords)
    np.savetxt(os.path.join(FRAME_DIR, "diag_masks.csv"), masks, fmt="%d", delimiter=",")
    balance = np.array([spatial_balance(pi, coords, np.flatnonzero(m)) for m in masks])
    np.savetxt(os.path.join(RESULTS, "python_diag_sb.csv"), balance,
               delimiter=",", header="balance", comments="")
    print("balance: %d cases (range %.5f to %.5f)" % (len(masks), balance.min(), balance.max()))

    # (b) pps probabilities
    sizes, targets = build_pps_cases()
    np.savetxt(os.path.join(FRAME_DIR, "diag_pps_sizes.csv"), sizes, delimiter=",")
    np.savetxt(os.path.join(FRAME_DIR, "diag_pps_n.csv"), targets, fmt="%d", delimiter=",")
    probabilities = np.column_stack(
        [pi_from_size(sizes[:, k], targets[k])[0] for k in range(sizes.shape[1])]
    )
    # A header row on both sides keeps every results file on one convention:
    # one header line, then data. compare.py's loader skips exactly one line,
    # and a headerless matrix would silently lose its first unit.
    pps_header = ",".join("case%d" % k for k in range(probabilities.shape[1]))
    np.savetxt(os.path.join(RESULTS, "python_diag_pps.csv"), probabilities,
               delimiter=",", header=pps_header, comments="")
    capped = int((probabilities >= 1.0 - 1e-12).sum())
    print("pps:     %d cases over %d units (%d certainty units capped in total)"
          % (sizes.shape[1], PPS_UNITS, capped))


if __name__ == "__main__":
    main()
