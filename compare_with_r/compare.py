"""Compare this package against BalancedSampling.

Two kinds of check, because two kinds are possible:

* Deterministic (test 0).  spatial_balance and pi_from_size are plain functions
  of their inputs, so given the same inputs both sides must return the same
  numbers.  These are equality checks.
* Distributional (tests 1-4).  The two samplers use different random number
  generators and break ties differently, so they cannot and should not produce
  identical samples.  What must agree is the DESIGN: the inclusion probabilities
  each one realises, the distribution of spatial balance it achieves, and the
  spatial dependence it induces between units.

Run, from this directory -- or just use ./run_comparison.sh, which does all of
it and finds the interpreters itself:

    python frames.py
    python diagnostic.py
    Rscript diagnostic.R
    python run_python.py --replicates 5000 --method lpm2
    Rscript run_r.R 5000 lpm2
    python compare.py
"""

from __future__ import annotations

import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")


def load(name):
    path = os.path.join(RESULTS, name)

    if not os.path.isfile(path):
        raise SystemExit("missing %s -- run the driver that produces it first" % path)

    return np.genfromtxt(path, delimiter=",", skip_header=1)


def ks_two_sample(a, b):
    """Kolmogorov-Smirnov statistic and its asymptotic p-value."""
    a = np.sort(a)
    b = np.sort(b)
    grid = np.concatenate([a, b])
    cdf_a = np.searchsorted(a, grid, side="right") / len(a)
    cdf_b = np.searchsorted(b, grid, side="right") / len(b)
    statistic = float(np.max(np.abs(cdf_a - cdf_b)))
    scale = math.sqrt(len(a) * len(b) / (len(a) + len(b)))
    lam = (scale + 0.12 + 0.11 / scale) * statistic
    p = 2.0 * sum((-1.0) ** (k - 1) * math.exp(-2.0 * k * k * lam * lam)
                  for k in range(1, 101))

    return statistic, min(1.0, max(0.0, p))


def report_diagnostic():
    """Tests 0a/0b: do both sides compute the same deterministic functions?

    Deterministic, so these are equality checks rather than statistical ones.
    0a has to pass before test 3 means anything: test 3 compares the two
    samplers using spatial balance as the yardstick, which is only a fair
    comparison if both sides agree on what the yardstick measures.
    """
    print("")
    print("0. Deterministic checks -- same inputs, so the outputs must match exactly")

    py_sb = load("python_diag_sb.csv")
    r_sb = load("r_diag_sb.csv")

    if py_sb.shape != r_sb.shape:
        raise SystemExit("balance case counts differ: python %s, R %s"
                         % (py_sb.shape, r_sb.shape))

    sb_rel = float(np.max(np.abs(py_sb - r_sb) / np.maximum(np.abs(r_sb), 1e-300)))
    print("   a. spatial_balance vs sb:      %d samples (balance %.4f to %.4f)"
          % (len(py_sb), min(py_sb.min(), r_sb.min()), max(py_sb.max(), r_sb.max())))
    print("      max relative difference = %.3e  (%s)"
          % (sb_rel, "identical to floating point" if sb_rel < 1e-12 else "DIFFERENT"))

    if sb_rel >= 1e-12:
        print("      WARNING: the yardsticks disagree, so test 3 below is not meaningful.")

    py_pps = load("python_diag_pps.csv")
    r_pps = load("r_diag_pps.csv")

    if py_pps.shape != r_pps.shape:
        raise SystemExit("pps case shapes differ: python %s, R %s"
                         % (py_pps.shape, r_pps.shape))

    # Absolute, not relative: these are probabilities on a fixed [0, 1] scale,
    # and a relative measure would be dominated by the near-zero entries.
    pps_abs = float(np.max(np.abs(py_pps - r_pps)))
    capped = int((r_pps >= 1.0 - 1e-12).sum())
    print("   b. pi_from_size vs getPips:    %d cases over %d units, %d certainty units"
          % (py_pps.shape[1], py_pps.shape[0], capped))
    print("      max absolute difference = %.3e  (%s)"
          % (pps_abs, "identical to floating point" if pps_abs < 1e-12 else "DIFFERENT"))

    return sb_rel, pps_abs


def morans_i(values, neighbours):
    """Moran's I of ``values`` under row-standardised k-nearest-neighbour weights."""
    centred = values - values.mean()
    lag = centred[neighbours].mean(axis=1)

    return float((centred * lag).sum() / (centred ** 2).sum())


def report_spatial_structure(pi, coords, py_freq, r_freq, reps, k=8):
    """Test 4: do the two induce the same spatial dependence between units?

    Tests 1-3 look at each unit's marginal frequency and at the balance
    distribution.  Neither looks at how selections are correlated *between*
    neighbouring units, which is the defining property of the method.

    Expect Moran's I to be NEGATIVE for both: neighbours compete for the same
    probability mass, so a unit over-selected in a finite run leaves its
    neighbours under-selected.  The comparison that matters is between the two
    columns, not against zero -- testing I against an exchangeable permutation
    null would reject for both implementations and prove nothing, because the
    correlation is a property of the design rather than a defect.
    """
    distance = ((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(distance, np.inf)
    neighbours = np.argsort(distance, axis=1)[:, :k]
    se = np.sqrt(np.maximum(pi * (1.0 - pi), 1e-12) / reps)

    py_i = morans_i((py_freq - pi) / se, neighbours)
    r_i = morans_i((r_freq - pi) / se, neighbours)

    print("")
    print("4. Do they induce the same spatial dependence between units?")
    print("   Moran's I of (freq - pi), %d-nearest-neighbour weights:" % k)
    print("   python %+.4f   R %+.4f   difference %+.4f" % (py_i, r_i, py_i - r_i))
    print("   Both negative is the expected signature of spatial balance; what")
    print("   matters here is that the two agree, not that they differ from zero.")

    return py_i, r_i


def report(method="lpm2"):
    frame = np.genfromtxt(os.path.join(HERE, "frames", "frame.csv"),
                          delimiter=",", names=True)
    pi = np.asarray(frame["pi"], dtype=float)

    py_freq = load("python_%s_freq.csv" % method)
    r_freq = load("r_%s_freq.csv" % method)
    py_bal = load("python_%s_balance.csv" % method)
    r_bal = load("r_%s_balance.csv" % method)

    print("")
    print("Comparison against BalancedSampling::%s" % method)
    print("=" * (36 + len(method)))
    print("frame: N = %d, n = %.0f, replicates: python %d, R %d"
          % (len(pi), pi.sum(), len(py_bal), len(r_bal)))

    print("")
    print("1. Does each implementation realise the prescribed inclusion probabilities?")

    for label, freq, reps in (("python", py_freq, len(py_bal)), ("R", r_freq, len(r_bal))):
        se = np.sqrt(np.maximum(pi * (1.0 - pi), 1e-12) / reps)
        z = (freq - pi) / se
        print("   %-7s max |freq - pi| = %.5f  (worst %.2f SE, mean z^2 = %.3f, expect 1.0)"
              % (label, np.max(np.abs(freq - pi)), np.max(np.abs(z)), float(np.mean(z ** 2))))

    print("")
    print("2. Are the two consistent with the SAME design?")
    pooled_se = np.sqrt(np.maximum(pi * (1.0 - pi), 1e-12)
                        * (1.0 / len(py_bal) + 1.0 / len(r_bal)))
    z_pair = (py_freq - r_freq) / pooled_se
    print("   paired z over %d units: max |z| = %.2f, mean z^2 = %.3f (expect ~1.0)"
          % (len(pi), np.max(np.abs(z_pair)), float(np.mean(z_pair ** 2))))
    print("   units beyond 3 SE: %d (expect ~%.1f by chance)"
          % (int((np.abs(z_pair) > 3.0).sum()), 0.0027 * len(pi)))

    print("")
    print("3. Do they spread equally well?")
    statistic, p_value = ks_two_sample(py_bal, r_bal)
    print("   mean spatial balance: python %.5f, R %.5f (difference %+.5f, %.1f%%)"
          % (py_bal.mean(), r_bal.mean(), py_bal.mean() - r_bal.mean(),
             100.0 * (py_bal.mean() - r_bal.mean()) / r_bal.mean()))
    print("   Kolmogorov-Smirnov D = %.4f, p = %.3f" % (statistic, p_value))
    print("")
    print("   Read: p above ~0.01 means the balance distributions are indistinguishable.")
    print("   A large mean difference with a tiny p would be the real finding, and would")
    print("   mean the implementations differ in how they spread, not merely in RNG.")

    coords = np.stack([frame["x"], frame["y"]], axis=1)
    report_spatial_structure(pi, coords, py_freq, r_freq, len(py_bal))


if __name__ == "__main__":
    try:
        report_diagnostic()
    except SystemExit as error:
        print("\nskipping the deterministic check: %s" % error)

    for method in ("lpm2", "lpm1"):
        try:
            report(method)
        except SystemExit as error:
            print("\nskipping %s: %s" % (method, error))
