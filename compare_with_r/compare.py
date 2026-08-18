"""Compare this package against BalancedSampling, distributionally.

The two implementations use different random number generators and break ties
differently, so they cannot and should not produce identical samples. What must
agree is the DESIGN: the inclusion probabilities each one realises, and the
distribution of spatial balance it achieves.

Run, from this directory:

    python frames.py
    python run_python.py --replicates 5000
    Rscript run_r.R 5000
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


if __name__ == "__main__":
    for method in ("lpm2", "lpm1"):
        try:
            report(method)
        except SystemExit as error:
            print("\nskipping %s: %s" % (method, error))
