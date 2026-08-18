"""Draw many LPM samples with this package and record them for comparison."""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from lpm_sampling import lpm, spatial_balance  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--frame", default=os.path.join(HERE, "frames", "frame.csv"))
    parser.add_argument("--replicates", type=int, default=5000)
    parser.add_argument("--method", default="lpm2", choices=("lpm1", "lpm2"))
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--out", default=os.path.join(HERE, "results"))
    args = parser.parse_args()

    data = np.genfromtxt(args.frame, delimiter=",", names=True)
    coords = np.stack([data["x"], data["y"]], axis=1)
    pi = np.asarray(data["pi"], dtype=float)
    rng = np.random.default_rng(args.seed)

    counts = np.zeros(len(pi))
    balances = np.empty(args.replicates)
    started = time.time()

    for replicate in range(args.replicates):
        selected = lpm(pi, coords, method=args.method, rng=rng)
        counts[selected] += 1.0
        balances[replicate] = spatial_balance(pi, coords, selected)

    elapsed = time.time() - started
    os.makedirs(args.out, exist_ok=True)
    np.savetxt(os.path.join(args.out, "python_%s_freq.csv" % args.method),
               counts / args.replicates, delimiter=",", header="frequency", comments="")
    np.savetxt(os.path.join(args.out, "python_%s_balance.csv" % args.method),
               balances, delimiter=",", header="balance", comments="")

    print("python %s: %d replicates in %.1f s (%.2f ms per draw)"
          % (args.method, args.replicates, elapsed, 1000.0 * elapsed / args.replicates))
    print("mean spatial balance %.5f" % balances.mean())


if __name__ == "__main__":
    main()
