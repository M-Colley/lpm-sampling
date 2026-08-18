"""Fixed frames shared by both implementations, so the comparison is paired.

Both sides must sample the SAME units with the SAME inclusion probabilities;
only the sampler differs. The frames are written to CSV once and read by both
the Python and the R driver.
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FRAME_DIR = os.path.join(HERE, "frames")

# Same as the other drivers here: make the package importable straight from a
# checkout, so the harness needs no install and no PYTHONPATH set by the caller.
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from lpm_sampling import pi_from_size  # noqa: E402


def clustered_frame(n=300, seed=1):
    """Clustered 2-D frame: clusters are what spatial balance is measured on."""
    rng = np.random.default_rng(seed)
    centres = rng.uniform(0.0, 100.0, size=(5, 2))
    points = [rng.normal(c, 4.0, size=(n // 8, 2)) for c in centres]
    points.append(rng.uniform(0.0, 100.0, size=(n - sum(len(p) for p in points), 2)))
    coords = np.vstack(points)
    sizes = rng.gamma(2.0, 50000.0, size=len(coords))

    return coords, sizes


def write_frames(n=300, seed=1, sample_size=30):
    os.makedirs(FRAME_DIR, exist_ok=True)
    coords, sizes = clustered_frame(n=n, seed=seed)
    pi, _certainty = pi_from_size(sizes, sample_size)
    path = os.path.join(FRAME_DIR, "frame.csv")

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("x,y,pi\n")

        for (x, y), p in zip(coords, pi):
            handle.write("%.10f,%.10f,%.12f\n" % (x, y, p))

    print("wrote %s (N=%d, n=%d, sum pi=%.9f)" % (path, len(coords), sample_size, pi.sum()))

    return path


if __name__ == "__main__":
    write_frames()
