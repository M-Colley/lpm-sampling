"""Draw a spatially balanced sample of cities weighted by population.

    python examples/select_cities.py cities.csv 150

The CSV needs columns: name, lat, lon, population.
"""

import csv
import sys

import numpy as np

from lpm_sampling import (
    lpm2,
    nearest_distance_km,
    pi_from_size,
    spatial_balance,
    systematic_pps,
    unit_sphere_xyz,
)


def main():
    path = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 100

    with open(path, encoding="utf-8-sig", newline="") as handle:
        rows = [r for r in csv.DictReader(handle) if r.get("lat") and r.get("lon")]

    lat = np.array([float(r["lat"]) for r in rows])
    lon = np.array([float(r["lon"]) for r in rows])
    population = np.array([max(float(r.get("population") or 0.0), 1.0) for r in rows])

    xyz = unit_sphere_xyz(lat, lon)
    pi, certainty = pi_from_size(population, n)
    rng = np.random.default_rng(20260818)
    sample = lpm2(pi, xyz, rng=rng)

    distance = nearest_distance_km(xyz, xyz[sample])
    baseline = systematic_pps(pi, rng=rng)

    print("frame %d units, drew %d, %d certainty units" % (len(rows), len(sample), len(certainty)))
    print("spatial balance  : %.4f  (spatially unaware baseline %.4f, lower is better)"
          % (spatial_balance(pi, xyz, sample), spatial_balance(pi, xyz, baseline)))
    print("population-weighted distance to nearest selected unit: %.0f km"
          % float((population * distance).sum() / population.sum()))
    print("")

    for rank, index in enumerate(sorted(sample, key=lambda i: -population[i])[:15], start=1):
        print("  %2d. %-24s pi=%.4f  weight=%8.1f" % (
            rank, rows[index].get("name", "?"), pi[index], 1.0 / pi[index]))


if __name__ == "__main__":
    main()
