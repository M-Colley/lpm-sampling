# Changelog

## Unreleased

- Verified against BalancedSampling 2.1.1 (R 4.6.0). `spatial_balance` and
  `pi_from_size` match `sb` and `getPips` to floating point; both LPM variants
  realise the prescribed inclusion probabilities and are statistically
  indistinguishable from the reference in spread. Results are in the README.
- `compare_with_r/`: added the deterministic checks (`diagnostic.py`,
  `diagnostic.R`) and a Moran's I check on the spatial dependence each
  implementation induces between units.
- `compare_with_r/run_comparison.sh`: finds Python and Rscript itself, including
  R's default Windows install location, and probes candidate interpreters by
  running them — `python3` on Windows resolves to a Store stub that `command -v`
  reports as success. Override with `PYTHON=` / `RSCRIPT=`.
- `compare_with_r/frames.py` no longer needs `PYTHONPATH` to be set by the caller.

## 0.1.0

First release.

- `lpm1` / `lpm2`: the Local Pivotal Method (Grafström, Lundström & Schelin 2012).
- `pi_from_size`: probabilities proportional to size, with the certainty-unit
  iteration (clipping once silently breaks `sum(pi) == n`).
- `spatial_balance`: the Voronoi balance index; `systematic_pps` as the
  spatially unaware baseline to compare against.
- `unit_sphere_xyz`, `great_circle_km`, `nearest_distance_km`: geographic
  helpers, so a global frame is sampled on a sphere rather than on a rectangle.
- `compare_with_r/`: harness comparing this implementation against
  BalancedSampling on a shared frame, distributionally.
