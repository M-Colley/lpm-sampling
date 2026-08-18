# Changelog

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
