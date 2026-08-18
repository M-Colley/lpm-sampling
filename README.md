# lpm-sampling

Spatially balanced probability sampling in Python: the **Local Pivotal Method**
(LPM1 and LPM2), probabilities proportional to size, and the diagnostics needed
to show a reviewer that the sample is what you claim it is.

```python
import numpy as np
from lpm_sampling import lpm2, pi_from_size, unit_sphere_xyz, spatial_balance

lat, lon, population = load_your_frame()          # N candidate locations

pi, certainty = pi_from_size(population, n=150)   # probability proportional to size
xyz = unit_sphere_xyz(lat, lon)                   # sample on the sphere, not a rectangle
sample = lpm2(pi, xyz, rng=np.random.default_rng(20260818))

weights = 1.0 / pi[sample]                        # Horvitz-Thompson design weights
print(spatial_balance(pi, xyz, sample))           # lower is better spread
```

## Why this exists

Hand-picking sites — or taking one per country, one per region — gives you no
inclusion probabilities, and therefore no design weights and no defensible
standard errors. Simple random sampling gives you the probabilities but clumps:
you get three sites in one city and nothing for 400 km.

The Local Pivotal Method gives you both. Every unit carries an inclusion
probability (proportional to population, area, or whatever size measure you
choose). The sampler repeatedly picks a unit and its **nearest neighbour** and
lets them compete for their combined probability: one moves toward selection,
the other toward rejection, in a way that leaves each unit's probability
unchanged in expectation. Because the competitors are always neighbours,
selecting one pushes probability *out of its neighbourhood*, so the realised
sample spreads itself across the space. When every probability has resolved to
0 or 1, the 1s are your sample.

The result is a strict probability design: `pi` comes out exactly as assigned,
so `1/pi` is a valid design weight.

**LPM1** requires the two competitors to be *mutual* nearest neighbours — a
slightly better spread for more work. **LPM2** takes the chosen unit's nearest
neighbour and is the practical default.

## Install

```bash
pip install lpm-sampling
```

Only dependency is numpy.

## What is in the box

| Function | What it does |
|---|---|
| `lpm2(pi, coords)` / `lpm1(pi, coords)` | Draw a spatially balanced sample |
| `lpm(pi, coords, method=...)` | Either variant, by name |
| `pi_from_size(sizes, n)` | Inclusion probabilities proportional to size, with the certainty-unit iteration |
| `spatial_balance(pi, coords, sample)` | Voronoi balance index; 0 is perfect, lower is better |
| `systematic_pps(pi)` | Randomised systematic pps — the spatially *unaware* baseline to compare against |
| `unit_sphere_xyz(lat, lon)` | Degrees to unit-sphere coordinates |
| `great_circle_km`, `nearest_distance_km` | Distances for diagnostics |

## Two traps this package exists to avoid

**Never hand raw latitude/longitude to a spatial sampler.** Longitude is
circular but cut at ±180°, and Euclidean distance treats that cut as a real
gap: Suva (+178.4°) and Apia (−171.8°) are about 1,200 km apart but look 350°
apart, i.e. maximally distant. A degree-space sampler therefore believes the
Pacific rim is empty on both sides and cheerfully selects both. Degrees are also
anisotropic — one degree of longitude is 111 km at the equator and 48 km at
Reykjavik. `unit_sphere_xyz` fixes both: chord distance on the sphere is
monotone in great-circle distance, with no seam and no pole degeneracy. Both
traps have tests.

**Apply eligibility to the frame, before the draw.** If you sample first and
then drop units that turn out to be unusable — no data, no access, no footage —
the realised inclusion probabilities are no longer the ones you assigned, and
every design weight becomes fiction. Filter the frame first; your sample then
represents the filtered frame, and you say so.

## Reference implementation and how to check this one

The method's authors maintain the R package
[BalancedSampling](https://github.com/envisim/BalancedSampling) (AGPL-3), which
is the reference for `lpm1`/`lpm2`. This package is an **independent
implementation written from the published algorithm**, not a translation of
that source, which is why it can be MIT licensed.

`compare_with_r/` contains a harness that checks this implementation against
that one on a shared frame. Run it on a machine with R:

```bash
cd compare_with_r
./run_comparison.sh 5000
```

The two implementations use different random number generators and break ties
differently, so identical samples are impossible — identical output would in
fact be suspicious. What must agree is the *design*, so the harness compares:

1. whether each implementation realises the prescribed inclusion probabilities
   (per-unit z against Monte Carlo standard error);
2. whether the two are consistent with the same design (paired z over all
   units);
3. whether they spread equally well (Kolmogorov–Smirnov on the distribution of
   the spatial balance index).

## Reference

Grafström, A., Lundström, N.L.P. & Schelin, L. (2012). Spatially balanced
sampling through the pivotal method. *Biometrics* 68(2), 514–520.
[doi:10.1111/j.1541-0420.2011.01699.x](https://doi.org/10.1111/j.1541-0420.2011.01699.x)

Related, for the diagnostics and the alternative design:

- Stevens, D.L. Jr. & Olsen, A.R. (2004). Spatially balanced sampling of natural
  resources. *JASA* 99(465), 262–278. — GRTS, and the Voronoi balance measure
  used here.
- Grafström, A. & Lundström, N.L.P. (2013). Why well spread probability samples
  are balanced. *Open Journal of Statistics* 3(1), 36–41.
- Deville, J.-C. & Tillé, Y. (2004). Efficient balanced sampling: the cube
  method. *Biometrika* 91(4), 893–912.

## Performance

Nearest neighbours are found by brute force over the units still in play, which
shrinks by at least one per iteration. A frame of 8,856 units drawing 150 takes
about 0.2 s; 300 units drawing 30 takes about 2 ms. Beyond ~10⁴ units a k-d tree
over the active set would be worth adding — the bottleneck is Python loop
overhead per pivot, not the vectorised distance computation.

## Licence

MIT.
