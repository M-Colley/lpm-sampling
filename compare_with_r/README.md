# Comparing against the reference R implementation

The authors of the method maintain
[BalancedSampling](https://github.com/envisim/BalancedSampling) in R. This
directory checks that this Python package realises the same *design*.

Nothing from that package is included or redistributed here. You install it
yourself; these scripts only run it and compare summary statistics. This
directory is not part of the built wheel.

## Run it

On a machine with R:

```bash
cd compare_with_r
./run_comparison.sh 5000        # replicates per implementation per variant
```

That builds a shared frame, runs the deterministic checks, draws with both
implementations, and prints the comparison. It installs the R package itself if
it is missing:

```r
install.packages("BalancedSampling")
```

The script finds the interpreters itself, including R's default Windows install
location. Override either if they live somewhere unusual:

```bash
PYTHON=/path/to/python RSCRIPT=/path/to/Rscript ./run_comparison.sh 5000
```

## What is compared, and why not the samples themselves

The two implementations use different random number generators (numpy's PCG64
versus R's Mersenne-Twister), consume randomness in a different order, and break
nearest-neighbour ties differently. **No seed can align them**, and identical
output would be evidence of a broken harness rather than of correctness.

Two things *can* be compared exactly, and are checked first, because they are
deterministic functions of their inputs rather than draws:

0. **The shared functions.** `spatial_balance` against `sb`, over 200 fixed
   samples ranging from well spread to deliberately clumped; and `pi_from_size`
   against `getPips`, over 12 size distributions whose skew forces the
   certainty-unit capping iteration to run zero, one and many rounds. Both must
   agree to floating point. The first is a precondition for test 3 below: that
   test compares the two samplers *using spatial balance as the yardstick*,
   which is only fair if both sides agree on what the yardstick measures.
   Samples cross the boundary as 0/1 masks, so the 0-based / 1-based index
   convention gap cannot silently shift one.

Beyond those, what must agree is the design each one realises:

1. **First-order inclusion probabilities.** Each implementation's empirical
   selection frequency per unit, against the prescribed `pi`, in units of Monte
   Carlo standard error. This is the property the whole method rests on: if it
   drifts, design weights are fiction. It is also the sensitive test — a
   plausible transcription slip in the pair update leaves sample size and total
   mass correct while silently biasing every unit, and only this test sees it.
2. **Cross-implementation agreement.** Paired comparison of the two frequency
   vectors, to check they are consistent with the *same* underlying design
   rather than merely each being self-consistent.
3. **Spread.** The distribution of the Voronoi spatial balance index across
   replicates, compared with a two-sample Kolmogorov–Smirnov test. Equal means
   with a large p-value is the expected result; a large mean difference with a
   tiny p-value would be the real finding — the implementations would then
   differ in how they spread, not merely in their RNG.
4. **Spatial dependence between units.** Moran's I of each implementation's
   per-unit deviation from `pi`, under 8-nearest-neighbour weights. Tests 1-3
   look at marginal frequencies and at the balance distribution; neither looks
   at how selections are correlated *between* neighbours, which is the defining
   property of the method. Both implementations should come out negative —
   neighbours compete for the same mass, so a unit over-selected in a finite run
   leaves its neighbours under-selected — and, more to the point, they should
   agree with each other. Note that testing I against an exchangeable
   permutation null is the wrong test here: it rejects for both implementations,
   because the correlation is a property of the design rather than a defect.

With 5,000 replicates the Monte Carlo standard error on a unit with `pi = 0.1`
is about 0.004, so differences below roughly 1% of `pi` are not resolvable. Use
more replicates if you want a tighter bound.

## Getting R without root

If the target machine has no R and you cannot use `sudo`, conda-forge ships
both R and the package for linux-64 and aarch64, installable entirely under
`$HOME` with micromamba — no administrator rights at any step.
