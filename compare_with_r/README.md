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

That builds a shared frame, draws with both implementations, and prints the
comparison. Installing the R package, if needed:

```r
install.packages("BalancedSampling")
```

## What is compared, and why not the samples themselves

The two implementations use different random number generators (numpy's PCG64
versus R's Mersenne-Twister), consume randomness in a different order, and break
nearest-neighbour ties differently. **No seed can align them**, and identical
output would be evidence of a broken harness rather than of correctness.

What must agree is the design each one realises:

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

With 5,000 replicates the Monte Carlo standard error on a unit with `pi = 0.1`
is about 0.004, so differences below roughly 1% of `pi` are not resolvable. Use
more replicates if you want a tighter bound.

## Getting R without root

If the target machine has no R and you cannot use `sudo`, conda-forge ships
both R and the package for linux-64 and aarch64, installable entirely under
`$HOME` with micromamba — no administrator rights at any step.
