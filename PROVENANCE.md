# Provenance and licensing

This package implements the **Local Pivotal Method** as published in:

> Grafström, A., Lundström, N.L.P. & Schelin, L. (2012). Spatially balanced
> sampling through the pivotal method. *Biometrics* 68(2), 514–520.
> doi:10.1111/j.1541-0420.2011.01699.x

## Why this can be MIT licensed

The method's authors maintain the R package
[BalancedSampling](https://github.com/envisim/BalancedSampling), which is
licensed **AGPL-3**. That licence binds code copied from or derived from *that
codebase*. It does not, and cannot, restrict independent implementations of the
published method: copyright protects the expression of a program, not the
algorithm or the mathematics it implements.

This package is an **independent implementation written from the published
description of the algorithm**. It is not a translation or port of the R or C++
source. Specifically:

- No source code was copied, translated, or paraphrased from BalancedSampling.
- No test fixtures, saved seeds, expected-output arrays, or example datasets
  were taken from it. The tests here are *property* tests — Monte Carlo checks
  that inclusion probabilities are realised, that sample size is exact, that
  certainty units are respected — rather than numeric arrays traceable to
  another package.
- No documentation text, argument descriptions, or vignette prose was copied.
- The API is deliberately *not* a mirror of theirs. Where a name is shared it
  is because the paper itself uses it: LPM1 and LPM2 are the authors' names for
  the two variants. Their `getPips` and `sb` are `pi_from_size` and
  `spatial_balance` here.
- Nothing from the R package is vendored, bundled, or shipped. The comparison
  harness in `compare_with_r/` runs against an installation the *user* provides
  on their own machine, and is excluded from the built wheel.

The Biometrics paper itself is under Wiley's copyright: it is cited, never
quoted at length or reproduced.

## If you extend this package

Keep the separation intact. Read the paper, not their source. Checking your
results against BalancedSampling numerically — which is exactly what
`compare_with_r/` is for — is fine and is good practice; adjusting your code to
match their *implementation* after reading it is not.
