# Primary run results

The committed primary run used seed `20260902`, 80 dimensions, 300 replications in each of 12 predeclared cells (3 values of \(q\) × 4 values of \(\rho\)). It generated 3,600 method-level replication records. The primary error is off-diagonal Frobenius error; each uncertainty interval in `data/summary.json` is a normal 95% Monte Carlo interval for the cell mean.

## Observed pattern

- MP clipping beat raw sample correlation in all 12 cells in paired mean error, by 0.017–0.099.
- Its largest absolute improvement was 0.099 at \(q=.8,\rho=0\); its smallest was 0.017 at \(q=.25,\rho=.6\).
- The design-only shrinker beat raw in the six cells with \(\rho\le .1\), and lost to raw in all six with \(\rho\ge .3\). Its largest loss was 0.205 at \(q=.8,\rho=.6\).
- The separate null-spectrum diagnostic pools 12,800 eigenvalues from 160 \(N=80,T=160\) Gaussian draws. It is shown with, rather than substituted for, the MP density in `figures/mp-null-sanity.png`.

## Interpretation boundary

This supports the narrow statement that this implementation of bulk-mean clipping recovered this known one-factor Gaussian population matrix more closely than the two named baselines in the named grid. It does not establish an estimator ranking for empirical correlation matrices.

## Plain-language model

In these synthetic matrices, the true factor creates one strong shared direction. The MP rule drops the noisy directions into a common bulk level, which happens to help under this generator. The simple shrinker also reduces noise but suppresses strong genuine correlations when \(\rho\) is large. If real data have multiple factors, tails, serial dependence, or a different spectrum, this pattern could change.
