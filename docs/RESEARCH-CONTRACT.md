# Research contract — Spectral Recovery

## Question

For a known one-factor population correlation matrix, across predeclared sample-to-dimension ratios \(q=N/T\) and factor correlations \(\rho\), when does Marchenko–Pastur (MP) bulk-edge eigenvalue clipping improve recovery of the population correlation relative to raw sample correlation and a fixed, design-only constant-shrinkage baseline?

## Controlled DGP

For each replication, draw \(T\) independent observations from \(N(0,\Sigma_\rho)\), where

\[
\Sigma_\rho=(1-\rho)I_N+\rho\mathbf{1}\mathbf{1}^{\mathsf T}.
\]

This is already a correlation matrix for \(0\le\rho<1\). Sample columns are centred and scaled, then the sample correlation \(R\) is calculated. The grid fixes \(N=80\), uses \(q\in\{0.25,0.5,0.8\}\) with \(T=\mathrm{round}(N/q)\), and \(\rho\in\{0,0.1,0.3,0.6\}\); there are 300 seeded replications per cell.

## Methods and estimand

The recovery target is the known \(\Sigma_\rho\). The primary metric is off-diagonal Frobenius error, \(\|\operatorname{offdiag}(\hat R-\Sigma_\rho)\|_F/\sqrt{N(N-1)}\), which avoids double-counting the forced unit diagonal.

1. **Raw:** the unregularized sample correlation.
2. **MP bulk-edge clipping:** diagonalize \(R\), retain eigenvalues above \(\lambda_+=(1+\sqrt{q})^2\), and replace all remaining bulk eigenvalues by their bulk mean. Reconstruct and rescale to a unit diagonal. This uses the MP null edge as a heuristic signal/bulk separator; it does not estimate a universal latent-factor model.
3. **Constant shrinkage:** \(\hat R_{\rm shrink}=aR+(1-a)I\), where \(a=1/(1+q)\). The coefficient is a predeclared function of the observed matrix shape only; it is never fitted using \(\rho\), held-out truth, or simulation results.

For each method/cell, report the Monte Carlo mean error, standard deviation, and normal 95% interval for the mean error. The paired mean error difference versus raw is reported for decision relevance. A seeded null diagnostic independently compares pooled noise-only sample-covariance eigenvalues against the MP density at \(q=.5\).

## Hypotheses

- H1: Raw sample correlation has higher error as \(q\) rises, especially for weak signals.
- H2: MP clipping improves mean recovery in at least some moderate/high-\(q\), nonzero-factor cells; it need not win under pure noise or very strong factors.
- H3: Design-only constant shrinkage is a meaningful honest comparator and can outperform clipping in some cells.

## Scope and limitations

This is a synthetic one-factor Gaussian correlation experiment. It varies 3 ratios, 4 factor correlations and 300 random replications per cell. It does **not** vary factor rank, non-Gaussianity, heavy tails, missing data, time dependence, feature-selection, covariance scaling, unknown population structure, real datasets, or any downstream portfolio/trading decision. The MP histogram is a finite-sample sanity diagnostic, not a goodness-of-fit proof. Results establish only what happened under this generator and these implementations.

## Plain-language model

When there are not many observations per variable, sample correlations contain a noisy cloud of eigenvalues. MP clipping treats the cloud as bulk noise and preserves unusually large directions; constant shrinkage simply pulls every correlation toward zero by a fixed amount. The known population matrix lets this study check which reconstruction is closer. If the real data are not like the one-factor Gaussian generator, the apparent advantage may disappear or reverse. The study does not resolve which method is best for any real dataset.

