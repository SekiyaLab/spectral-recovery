# ADR-001: Use MP bulk-mean clipping and design-only constant shrinkage

## Status

Accepted — 2026-09-02

## Context

The study requires a spectral method and a comparator that does not use knowledge of the population correlation. Choosing a shrinkage intensity after inspecting simulated recovery would invalidate the comparison.

## Decision

Use the MP upper bulk edge \((1+\sqrt q)^2\). Replace eigenvalues at or below it with their empirical bulk mean, preserve eigenvalues above it, then renormalize the reconstructed matrix to a correlation matrix. Compare with raw sample correlation and \(aR+(1-a)I\), \(a=1/(1+q)\).

## Alternatives considered

- Oracle shrinkage: rejected because it tunes against known truth.
- Ledoit–Wolf/OAS: useful but a different covariance-estimation study with additional implementation and target choices.
- Eigenvalue deletion: rejected because it destroys trace information and is less defensible as a correlation reconstruction.

## Consequences

MP clipping is evaluated as a declared heuristic, not claimed as an optimal estimator. The constant shrinker is intentionally simple and honest, not presented as a state-of-the-art optimum.

