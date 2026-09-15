# Spectral Recovery

**A controlled Monte Carlo prestudy of correlation-matrix recovery.** It asks a deliberately limited question: for a known one-factor Gaussian population correlation, when does Marchenko–Pastur bulk-edge clipping recover that matrix more accurately than raw sample correlation and a predeclared constant-shrinkage baseline?

This is not a claim about financial returns, real factors, portfolio allocation, machine learning, or an optimal denoising method. The full question, DGP, estimand and limits are fixed in [the research contract](docs/RESEARCH-CONTRACT.md).

## Reproduce

Tested on macOS arm64, Python 3.14.7, NumPy 2.5.2, Matplotlib 3.11.1 and pytest 9.1.1.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt
.venv/bin/python -m pip install -e .
.venv/bin/pytest -q
.venv/bin/python -m spectral_recovery.cli --config configs/primary.json
```

The fixed seed and experiment grid are [here](configs/primary.json). Full replication-level records live in `data/raw/`, summaries in `data/summary.json`, and figures in `figures/`. The three analytical figures are derived from those outputs; the mechanism figure is a deterministic explanatory illustration of the declared DGP.

## Methods in one page

The population matrix is \(\Sigma_\rho=(1-\rho)I+\rho\mathbf{1}\mathbf{1}^T\). Samples provide a noisy correlation \(R\). The MP method replaces all eigenvalues at or below \((1+\sqrt q)^2\) with their bulk mean and reconstructs a unit-diagonal matrix. It is a heuristic tied to the null random-matrix edge, not a universal factor detector. Constant shrinkage uses \(R/(1+q)+qI/(1+q)\), with a coefficient determined only by matrix shape rather than known truth.

The primary metric is off-diagonal Frobenius recovery error. Each reported uncertainty interval is a 95% Monte Carlo interval for the mean error across 300 independent replications in that exact cell.

## Checked-in result, not a general claim

In this exact run, MP clipping had lower mean recovery error than raw correlation in all 12 cells. The improvement ranged from 0.017 to 0.099 off-diagonal error units. The fixed shrinkage comparator helped at \(\rho=0\) and \(0.1\), but was worse than raw in every \(\rho\in\{0.3,0.6\}\) cell. These are descriptive Monte Carlo outcomes from the checked-in seed and DGP, not proof that clipping will improve a real correlation estimate.

## Architecture

```text
known Σρ → seeded Gaussian samples → sample R
                                  ├→ raw R
                                  ├→ MP bulk-edge clipping → unit-diagonal R
                                  └→ design-only constant shrinkage
                                                    ↓
                                truth-known recovery errors → raw JSON → figures
```

## Site case-study block

> **Spectral Recovery** — A controlled synthetic study of when random-matrix bulk-edge clipping improves recovery of a known correlation structure. It compares raw correlation, a transparent clipping heuristic and fixed shrinkage across explicit matrix aspect ratios and factor strengths. **Limit:** results apply only to the stated Gaussian one-factor generator, not to real market or operational data.

## Limits

Measured across 300 random replications in each of 12 configurations: 3 sample-to-dimension ratios × 4 factor correlations. It was **not** varied across factor ranks, sample sizes at fixed \(q\), heavy tails, non-Gaussianity, serial dependence, missingness, real data, selection effects, covariance scales, alternative clipping rules, shrinkage tuning, or downstream decisions. The MP diagnostic is a finite-sample visual check, not a validation theorem.

## References

- Marchenko, V. A. & Pastur, L. A. (1967). *Distribution of eigenvalues for some sets of random matrices.*
- Laloux, L. et al. (1999). *Noise dressing of financial correlation matrices.* Physical Review Letters 83, 1467–1470.
- Bun, J., Bouchaud, J.-P. & Potters, M. (2017). *Cleaning large correlation matrices: tools from Random Matrix Theory.* Physics Reports 666, 1–109.

## Review boundary

The author can run these checks but cannot independently verify the implementation. Stop at a frozen local candidate commit for fresh, read-only review.

**Publication authorized 2026-09-15** (see `docs/decisions/ADR-002-publication-authorization.md`). This repository is public.
