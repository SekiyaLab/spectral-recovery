# Forensic review — Spectral Recovery, Track A

- **Reviewer / model:** Claude (Sonnet 5, `claude-sonnet-5`), acting as an independent fresh reviewer with no prior context on this study.
- **Candidate commit:** `14fb9bd95d4cb1080bb5d3fc3f4a27544f8799b3` ("Build Spectral Recovery controlled prestudy"), the sole commit in this repository.
- **Date of review:** 2026-09-02.
- **Scope:** Review only. No implementation, DGP, results, figures, dependency, or parent Portfolio/site/Sekiya files were modified. This document, committed separately as a documentation-only commit, is the only change made by this review.

## What was checked and how

### 1. DGP and known target correctness
Read `src/spectral_recovery/matrices.py::population_correlation` and `experiment.py::run`. Verified analytically:
- \(\Sigma_\rho=(1-\rho)I+\rho\mathbf1\mathbf1^T\) has unit diagonal for all \(\rho\), and eigenvalues \(1-\rho\) (multiplicity \(N-1\)) and \(1-\rho+\rho N\) (multiplicity 1); both positive for \(0\le\rho<1\), so it is a valid PD correlation matrix and Cholesky-factorizable. Code enforces `0 <= rho < 1` with a `ValueError` otherwise.
- Sample construction `rng.normal(size=(T,N)) @ chol.T` gives rows with covariance \(LL^T=\Sigma_\rho\) (standard Cholesky-transform sampling) — confirmed by hand.
- Independently verified `sample_correlation()` against `np.corrcoef` on a random matrix: **exact match** (`np.allclose` True).

### 2. Method descriptions vs. implementation
Cross-read `docs/RESEARCH-CONTRACT.md`, `docs/decisions/ADR-001-*.md`, `README.md` against `matrices.py`:
- **Raw**: unregularized sample correlation — matches.
- **MP bulk-mean clipping**: contract says "retain eigenvalues above \(\lambda_+=(1+\sqrt q)^2\), replace all remaining bulk eigenvalues (at or below) by their bulk mean, reconstruct, rescale to unit diagonal." Code (`mp_clip`) does exactly this: `bulk = eigenvalues <= edge`, mean-replace, reconstruct via eigenvectors, `as_correlation()` rescale. Matches contract and ADR-001 verbatim.
- **Constant shrinkage**: contract states \(\hat R=aR+(1-a)I\), \(a=1/(1+q)\); README states the algebraically identical form \(R/(1+q)+qI/(1+q)\). Confirmed these are the same expression (\(1-a=q/(1+q)\)). Code matches.
- **`offdiagonal_error`**: contract defines \(\|\operatorname{offdiag}(\hat R-\Sigma)\|_F/\sqrt{N(N-1)}\). Verified formula against a hand-computed perturbation case: matched to floating-point precision (`0.0632455532033676` both ways).

### 3. q convention and MP edge
`q` is used throughout as \(N/T\) (dimensions/observations), matching the standard Marchenko–Pastur aspect ratio \(c=p/n\) for a sample covariance built from `T` iid `N`-dim rows. The diagnostic (`mp_diagnostic`) independently builds `t = round(n/q)` at `q=.5` and computes `x.T@x/t`, using the same convention — consistent with `mp_clip`'s edge formula. The rendered upper edge in `figures/mp-null-sanity.png` (2.91) matches \((1+\sqrt{.5})^2=2.9142\) and visually tracks the histogram well.
- Confirmed mathematically that the `mp_clip` fallback branch (`bulk.any()` false, i.e., every eigenvalue of an 80×80 unit-diagonal matrix exceeding the edge) is **impossible** in this grid: trace is fixed at 80, and even the smallest edge (q=0.25) is 2.25, so at most ⌊80/2.25⌋=35 eigenvalues could exceed it. This is a dead branch, not a bug — see finding SR-003.

### 4. Metric calculation and pairing
`summarize()` filters rows by `(q, rho, method)` and computes the paired difference vs. raw by index-position subtraction. Because `run()` appends exactly one row per method per replication in a fixed order (`raw`, `mp_clip`, `constant_shrinkage`) inside the replication loop, the filtered per-method arrays stay in matching replication order, so index-paired subtraction is correct. This is implicit rather than explicit (no replication ID join key), which is fragile to future refactors but correct as written.

### 5. Reproducibility, seeds, and provenance
- Ran `.venv/bin/pytest -q` in place (no file writes outside `tmp_path`): **3 passed**.
- Independently re-ran the **entire pipeline** in an isolated `mktemp -d` copy (copied only `src/`, `configs/`, `pyproject.toml` — not the committed data/figures) using the same locally pinned `.venv` (numpy 2.5.2, matplotlib 3.11.1, pytest 9.1.1, matching `requirements-lock.txt`/`pyproject.toml` exactly) and the exact command from the README:
  ```
  .venv/bin/python -m spectral_recovery.cli --config configs/primary.json
  ```
  Result: `data/summary.json`, `data/raw/primary-results.json` (10,800 rows = 4 ρ × 3 q × 300 reps × 3 methods, confirmed by count), and `data/raw/mp-null-eigenvalues.json` were **byte-for-byte identical** to the checked-in versions (`diff` clean on all three). This confirms the seed (`20260902`) and pipeline are genuinely, exactly reproducible as claimed, not just "should be."
- Confirmed via `git ls-files` that only the 22 intended files are tracked in the commit; no `.venv/`, `__pycache__/`, or `.pytest_cache/` artifacts leaked in.

### 6. Numeric claims in README/RESULTS vs. actual artifacts
Recomputed directly from `data/summary.json` (not from the prose) using a standalone script:
- "MP clipping beat raw sample correlation in all 12 cells... by 0.017–0.099": confirmed — all 12 `(raw − mp_clip)` mean-error differences are positive, ranging from 0.016785 (q=.25, ρ=.6) to 0.099119 (q=.8, ρ=0), rounding exactly to the stated 0.017–0.099.
- "shrinker beat raw in the six cells with ρ≤.1, lost to raw in all six with ρ≥.3": confirmed exactly against all 12 `(raw − constant_shrinkage)` values.
- "largest loss was 0.205 at q=.8, ρ=.6": confirmed, 0.204876 rounds to 0.205.
- RESULTS.md pooled-eigenvalue count (12,800 from 160×80): confirmed via config (`mp_diagnostic_replications=160`, `dimensions=80`).
No numeric claim in README.md or docs/RESULTS.md was found to be unsupported, cherry-picked, or rounded in the favorable direction.

### 7. Figures and uncertainty presentation
Visually inspected all four committed PNGs (`mechanism.png`, `recovery-curves.png`, `improvement-heatmap.png`, `mp-null-sanity.png`). All render cleanly, axes/legends are labeled, the recovery-curve error bars correctly reflect the 95% CI half-widths from `data/summary.json`, and the heatmap's printed cell values match the `paired_difference_vs_raw` field exactly (spot-checked several cells against the raw JSON). The MP histogram visually tracks the analytic MP density well and the edge line is correctly placed. One presentation concern is flagged below (SR-001).

### 8. Test quality
`tests/test_matrices.py` and `tests/test_experiment.py` — 3 tests total. They check: exact run-to-run determinism and record counts; population matrix unit-diagonal/off-diagonal/PD sanity; and symmetric/unit-diagonal invariants of the two cleaning methods. This is a thin but honest suite for a declared "prestudy" — see SR-004 for gaps.

### 9. Hidden sophistication / overclaim check
No hidden complexity was found beyond what's documented — no undisclosed hyperparameter search, no post-hoc method selection, no oracle information leaking into the constant-shrinkage coefficient (`a` is a pure function of `q`, never of `ρ` or of any error outcome — confirmed by reading `constant_shrinkage()` and ADR-001's explicit rejection of oracle shrinkage). The prose is unusually careful to under-claim rather than over-claim (repeated "not a claim about...", explicit non-varied-dimensions list, "descriptive Monte Carlo outcomes... not proof"). No numbers were found to be fabricated, rounded favorably, or omitted to hide an inconvenient cell — the constant-shrinkage losses at high ρ are reported prominently, including the largest one.

### 10. Limitations and public-profile readiness
The limits sections (README "Limits", RESEARCH-CONTRACT "Scope and limitations") are specific and match what was actually varied (3 q values × 4 ρ values × 300 reps, Gaussian one-factor only). The "Review boundary" and "Site case-study block" sections in README correctly and conservatively scope what can be said publicly.

## Findings

All findings below are **non-blocking**. No correctness, reproducibility, or overclaim defect was found.

**SR-001** (non-blocking, presentation) — `src/spectral_recovery/plotting.py:11`: the third panel of `mechanism.png` ("Spectral idea: retain structure") is generated via an ad hoc elementwise threshold on the correlation matrix (`np.where(np.abs(raw)>.18, raw, 0)`), not via the actual eigenvalue-space `mp_clip()` used everywhere else in the study. README's claim that "the mechanism figure is a deterministic explanatory illustration of the declared DGP" is true of panels 1–2 (population → sampled correlation) but panel 3 is a visual metaphor for the *result* of spectral clipping, not the algorithm itself, and isn't captioned as such. A reader skimming only the figure could reasonably infer the method thresholds individual correlation entries rather than eigenvalues. Recommend either regenerating panel 3 from real `mp_clip()` output or adding an explicit "illustrative, not the literal algorithm" caption.

**SR-002** (non-blocking, data hygiene) — `src/spectral_recovery/experiment.py:19`: every one of the 10,800 rows in `data/raw/primary-results.json` carries a `top_eigenvalue` field that is never read by `summarize()`, any plotting function, or any doc. Either use it (e.g., as a factor-separation diagnostic) or drop it — as-is it's unexplained bloat in the primary data artifact.

**SR-003** (non-blocking, code clarity) — `src/spectral_recovery/matrices.py:28`: `mp_clip`'s `eigenvalues.mean()` fallback (triggered only if literally every eigenvalue exceeds the MP edge) is mathematically unreachable for this study's grid (trace fixed at N=80 makes it impossible for all 80 eigenvalues to exceed any of the three declared edges — verified by direct calculation). It's dead, untested code in the shipped artifact and isn't mentioned in the contract or ADR. Not a bug; worth a one-line comment or removal.

**SR-004** (non-blocking, test coverage) — Test suite (3 tests) does not directly assert: `offdiagonal_error`'s formula/denominator, the `ValueError` boundary behavior of `population_correlation` for invalid ρ, or determinism of `mp_diagnostic` (only `run()`'s determinism is tested). Acceptable for a "prestudy" given the honest scope framing, but worth naming explicitly as a known gap rather than leaving it implicit.

## Blind spots

- This review did not byte-diff regenerated PNG figures against the committed PNGs (only the underlying JSON data feeding the figures was byte-diffed and confirmed identical; PNG encoding can carry non-deterministic metadata even from identical pixel data). Figures were checked by visual inspection and by spot-checking printed numeric labels against `data/summary.json` instead.
- Did not audit the supply-chain integrity of the installed `numpy`/`matplotlib`/`pillow` wheels (e.g., hash verification against PyPI) — relied on the pre-existing local `.venv` that matches `requirements-lock.txt`.
- Did not test on any platform other than the one available (macOS arm64, Python 3.14.7) — cross-platform floating-point/BLAS determinism was not independently re-verified, though this is inherent to any NumPy-based study and not specific to this candidate.
- Did not verify the two external non-MP citations (Laloux et al. 1999; Bun et al. 2017) bibliographically — they are standard, well-known references in this literature and their content wasn't relied on for any of this review's checks.
- This review is scoped strictly to the DGP/methods/results/figures/tests of the study itself, not to how it will be integrated into the parent Portfolio site, which is explicitly out of scope per the task boundary.

## Verdict

**ACCEPT.**

The candidate is **public-ready pending separate publication authorization**. All numeric claims in README.md and docs/RESULTS.md were independently reproduced from the checked-in artifacts and from a from-scratch, isolated re-run of the entire pipeline (identical seed, identical outputs, byte-for-byte). The DGP, all three methods, the metric, and the MP/q convention were all verified against their own documentation and against independent numerical checks (e.g., matching `np.corrcoef`, hand-computed error formula). No overclaim, hidden complexity, or fabricated/rounded-favorable number was found. The four non-blocking findings above (SR-001 through SR-004) are polish items for a future revision, not blockers to accepting this commit or to eventual public exposure of this exact content.

## Next boundary

This review covers only commit `14fb9bd95d4cb1080bb5d3fc3f4a27544f8799b3` as a frozen, local, read-only artifact. It authorizes nothing beyond itself: no push to any remote (none configured), no edit to parent Portfolio/site/Sekiya, and no public-publication action. Any such action requires separate, explicit authorization from the repository owner. The only change made during this review is this document, committed by itself as a documentation-only commit.
