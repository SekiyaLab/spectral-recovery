from __future__ import annotations

import numpy as np


def population_correlation(n: int, rho: float) -> np.ndarray:
    if not 0 <= rho < 1: raise ValueError("rho must be in [0, 1)")
    return (1.0 - rho) * np.eye(n) + rho * np.ones((n, n))


def sample_correlation(samples: np.ndarray) -> np.ndarray:
    centered = samples - samples.mean(axis=0, keepdims=True)
    scaled = centered / centered.std(axis=0, ddof=1, keepdims=True)
    return (scaled.T @ scaled) / (len(samples) - 1)


def as_correlation(matrix: np.ndarray) -> np.ndarray:
    scales = np.sqrt(np.clip(np.diag(matrix), 1e-12, None))
    result = matrix / np.outer(scales, scales)
    return (result + result.T) / 2


def mp_clip(correlation: np.ndarray, q: float) -> np.ndarray:
    """Bulk-mean MP clipping, followed by diagonal re-normalization."""
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    edge = (1.0 + np.sqrt(q)) ** 2
    bulk = eigenvalues <= edge
    replacement = eigenvalues[bulk].mean() if bulk.any() else eigenvalues.mean()
    cleaned = np.where(bulk, replacement, eigenvalues)
    return as_correlation((eigenvectors * cleaned) @ eigenvectors.T)


def constant_shrinkage(correlation: np.ndarray, q: float) -> np.ndarray:
    alpha = 1.0 / (1.0 + q)
    return alpha * correlation + (1.0 - alpha) * np.eye(len(correlation))


def offdiagonal_error(estimate: np.ndarray, truth: np.ndarray) -> float:
    n = len(estimate); difference = estimate - truth
    return float(np.linalg.norm(difference - np.diag(np.diag(difference)), ord="fro") / np.sqrt(n * (n - 1)))

