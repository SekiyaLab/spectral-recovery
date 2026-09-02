import numpy as np
from spectral_recovery.matrices import as_correlation,constant_shrinkage,mp_clip,population_correlation,sample_correlation

def test_population_correlation_has_unit_diagonal_and_expected_offdiagonal():
    matrix=population_correlation(5,.3); assert np.allclose(np.diag(matrix),1); assert matrix[0,1]==.3; assert np.linalg.eigvalsh(matrix)[0]>0

def test_cleaners_return_symmetric_unit_diagonal_correlations():
    raw=sample_correlation(np.random.default_rng(2).normal(size=(40,8)))
    for cleaned in (mp_clip(raw,.2),constant_shrinkage(raw,.2),as_correlation(raw)):
        assert np.allclose(cleaned,cleaned.T); assert np.allclose(np.diag(cleaned),1)

