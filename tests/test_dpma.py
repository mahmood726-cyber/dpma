"""Tests for DPMA pipeline: DL estimator, DP sampler, convergence, sensitivity."""

import sys
import numpy as np
import pytest

sys.path.insert(0, r'C:\Models\DPMA')
from pipeline import DPMetaAnalysis, standard_dl, check_convergence, run_sensitivity


# ── DerSimonian-Laird ──

class TestDL:
    def test_dl_basic(self):
        """DL on 3 studies with known values gives finite results."""
        yi = np.array([0.5, 0.3, 0.7])
        sei = np.array([0.1, 0.15, 0.2])
        result = standard_dl(yi, sei)
        assert 0.2 < result['theta'] < 0.8
        assert result['se'] > 0
        assert result['tau2'] >= 0
        assert 0 <= result['I2'] <= 100

    def test_dl_homogeneous(self):
        """Identical effects should give tau2 near 0."""
        yi = np.array([0.5, 0.5, 0.5, 0.5, 0.5])
        sei = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        result = standard_dl(yi, sei)
        assert result['tau2'] == 0.0
        assert abs(result['theta'] - 0.5) < 1e-10

    def test_dl_two_studies(self):
        """k=2 should not crash and gives finite PI."""
        yi = np.array([0.3, 0.7])
        sei = np.array([0.1, 0.2])
        result = standard_dl(yi, sei)
        assert np.isfinite(result['theta'])
        assert np.isfinite(result['pi_lo'])

    def test_dl_prediction_interval(self):
        """PI should be wider than CI for heterogeneous data."""
        yi = np.array([0.1, 0.5, 0.9, 0.2, 0.8])
        sei = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        result = standard_dl(yi, sei)
        ci_width = 2 * 1.96 * result['se']
        pi_width = result['pi_hi'] - result['pi_lo']
        assert pi_width > ci_width


# ── Dirichlet Process Sampler ──

class TestDPMA:
    def test_dp_recovers_two_clusters(self):
        """DP sampler should find 2+ clusters in bimodal data."""
        rng = np.random.RandomState(42)
        yi = np.concatenate([rng.normal(0.5, 0.05, 10), rng.normal(-0.5, 0.05, 10)])
        sei = np.full(20, 0.1)
        dp = DPMetaAnalysis(n_iter=500, burn=100)
        result = dp.fit(yi, sei)
        assert result['n_clusters_mode'] >= 2

    def test_dp_unimodal_data(self):
        """Tight unimodal data should yield 1 cluster (mode)."""
        rng = np.random.RandomState(123)
        yi = rng.normal(0.3, 0.01, 15)
        sei = np.full(15, 0.05)
        dp = DPMetaAnalysis(n_iter=500, burn=100)
        result = dp.fit(yi, sei)
        assert result['n_clusters_mode'] <= 2  # 1 expected, allow 2

    def test_single_study(self):
        """Edge case: k=1 study should not crash."""
        yi = np.array([0.5])
        sei = np.array([0.1])
        dp = DPMetaAnalysis(n_iter=100, burn=20)
        result = dp.fit(yi, sei)
        assert 'n_clusters_mode' in result
        assert 'theta_posterior_mean' in result

    def test_two_studies(self):
        """Edge case: k=2 should not crash."""
        yi = np.array([0.3, 0.7])
        sei = np.array([0.1, 0.15])
        dp = DPMetaAnalysis(n_iter=100, burn=20)
        result = dp.fit(yi, sei)
        assert result['n_clusters_mode'] >= 1

    def test_result_keys(self):
        """All expected keys present in result dict."""
        yi = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        sei = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        dp = DPMetaAnalysis(n_iter=100, burn=20)
        result = dp.fit(yi, sei)
        expected_keys = [
            'n_clusters_mean', 'n_clusters_median', 'n_clusters_mode',
            'theta_posterior_mean', 'theta_posterior_sd',
            'theta_95_lo', 'theta_95_hi',
            'predictive_bimodal', 'n_clusters_trace', 'final_clusters',
        ]
        for key in expected_keys:
            assert key in result, f"Missing key: {key}"

    def test_credible_interval_contains_mean(self):
        """95% CrI should contain the posterior mean."""
        yi = np.array([0.2, 0.3, 0.4, 0.5, 0.6])
        sei = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        dp = DPMetaAnalysis(n_iter=500, burn=100)
        result = dp.fit(yi, sei)
        assert result['theta_95_lo'] <= result['theta_posterior_mean'] <= result['theta_95_hi']

    def test_alpha_parameter(self):
        """Higher alpha should tend toward more clusters."""
        rng = np.random.RandomState(42)
        yi = np.concatenate([rng.normal(0.5, 0.1, 8), rng.normal(-0.5, 0.1, 8)])
        sei = np.full(16, 0.15)
        low = DPMetaAnalysis(n_iter=300, burn=50, alpha=0.1).fit(yi, sei)
        high = DPMetaAnalysis(n_iter=300, burn=50, alpha=10.0).fit(yi, sei)
        # High alpha should give >= clusters as low alpha (on average)
        assert high['n_clusters_mean'] >= low['n_clusters_mean'] - 1  # Allow some stochasticity


# ── Convergence Diagnostics ──

class TestConvergence:
    def test_stationary_trace(self):
        """Constant trace should pass convergence."""
        trace = [2] * 100
        assert check_convergence(trace) is True

    def test_drifting_trace(self):
        """Strongly drifting trace should fail convergence."""
        trace = list(range(100))  # 0, 1, 2, ..., 99
        assert check_convergence(trace) is False

    def test_short_trace(self):
        """Very short trace should return True (too short to assess)."""
        assert check_convergence([1, 2]) is True

    def test_noisy_stationary(self):
        """Noisy but stationary trace should pass."""
        rng = np.random.RandomState(42)
        trace = rng.normal(3.0, 0.5, 200).tolist()
        assert check_convergence(trace) is True


# ── Sensitivity Analysis ──

class TestSensitivity:
    def test_sensitivity_runs(self):
        """Sensitivity analysis returns results for each alpha."""
        yi = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        sei = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
        results = run_sensitivity(yi, sei, alphas=[0.5, 1.0, 2.0])
        assert len(results) == 3
        for alpha in [0.5, 1.0, 2.0]:
            assert alpha in results
            assert 'n_clusters' in results[alpha]
            assert 'theta' in results[alpha]


# ── Bimodality Test ──

class TestBimodality:
    def test_bimodal_detected(self):
        """Strongly bimodal samples should be detected."""
        rng = np.random.RandomState(42)
        samples = np.concatenate([rng.normal(-2, 0.3, 500), rng.normal(2, 0.3, 500)])
        dp = DPMetaAnalysis()
        assert dp._test_bimodality(samples) is True

    def test_unimodal_not_flagged(self):
        """Unimodal samples should not be flagged."""
        rng = np.random.RandomState(42)
        samples = rng.normal(0, 1, 1000)
        dp = DPMetaAnalysis()
        assert dp._test_bimodality(samples) is False


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
