"""Dirichlet Process Meta-Analysis (DPMA): Bayesian Nonparametric Evidence Synthesis.

Standard random-effects: yi ~ N(theta, tau2 + sei2)  → ONE subgroup, normal distribution
DPMA: yi ~ DP-mixture of Gaussians → AUTOMATIC discovery of subgroups

For each of 307 Cochrane reviews:
  1. Fit standard DL random-effects (1 component)
  2. Fit DPMA via Gibbs sampling (discovers K components automatically)
  3. Compare: does DPMA find subgroups? How many? How does the pooled estimate change?
  4. Compute predictive distribution (full posterior predictive, not just PI)

Key innovation: the number of subgroups is LEARNED from data, not pre-specified.
"""

import csv
import json
import math
import time
import numpy as np
import pyreadr
from pathlib import Path
from scipy import stats as sp_stats
from collections import Counter

PAIRWISE_DIR = Path(r'C:\Models\Pairwise70\data')
OUTPUT_DIR = Path(r'C:\Models\DPMA\data\output')


# ═══════════════════════════════════════════════════
# DIRICHLET PROCESS MIXTURE MODEL (Gibbs Sampler)
# ═══════════════════════════════════════════════════

class DPMetaAnalysis:
    """Dirichlet Process Gaussian Mixture for meta-analysis.

    Model:
      yi | mu_i, sei ~ N(mu_i, sei^2)       [known sampling variance]
      mu_i | G ~ G                            [true effect for study i]
      G ~ DP(alpha, G0)                       [Dirichlet Process prior]
      G0 = N(mu0, sigma0^2)                   [base distribution]

    Gibbs sampler via Chinese Restaurant Process (CRP) representation.
    """

    def __init__(self, alpha=1.0, mu0=0.0, sigma0=1.0, n_iter=500, burn=200):
        self.alpha = alpha
        self.mu0 = mu0
        self.sigma0 = sigma0
        self.n_iter = n_iter
        self.burn = burn

    def fit(self, yi, sei):
        k = len(yi)
        vi = sei**2

        # Initialize: all in one cluster
        labels = np.zeros(k, dtype=int)
        cluster_means = {0: np.mean(yi)}

        # Storage
        n_clusters_trace = []
        theta_samples = []  # posterior predictive samples

        rng = np.random.RandomState(42)

        for iteration in range(self.n_iter):
            # For each study, resample its cluster assignment (CRP)
            for i in range(k):
                old_label = labels[i]

                # Remove study i from its cluster
                members = [j for j in range(k) if j != i and labels[j] == old_label]
                if len(members) == 0:
                    del cluster_means[old_label]

                # Compute CRP probabilities for existing clusters
                probs = []
                cluster_ids = []
                for c, mu_c in cluster_means.items():
                    n_c = np.sum(labels == c) - (1 if labels[i] == c else 0)
                    if n_c <= 0:
                        continue
                    # Likelihood: N(yi | mu_c, vi)
                    ll = sp_stats.norm.logpdf(yi[i], loc=mu_c, scale=math.sqrt(vi[i]))
                    probs.append(math.log(n_c) + ll)
                    cluster_ids.append(c)

                # Probability of new cluster
                # Marginal likelihood under base distribution G0
                sigma_new = math.sqrt(vi[i] + self.sigma0**2)
                ll_new = sp_stats.norm.logpdf(yi[i], loc=self.mu0, scale=sigma_new)
                probs.append(math.log(self.alpha) + ll_new)
                cluster_ids.append(-1)  # new cluster

                # Normalize (log-sum-exp)
                max_p = max(probs)
                probs_exp = [math.exp(p - max_p) for p in probs]
                total = sum(probs_exp)
                probs_norm = [p / total for p in probs_exp]

                # Sample
                choice = rng.choice(len(cluster_ids), p=probs_norm)
                chosen = cluster_ids[choice]

                if chosen == -1:
                    # New cluster: sample mu from posterior given yi
                    prec_prior = 1.0 / self.sigma0**2
                    prec_data = 1.0 / vi[i]
                    prec_post = prec_prior + prec_data
                    mu_post = (prec_prior * self.mu0 + prec_data * yi[i]) / prec_post
                    sigma_post = math.sqrt(1.0 / prec_post)
                    new_label = max(cluster_means.keys(), default=-1) + 1
                    cluster_means[new_label] = rng.normal(mu_post, sigma_post)
                    labels[i] = new_label
                else:
                    labels[i] = chosen
                    # Re-add to existing cluster
                    if chosen not in cluster_means:
                        cluster_means[chosen] = yi[i]

            # Update cluster means given current assignments
            for c in list(cluster_means.keys()):
                members = np.where(labels == c)[0]
                if len(members) == 0:
                    del cluster_means[c]
                    continue

                # Posterior: conjugate normal update
                prec_prior = 1.0 / self.sigma0**2
                prec_data = np.sum(1.0 / vi[members])
                prec_post = prec_prior + prec_data
                mu_post = (prec_prior * self.mu0 + np.sum(yi[members] / vi[members])) / prec_post
                sigma_post = math.sqrt(1.0 / prec_post)
                cluster_means[c] = rng.normal(mu_post, sigma_post)

            # Record after burn-in
            if iteration >= self.burn:
                n_clusters_trace.append(len(cluster_means))

                # Posterior predictive sample: pick a cluster from CRP, then sample from it
                counts = Counter(labels)
                probs_pred = []
                means_pred = []
                for c, mu_c in cluster_means.items():
                    probs_pred.append(counts.get(c, 0))
                    means_pred.append(mu_c)
                probs_pred.append(self.alpha)
                means_pred.append(rng.normal(self.mu0, self.sigma0))

                probs_pred = np.array(probs_pred, dtype=float)
                probs_pred /= probs_pred.sum()
                chosen_idx = rng.choice(len(means_pred), p=probs_pred)
                # Effect-level predictive: sample from posterior mixture (no observation noise)
                theta_samples.append(means_pred[chosen_idx])

        # Summarize
        n_clusters_post = np.array(n_clusters_trace)
        theta_post = np.array(theta_samples)

        # Final cluster assignments (last iteration)
        final_labels = labels.copy()
        final_clusters = {}
        for c in set(final_labels):
            members = np.where(final_labels == c)[0]
            final_clusters[c] = {
                'k': len(members),
                'mean': float(cluster_means.get(c, 0)),
                'studies': members.tolist(),
            }

        return {
            'n_clusters_mean': float(np.mean(n_clusters_post)),
            'n_clusters_median': int(np.median(n_clusters_post)),
            'n_clusters_mode': int(Counter(n_clusters_post.tolist()).most_common(1)[0][0]),
            'theta_posterior_mean': float(np.mean(theta_post)),
            'theta_posterior_sd': float(np.std(theta_post)),
            'theta_95_lo': float(np.percentile(theta_post, 2.5)),
            'theta_95_hi': float(np.percentile(theta_post, 97.5)),
            'predictive_bimodal': self._test_bimodality(theta_post),
            'n_clusters_trace': n_clusters_post.tolist(),
            'final_clusters': final_clusters,
        }

    def _test_bimodality(self, samples):
        """Check if posterior predictive is bimodal via Hartigan's dip test."""
        try:
            import diptest
            dip_stat, p_value = diptest.diptest(np.array(samples))
            return p_value < 0.05  # Reject unimodality
        except ImportError:
            # Fallback to KDE mode counting
            from scipy.stats import gaussian_kde
            try:
                kde = gaussian_kde(samples)
                x = np.linspace(np.min(samples), np.max(samples), 200)
                density = kde(x)
                modes = 0
                for i in range(1, len(density) - 1):
                    if density[i] > density[i-1] and density[i] > density[i+1]:
                        modes += 1
                return modes >= 2
            except Exception:
                return False


def run_sensitivity(yi, sei, alphas=None):
    """Run DPMA with different concentration parameters to test robustness."""
    if alphas is None:
        alphas = [0.1, 0.5, 1.0, 2.0, 5.0]
    results = {}
    for alpha in alphas:
        dp = DPMetaAnalysis(n_iter=1000, burn=200, alpha=alpha)
        result = dp.fit(yi, sei)
        results[alpha] = {
            'n_clusters': result['n_clusters_mode'],
            'theta': result['theta_posterior_mean'],
        }
    return results


# ═══════════════════════════════════════════════════
# CONVERGENCE DIAGNOSTICS
# ═══════════════════════════════════════════════════

def check_convergence(traces):
    """Basic convergence check: split-chain stationarity test.

    Compares first-half vs second-half means of the trace.
    Returns True if trace appears stationary (split difference < 0.5 SD).
    """
    traces = np.asarray(traces, dtype=float)
    n = len(traces)
    if n < 4:
        return True  # Too short to assess
    half = n // 2
    first_half_mean = np.mean(traces[:half])
    second_half_mean = np.mean(traces[half:])
    sd = np.std(traces)
    split_diff = abs(first_half_mean - second_half_mean) / (sd + 1e-10)
    return bool(split_diff < 0.5)


# ═══════════════════════════════════════════════════
# STANDARD DL FOR COMPARISON
# ═══════════════════════════════════════════════════

def standard_dl(yi, sei):
    k = len(yi)
    wi = 1.0 / sei**2
    sw = np.sum(wi)
    theta_fe = np.sum(wi * yi) / sw
    Q = float(np.sum(wi * (yi - theta_fe)**2))
    C = float(sw - np.sum(wi**2) / sw)
    tau2 = max(0, (Q - (k-1)) / C) if C > 0 else 0
    ws = 1.0 / (sei**2 + tau2)
    sws = np.sum(ws)
    theta = float(np.sum(ws * yi) / sws)
    se = float(1.0 / math.sqrt(sws))
    I2 = max(0, (Q - (k-1)) / Q * 100) if Q > 0 else 0

    # Standard PI
    if k >= 3:
        t_crit = sp_stats.t.ppf(0.975, k - 2)
        pi_se = math.sqrt(tau2 + se**2)
        pi_lo = theta - t_crit * pi_se
        pi_hi = theta + t_crit * pi_se
    else:
        pi_lo, pi_hi = -999, 999

    return {
        'theta': theta, 'se': se, 'tau2': tau2, 'I2': I2,
        'pi_lo': pi_lo, 'pi_hi': pi_hi,
    }


# ═══════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════

def load_review(rda_path):
    result = pyreadr.read_r(str(rda_path))
    df = list(result.values())[0].copy()
    df.columns = df.columns.str.replace(' ', '.', regex=False)
    review_id = rda_path.stem.split('_')[0]

    import pandas as pd
    groups = []
    for (grp, num), sub in df.groupby(['Analysis.group', 'Analysis.number']):
        has_binary = (sub['Experimental.cases'].notna() & (sub['Experimental.cases'] > 0)).any()
        groups.append({'grp': grp, 'num': num, 'k': len(sub), 'binary': has_binary})
    if not groups: return None
    gdf = pd.DataFrame(groups)
    binary = gdf[gdf['binary']]
    best = binary.loc[binary['k'].idxmax()] if len(binary) > 0 else gdf.loc[gdf['k'].idxmax()]
    primary = df[(df['Analysis.group'] == best['grp']) & (df['Analysis.number'] == best['num'])]

    has_binary = (primary['Experimental.cases'].notna() & (primary['Experimental.cases'] > 0)).any()
    scale = 'ratio' if has_binary else ('ratio' if (primary['Mean'].dropna() > 0).all() else 'difference')

    if scale == 'ratio':
        v = (primary['Mean'].notna() & (primary['Mean'] > 0) & primary['CI.start'].notna() & (primary['CI.start'] > 0) & primary['CI.end'].notna() & (primary['CI.end'] > 0))
        sub = primary[v]
        if len(sub) < 5: return None
        yi = np.log(sub['Mean'].values.astype(float))
        sei = (np.log(sub['CI.end'].values.astype(float)) - np.log(sub['CI.start'].values.astype(float))) / (2 * 1.96)
    else:
        v = primary['Mean'].notna() & primary['CI.start'].notna() & primary['CI.end'].notna()
        sub = primary[v]
        if len(sub) < 5: return None
        yi = sub['Mean'].values.astype(float)
        sei = (sub['CI.end'].values.astype(float) - sub['CI.start'].values.astype(float)) / (2 * 1.96)

    ok = (sei > 0) & np.isfinite(yi) & np.isfinite(sei)
    yi, sei = yi[ok], sei[ok]
    if len(yi) < 5: return None
    return {'review_id': review_id, 'yi': yi, 'sei': sei, 'k': len(yi), 'scale': scale}


# ═══════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Dirichlet Process Meta-Analysis (DPMA)")
    print("=" * 45)

    t0 = time.time()
    rda_files = sorted(PAIRWISE_DIR.glob('*.rda'))

    dpma = DPMetaAnalysis(alpha=1.0, mu0=0.0, sigma0=1.0, n_iter=5000, burn=1000)
    results = []
    n_processed = 0

    for rda in rda_files:
        review = load_review(rda)
        if review is None:
            continue

        yi, sei, k = review['yi'], review['sei'], review['k']

        # Standard DL
        dl = standard_dl(yi, sei)

        # DPMA
        dp = dpma.fit(yi, sei)

        # Convergence check
        if not check_convergence(dp['n_clusters_trace']):
            print(f"  WARNING: convergence concern for {review['review_id']} (split-chain test)")

        # Comparison
        theta_shift = abs(dl['theta'] - dp['theta_posterior_mean'])
        width_dl = dl['pi_hi'] - dl['pi_lo']
        width_dp = dp['theta_95_hi'] - dp['theta_95_lo']
        width_ratio = width_dp / width_dl if width_dl > 0 else 1

        row = {
            'review_id': review['review_id'],
            'k': k,
            'scale': review['scale'],
            # DL
            'dl_theta': round(dl['theta'], 4),
            'dl_se': round(dl['se'], 4),
            'dl_tau2': round(dl['tau2'], 4),
            'dl_I2': round(dl['I2'], 1),
            'dl_pi_lo': round(dl['pi_lo'], 4),
            'dl_pi_hi': round(dl['pi_hi'], 4),
            'dl_pi_width': round(width_dl, 4),
            # DPMA
            'dp_theta': round(dp['theta_posterior_mean'], 4),
            'dp_sd': round(dp['theta_posterior_sd'], 4),
            'dp_ci_lo': round(dp['theta_95_lo'], 4),
            'dp_ci_hi': round(dp['theta_95_hi'], 4),
            'dp_width': round(width_dp, 4),
            'dp_n_clusters': dp['n_clusters_mode'],
            'dp_bimodal': dp['predictive_bimodal'],
            # Comparison
            'theta_shift': round(theta_shift, 4),
            'width_ratio': round(width_ratio, 3),
            'dp_finds_subgroups': dp['n_clusters_mode'] >= 2,
        }
        results.append(row)
        n_processed += 1

        if n_processed % 50 == 0:
            print(f"  {n_processed} reviews processed...")

    elapsed = time.time() - t0
    n = len(results)
    print(f"  Total: {n} reviews in {elapsed:.1f}s")

    # HEADLINES
    n_subgroups = sum(1 for r in results if r['dp_finds_subgroups'])
    n_bimodal = sum(1 for r in results if r['dp_bimodal'])
    shifts = [r['theta_shift'] for r in results]
    widths = [r['width_ratio'] for r in results]
    clusters = [r['dp_n_clusters'] for r in results]

    print(f"\n{'='*55}")
    print("HEADLINE FINDINGS")
    print(f"{'='*55}")
    print(f"  DPMA finds subgroups: {n_subgroups}/{n} ({100*n_subgroups/n:.1f}%)")
    print(f"  Bimodal predictive: {n_bimodal}/{n} ({100*n_bimodal/n:.1f}%)")
    print(f"  Mean clusters: {np.mean(clusters):.2f} (median: {np.median(clusters):.0f})")
    print(f"  Cluster distribution: {Counter(clusters)}")
    print(f"  Mean theta shift (DL vs DP): {np.mean(shifts):.4f}")
    print(f"  Mean width ratio (DP/DL PI): {np.mean(widths):.2f}")
    print(f"  DP wider in {sum(1 for w in widths if w > 1)}/{n} ({100*sum(1 for w in widths if w > 1)/n:.1f}%)")

    # EXPORT
    fields = list(results[0].keys())
    with open(OUTPUT_DIR / 'dpma_results.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)

    summary = {
        'n_reviews': n,
        'n_subgroups_found': n_subgroups,
        'pct_subgroups': round(100 * n_subgroups / n, 1),
        'n_bimodal_predictive': n_bimodal,
        'cluster_distribution': dict(Counter(clusters)),
        'mean_theta_shift': round(float(np.mean(shifts)), 4),
        'mean_width_ratio': round(float(np.mean(widths)), 2),
        'elapsed_seconds': round(elapsed, 1),
    }
    with open(OUTPUT_DIR / 'dpma_summary.json', 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)

    print(f"\n  Saved to {OUTPUT_DIR}/")


if __name__ == '__main__':
    main()
