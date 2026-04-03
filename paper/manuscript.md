# Dirichlet Process Meta-Analysis: Bayesian Nonparametric Discovery of Latent Subgroups Across 307 Cochrane Reviews

**Mahmood Ahmad**^1

1. Royal Free Hospital, London, United Kingdom

**Correspondence:** Mahmood Ahmad, mahmood.ahmad2@nhs.net | **ORCID:** 0009-0003-7781-4478

---

## Abstract

**Background:** Standard random-effects meta-analysis assumes study effects arise from a single normal distribution. If the true effect distribution is multimodal — for example, a treatment benefits one subpopulation but not another — this assumption masks clinically important heterogeneity within a single variance component.

**Methods:** We applied Dirichlet process mixture (DPM) modeling to 307 Cochrane systematic reviews from the Pairwise70 dataset (k >= 5). For each review, a collapsed Gibbs sampler with 5,000 posterior iterations discovered the number of latent effect clusters nonparametrically, without pre-specifying subgroup count or membership. Results were compared against DerSimonian-Laird random-effects estimates. Bimodality was confirmed via Hartigan's dip test on the posterior predictive distribution.

**Results:** Among 307 reviews, 282 (91.9%) revealed two or more distinct effect clusters. The most common cluster counts were 2 (158 reviews, 51.5%), 3 (84 reviews, 27.4%), and 4-11 (40 reviews, 13.0%). The mean absolute shift in the pooled estimate relative to DerSimonian-Laird was 0.22 (95% CI 0.18-0.26). Bimodal posterior predictive distributions were confirmed by Hartigan's dip test in 161 reviews (52.4%). The Dirichlet process predictive distribution was 3.08 times wider than the standard prediction interval on average, reflecting the multimodal uncertainty that normal-theory prediction intervals understate. Results were robust across concentration parameters (alpha = 0.1 to 5.0).

**Conclusions:** Over 90% of Cochrane meta-analyses contain evidence of latent subgroups that the standard random-effects model absorbs into a single variance component. Bayesian nonparametric meta-analysis reveals this hidden structure and provides more honest predictive distributions for clinical decision-making.

**Keywords:** Dirichlet process, Bayesian nonparametric, meta-analysis, heterogeneity, latent subgroups, predictive distribution

---

## 1. Introduction

The random-effects model, introduced by DerSimonian and Laird in 1986,^1 assumes that study-specific true effects theta_i are drawn from a normal distribution N(mu, tau^2). This model has a single mode: all studies are exchangeable draws from a unimodal distribution, and heterogeneity is characterised entirely by the variance tau^2.

In practice, many meta-analyses combine studies that differ qualitatively — different patient populations, different drug doses, different control conditions — and the true effect distribution may have multiple modes. A treatment that benefits younger patients (theta ~ -0.5) but harms older patients (theta ~ +0.2) would produce a bimodal effect distribution, but the standard model would report a single pooled estimate near zero with high I-squared, obscuring both the benefit and the harm.

The Dirichlet process^2 provides a Bayesian nonparametric alternative. Instead of assuming theta_i ~ N(mu, tau^2), we assume theta_i ~ G, where G itself is drawn from a Dirichlet process: G ~ DP(alpha, G_0). This allows the effect distribution to have an unbounded number of components, with the actual number of clusters learned from the data. The concentration parameter alpha controls the expected number of clusters: small alpha favours few clusters; large alpha favours many.

Cao et al. (2024) applied DPM to reference interval estimation;^3 we extend the approach to meta-analytic effect estimation at scale, applying it to 307 Cochrane reviews.

---

## 2. Methods

### 2.1 Data Source

The Pairwise70 dataset provides study-level data from 501 Cochrane reviews. We included reviews with k >= 5 studies to ensure sufficient data for cluster discovery, yielding 307 reviews.

### 2.2 Dirichlet Process Mixture Model

For study i in review j:

    y_i | theta_i, s_i ~ N(theta_i, s_i^2)
    theta_i | G ~ G
    G ~ DP(alpha, G_0)
    G_0 = N(mu_0, sigma_0^2)

where y_i is the observed effect, s_i is the known standard error, G is the random effect distribution (a discrete distribution with probability 1 under the Dirichlet process), alpha is the concentration parameter, and G_0 is the base measure.

### 2.3 Collapsed Gibbs Sampler

We integrate out G analytically (Blackwell-MacQueen Polya urn scheme)^4 and sample cluster assignments directly. At each iteration, study i is assigned to:

- An existing cluster c with probability proportional to n_c * f(y_i | mu_c, s_i^2 + sigma_c^2)
- A new cluster with probability proportional to alpha * integral f(y_i | theta, s_i^2) G_0(theta) d(theta)

Cluster parameters are updated via conjugate normal-normal updates. We run 5,000 iterations with a 1,000-iteration burn-in.

### 2.4 Outputs

For each review:
- Number of clusters (posterior mode)
- Cluster-specific pooled estimates and sizes
- Overall DPM pooled estimate (posterior mean of the mixture)
- Posterior predictive distribution (full distribution, not just mean + SD)
- Hartigan's dip test for bimodality of the predictive distribution

### 2.5 Comparison with Standard Methods

Each review was also analysed with DerSimonian-Laird (DL) random-effects. We computed:
- **Theta shift:** |DPM estimate - DL estimate|
- **Width ratio:** DPM predictive distribution width / DL prediction interval width
- **Subgroup discovery:** Whether DPM found >= 2 clusters

---

## 3. Results

### 3.1 Prevalence of Latent Subgroups

Of 307 reviews, 282 (91.9%) showed two or more latent clusters under the DPM model (Table 1).

**Table 1. Distribution of discovered cluster counts**

| Clusters | Reviews | % |
|----------|---------|---|
| 1 | 25 | 8.1 |
| 2 | 158 | 51.5 |
| 3 | 84 | 27.4 |
| 4 | 20 | 6.5 |
| 5 | 11 | 3.6 |
| 6-11 | 9 | 2.9 |
| **Total** | **307** | **100.0** |

### 3.2 Impact on Pooled Estimates

The mean absolute shift between the DPM and DL pooled estimates was 0.22 (95% CI 0.18-0.26) on the log scale. This shift was not random noise: it was systematically associated with the number of discovered clusters (Spearman rho = 0.41 between cluster count and theta shift).

### 3.3 Predictive Distributions

The DPM posterior predictive distribution was 3.08 times wider than the standard DL prediction interval on average. This reflects the multimodal uncertainty: when the effect distribution has two modes, the predictive distribution is bimodal, and its width reflects the separation between modes rather than just the variance within a single mode.

Hartigan's dip test confirmed bimodality in 161 reviews (52.4%). Among reviews with I-squared < 50% under the standard model, 38% still showed DPM-detected multimodality — indicating that moderate I-squared can mask bimodal structure.

### 3.4 Sensitivity to Concentration Parameter

Results were robust across concentration parameters alpha = 0.1, 0.5, 1.0, 2.0, and 5.0. The subgroup discovery rate ranged from 87% (alpha = 0.1, favouring fewer clusters) to 96% (alpha = 5.0, favouring more clusters). The core finding — that the vast majority of reviews show multi-cluster structure — was insensitive to the prior choice.

---

## 4. Discussion

### 4.1 Principal Findings

Over 90% of Cochrane meta-analyses contain latent subgroup structure detectable by Bayesian nonparametric methods. The standard random-effects model absorbs this structure into a single variance component (tau-squared), producing a pooled estimate that may not represent any actual subgroup and a prediction interval that understates the true spread of treatment effects.

### 4.2 Clinical Implications

The practical consequence is that clinical decisions based on the standard pooled estimate may be misleading. If a drug has a large benefit in one subgroup and no benefit in another, the pooled estimate shows a moderate benefit — which applies to neither subgroup. The DPM predictive distribution reveals this bimodality, showing that the next patient might experience either a large benefit or no benefit, rather than a moderate benefit.

### 4.3 Relationship to Existing Methods

Subgroup analysis in meta-analysis typically requires pre-specified subgroup variables (e.g., age, sex, dose). The Dirichlet process discovers subgroups from the effect size distribution itself, without any covariate information. This is both a strength (no covariates needed) and a limitation (discovered clusters cannot be labelled without additional investigation).

### 4.4 Limitations

The exchangeability assumption within discovered clusters may not hold if studies within a cluster differ in important ways. The method requires k >= 5 studies for stable cluster discovery. The Gibbs sampler convergence was assessed by trace plots but not formally diagnosed for all 307 reviews. Runtime was approximately 25 minutes for the full pipeline (5,000 iterations per review, 307 reviews).

---

## 5. Conclusion

Bayesian nonparametric meta-analysis via Dirichlet process mixtures reveals that over 90% of Cochrane meta-analyses contain latent subgroup structure invisible to standard methods. This hidden heterogeneity produces pooled estimates that may not represent any actual patient subgroup and prediction intervals that understate the true uncertainty. We recommend DPM analysis as a complement to standard random-effects meta-analysis for any review where clinical heterogeneity is suspected.

---

## References

1. DerSimonian R, Laird N. Meta-analysis in clinical trials. *Control Clin Trials*. 1986;7:177-188.
2. Ferguson TS. A Bayesian analysis of some nonparametric problems. *Ann Stat*. 1973;1:209-230.
3. Cao J, et al. A Bayesian nonparametric meta-analysis model for estimating the reference interval. *Stat Med*. 2024;43(5).
4. Blackwell D, MacQueen JB. Ferguson distributions via Polya urn schemes. *Ann Stat*. 1973;1:353-355.
5. Neal RM. Markov chain sampling methods for Dirichlet process mixture models. *J Comput Graph Stat*. 2000;9:249-265.
6. Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. *Stat Med*. 2002;21:1539-1558.

---

## Data Availability

Code and results at https://github.com/mahmood726-cyber/dpma (MIT licence). Pipeline requires NumPy; runs on the Pairwise70 dataset.
