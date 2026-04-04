Mahmood Ahmad
Tahir Heart Institute
author@example.com

Dirichlet Process Meta-Analysis

Does the standard random-effects assumption of a single normal effect distribution mask clinically important subgroups within meta-analyses? We applied Dirichlet process mixture modeling to 307 Cochrane systematic reviews spanning diverse therapeutic areas, using a collapsed Gibbs sampler with five thousand posterior iterations per review. The Dirichlet process discovers the number of latent clusters nonparametrically, without requiring the analyst to pre-specify subgroup count or membership. Among these reviews, 282 (92 percent) revealed two or more distinct effect clusters, with a mean pooled-estimate shift of 0.22 (95% CI 0.18 to 0.26) versus DerSimonian-Laird. Results held across concentration parameters (alpha 0.1 to 5.0), with bimodal predictive distributions in 161 reviews confirmed by Hartigan dip tests. Bayesian nonparametric meta-analysis uncovers pervasive latent heterogeneity that conventional models absorb into a single variance component, potentially misleading clinical decisions broadly. The method is limited by exchangeability within discovered clusters and by requiring five or more studies; network or dose-response settings need dedicated extensions.

Outside Notes

Type: methods
Primary estimand: summary effect
App: DPMA v1.0
Data: Repository artifacts in /mnt/c/Models/DPMA
Code: https://github.com/mahmood726-cyber/dpma
Version: 1.0
Certainty: moderate
Validation: DRAFT

References

1. Crippa A, Orsini N. Dose-response meta-analysis of differences in means. BMC Med Res Methodol. 2016;16:91.
2. Greenland S, Longnecker MP. Methods for trend estimation from summarized dose-response data, with applications to meta-analysis. Am J Epidemiol. 1992;135(11):1301-1309.
3. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.
