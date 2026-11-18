# Pilot report (Phase 3, deliverable D5)

**Date:** 17 Nov 2026. **Data:** 400 held-out species (272 land, 128 marine) drawn from the lowest and
highest 20 % of the primary predictor (seed 20261117). These species are excluded from every
confirmatory test.

## Purpose
The pilot checks that the pipeline runs end to end and estimates the residual structure needed for a
power analysis. It is not used to judge any hypothesis.

## Pipeline check: passed
GTDB r232 metadata, then habitat classification, then OISST/GHCN_CAMS/PREC/L variability (1991–2020),
then NCBI genome download and pyrodigal gene calls, then pyhmmer with the 174-model dictionary, then
module rules (`scripts/08_aggregate_modules.py`), then phylolm/phyloglm on the pruned GTDB tree.
- All 400 genomes were downloaded and annotated. No genome failed QC (≥ 500 proteins).
- Module prevalence in the pilot:
  - 11 sporulators (2.75 %; all Bacillota).
  - Rpf in 20 %.
  - Mean of 1.07 hibernation systems.
  - Mean of 1.48 stringent-response proteins.
  - Mean of 3.5 paired TA loci.
- The phylogenetic signal of the residuals is high (Pagel's λ = 0.89 for DORM_z). Phylogenetic
  control is therefore essential.

## Pilot estimates (all models; not inferential)
| id               |   n |   estimate |     se |      p |   lambda_or_alpha |
|:-----------------|----:|-----------:|-------:|-------:|------------------:|
| H1_DORM_P0       | 400 |     0.009  | 0.0168 | 0.5919 |            0.8932 |
| H1_DORM_P0_total | 400 |    -0.0029 | 0.0174 | 0.866  |            0.9194 |
| PATH_gsize_P0    | 400 |    -0.0128 | 0.006  | 0.034  |            0.9823 |
| H1b_rpf          | 400 |    -0.0017 | 0.0054 | 0.751  |            1      |
| H1c_hib          | 400 |     0.0224 | 0.015  | 0.1356 |            0.5944 |
| H1d_sr           | 400 |    -0.0002 | 0.0139 | 0.9888 |            0.7298 |
| H6_seasonal      | 400 |     0.009  | 0.0168 | 0.5928 |            0.8934 |
| H6_resid         | 400 |     0.007  | 0.0165 | 0.673  |            0.8929 |
| H6_ac1           | 400 |    -0.0254 | 0.0164 | 0.1211 |            0.8937 |
| H6_colwellP      | 400 |    -0.0044 | 0.0165 | 0.7911 |            0.8935 |
| H6_colwellM      | 400 |     0.0131 | 0.0168 | 0.436  |            0.8934 |
| H6_colwellC      | 400 |    -0.0102 | 0.0168 | 0.5441 |            0.8932 |
| H3_TA_marine     | 400 |    -0.2716 | 0.1013 | 0.0077 |            0.7213 |
| H3_TA_marine_mob | 400 |    -0.243  | 0.0961 | 0.0118 |            0.7399 |
| H3_TA_P0         | 400 |    -0.0803 | 0.039  | 0.04   |            0.7377 |
| S_isolates       | 107 |     0.023  | 0.0278 | 0.4106 |            0.9658 |
| S_comp95         | 273 |     0.0121 | 0.0185 | 0.5116 |            0.9303 |
| S_spmed          | 400 |     0.0099 | 0.0169 | 0.559  |            0.8934 |
| S_marine         | 128 |     0.0029 | 0.0316 | 0.927  |            0.9604 |
| S_land           | 272 |     0.0167 | 0.0213 | 0.4344 |            0.8637 |
| S_noBacCPR       | 373 |     0.0124 | 0.0169 | 0.4622 |            0.8095 |
| S_mobile         | 400 |     0.0111 | 0.0169 | 0.5122 |            0.8936 |
| S_sum            | 400 |     0.0024 | 0.0161 | 0.8836 |            0.9397 |
| S_Tmean          | 400 |     0.0162 | 0.0201 | 0.4207 |            0.8936 |
| X_chl_total      | 119 |    -0.0022 | 0.0296 | 0.9396 |            0.9617 |
| X_chl_resid      | 117 |     0      | 0.0286 | 0.9996 |            0.9604 |
| X_prec_cv        | 272 |    -0.0041 | 0.0232 | 0.8605 |            0.8684 |
| X_prec_cvres     | 272 |     0.0103 | 0.0241 | 0.6681 |            0.8674 |

Sporulation models did not run in the pilot (fewer than 30 Bacillota). In the confirmatory set they
are fitted within Bacillota only (about 330 species).

## Power (simulation, `scripts/11_power.R`, `results/power_simulation.tsv`)
- **Method.** Residuals were simulated under Pagel's λ = 0.893 and σ² = 0.388, taken from the pilot
  H1 fit, on the confirmatory tree (n = 7,467) using the real confirmatory predictor values. Each
  scenario had 100 replicates and used a two-sided α of 0.05.

| true β (DORM_z units per 1 SD of P0) | power | SD of estimate |
|---|---|---|
| 0 | 0.01 (type I) | 0.0027 |
| 0.005 | 0.37 | 0.0029 |
| 0.010 | 0.89 | 0.0031 |
| 0.015 | 1.00 | 0.0028 |
| 0.020 | 1.00 | 0.0032 |

- **Effect size.** The pilot point estimate for H1 is +0.006 to +0.009 (SE 0.012–0.017). That is
  smaller than the SESOI (below), so the pilot gives no hint of a large effect.
- **SESOI.** The smallest effect size of interest is **β = 0.015 DORM_z units per SD of P0**. This is
  about 0.027 SD of the composite per SD of variability, an explained variance of roughly 0.07 %.
  Anything smaller is biologically negligible at the genome scale.
- **What the design can detect.** It has essentially 100 % power at the SESOI. A null result can
  therefore be read as evidence of absence through an equivalence test: the 95 % CI upper bound is
  below 0.015.

## Problems found and fixed during the pilot
1. **Wrong reference accession.** My guessed "Synechococcus CC9311" accession was in fact the spore
   former *Syntrophomonas wolfei*; the dictionary called it correctly. The reference table was fixed.
2. **KOfam K06381 (SpoIID) was non-specific.** It hit amidases outside Bacillota, so it was removed
   from the scoring.
3. **The 3D (Sps) domain also occurs in MltA.** Sps is therefore counted only in Bacillota.
4. **Dependencies.** `data.table::fread` on gz needed R.utils, and the chlorophyll stack needed a
   non-dask loader.
5. **Literature audit changes to the design.**
   - Genome size is treated as a possible mediator: both the total and the direct effect are
     reported, plus the path variability → genome size.
   - Rpf is also modelled within Actinomycetota.
   - Precipitation for land sites and the niche-breadth predictors are added as secondary.
