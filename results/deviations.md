# Deviations from the pre-registration

The pre-registration was frozen on 17 Nov 2026 at 10:16 UTC (`prereg/PREREGISTRATION.md`). Every
departure or addition made after that point is listed here, together with the reason for it.

## Procedural

1. **The DefenseFinder TA annotation was not used.** The TA module uses only the pre-registered
   HMM dictionary rules: toxin and antitoxin hits paired within 3 genes. DefenseFinder had been
   considered as a cross-check, but it cost about 70 CPU-seconds per genome, which was too much
   for this compute budget. This changes no pre-registered model.
2. **gRodon2 was re-run (and resumed after a container restart).** The pilot gRodon run was killed by the session time limit after 200
   genomes. gRodon2 was then run on all analysis genomes with a restart-safe chunk runner
   (`scripts/run_growth_chunks.sh`). The method and code are unchanged.
3. **HMM annotation was re-run for some batches.** Batch 2 was interrupted by a worker restart and
   repeated from the start (`scripts/run_chunks.sh`). A batch's output is written only once the
   batch completes, so no partial outputs were used.
4. **A preliminary model run was made to debug the code.** Before gRodon2 had finished, the
   confirmatory models were run once on a table whose growth-rate columns were incomplete. This
   checked that the code works at full scale. Its estimates are not reported. The reported run
   uses the final table, built after gRodon2 had finished, because the QC rule on ribosomal genes
   (n_HE < 20) depends on gRodon2 output.

## Additional analyses (not pre-registered; labelled exploratory wherever reported)

5. **Sporulator threshold ≥ 27.** An independent phenotype calibration (`validation/`) gave a
   Youden-optimal threshold of ≥ 27 core genes. The pre-registered threshold is ≥ 25, and the
   pre-registration does **not** list a ≥ 27 sensitivity analysis. (`validation/README.md` says
   it does; that statement is wrong.) H1a with a ≥ 27 threshold is reported as exploratory only.
   The confirmatory decision uses ≥ 25.
6. **Genus-level sanity check** (`scripts/15_genus_sanity.py`). This is descriptive only.
7. **Exploratory model set** (`scripts/exploratory_config.tsv`, run through `09_models.R` with an optional
   5th argument). It adds the MAG × P0 interaction, MAG-only and isolate-only fits, per-phylum fits,
   per-domain module fits, single hibernation/SR factor presence, assembly-quality covariates (N50,
   contamination), absolute latitude, and TA checks. All are labelled exploratory.
8. **Site pseudoreplication check** (`scripts/16_cluster_robust.R`). This problem was not anticipated
   in the pre-registration. The confirmatory MAGs come from only about 370 sites, and the 10 largest
   sites hold about 41 % of them. Species from the same site share an identical predictor value,
   so phylolm SEs overstate the amount of environmental replication. We therefore added:
   - a site-grouped jackknife SE;
   - site-capped resampling (≤ 1 or ≤ 3 species per 0.5° site);
   - MCMCglmm with phylogenetic and site random effects.
   These do not replace the pre-registered decision. They are reported next to it, and they decide
   how far that decision can be generalised.
9. **Re-parameterised variability (exploratory).** Following the literature update
   (`docs/phase0/novelty_update_2026-11.md`; Marshall & Burgess 2015; Barneche et al. 2018), we
   added joint models of:
   - log total SD;
   - the seasonal fraction of variance (`T_frac_seasonal`);
   - residual lag-1 autocorrelation ("colour").
   We also added marine-only seasonal + residual models. The reason is that on land, log total SD
   and log seasonal SD correlate at r ≈ 1.00, so the predictable/unpredictable split can only be
   identified in the ocean. These models do not replace the pre-registered H6 tests; they sit
   beside them.
10. **Follow-up tests of the post-hoc explanations** (`scripts/followup_config.tsv`,
    `results/followup_models.tsv`) and a pilot-set replication
    (`scripts/pilot_replication_config.tsv`, `results/pilot_replication_models.tsv`). These were
    designed after the confirmatory results were seen. All are exploratory.
    - New columns in `10_build_table.py`: coordinate precision, hibernation and RSH paralogue
      copies, and an isolate flag.
    - **The collider-style model** (is_isolate ~ P0 × DORM_z) is the mathematical mirror of the
      MAG × P0 interaction. It is a reframing, not independent evidence.
11. **MCMCglmm implementation.** With the inverse relationship matrix built on tips only
    (`inverseA(nodes = "TIPS")`), 300 MCMC iterations took more than 10 minutes on the
    7,327-species tree, so the pre-registered sensitivity analysis 10 could not finish. We switched
    to the standard sparse formulation (`nodes = "ALL"`), which includes the internal nodes. It
    implies the same phylogenetic covariance among tips and ran about 200 times faster. Chain
    lengths are as pre-registered.
    - We also added an exploratory Gaussian MCMCglmm with phylogenetic and site random effects for
      the MAG × P0 interaction (`scripts/21_mcmc_site_interaction.R`).
    - The `run_final.sh` / `run_final2.sh` runners were stopped at that point. Their remaining steps
      (MCMCglmm and the second-stage MAG/isolate script) were run directly. The second-stage script
      had already been run on the final table.

## Corrections after the independent audit (`audit/README.md`) and the within-species replication (`replication_within_species/`)
12. **MCMCglmm branch lengths (audit M1).** The MCMCglmm models (scripts 14, 16, 21) had used
    Grafen-rescaled branch lengths, unlike every phylolm model. All three were refit on the GTDB
    bac120 branch lengths (non-positive edges set to 1e-8; unscaled inverse). The reported values
    are the refit ones; the Grafen versions are superseded.
13. **QC rule (audit M7).** `10_build_table.py` let through 6 genomes whose gRodon run failed (no
    ribosomal-gene count), contrary to the pre-registered rule. The code is now fixed. The reported
    models keep the as-run table (n = 7,327). The corrected rerun (n = 7,321;
    `audit/qc_fix_rerun.tsv`) gives H1 +0.0064 [0.0016, 0.0113] and interaction +0.024. No decision
    changes.
14. **Pilot wording (audit M2).** The pilot is the *discovery* sample for the marine TA deficit, so
    it is not a replication. The confirmatory set replicated the pilot estimate (−0.27). The
    earlier "held-out pilot replication" wording has been withdrawn.
15. **SAGs (audit M4).** 92 single-amplified genomes were coded as MAGs. A SAG-only slope is
    uninformative (+0.002 ± 0.024), and removing SAGs changes nothing in the within-species
    replication.
16. **Reinterpretation of the MAG/isolate reversal (audit M3–M4; within-species replication).**
    - Only 3 of the 4,919 MAG-represented species also have an in-domain isolate genome. "MAG vs
      isolate" therefore means "never-cultured vs cultured lineages".
    - Within species (7,098 genomes, 1,473 species, species fixed effects, site-clustered CIs), the
      composite reversal is absent.
    - The audit's strain-heterogeneity tests do not support the "population consensus"
      explanation, so it has been dropped.
    - The extra-hibernation-copy MAG slope remains positive within species. It is not robust to
      correction across the ten responses tested.
17. **Balon numbers (audit M5)** have been corrected to the logit slope plus the full octile curve.
    **Marine TA citation (audit M6):** Lennon & Jones (2011) has been replaced by Meaden et al. (2025)
    and Beavogui et al. (2024).
18. **Post-hoc checks, set 1 (18 Nov 2026).** All are exploratory additions; no
    pre-registered decision changed.
    - Formal TOST equivalence tests (`scripts/22_equivalence_tests.py`).
    - TA-deficit estimates by genome-quality tier (`scripts/review_config.tsv`). The paper's headline
      TA figure is now qualified as 21 % overall and 10–15 % in the highest-quality genomes.
    - BioSample collection dates and a within-climate-window sensitivity analysis
      (`scripts/23_collection_dates.py`, `scripts/dates_config.tsv`).
19. **Post-hoc checks, set 2 (18 Nov 2026).** All additions are exploratory; no
    pre-registered decision changed.
    - TOST operating-characteristic simulation (`scripts/24_tost_power.R`).
    - Spatial-block sensitivity at 0.5–5° (`scripts/25_spatial_blocks.R`).
    - Clade-aware dictionary validation (`validation/06_clade_holdout.py`).
    - Variability-component analysis, including interannual SD (`scripts/variability_review_config.tsv`).
    - TA tier-interaction and single-component tests (`scripts/ta_review_config.tsv`).
    - SESOI translation into gene counts (`results/sesoi_translation.tsv`).
    - Collection-date parsing now records the format and requires the whole stated period to lie
      within 1991–2020.
    - **Wording changes:**
      - the TA claim was downgraded;
      - P0 is described as seasonal amplitude;
      - the pre-registration is described as an internal dated plan, not a registry entry;
      - the audit is described as a computational check.
20. **Post-hoc checks, set 3 (18 Nov 2026).** All additions are exploratory; no
    pre-registered decision changed.
    - **Primary equivalence procedure P** defined and calibrated by simulation
      (`scripts/26_tost_calibration.R`): 5,000 replicates at each boundary, 2,000 at 0, 1,000 at
      0.005 and 0.010. Boundary false-equivalence rate 4.75 % (4.3–5.2 %), so P is calibrated.
    - **Slopes by genome source**, with interaction tests and within-group equivalence, plus margin
      sensitivity at ±0.010 and ±0.020 (`scripts/27_group_equivalence.py`).
    - **BioSample retrieval bug fixed.** 121 accessions were resolved to UIDs via esearch and
      verified. All date-based numbers were recomputed.
    - **Title and conclusions reworded** to describe the pooled average, not an absence of any
      relationship.
21. **Post-hoc checks, set 4 (18 Nov 2026).** All additions are exploratory; no
    pre-registered decision changed.
    - **Noise scale of the calibration** (`scripts/29_calibration_noise_scale.py`,
      `scripts/29b_linearity_check.R`). **Withdrawn in item 22.** The claim was that the simulated
      data were noisier than the real data (slope SD 0.0035 vs real-data jackknife SE 0.0022), and
      that power at ±0.010 at the real noise level was 69–70 %. The linearity argument holds, but
      its premise was wrong.
    - **One-sided tests** for the genome-source slopes, against 0 and against the SESOI
      (`scripts/27_group_equivalence.py`). MAG slope: p = 0.003 for β > 0, p = 0.33 for β > 0.015.
    - **Composition-adjusted interaction** (`scripts/28_interaction_composition.R`): habitat,
      phylum, absolute latitude and mean temperature, each with its own P0 slope, plus a
      common-support model over 18 phylum × habitat strata. Interaction +0.024 to +0.028 throughout.
    - **BioSample provenance log** committed (`scripts/23b_biosample_provenance.py`,
      `data_small/final/biosample_retrieval_log.tsv`).
    - **Calibration file header corrected.** The first three column names of
      `results/tost_calibration_sims.tsv` are now true_beta, seed, estimate. The values are unchanged.
    - **Wording:**
      - the within-species result is described as "not detected", not as absent;
      - the MAG slope is a positive association, not a demonstrated effect above the SESOI;
      - Simpson's paradox is not invoked;
      - new title: "Dormancy-associated genomic repertoire across bacterial taxa and the seasonal
        temperature amplitude of their sampling environments".
22. **Post-hoc checks, set 5 (18 Nov 2026).** All additions are exploratory; no
    pre-registered decision changed.
    - **Grouping sensitivity** (`scripts/30_grouping_sensitivity.R`).
      - The single random grouping used for procedure P gave the 4th-smallest of 200 SEs.
      - The equivalence decision holds in 199 of 200 groupings, under delete-one-cell, and with
        contiguous groups.
      - The reported precision changes: p_TOST goes from 0.0005 to 0.012.
    - **Final procedure P̄.** The jackknife variance is averaged over 200 random groupings.
      - This is a change in how the SE of P is computed, made after seeing that one grouping was
        unrepresentative. It is reported as such, with the single-grouping result,
        delete-one-cell and contiguous groups alongside.
      - P̄ was calibrated by simulation before being adopted (scripts 31 and 31c; averaged over
        5 groupings in the simulations).
    - **λ re-estimated in every simulation**, with three data-generating models
      (`scripts/31_calibration_lambda_free.R`, `scripts/31c_calibration_real_structure.R`).
      - DGM D keeps most of the real residual clustering by flipping the signs of the real
        residual parts per order and per cell.
      - Boundary false-equivalence rates: 3.9–5.9 % for P; 3.9–4.9 % for P̄.
      - A naive cell-level wild bootstrap (`scripts/31b`) destroyed the phylogenetic structure
        (λ̂ ≈ 0.04). It was stopped after 80 runs and is not used.
    - **Earlier noise-scale claim withdrawn** (item 21): the simulations were not noisier than the
      data.
    - **Final genome-source table and reconciliation of superseded numbers**
      (`scripts/32_final_genome_source_table.py`, `results/genome_source_final_table.tsv`,
      `results/superseded_numbers.tsv`).
    - **Composition-adjusted interaction re-run** with P̄ (50 groupings). The models are described as
      robustness checks of one result.
    - **Within-species inference scope** (`replication_within_species/11_inference_scope.py`).
    - **BioSample retrieval rebuilt from scratch without caches**
      (`scripts/23c_biosample_reconstruct.py`). The per-request log replaces the cache-derived log
      of item 21; `scripts/23b_biosample_provenance.py` is superseded.
23. **Post-hoc checks, set 6 (18 Nov 2026).** No new estimates of the hypotheses. All additions
    are audits of the inference.
    - **Formal definition of P̄; SE distribution; estimator family with an intersection–union
      rule; calibration by boundary with a large-K (oracle SE) bracket** (`scripts/34_pbar_audit.py`).
    - **Multi-metric comparison of real and simulated dependence structure**
      (`scripts/33_structure_diagnostics.R`).
    - **Within-species coverage by species type**
      (`replication_within_species/inference_scope_by_type.tsv`).
    - **Wording:**
      - BioSample table rows clarified;
      - stale "89 %" in Limitations corrected to 90.7 %;
      - limitations added on TA versus sequencing depth and on the residual skewness that no
        simulation reproduces.
24. **Wording revision (18 Nov 2026).** Wording only; no new analyses.
    - The intersection–union argument is restated. It rests on P with the original assignment,
      which was fixed before it was run. P̄ is the stricter procedure. The other estimators are
      sensitivity analyses, not calibrated components.
    - Calibration is described as "broadly compatible with 5 %", not as proven.
    - The limited leverage of the within-species analysis is stated explicitly.
    - The abstract conclusion separates the confirmatory result, the exploratory finding and the
      limits of inference.
    - The effect of residual skewness on boundary calibration, and the validity of the decision
      rule, are referred to independent statistical review.
25. **Final wording revision (18 Nov 2026). Analyses frozen.** Wording and documentation only.
    - **Decision rule stated in Methods, Abstract, Results and README:** equivalence requires both
      P and P̄.
    - **Provenance of P documented** from git: defined, seeded and run once in `f8b7719`.
    - **"P̄ is more conservative" replaced** by a statement about the observed data only.
    - **"Dominated by ancestry" replaced** by a quantitative comparison: the phylogenetic share of
      residual variance is 95 %, and P0 accounts for about 0.01 % of the variance in DORM_z.
    - **Student report:** "sleep genes" flagged as a nickname; "publicly" pre-registered corrected
      to "a dated git commit".
