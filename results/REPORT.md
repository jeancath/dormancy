# The Dormancy Rule: results report

**Status.** Confirmatory analysis complete; analysis date 17 Nov 2026. The pre-registration was
frozen at commit 9fae630 (10:16 UTC), before any confirmatory model was fitted (annotation was running; no predictor was involved).

**Supporting files.**
- `results/decisions.md`: every model with its pre-registered decision.
- `results/deviations.md`: every departure from the plan.
- `results/figures/`: the figures.

## 1. Data

**Species.**
- 7,327 bacterial species, one genome per GTDB r232 species, in the confirmatory set; the 400
  pilot species were held out. 152 genomes failed the pre-registered QC: fewer than 500 proteins,
  or fewer than 20 ribosomal-protein genes.
- 4,994 land species (soil and freshwater) and 2,333 marine species, from 81 phyla.
- 2,835 isolate genomes and 4,492 MAGs, from 1,680 distinct 0.5° sites (Fig. 1).

**Predictor (P0).** P0 is the log SD of 30 years of monthly temperature at the isolation site,
standardised within each domain.
- Marine sites use SST; the median SD is 3.0 °C.
- Land sites use air temperature; the median SD is 6.9 °C.

**Responses.** The dormancy modules come from the validated profile-HMM dictionary (`dictionary/`).
Independent phenotype calibration is in `validation/`:
- sporulation AUC = 0.92;
- the pre-registered rule has sensitivity 0.90 and specificity 0.78, with no positives outside
  Bacillota;
- Rpf counts match the literature exactly in 23 of 24 Actinomycetota.

On the full set the sanity checks agree. All 15 spore-forming genera tested are called sporulators,
0 of 10 non-sporulating genera are, and Spo0A is never called outside Bacillota
(`results/genus_sanity.tsv`).

## 2. Pre-registered (confirmatory) results

| Hypothesis | Estimate [95 % CI] | p | Pre-registered decision |
|---|---|---|---|
| **H1** DORM_z ~ P0, direct effect | **+0.0065 [+0.0016, +0.0113]** per SD | 0.009 | **Practically null.** The effect is statistically non-zero, but the whole CI lies inside the SESOI of ±0.015. |
| H1-total (genome size not adjusted) | +0.0057 [+0.0007, +0.0106] | 0.025 | Practically null |
| PATH: genome size ~ P0 | +0.0004 [−0.0010, +0.0019] | 0.56 | No path. Within a domain, variability does not select smaller genomes. |
| H1a sporulator (Bacillota, n = 394; logit) | +0.024 [−0.16, +0.21] | 0.80 | Not supported (q = 0.80) |
| H1a sporulation genes (log1p) | −0.0037 [−0.0084, +0.0010] | 0.12 | Not supported (q = 0.24) |
| H1b Rpf, all species | −0.0014 [−0.0033, +0.0005] | 0.16 | Not supported (q = 0.24) |
| H1b Rpf, within Actinomycetota (n = 1,029) | −0.0053 [−0.0174, +0.0067] | 0.39 | Not supported (q = 0.46) |
| H1c hibernation | +0.0051 [+0.0003, +0.0099] | 0.037 | Not significant after BH (q = 0.12) |
| H1d stringent response | +0.0049 [+0.0002, +0.0095] | 0.040 | Not significant after BH (q = 0.12) |
| H6 seasonal SD vs residual SD | +0.0066 vs +0.0061; Wald difference (resid − seasonal) −0.002 ± 0.008 | 0.76 | **Not supported.** The two components correlate at r = 0.82 and cannot be separated. No Colwell index survives BH. |
| **H3** marine TA deficit after genome-size correction | **−0.231 log1p units (≈ 21 % fewer TA pairs) [−0.28, −0.18]** | 2 × 10⁻¹⁹ | **Not supported.** Marine genomes *are* TA-depleted. The deficit is unchanged after adding mobile-element load (−0.224). |
| H3b TA ~ P0 + mobile elements | +0.021 [+0.006, +0.037] | 0.007 | Reported. TA content rises weakly with variability. |
| H4 streamlined marine clades | −0.10 [−0.29, +0.09] | 0.29 | Exploratory; null |
| H5 16S copy number → DORM_z (complete isolates, n = 864) | +0.035 [−0.003, +0.073] | 0.068 | Replication failed (p > 0.05). The direction matches Kearns & Shade 2018. |

**Sensitivity analyses (pre-registered).** The estimates lie between +0.0046 and +0.0072 in every
variant except one:
- species-median environment;
- completeness ≥ 95 %;
- marine only and land only;
- excluding Bacillota and CPR;
- adding growth rate, mobile elements or mean temperature as covariates;
- DORM_sum instead of DORM_z.

All of these CIs fall inside ±SESOI except marine only, whose upper limit of +0.0155 makes it
"inconclusive". DORM_sum is on a log1p scale, where the SESOI does not apply. The one variant
that disagrees is **isolates only, where the sign reverses**: −0.0079 [−0.0151, −0.0008],
p = 0.03.

**Secondary predictors.**
- Chlorophyll variability (marine): null (+0.0003).
- Precipitation CV (land): +0.0078 [+0.0016, +0.014], p = 0.014. This is also below the SESOI.
- Habitat niche breadth and social niche breadth: null.

**Confirmatory conclusion.** Across 7,327 species, phylogeny explains most of the variation in
dormancy gene content (Pagel's λ = 0.94). After phylogeny is controlled, isolation-site temperature
variability predicts at most a trivially small amount: under 0.015 SD of the composite per SD of
variability, about 1/40 of an SD of DORM_z. This holds for every dormancy mechanism. **The genome-level
seed-bank prediction is not supported at a biologically meaningful effect size.** That is the
pre-specified falsification outcome (§8 of the pre-registration).

## 3. Exploratory results (not pre-registered; `results/deviations.md` §5–9)

### 3.1 The small pooled effect hides opposite slopes in MAGs and isolates
- **MAGs (n = 4,492):** +0.018 [+0.011, +0.025], p = 1 × 10⁻⁷.
- **Isolates (n = 2,835):** −0.008 [−0.015, −0.001], p = 0.03.
- **Interaction** (P0 × MAG): +0.024, p = 2 × 10⁻⁶ (Fig. 5).

**Robustness of the split.** It survives every check we could devise:

| Check | MAG slope | Isolate slope |
|---|---|---|
| Phylolm, all genomes of that type | +0.018 | −0.008 |
| Site-grouped jackknife SE (100 groups of sites) | p = 0.003 | p = 0.054 |
| One random species per 0.5° site (200 resamples) | median +0.023; 99 % of resamples > 0 | median −0.008; 0 % > 0 |
| Leave out each of the 15 largest MAG sites | +0.017 to +0.023 | — |
| Drop all 15 largest MAG sites | +0.021 (p = 2 × 10⁻⁴) | — |
| By domain | marine +0.024, land +0.013 | marine −0.016, land −0.004 |
| Within the 39 families that contain both genome types | interaction +0.023, p = 0.002 | |
| With assembly N50 and contamination as covariates | +0.018 | — |
| With soil/sediment covariates | +0.018 | — |
| With gRodon doubling time as a covariate | interaction unchanged, +0.024 | |

**Which modules carry it.** In MAGs the effect comes from ribosome hibernation (+0.017,
p = 1 × 10⁻⁶) and the stringent response (+0.009, p = 0.006). Sporulation and Rpf contribute
nothing. In isolates hibernation runs the other way (−0.008, p = 0.017).

**Ruling out a gene-recovery artifact.** If high-variability MAG studies simply recovered genes
better than CheckM2 completeness reflects, near-universal genes should also rise with P0 in MAGs.
They do not:
- ribosomal-protein HE genes fall slightly (−0.005, p = 1 × 10⁻⁴);
- long RSH is null (+0.003, p = 0.12);
- coding density is null.

Growth rate does not explain the split either. Predicted doubling time does not track P0 in MAGs
(p = 0.14). It is marginal in isolates (+0.013, p = 0.06, slower growth at variable sites).

### 3.2 Site pseudoreplication
The 4,492 MAGs come from 368 sites, and the 10 largest sites hold 42 % of them. The pooled H1
estimate is therefore fragile:
- The site-grouped jackknife SE (0.0038) is 1.5 × the phylolm SE, which gives p = 0.089.
- With one species per site, the pooled slope is **negative**: median −0.005, positive in 1 % of
  200 resamples. Capping a site at one species makes isolates (1,399 sites) dominate.
- With up to three species per site the median is −0.003. In that resampling MAGs were positive in
  100 % of resamples and isolates in 0 %.

The confirmatory decision (practically null) does not change, but the small positive estimate
should not be read as evidence for the hypothesis.

### 3.3 Predictable vs unpredictable variability
**Re-parameterisation.** Following Marshall & Burgess (2015), we used three predictors: log total
SD, the seasonal fraction of variance, and the residual lag-1 autocorrelation. Neither the seasonal
fraction (+0.0009, p = 0.84) nor the autocorrelation (+0.0012, p = 0.68) adds anything once total
SD is in the model.

**Land.** Total and seasonal SD correlate at r = 1.00, so the decomposition is not identifiable.

**Ocean.** Residual and seasonal variability are separable there (r = 0.24 with the seasonal
fraction). The residual slope, +0.006 (p = 0.34), exceeds the seasonal slope, +0.003 (p = 0.63),
but neither is significant.

The predictable/unpredictable question cannot be answered with surface temperature at this sample
size.

### 3.4 The marine TA paradox (H3)
- **Overall.** The marine deficit of ≈ 21 % is robust to mobile-element load, assembly N50 and
  contamination, and to excluding plasmid-borne loci (−0.232).
- **MAGs:** −0.33, p = 7 × 10⁻²⁶.
- **Isolates:** −0.016, p = 0.70.
- **Complete isolate genomes (n = 865):** −0.155, p = 0.03.

**Why it is concentrated in MAGs.** Short-read binning recovers only 38–44 % of genomic islands and
1–29 % of plasmid DNA (Maguire et al. 2020). Marine communities are typically rich in co-occurring
strains, which makes binning harder. Part of the marine TA deficit reported in metagenomes is
therefore probably a recovery artifact.

**But it is not all artifact.** A real deficit of about 14 % remains in complete marine isolates.
That is consistent with Fucich & Chen (2020) for picocyanobacteria.

TA content also rises slightly with P0 (H3b), while it remains a separate module from dormancy.

### 3.5 Other exploratory signals
- **Cyanobacteriota:** the one phylum with a positive within-clade slope (+0.036, n = 311,
  p = 0.004, uncorrected).
- **Balon:** presence rises with P0 (logit +0.071, p = 0.006). Balon was discovered in a polar
  *Psychrobacter* (Helena-Bueno et al. 2024).
- Both are hypotheses for follow-up, not findings.


### 3.6 Follow-up tests of the post-hoc explanations (`results/followup_models.tsv`)
- **Coordinate precision does not explain the negative isolate slope.** The slope is more negative
  in isolates with precise (≥ 3-decimal) coordinates (−0.010, p = 0.05) than in those with coarse
  coordinates (−0.001). MAGs are positive in both groups (+0.014 and +0.042).
- **The reversal sits in hibernation paralogues.**
  - Extra hibernation-factor copies: MAGs +0.013 (p = 2 × 10⁻⁵), isolates −0.011 (p = 0.0015).
  - Presence of hibernation systems: MAGs +0.009, isolates −0.003 (ns).
  - The genes involved are short HPF/RaiA (MAGs +0.013, isolates −0.006) and Balon (MAGs +0.008).
    Small alarmone synthetases also rise in MAGs (+0.0075).
- **Contamination does not explain the MAG effect:**
  - MAGs with < 1 % contamination: DORM_z +0.022, hibernation paralogues +0.017.
  - Contamination as a covariate leaves the paralogue slope unchanged (+0.0135).
  - In near-complete, clean MAGs (≥ 98 % completeness, < 1 % contamination; n = 906) the slope is
    +0.007 (p = 0.35). It is attenuated, but the CI still includes +0.018.
- **The "collider" model is not independent evidence.** The isolate-status model
  (is_isolate ~ P0 × DORM_z; logit −0.56, p = 2 × 10⁻¹¹) is the mathematical mirror of the
  MAG × P0 interaction, so it is a reframing only.
- **Balon and cold seas:**
  - In marine *isolates*, Balon presence falls with mean SST (logit −0.34 per SD, p = 4 × 10⁻⁴):
    - by threshold: 21 % of isolates below 5 °C (n = 80) and 8 % above 20 °C (n = 379);
    - across SST octiles the curve is not monotonic: 17, 12, 10, 9, 4, 12, 5 and 12 %;
    - Helena-Bueno et al. (2024) ran no temperature test, so this association is new.
  - The association is absent in marine MAGs (p = 0.13) and reversed on land (+0.16).
  - Status: a hypothesis, not a result.
- **Cyanobacteriota:** the positive within-phylum slope is carried by hibernation genes
  (+0.042, p = 0.003) and comes from MAGs only (MAGs +0.039, isolates −0.024).
- **Pilot re-analysis** (n = 388; `results/pilot_replication_models.tsv`):
  - The pilot was the *discovery* sample for the marine TA deficit (−0.27 before the freeze), so
    it does not replicate the confirmatory result; the confirmatory set replicated the pilot.
  - The MAG × P0 interaction cannot be assessed at this size: +0.013 ± 0.034, which neither
    confirms nor refutes it.

### 3.7 Independent audit (`audit/README.md`)
- **Overall.** No critical issue was found. The pre-registration is intact, and the confirmatory
  conclusion holds.
- **Fixes applied:**
  - MCMCglmm refit on GTDB branch lengths (§6).
  - QC rule for genomes without a ribosomal-gene count; the rerun changes nothing.
  - Pilot wording.
  - Balon numbers.
  - Marine TA citation.
  - Minor wording and rounding errors.
- **Alternative refit:** nested taxonomy plus site random effects (lme4) gives pooled +0.0049
  [−0.0015, +0.0114] and interaction +0.028.
- **Strain-heterogeneity tests** do not support the "MAG population consensus" explanation:
  - CheckM strain heterogeneity falls with P0 in MAGs.
  - The MAG slope is equal in MAGs with zero heterogeneity (+0.024).
  - In the ocean, the "reversal" comes from isolates losing hibernation copies, not MAGs gaining
    them.
- **The lineage confound:** only 3 of 4,919 MAG-represented species have an in-domain isolate
  genome. "MAG vs isolate" therefore means never-cultured vs cultured lineages.

### 3.8 Within-species replication (`replication_within_species/README.md`)
**Data.** 7,098 high-quality genomes from 1,473 species with ≥ 2 genomes, 71 of them with both
genome types. Species fixed effects, site-clustered 95 % CIs.

| DORM_z ~ P0 | Within species | Within, ≥ 95 % complete | Between species |
|---|---|---|---|
| All genomes | +0.001 [−0.006, +0.009] | — | +0.0065 |
| MAG genomes | −0.008 [−0.023, +0.006] | +0.009 [−0.001, +0.019] | +0.018 |
| Isolate genomes | +0.003 [−0.006, +0.011] | +0.004 [−0.005, +0.012] | −0.008 |
| MAG × P0 interaction | −0.011 [−0.028, +0.005] | +0.005 [−0.007, +0.018] | +0.024 |

- **Composite reversal.** Within species it is absent. The interaction and isolate CIs exclude the
  between-species values at every completeness threshold, so **the reversal is a between-lineage
  pattern, not a genome-source effect**.
- **MAG slope.** It depends on completeness: missing single-copy RSH genes in less complete MAGs
  pull it negative. At ≥ 95 % it does not cleanly exclude +0.018.
- **Extra hibernation copies.** This is the one MAG-side within-species signal: +0.006 to +0.017
  across thresholds (between species +0.013). It is not robust to multiplicity (q = 0.26).
- **Medium-quality extension** (11,050 genomes, 2,515 species): same interaction result. Its MAG
  slope is negative (−0.010), a completeness artefact that vanishes at ≥ 95 %.


### 3.9 Checks added later (18 Nov 2026)
- **Formal equivalence (TOST against ±0.015; `results/equivalence_tests.tsv`).**
  - H1 direct: 90 % CI [0.0024, 0.0105], p_TOST = 0.0003.
  - H1 total: p_TOST = 0.0001.
  - With the site-grouped jackknife SE: p_TOST = 0.012.
  - Within species: p_TOST = 0.0002.
  - MCMCglmm phylogeny + site: the 95 % HPD lies inside ±0.015.
  - All pre-registered sensitivity variants are equivalent. The exceptions are the two
    niche-breadth predictors, which are inconclusive at their small n.
  - **Conclusion:** the null is evidence of absence of a meaningful effect, not merely a failure to
    detect one.
- **Marine TA deficit by genome quality** (`results/review_models.tsv`, `results/ta_review_models.tsv`).
  - **The tiers differ.** Marine × high-quality interaction +0.141 (p = 1 × 10⁻⁴); marine × MAG
    interaction −0.296 (p = 1 × 10⁻¹²).
  - **Clean genomes** (≥ 98 % complete, < 1 % contaminated): MAGs −26 % (p = 2 × 10⁻⁶); isolates
    −4 % (p = 0.37).
  - **Not a pairing artefact.** Single toxins, single antitoxins and orphan toxins show the same
    deficit as pairs.
  - **Conclusion.** The 21 % is not a stable estimate. The deficit is concentrated in MAG-represented
    lineages, and among cultured isolates the evidence is weak (−4 % to −14 % depending on the
    subset).
- **Collection dates (`data_small/final/collection_dates.tsv`, `results/collection_date_models.tsv`).**
  - **Retrieval bug, now fixed.** 121 DDBJ/EBI BioSamples had been fetched by accession string, and
    NCBI resolved them to unrelated records. Those records were discarded unused. The 121 were then
    fetched by UID, with the accession verified.
  - All 7,714 BioSamples are now retrieved.
  - 97.4 % of confirmatory genomes are dated (7,140 / 7,327): median 2016, IQR 2013–2018.
  - 90.7 % of all confirmatory genomes, or 93.1 % of the dated ones, fall entirely within the
    1991–2020 climate window.
  - Restricting to those samples:
    - H1 +0.0044 [−0.0010, +0.0098];
    - interaction +0.024;
    - TA −0.24.
  - Adding collection year as a covariate leaves H1 at +0.0054.

- **What P0 measures** (`results/variability_component_correlations.tsv`).
  - About 94–95 % of monthly variance is seasonal, and P0 correlates r = 1.00 with seasonal SD.
    P0 is therefore regional seasonal amplitude.
  - Interannual SD is a separate axis (r = 0.55–0.61). Given seasonality it is null: +0.0029,
    p = 0.30 (`results/variability_review_models.tsv`).
- **SESOI in genes** (`results/sesoi_translation.tsv`). Across the 5th–95th percentile of sites:
  - about 2–3.5 % more genes per module if all modules shift together;
  - about 7–15 % if one module carries the effect;
  - for scale, the SD of phylum means is about 30 × the SESOI.
- **Dictionary, clade-aware validation** (`validation/clade_holdout*.tsv`).
  - Genus-clustered bootstrap: AUC 0.92 (0.88–0.95).
  - Leave-one-order-out accuracy: 0.87 for a threshold chosen on the other orders, 0.82 for the
    fixed rule.
- **TOST operating characteristics** (`results/tost_power_simulation.tsv`) and **spatial blocks of
  0.5–5°** (`results/spatial_blocks_*.tsv`): see the paper's Methods and Results.

### 3.10 calibration, genome-source slopes, margins (`results/group_equivalence.tsv`)
- **Primary equivalence procedure P.**
  - The procedure: λ fixed; jackknife over 20 groups of 0.5° cells; TOST on t₁₉.
  - On the real data, with one random grouping: 90 % CI [0.0026, 0.0103], p_TOST = 0.0005.
    **Superseded in §3.12:** that grouping was unusually favourable.
  - Calibration by simulation (`results/tost_calibration_summary.tsv`): false-equivalence rate
    5.0 % at −0.015 and 4.5 % at +0.015 (5,000 replicates each); 4.75 % combined, Monte Carlo
    95 % CI 4.3–5.2 %. P(equivalence) is 97.5 % at 0, 83.8 % at 0.005 and 37.1 % at 0.010. The
    jackknife SE (0.0036–0.0037) matches the simulated slope SD (0.0035). P is calibrated.
  - Note: the R summary written by the script was wrong (two columns named `beta`); the summary
    was recomputed by column position in Python, with Clopper–Pearson CIs.
- **The earlier p_TOST values differ only by SE estimator.**
  - 0.012: jackknife over 100 groups, centred 0.5° cells.
  - 0.007: jackknife over 100 groups, grid-aligned 0.5° cells.
  - 0.0003: model-based SE.
- **Margins.**
  - The pooled slope is equivalent at ±0.015 and ±0.020.
  - At ±0.010 it is equivalent only under MCMCglmm and within species.
- **Slopes by genome source.**
  - MAG-represented species: +0.018, 90 % CI 0.008–0.028 (jackknife). Not equivalent; at or above
    the SESOI.
  - Isolate-represented species: −0.008; equivalence is method-dependent at ±0.015.
  - Interaction +0.024 (phylolm) and +0.026 (MCMCglmm); within species −0.011 [−0.028, +0.005].
  - **Reading:** the pooled near-zero slope averages over opposite between-species signals.
- **Wording now in the paper.** Title, abstract, conclusions and the Discussion heading were changed
  from "tracks ancestry, not variability" to statements about the pooled average below the SESOI.

### 3.11 SE levels, one-sided tests, composition-adjusted interaction
> **Partly superseded by §3.12.** The "noisier simulation" explanation and the real-noise power
> figures below are withdrawn. The single-grouping MAG and isolate tests are replaced by P̄.
- **Three different SEs.**
  - Procedure P on the real data: SE 0.0022. This SE gives the reported 90 % CI
    [0.0026, 0.0103] and p_TOST = 0.0005.
  - Mean jackknife SE across simulations: 0.0036–0.0037.
  - SD of the simulated slopes: 0.0035.
  - The simulation was noisier than the real data, because a site effect was added on top of a
    residual variance that already included site variation.
- **The noise level does not change the boundary Type I rate** (linearity; `scripts/29b_linearity_check.R`).
  *Withdrawn in §3.12:* the real-noise power figures derived from the 0.0022 SE.
- **MAG-represented species, procedure P** (`results/group_equivalence.tsv`):
  - slope +0.0175, 95 % CI [0.006, 0.029];
  - p = 0.005 two-sided, p = 0.003 for β > 0, and p = 0.33 for β > 0.015.
  - So this is a positive association, but the data do not show that the slope exceeds the SESOI.
- **Isolate-represented species, procedure P:** −0.009, 90 % CI [−0.016, −0.001],
  p_TOST = 0.090, two-sided p = 0.055.
- **Interaction adjusted for composition** (`results/interaction_composition.tsv`). Each variable
  enters as a main effect and × P0.

  | Model | Interaction | 95 % CI, P̄ (50 groupings) |
  |---|---|---|
  | Unadjusted | +0.024 | [0.009, 0.039] |
  | Habitat | +0.027 | [0.015, 0.040] |
  | Phylum | +0.026 | [0.010, 0.041] |
  | Latitude and mean temperature | +0.024 | [0.008, 0.039] |
  | All of these | +0.028 | [0.010, 0.045] |
  | Common support, 18 strata | +0.025 | [0.007, 0.043] |

- **Within species:** −0.011 [−0.028, +0.005]. This is not evidence of absence: p_TOST = 0.33 at
  ±0.015.
- **Provenance.** `data_small/final/biosample_retrieval_log.tsv` has one row per BioSample record:
  - 7,593 records fetched by accession;
  - 121 fetched by UID with the accession verified;
  - 115 misresolved returns, discarded.

### 3.12 grouping sensitivity, λ re-estimation, final genome-source table
- **How the 20 jackknife groups are built.**
  - The 1,642 occupied 0.5° cells are assigned at random to 20 groups of about 82 cells, which are
    scattered, not contiguous.
  - The deletion unit is the whole cell.
- **Grouping sensitivity of P on the real data** (`results/grouping_sensitivity_summary.tsv`):

  | Variant | SE | p_TOST | Equivalent |
  |---|---|---|---|
  | Original grouping | 0.0022 (rank 4 of 200) | 0.0005 | yes |
  | 200 random groupings | median 0.0034 (5–95 %: 0.0025–0.0043) | median 0.011 | 199/200 |
  | Delete-one-cell | 0.0035 | 0.007 | yes |
  | 20 contiguous (k-means) groupings | median 0.0035 | median 0.012 | 20/20 |
  | λ re-estimated per replicate | 0.0023 | 0.0007 | yes |

- **Correction of an earlier statement.**
  - The typical real-data SE (0.0035) equals the simulated slope SD (0.0035), so the simulations
    were not noisier than the data.
  - The earlier "noisier simulation" explanation and the 69–70 % power figure are withdrawn.
  - Power at ±0.010 is 37 %.
- **Final procedure P̄.** The jackknife variance is averaged over 200 random groupings, with t₁₉.
  - Pooled slope: 90 % CI [0.0004, 0.0125], p_TOST = 0.012.
  - Two-sided p against zero: 0.08.
- **Calibration with λ re-estimated** (`results/calibration_lambda_free_summary.tsv`,
  `results/calibration_real_structure_summary.tsv`). False-equivalence rate at the boundaries:

  | DGM | P, one grouping | P̄, averaged |
  |---|---|---|
  | A, original (λ fixed; 10,000 runs, first calibration) | 4.8 % (4.3–5.2 %; n = 10,000) | — |
  | A, original, λ re-estimated | 3.9 % (2.8–5.3 %; n = 1,000) | — |
  | B, single fitted model, λ re-estimated | 5.9 % (4.5–7.5 %; n = 1,000) | 4.9 % (3.6–6.4 %) |
  | D, real residual structure, λ re-estimated | 4.2 % (3.0–5.6 %; n = 1,000) | 3.9 % (2.8–5.3 %) |

- **Variance structure.**
  - DGM B matches λ̂ (0.93 vs 0.94) and the SE (0.0037 vs 0.0035), but not the between-cell share of
    residual variance (0.12 vs 0.50).
  - The real excess clustering is phylogenetic: share 0.45 in the phylogenetic part of the
    residuals, against 0.06 in the non-phylogenetic part.
  - DGM D (order- and cell-level sign flips of the real residual parts) reaches 0.36, against 0.41
    for the same model form.
  - A naive cell-level wild bootstrap destroyed the phylogenetic structure (λ̂ ≈ 0.04) and was not
    used.

- **Final genome-source table:** `results/genome_source_final_table.tsv`. Superseded numbers are
  reconciled in `results/superseded_numbers.tsv`.
  - MAG slope: +0.017, 95 % CI [0.004, 0.031]; p = 0.013 two-sided, p = 0.006 for β > 0,
    p = 0.35 for β > 0.015.
  - Isolate slope: −0.009, 90 % CI [−0.016, −0.002], p_TOST = 0.084.
  - Interaction: +0.024, 95 % CI [0.009, 0.039] (P̄, 50 groupings); robust to composition adjustment (+0.024 to +0.028).
- **Within-species scope** (`replication_within_species/inference_scope.tsv`).
  - 700 of 1,473 species vary in P0 within the species.
  - 238 species (1,214 MAGs) identify the MAG slope; 431 species (3,295 isolates) identify the
    isolate slope.
  - Only 2 of the 71 mixed species vary within both sources.
  - MAG-only species supply 11 % of the within-species P0 variation.
- **BioSample retrieval rebuilt from scratch** (`results/biosample_retrieval_flow.tsv`,
  `scripts/23c_biosample_reconstruct.py`):

  | Step | Count |
  |---|---|
  | Requests | 7,714 |
  | Step 1 records returned | 7,708 |
  | — matched | 7,593 |
  | — misresolved and discarded (each attributed to its request) | 115 |
  | Requests unmatched after step 1 (115 got a wrong record, 6 got nothing) | 121 |
  | Step 2 verified | 121 |
  | Linked | 7,714 |
  | Missing | 0 |

  The raw dates are identical to the original run for all 7,714.

### 3.13 audit of P̄, dependence diagnostics, coverage, provenance wording
- **P̄ defined formally** (`scripts/34_pbar_audit.py`).
  - SĒ = √((1/K) Σₖ Vₖ): variances are averaged, then the root is taken. K = 200.
  - Averaging SEs instead gives 0.00344 vs 0.00348 (p_TOST 0.011 vs 0.012).
  - SE over 200 groupings: 0.0020–0.0051, median 0.0034, CV 16 % (`results/pbar_se_distribution.tsv`).
- **Estimator family and intersection–union rule** (`results/pbar_estimator_family.tsv`).
  - Every non-superseded estimator declares equivalence; the most conservative (5° blocks) gives
    p_TOST = 0.015.
  - Revised. The calibrated bound comes from P with the original assignment, which was
    fixed before it was run and declared equivalence (p_TOST = 0.0005).
  - Requiring both P and the stricter P̄ (p_TOST = 0.012) keeps the Type I rate at or below P's.
  - The other estimators were assembled during review. Some are uncalibrated (phylolm) or give no
    TOST p-value (MCMCglmm), so they are sensitivity analyses only.
  - Whether the rule is appropriate is left for independent statistical review.
- **Calibration by boundary** (`results/pbar_calibration_by_boundary.tsv`).
  - Rates are 3.4–6.2 %. No Monte Carlo CI lies above 5 %, and the oracle-SE limit lies at or
    below 5.2 %.
  - The two boundaries differ by at most 1.2 points.
- **Dependence diagnostics** (`results/structure_diagnostics_summary.tsv`). B and D bracket the
  real between-cell and between-order clustering.
  - Neighbouring-cell correlation: real 0.16, within the simulated ranges but above both means.
  - Real residual skewness is 0.88. Neither DGM reproduces it, and D is symmetric by construction.
    This is listed as a limitation.
- **Within-species coverage by species type** (`replication_within_species/inference_scope_by_type.tsv`).
  This is now in Methods and Results.
- **BioSample wording.**
  - The 6 requests with nothing returned in step 1 were all recovered in step 2. "Not retrieved"
    (0) is distinguished from "no usable date" (187).
  - A stale "89 %" in the paper's Limitations was corrected to 90.7 %.
  - No "97.2 %" appears in any version of the documents.

## 4. Interpretation

1. **The comparative seed-bank prediction fails at genome level.**
   - **What was tested.** Dormancy gene repertoires are phylogenetically conserved (λ ≈ 0.94).
     Across species, how many dormancy genes a lineage carries is not tuned to the temperature
     variability of the place its genome was sampled.
   - **What this does not rule out.** Dormancy may still be *expressed* more in variable
     environments, because regulation and expression are invisible to gene counts. Selection may
     also act at the level of whole clades over long timescales.
2. **The MAG–isolate split is a between-lineage pattern (never-cultured vs cultured species).**
   - **Between species** it survives every check: site clustering, shared families, assembly
     quality, gene recovery, coordinate precision, growth rate, and phylogenetic and site mixed
     models.
   - **Within species** it vanishes (§3.8). So it is not produced by the assembly method within a
     species. Between species, this design cannot separate a biological difference from
     differences in taxonomic, habitat and sampling composition (§3.10).
   - **"Population consensus" is dropped.** It was not supported (§3.7).
   - **Coordinate quality is excluded.** The isolate decline is strongest in precisely
     georeferenced isolates.
   - **Remaining explanations:**
     - **(a) Culturing selection.** Cultured lineages from variable sites may be a biased subset
       that relies less on hibernation.
     - **(b) Lineage-level adaptation fixed per species.**
   - **One open MAG-side signal:** extra hibernation copies within species. It needs long-read or
     single-cell genomes to resolve.
   - **Novelty.** Systematic MAG–isolate functional differences are known (Gupta & Almeida 2025).
     What is new is the sign reversal of an *environment–trait slope* between never-cultured and
     cultured lineages.
3. **Methodological lesson.** Global genome–environment analyses that pool MAGs with isolates, and
   that treat many species from one sample as independent environmental replicates, can produce
   small, "significant" pooled effects that do not generalise. We recommend:
   - fitting genome source as a moderator;
   - site-blocked inference (jackknife, site caps, or phylogenetic mixed models with a site
     random effect).

## 5. Limitations
- **Isolation coordinates** are taken from BioSample metadata and are not verified.
- **Surface temperature** is not microhabitat temperature, especially in soil and sediment.
- **Gene presence is not expression**, and the dictionary misses divergent sporulation genes in
  Erysipelotrichia.
- **Module tests within clades** (Bacillota n = 394) have low power.
- **The MAG/isolate analyses are exploratory** and were designed after seeing the preliminary
  pooled result. The within-species replication used independent genomes, but only 71 species
  carried both genome types.

## 6. MCMCglmm models
These use the sparse all-node inverse on **GTDB branch lengths** (audit M1; `results/deviations.md`
items 11–12). They supersede the earlier Grafen-length fits, and all chains have ESS ≥ 800.

**Pre-registered sensitivity analysis 10: Poisson model per module** (`results/confirmatory_mcmcglmm.tsv`;
log-rate per SD of P0; the sporulation module was not included):

| Module | Posterior mean [95 % HPD] | pMCMC |
|---|---|---|
| Rpf | −0.004 [−0.035, +0.023] | 0.848 |
| Hibernation | +0.017 [+0.003, +0.031] | 0.001 |
| Stringent response | +0.010 [+0.000, +0.021] | 0.058 |

**Gaussian DORM_z with phylogenetic and site random effects** (`results/siteMCMC_mcmcglmm_site.tsv`):
P0 slope +0.0039 [−0.0027, +0.0098] (pMCMC = 0.20). With phylogeny and site modelled together, the pooled slope is not
distinguishable from zero.

**Exploratory: MAG × P0 with phylogenetic and site random effects**
(`results/exploratory_mcmc_site_interaction.tsv`):

| Term | Posterior mean [95 % HPD] |
|---|---|
| Interaction | +0.026 [+0.015, +0.038] (pMCMC = 0.001) |
| MAGs | +0.022 [+0.012, +0.032] |
| Isolates | −0.009 [−0.017, −0.002] |

`results/siteMCMC_smoketest_G2R2_*` are smoke-test outputs (G = 2, R = 2) of the jackknife part of
script 16, not results.
