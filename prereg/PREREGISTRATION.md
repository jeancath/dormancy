# Pre-registration — The Dormancy Rule

**Status: FROZEN.** The frozen version is the git commit whose message starts with `PREREG FROZEN`
on branch `ccr-70a28546-5vpeux` of `genmnz/project-melad`. The GitHub push time of that commit is the
registration time. Until then, no association between a predictor and a response had been computed
for confirmatory genomes (`pilot == 0`). Confirmatory genome download, gene calling and HMM annotation
were running in the background when this was frozen; none of that involves a predictor. The only
association results seen were from the 400 held-out pilot species (`results/pilot_report.md`).

The binding list of models is `scripts/models_config.tsv` as committed in this commit. The analysis
code is `scripts/09_models.R`, and the module rules are `scripts/08_aggregate_modules.py`, both also
as committed here. The pipeline may later be changed to fix a bug. Any such change, and anything else
that departs from this document, will be listed in `results/deviations.md`.

## 1. Question and theory
Seed-bank theory (Lennon & Jones 2011, *Nat Rev Microbiol* 9:119) predicts that dormancy pays off
more when favourable periods are rarer or less predictable. Bet-hedging theory (Cohen 1966;
Kussell & Leibler 2005) separates two cases:
- **Predictable** (seasonal) variation favours responsive, sensing-based strategies.
- **Unpredictable** variation favours bet hedging.

If dormancy capacity is under selection at the genome level, lineages from more temporally variable
environments should carry more dormancy machinery than related lineages from stable environments.

## 2. Data
- **Genomes.** GTDB R11-RS232 bacterial genomes with NCBI isolation coordinates, CheckM2
  completeness ≥ 90 % and contamination < 5 %.
- **Habitat.** Rule-based classification of NCBI `isolation_source` (`scripts/01_habitat_classify.py`).
  - **marine**: marine keywords, no deep-water keywords, and a valid OISST cell within ≤ 50 km.
    Also unspecified "water" samples that fall on an ocean cell.
  - **land**: soil and freshwater, on a valid GHCN_CAMS land cell.
  - **Excluded by design**: host-associated, engineered, subsurface, hydrothermal, hypersaline and
    unclassifiable samples.
- **Unit.** One genome per GTDB species, chosen by quality: higher completeness, lower contamination,
  isolate, complete assembly. 7,867 species in total.
- **Confirmatory set.** n = 7,467 (pilot excluded); about 2,392 marine and 5,075 land, before QC.
- **Phylogeny.** The GTDB r232 bac120 tree, pruned to the analysis species.
- **QC exclusions.** Applied before analysis, and only these: < 500 predicted proteins, or < 20
  ribosomal-protein genes detected (gRodon input). No outlier removal on any variable.

## 3. Predictors
All predictors are organism-independent and are computed from 1991–2020 monthly series
(`scripts/02_env_variability.py`).
- **Marine sites:** NOAA OISST v2.1 sea-surface temperature, 0.25°.
- **Land sites:** GHCN_CAMS 2 m air temperature, 0.5°.

The predictors are standardised within domain, so that marine and land values share one scale:
- **P0 (primary):** log of the SD of the 360 monthly temperature values (total variability).
- **Decomposition (H6):**
  - log seasonal SD: the SD of the 12-month climatology (the predictable part).
  - log residual SD: the SD of de-seasonalised, de-trended anomalies (the unpredictable part).
  - lag-1 autocorrelation of the anomalies.
  - Colwell (1974) predictability P, constancy C and contingency M, using 2 °C classes.
- **Secondary predictors:**
  - Marine: ESA OC-CCI v6 log10 chlorophyll-a, 1998–2022 (log total SD and log residual SD).
  - Land: PREC/L precipitation (log CV of the total and residual series).
  - Habitat niche breadth (Hao et al. 2026 *PNAS*, mapped from r207 to r232).
  - Genus-level social niche breadth (von Meijenfeldt et al. 2023).

## 4. Responses
The responses use the dictionary in `dictionary/`, which was validated on 14 reference genomes.
- **SPO:** the number of distinct scored core sporulation genes (0–33). *Sporulator* means Spo0A is
  present and the strict score is ≥ 25.
- **RPF:** the number of Rpf-domain (PF06737) proteins.
- **HIB:** the number of hibernation-factor proteins (long HPF, short HPF/RaiA, RMF, Balon, MRF).
- **SR:** stringent response, counted as long RSH + SAS + SAH proteins.
- **DORM_z (primary composite):** the mean of the z-scores of log1p(SPO), log1p(RPF), log1p(HIB)
  and log1p(SR). Each mechanism gets equal weight.
- **TA:** the number of paired type II toxin–antitoxin loci. This is a separate module and is
  **never** part of DORM.

## 5. Models
Gaussian responses use `phylolm(model = "lambda")`. Binary responses use
`phyloglm(method = "logistic_MPLE")`. All tests are two-sided at α = 0.05. The standard covariates
are domain, log10 genome size, CheckM2 completeness and MAG/isolate status.

| ID | Prediction | Model | Decision rule |
|---|---|---|---|
| **H1 (primary, direct effect)** | DORM_z increases with P0 | DORM_z ~ P0 + covariates | **Supported** if β > 0 and p < 0.05 **and** the 95 % CI overlaps or exceeds the SESOI (0.015). **Practically null** if the 95 % CI lies inside (−0.015, +0.015). **Contradicted** if β < 0 and p < 0.05. |
| H1-total (co-primary) | Same as H1, as a total effect (genome size treated as a mediator) | as H1 without log genome size | same rules |
| PATH | Variability selects smaller genomes | log10 genome size ~ P0 | reported; it is used to interpret H1 against H1-total |
| H1a | Sporulation increases with P0 (within Bacillota) | phyloglm sporulator; log1p(SPO) | BH-FDR across the H1a–d family (q < 0.05) |
| H1b | Rpf increases with P0 (all, and within Actinomycetota) | log1p(RPF) | BH within the family |
| H1c | Hibernation increases with P0 | log1p(HIB) | BH within the family |
| H1d | Stringent response increases with P0 | log1p(SR) | BH within the family |
| **H6 (bet hedging)** | Unpredictable variability (residual SD, contingency M, low Colwell P) predicts DORM_z more strongly than seasonal SD | one model per component | BH across the decomposition family. A difference between the residual and seasonal slopes is assessed with a joint model (DORM_z ~ seasonal + residual + covariates) and a Wald test of the difference. |
| H3 (marine TA) | After genome-size correction, marine genomes are not TA-depleted | log1p(TA) ~ marine + covariates (± mobile-element load) | Supported if the marine coefficient is not significantly negative |
| H3b | TA tracks mobile-element load rather than variability | log1p(TA) ~ P0 + mobile + covariates | reported |
| H4 (streamlining) | Pelagibacterales, Prochlorococcus, Actinomarinales and SAR86 carry fewer dormancy genes than genome size predicts | phylolm DORM_z ~ streamlined + log genome size + covariates, within marine | exploratory |
| H5 (replication of Kearns & Shade 2018) | 16S copy number is positively associated with DORM_z | complete-genome isolates | replication succeeds if β > 0 and p < 0.05 |
| H2 (responsive vs stochastic) | dropped: gene presence cannot identify the switching mode | — | — |

## 6. Sensitivity analyses
These are pre-specified and run on H1:
1. Isolates only.
2. Completeness ≥ 95 %.
3. Species-median environment.
4. Marine only and land only.
5. Excluding Bacillota and Patescibacteriota.
6. With gRodon2 doubling time as a covariate.
7. With mobile-element load as a covariate.
8. DORM_sum instead of DORM_z.
9. With mean temperature as a covariate.
10. MCMCglmm Poisson for each module count, if time allows.

## 7. Power and SESOI
These come from the simulation on the confirmatory tree, using pilot residual parameters (λ = 0.89,
σ² = 0.39).
- **SE(β).** About 0.003 for H1.
- **Power.** 0.89 at β = 0.010 and about 1.00 at β = 0.015.
- **SESOI.** 0.015 DORM_z units per SD of P0, about 0.027 SD of the composite.
- **What a null means.** At this power, a null result can be read as equivalence to zero within ±SESOI.
- **Module tests (H1a–d).** These have lower power, mostly because of clade restriction (n ≈ 330
  Bacillota and ≈ 1,000 Actinomycetota). They are reported with CIs.

## 8. What would falsify the rule
H1 and H1-total are both practically null (CI inside ±0.015) or negative. If that happens, the
conclusion is that, across 7,000+ free-living bacterial species, dormancy gene content does not track
the temperature variability of the isolation site, after phylogeny is controlled. That null result
will be reported as the main finding.
