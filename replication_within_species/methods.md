# Methods: within-species replication 

**Status.** This is an **exploratory replication of an exploratory finding**: the MAG/isolate
reversal of the DORM_z ~ P0 slope, which was not pre-registered. Nothing here is confirmatory, and
nothing here changes the pre-registered H1 decision.

## Question

The main analysis found two opposite slopes between species:
- MAG-represented species: DORM_z rises with climate variability (P0);
- isolate-represented species: DORM_z falls with P0.

The audit (`audit/code_review.md` M3) showed what that comparison really is. Only 3 of the 4,919
MAG-represented species have any in-domain isolate genome, so it is "never-cultured vs cultured
species".

The within-species design separates the two explanations:

| Explanation | What it predicts within a species |
|---|---|
| Genome source (assembly or binning) | The reversal reappears between MAG and isolate genomes of the *same* species |
| Lineage (culturability) | Within-species slopes are flat or similar in both genome types |

## Genome selection (`01_select.py`, `05_select_mq.py`)

- **Source.** GTDB r232 bacterial metadata (the same file as `scripts/01`–`04`).
- **Eligibility.** Exactly the main-analysis rules of `04_select_dataset.py`:
  - habitat class and domain, with the majority-domain rule per species;
  - valid coordinates;
  - CheckM2 completeness ≥ 90 % and contamination < 5 % (high quality, HQ).
- **Unit.** Every eligible genome of every species with ≥ 2 eligible genomes in its domain. There is
  **no** one-genome-per-species pick and no isolate bonus.
- **HQ set.** 7,206 genomes, 1,503 species, 74 species with both MAG and isolate genomes. This is
  the complete eligible pool, not a sample. It equals the main `analysis_set` rows with
  `n_genomes_species_domain ≥ 2`.
- **Extension, reported separately.** The same rules with the MIMAG medium-quality threshold
  (completeness ≥ 50 %, contamination < 10 %): 11,788 genomes, 2,716 species, 111 mixed species.
  Site variability for the 4,592 new genomes was computed with the unchanged
  `scripts/02_env_variability.py` (`data/env_variability_mq.tsv.gz`).

## Annotation (`run_annotation.sh`, `run_annotation_list.sh`, `02_count_he.py`)

- The repository's own scripts were run unchanged:
  - `scripts/03_fetch_and_call_genes.py` (NCBI download; pyrodigal 3.7.1; translation table 11, or 25
    for Gracilibacteria);
  - `scripts/05_hmm_search.py` (pyhmmer 0.12.3, the 174-model dormancy library, `dictionary/thresholds.tsv`);
  - `scripts/08_aggregate_modules.py` (module counts).
- `02_count_he.py` re-implements only the ribosomal-protein count of `scripts/06_growth_rate.py`:
  Pfam GA cutoffs, CDS ≥ 240 nt and a multiple of 3. This count is needed for the QC rule.
- **Check.** 131 genomes are shared with the main analysis. All 30 module columns and n_HE are
  identical to the committed values.

## Variables (`03_build_table.py`)

| Variable | Definition |
|---|---|
| DORM_z | Mean of z-scored log1p(spo_strict, n_rpf, n_hib, n_sr), using the **main-analysis** means and SDs (`dormz_params_main.tsv`; reproduces the committed DORM_z to 5e-6) |
| P0 | log T_sd_total, z-scored within domain with the main confirmatory parameters (`p0_params_main.tsv`), so a slope of 0.01 means the same thing as in the paper |
| QC | n_proteins ≥ 500 and n_HE ≥ 20, with NA counted as a failure (the corrected rule; see `audit/` M2) |
| Site | 0.5° cell, `round(lat*2)/2`, as in scripts 16/18/21 |
| Covariates | log10 genome size, CheckM2 completeness |

Single-cell genomes are coded as MAGs (as in the main analysis). `is_sag` flags them for the
sensitivity run.

## Models (`04_models.R`)

- **Responses.** Each model is fitted to DORM_z and to nine module responses (log1p counts:
  hibernation, extra hibernation copies, short HPF, Balon, stringent response, SAS, sporulation,
  Rpf, TA pairs).
- **Primary estimate.** Species fixed effects, standard errors clustered by site (HC1).
  - Species fixed effects remove every between-species difference, including phylogeny.
  - Only genomes of one species sampled at different sites identify the slope.

| Label | Model |
|---|---|
| a_speciesFE_cl{species,site} | y ~ P0 + lgs + comp + is_mag + species |
| a_lmer_randslope | Mundlak: y ~ P0_within + P0_species-mean + covariates + (1 + P0_within ‖ species) |
| a_lmer_spec | As written in the replication plan: y ~ P0 + covariates + (1 + P0 ‖ species) |
| e_lmer_randslope_site | a_lmer_randslope plus a crossed (1 ‖ site) |
| **b_speciesFE_cl{species,site}** | **y ~ P0 × is_mag + lgs + comp + species**: within-species slope in isolates, the interaction, and the MAG slope (delta method) |
| b_e_lmer_randslope_site | Random-slope version of b with a site random intercept |
| c_speciesFE_within_type_clsite | Species fixed effects fitted separately to MAG genomes and to isolate genomes (species with ≥ 2 genomes of that type) |
| d_paired_* | Mixed species only: species-mean MAG − isolate difference (t-test and sign count), and a regression of that difference on ΔP0 + Δcompleteness + Δlog genome size (HC1) |
| d_mixed_speciesFE_clsite | Mixed species only, genome level: y ~ is_mag × P0 + covariates + species; is_mag is the within-species MAG − isolate gap |

**Not estimable.** If a response does not vary within any mixed species (e.g. Balon), the model is
reported as NA.

## Sensitivity (`08_sensitivity_tables.py`)

- **cap20.** At most 20 genomes per species, a random draw with seed 1. The species with the most
  genomes, *Burkholderia ubonensis*, has 250. Capping removes some mixed species, because a draw can
  lose a species' only MAG or isolate.
- **nosag.** SAGs removed.
- **Completeness strata** (`10_completeness_strata.R`). Model b is refitted at completeness ≥ 0,
  90, 95 and 98 %, for DORM_z, extra hibernation copies and the stringent response.
  - This was added after the MQ extension showed that missing single-copy genes in incomplete MAGs
    can create a slope.
  - The strata are nested subsets, not independent replications.

## Limits

- **Identification.** The within-species slope is identified only by species sampled at more than
  one site: 677 of the 1,473 analysed HQ species.
- **Precision.** Within-species variation in P0 (SD 0.42) is about half the total.
- **Mixed species.** Few species have both genome types (71 after QC at HQ), and many of them have
  only 1–2 genomes of one type.
- **Interpretation.** A null within-species interaction is therefore evidence against a *large*
  genome-source effect, not proof of none. The README reports each estimate with its confidence
  interval, compared against the between-species estimate.
