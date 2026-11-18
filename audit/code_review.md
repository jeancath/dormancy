# Code review (audit task 3)

**Scope.** The review covers:
- `scripts/08_aggregate_modules.py`, `10_build_table.py`, `09_models.R`, `16_cluster_robust.R`
  and `18_mag_isolate.R`, as requested;
- `01`, `04`, `14`, `21` and the config files, as far as they feed into those scripts.

**Method.** I read every line and checked the behaviour against the data in
`data_small/final/analysis_table.tsv.gz`.

**Severity scale:**
- **critical:** changes a reported conclusion;
- **major:** changes a reported number materially, or makes an analysis something other than
  what is described;
- **minor:** small numerical effect, or misleading documentation.

**No critical bug found.**
- The confirmatory pipeline does what the pre-registration says.
- An independent model family (nested taxonomic random effects, `audit/refit/`) reproduces the
  confirmatory H1 to the third decimal (+0.0066 vs +0.0065).

## Requested checks

| Check | Result |
|---|---|
| P0 z-scored within domain? | **Yes.** `09_models.R:27` (`zwithin(…, tab$domain)`), and identically in `16:20`, `18:11`, `14:10`, `21:11`. It is computed on the confirmatory subset after the pilot is dropped. |
| Pilot excluded? | **Yes**, in every model script: `09:12` (`pilot == 0`), and `16`, `18`, `14`, `21` (`tab$pilot == 0`). Pilot genomes *are* included in the DORM_z standardisation (see m3). |
| Filters in `models_config.tsv` / `exploratory_config.tsv` / `followup_config.tsv` | All 39 + 54 + 34 filters were checked against their ids. Each sample size matches the filter (e.g. H1a n = 394 Bacillota, H5 n = 864 complete isolates, F_MAG_DORM_comp98 n = 906). None is mis-specified. |
| Coefficient extraction | Correct. `09:53` uses `co[r$predictor, ·]` with the term name, including the interaction term `T_logsd_total:is_mag`. The p column is the last column (`p.value`) for both phylolm and phyloglm. `18:28` uses column 4 = p.value. |
| TA pairing | It matches the dictionary rule, with caveats (m5). |
| QC rules | One implementation error (M2). |

## Findings

### Major

**M1. Every MCMCglmm model uses Grafen branch lengths, not the GTDB tree.**
- Where: `scripts/14_extra_tests.R:29`, `16_cluster_robust.R:66`, `21_mcmc_site_interaction.R:22`
  (`compute.brlen(ph, method = "Grafen")`).
- What it does: it throws away the GTDB bac120 branch lengths and imposes Grafen's arbitrary
  ultrametric lengths.
- Why it matters:
  - The phylolm models, including the confirmatory H1, use the real branch lengths. The MCMCglmm
    results, including pre-registered sensitivity 10 and the "phylogeny + site" models now in the
    paper, therefore assume a different covariance structure.
  - The comment ("ultrametric version required by inverseA") is wrong: `inverseA` accepts
    non-ultrametric trees.
  - Deviation 11 calls the sparse inverse "the same model". It is the same as the tip-only inverse,
    but not the same as the phylolm models.
- Symptom: the posterior animal variance is about 3.95 against a residual of 0.04 (100:1). That is
  implausible on the DORM_z scale, whose total variance is 0.33, and is a sign of the rescaling.
- Fix: refit with `inverseA(ph, nodes = "ALL", scale = FALSE)` on the GTDB lengths. This is done in
  `audit/02_alt_refit.R` (method C); see `audit/README.md` for the comparison. Then report which
  tree each model used.

**M2. The QC rule lets through genomes with too few ribosomal genes.**
- Where: `10_build_table.py:73`, `qc = (n_proteins >= 500) & ~(n_HE < 20)`.
- What it does: gRodon2 returns `n_HE = NA` when a genome has fewer than 10 HE genes
  (`06b_grodon.R:12`). For those genomes `NA < 20` is False, so they **pass** QC, although they are
  exactly the genomes the rule is meant to exclude.
- Effect: six genomes in the final table have `n_HE = NA`. This is numerically negligible, but it
  is a rule violation.
- Fix: `qc = (n_proteins >= 500) & (n_HE >= 20)`, treating NA as a failure.

**M3. "MAG vs isolate" is confounded with culturability of the lineage by construction.**
- Where: `04_select_dataset.py:40`. The per-species pick gives isolates a +5 quality bonus.
- Effect: only 3 of the 4,919 MAG-represented species in the analysis set have *any* in-domain
  isolate genome in the pool. The main-analysis moderator therefore means "species never cultured
  vs species cultured", not "MAG assembly vs isolate assembly".
- This is not a coding bug, but the paper's two explanations ("population consensus in MAGs",
  "culturing selection") are presented as if the design could separate them, and it cannot. The
  within-species analysis (`replication_within_species/`) can.

**M4. SAGs are coded as MAGs.**
- Where: `01_habitat_classify.py:69`. `is_mag` matches "single cell", so 92 marine
  single-amplified genomes are treated as MAGs.
- Why it matters: a SAG comes from one cell, so the "population consensus" explanation does *not*
  apply to it. Pooling SAGs with MAGs blurs the very contrast the paper interprets.
- Fix: add a three-level `genome_source` (isolate / MAG / SAG) and report SAGs separately. A first
  test is in `audit/README.md`.

### Minor

**m1. Coordinate-precision parser** (`10_build_table.py:63`).
- What it does: it keeps the *minimum* number of decimals among coordinate parts that have a
  decimal point. Integer parts (e.g. "8 N 136.123 E") are skipped, so 2 such genomes are flagged
  "precise".
- Trailing zeros are lost in NCBI strings ("113.20"), so "2 decimals" partly means "3 decimals
  ending in 0".
- Impact: small, but `coord_precise` should be read as a rough proxy.

**m2. Site definition** (`16_cluster_robust.R:21`, `18`, `21`).
- What it does: `round(lat*2)/2` gives cells *centred* on half-degree points, not the 0.5° climate
  cells. R also rounds half to even.
- Impact: sites that differ by grid alignment are merged or split. Harmless as a clustering unit,
  but it is not the same cell as the predictor's grid cell.

**m3. DORM_z reference set** (`10_build_table.py:41–43`).
- What it does: the z-scores are computed over all 7,867 annotated genomes (both domains, pilot
  included) *before* QC.
- Not a bug, but it should be stated, because any replication must reuse these parameters
  (the within-species replication does).

**m4. Length thresholds include the stop codon** (`08_aggregate_modules.py:31`).
- `len_aa = (end − start + 1) // 3` counts the stop codon, so the SAS ≤ 350 aa and SAH ≤ 300 aa
  thresholds (`:85–86`) are effectively ≤ 349 and ≤ 299.
- Impact: negligible.

**m5. TA pairing** (`08_aggregate_modules.py:102–106`).
- What it does:
  - pairs are counted per toxin;
  - one antitoxin can pair with several toxins;
  - there is no strand or operon-orientation check;
  - "within 3 genes" uses pyrodigal's per-contig gene index, so it is correct across contigs.
- This matches the dictionary's stated rule. TA counts are an upper bound on true operons, which
  is the same for both domains, so H3 is not affected.

**m6. Sensitivity 10 omits sporulation** (`14_extra_tests.R:34`). The pre-registered "MCMCglmm
Poisson for each module count" runs only n_rpf, n_hib and n_sr.

**m7. Misleading comment** (`18_mag_isolate.R:6`). The comment says "with family fixed effects",
but the shared-family model (`:43–48`) only restricts to the 39 families and adds no family
terms. The paper's wording ("restriction to the 39 families") is accurate; the comment is not.

**m8. The interaction models force a common covariate slope** (`exploratory_config.tsv:2`,
`21:33`).
- What they do: genome size and completeness slopes are shared between MAGs and isolates.
- Why it matters: MAG completeness and genome size relate to gene counts differently from those
  of isolates.
- The separate MAG-only and isolate-only fits relax this and agree, so this is a robustness note,
  not an error.

**m9. Smoke-test output in `results/`.** `siteMCMC_cluster_robust.tsv` and
`siteMCMC_site_capped_betas.tsv` come from a G = 2, R = 2 test run (`run_final.sh`). Their
jackknife SE (0.0004) and p = 0 are meaningless.

**m10. The pilot was drawn from the P0 tails** (`04_select_dataset.py:66–72`).
- What it does: 400 species are taken from the bottom and top 20 % of P0.
- Effects:
  - the confirmatory P0 distribution is thinned at the tails;
  - the pilot "replication" models are fitted on an extreme-groups design;
  - the pilot estimates are not directly comparable to the confirmatory slopes.
