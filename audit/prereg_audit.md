# Pre-registration audit (audit task 2)

Auditor: separate AI-assisted audit session, 17 Nov 2026. The aim was to find problems, not to confirm.

## 1. Was the frozen file changed?

- `git log --follow -- prereg/PREREGISTRATION.md` lists two commits:
  - `73f9674` (09:37 UTC), the draft;
  - `9fae630` (10:16 UTC), "PREREG FROZEN".
- `git diff 9fae630 HEAD -- prereg/` is empty, so **the file has not changed since the freeze. Pass.**
- **Binding files.** The pre-registration names three binding files. Their history:
  - `scripts/models_config.tsv`: last changed at the freeze commit. Pass.
  - `scripts/08_aggregate_modules.py`: last changed at `c7f3cf0`, before the freeze. Pass.
  - `scripts/09_models.R`: changed after the freeze (`f8fedfe`, `752bcca`). The diff only adds an
    optional config argument and two exploratory predictors (`abslat_z`, `T_fracseas`). No
    pre-registered model is altered. Pass.
  - `scripts/10_build_table.py` (not named as binding): the post-freeze diffs only add columns. The
    pre-registered DORM_z, QC rule and predictors are byte-identical. Pass.
- **Dictionary.** `dictionary/` (thresholds, dictionary, sporulation core) was last changed at
  `c7f3cf0`, before the freeze. Pass.
- **Limitation.** Commit timestamps are set by the author. The pre-registration defines the
  registration time as the GitHub *push* time of `9fae630` on branch `ccr-70a28546-5vpeux`. This
  audit cannot see that push time from a clone. Someone with repository admin access should check
  it in the GitHub push/events log.

## 2. Decisions against the pre-registered rules

| Test | Rule (§5) | Result | Correct? |
|---|---|---|---|
| H1 direct | practically null if the 95 % CI ⊂ (−0.015, +0.015) | [0.0016, 0.0113] | yes, practically null |
| H1 total | same | [0.0007, 0.0106] | yes |
| H1a–d | BH q < 0.05 within the family | smallest q = 0.12 | yes, none supported |
| H6 | BH across the decomposition family; Wald test of residual − seasonal | q ≥ 0.052; Wald p = 0.76 | yes |
| H3 | supported if the marine coefficient is not significantly negative | −0.231, p = 2e-19 | yes, not supported |
| H5 | succeeds if β > 0 and p < 0.05 | p = 0.067 | yes, failed |
| §8 falsification | H1 and H1-total both practically null | both | yes, the falsification outcome applies |

**Two gaps in the rules:**
- **No "inconclusive" category.** The rules do not define it. `scripts/17_decisions.py` adds one
  ("anything else"), and S_marine falls into it. That is a reasonable completion of the rules, but
  `results/decisions.md` should say that the category was added after the freeze.
- **Sensitivity analyses have no decision rule.** The paper's statement that "all [sensitivity CIs]
  fall inside the SESOI" is false for S_marine, whose upper bound is 0.0155 (see
  `number_check.tsv`). This does not change the H1 decision, because the sensitivity analyses are
  not decisive under the pre-registration.

## 3. Deviations: disclosed and undisclosed

**Disclosed** in `results/deviations.md` items 1–11, and checked against the code:
- DefenseFinder not used;
- gRodon re-run;
- HMM batch re-runs;
- a debugging run of the models;
- sporulator ≥ 27;
- the genus sanity check;
- the exploratory config;
- the site-pseudoreplication analyses;
- the re-parameterised variability models;
- the follow-ups;
- the MCMCglmm sparse inverse.

**Undisclosed or misdescribed** (ranked):

1. **MCMCglmm branch lengths (major).**
   - What the code does: `scripts/14_extra_tests.R` (pre-registered sensitivity 10),
     `scripts/16_cluster_robust.R` and `scripts/21_mcmc_site_interaction.R` all call
     `compute.brlen(ph, method = "Grafen")`. This **discards the GTDB branch lengths** and replaces
     them with arbitrary Grafen lengths.
   - Why it matters: the phylolm models, including the confirmatory H1, use the real GTDB branch
     lengths. Deviation 11 says the sparse formulation "implies the same phylogenetic covariance".
     That is true for `nodes = "TIPS"` versus `nodes = "ALL"`. It is false as a comparison with
     the phylolm models, because both MCMCglmm variants use a different tree.
   - The pre-registration says nothing about branch lengths for sensitivity 10, so this is an
     undisclosed analytic choice. It is not a rule violation.
   - Independent re-fit: `audit/refit/` re-fits with the real GTDB branch lengths.
2. **Timeline misstatement (minor; it affects credibility, not validity).**
   - What the documents claim: `results/REPORT.md` ("frozen … before any confirmatory genome was
     annotated") and the paper's Methods ("before any confirmatory genome was annotated").
   - What the pre-registration says: "Confirmatory genome download, gene calling and HMM annotation
     were running in the background when this was frozen."
   - The pre-registration's own wording (no predictor–response association had been computed) is
     the correct claim, and the documents should use it.
3. **The pilot is described as a replication (major for interpretation).**
   - The pilot's H3 estimate (−0.272, p = 0.008, `results/pilot_report.md`) was seen *before* the
     freeze. H3 was then pre-registered as "marine genomes are **not** TA-depleted".
   - The confirmatory result (−0.231) therefore replicates the pilot, not the other way round.
   - Re-fitting the pilot after the confirmatory analysis (`pilot_replication_models.tsv`,
     −0.251) re-describes the discovery data. It is not independent replication. The paper's
     Results ("It replicated in the held-out pilot") and REPORT §3.6 should be corrected.
   - The same applies to describing the pilot as used "only to estimate the residual structure":
     the pilot association estimates were computed and reported (`results/pilot_report.md`).
4. **QC rule implementation (minor).**
   - The rule is: exclude if there are "< 20 ribosomal-protein genes detected".
   - What the code does: `10_build_table.py` tests `~(n_HE < 20)`. gRodon2 returns NA when there
     are fewer than 10 HE genes, so NA genomes pass the filter.
   - Effect: six genomes in the final table have `n_HE = NA` and were retained. The direction of
     the error is the opposite of the rule's intent. A rerun without them is recommended
     (negligible impact is expected at n = 7,327).
5. **DORM_z standardisation scope (minor, not specified in the pre-registration).**
   - DORM_z is z-scored over all 7,867 annotated genomes. That covers both domains and the pilot,
     and is done *before* QC.
   - The pre-registration does not specify the reference set. The choice should be stated in the
     Methods, because any future replication must reuse these exact parameters.
6. **A test run left in `results/` (minor).** `results/siteMCMC_cluster_robust.tsv` and
   `results/siteMCMC_site_capped_betas.tsv` come from a run with G = 2 jackknife groups and R = 2
   resamples (`run_final.sh`). Their jackknife SE (0.0004) and "p = 0" are meaningless. They
   should be labelled as smoke tests or deleted, so that nobody cites them.
7. **"Nine of ten" sensitivity analyses (minor).**
   - The pre-registration lists 10 items. Item 4 ("marine only and land only") is one item, and
     item 10 (MCMCglmm Poisson) has not finished.
   - The paper's "nine of ten gave slopes of +0.0046 to +0.0072" counts marine and land separately
     and includes DORM_sum, which is on a different scale (log1p of a sum). That scale makes "all
     inside the SESOI" meaningless for S_sum.
8. **Sensitivity 10 covers only three of the four modules (minor).** `14_extra_tests.R` fits the
   MCMCglmm Poisson to n_rpf, n_hib and n_sr only. The pre-registration says "for each module
   count", and spo_strict is omitted. This should be disclosed when the result is reported.

## 4. Analytic flexibility

- **Specification count.** The confirmatory set is fixed. Around it there are 54 exploratory,
  34 follow-up, 5 pilot-replication and about 10 site-robustness specifications. They are all
  labelled exploratory, which is good practice.
- **Problem 1: adjustment within families only.** The FDR "families" for the exploratory
  (`expl`) and follow-up (`follow`) sets are adjusted within each set only.
- **Problem 2: the follow-ups were chosen after seeing the result.** Every follow-up
  explanation was selected after the MAG/isolate split was seen. Some of them were then reported
  as "ruling out" alternatives (coordinate precision, contamination).
  - A non-rejection at these sample sizes is weak evidence.
  - Example: in near-complete clean MAGs the slope falls by 60 % (+0.007, p = 0.35).
  - Example: MAGs with coarse coordinates have a **3× larger** slope (+0.042) than MAGs with precise
    coordinates (+0.014). This was not discussed as possible evidence *for* a coordinate or site
    artefact in MAGs.
- **Problem 3: an unlisted moderator.** The "MAG vs isolate" moderator was not a pre-registered
  moderator. `04_select_dataset.py` gives a +5 quality bonus to isolates. A species is therefore
  represented by a MAG essentially only when no comparable isolate exists, which makes "MAG" largely
  a proxy for "species never cultured". The main analysis cannot separate a genome-source
  (assembly) effect from a lineage culturability effect. The within-species design can, and
  this framing should be added to the paper.
