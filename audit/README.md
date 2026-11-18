# Independent audit of "The Dormancy Rule" 

Auditor: separate AI-assisted audit session, 17 Nov 2026. The brief was to find errors, not to agree.

**Supporting files:**
- `number_check.tsv` (223 claims);
- `prereg_audit.md`;
- `code_review.md`;
- `refit/`, `refit_AC/` (alternative models);
- `strain_heterogeneity_tests.tsv`, `sag_test.tsv` (new tests);
- `literature_second_opinion.md`, `literature_queries.tsv` (67 queries),
  `study_scale_comparison.tsv`.

## Verdict

**The confirmatory conclusion (H1 practically null, phylogeny dominates) is sound, and it gets
*stronger* under scrutiny.**
- The pre-registration is intact, and no pre-registered model was changed after the freeze.
- No critical code bug was found.

**The interpretation of the exploratory MAG/isolate reversal needs to be rewritten.**
- The reversal itself is robust to every alternative model.
- But the evidence now points away from the paper's favoured "MAG population-consensus" reading.
- The main analysis cannot separate genome source from lineage culturability.

## Alternative re-fits (task 4; confirmatory n = 7,327)

| Model | Pooled H1 (DORM_z ~ P0) | MAG × P0 interaction | MAGs | Isolates |
|---|---|---|---|---|
| Original phylolm λ | +0.0065 [0.0016, 0.0113] | +0.024 | +0.018 | −0.008 |
| lmer, nested taxonomy (no site) | +0.0066 [0.0015, 0.0118] | +0.025 | +0.016 | −0.011 |
| **lmer, nested taxonomy + site RE** | **+0.0049 [−0.0015, 0.0114], p = 0.14** | **+0.028, p = 1e-5** | +0.022 | −0.011 (p = 0.012) |
| **MCMCglmm, GTDB branch lengths + site RE** | **+0.0039 [−0.0019, 0.0098], p = 0.21** | **+0.026 [0.016, 0.039]** | — | — |
| MCMCglmm, *Grafen* lengths + site (now in the paper) | +0.0079 [0.0003, 0.0143] | +0.029 [0.016, 0.044] | +0.024 | −0.009 |
| phylolm OU (ΔAIC ≈ +2,200 vs λ: a poor fit) | −0.0023, p = 0.06 | +0.021 | +0.014 | −0.005 |

**Conclusion of the re-fits.**
- **Pooled slope:** once site and phylogeny are modelled together on the real tree, the pooled
  slope is no longer different from zero.
- **Interaction:** it survives every method (+0.021 to +0.029) *between* species, but not within
  species (−0.011 [−0.028, +0.005]; see issue 3).
- **The paper's MCMCglmm:** it uses Grafen branch lengths. Its animal variance (3.95) is 12× the
  variance under the real tree (0.33), and it is the only phylogeny + site model that still calls
  the pooled slope "significant".

## Issues ranked by severity

### Critical: none

### Major

1. **MCMCglmm uses Grafen branch lengths** (`code_review.md` M1; scripts 14, 16, 21).
   - Effect: pre-registered sensitivity 10 and the paper's "phylogeny + site" models assume a
     different tree from every phylolm model, and this is not disclosed.
   - **Fix:** refit on the GTDB lengths (`audit/02_alt_refit.R` method C) and report those values.
     On the real tree the pooled estimate is +0.0039 [−0.0019, 0.0098].
2. **The pilot is presented as a "held-out replication" of the marine TA deficit.**
   - The pilot is the *discovery* sample for this effect: its H3 estimate, −0.27 with p = 0.008,
     was seen before the freeze. H3 was then registered in the opposite direction.
   - **Fix:** say that the confirmatory set replicated the pilot. Delete "replicated in the
     held-out pilot" from the paper and from REPORT §3.6.
3. **"MAG vs isolate" means "never-cultured vs cultured species".**
   - Only 3 of the 4,919 MAG-represented species have any in-domain isolate genome, because of the
     +5 isolate bonus in `04_select_dataset.py`. The two explanations in the paper therefore cannot
     be separated in this design.
   - **Fix:** state this explicitly, and lead with the within-species replication
     (`replication_within_species/`).
   - **Now done:** 7,098 HQ genomes in 1,473 species, species fixed effects, site-clustered.
     | Term | Within species, 95 % CI | Between species |
     |---|---|---|
     | MAG slope | −0.008 [−0.023, +0.006] | +0.018 |
     | Isolate slope | +0.003 [−0.006, +0.011] | −0.008 |
     | Interaction | −0.011 [−0.028, +0.005] | +0.024 |

     The interaction and isolate CIs exclude the between-species values at every completeness
     threshold, so **the composite reversal is a between-lineage pattern, not a genome-source
     effect**.
   - Two caveats:
     - The MAG DORM_z slope depends on completeness: +0.009 [−0.001, +0.019] at ≥ 95 %, so it does
       not cleanly exclude +0.018 there.
     - Extra hibernation copies keep a positive within-species slope in MAGs (+0.006 to +0.017
       across thresholds; q = 0.26 in the main test). This is the one MAG-side signal still open.
       See `replication_within_species/README.md`.
4. **The favoured explanation is not supported by three direct tests run in this audit:**
   - (a) CheckM strain heterogeneity *falls* with P0 in MAGs (p = 2e-7).
   - (b) The MAG slope is as large in MAGs with zero strain heterogeneity (+0.024, p = 8e-7) as in
     heterogeneous ones (+0.013). The interaction is null for DORM_z, hib_extra and n_hib.
   - (c) In the ocean, hibernation paralogue copies are *flat* in MAGs (+0.002) but *fall* in
     isolates (−0.023, p = 2e-4). Here the "reversal" is isolates losing copies, not bins gaining
     them.
   - The 92 SAGs, which are single cells and therefore not consensus genomes, are too few to be
     informative (+0.002 ± 0.024).
   - **Fix:** rewrite Discussion §"Genome source" as an open question that tilts towards
     lineage/culturing effects. Report SAGs as a separate source (they are currently coded as MAGs;
     code review M4).
5. **The Balon numbers are cherry-picked.**
   - "About 17 % below 5 °C vs about 5 % above 20 °C": the octile curve is non-monotonic (17, 12,
     10, 9, 4, 12, 5, 12 %). The coldest octile spans −1.6 to 9.7 °C, and the warmest is 12 %. At
     the stated thresholds the values are 21 % (n = 80) and 8 % (n = 379).
   - **Fix:** report the logit slope and the full curve, and soften the student report's "bonus
     finding".
   - Also: Helena-Bueno et al. (2024) ran no temperature test (`literature_second_opinion.md` §C),
     so the result is new, not "fitting" a known pattern.
6. **The Lennon & Jones (2011) citation for a metagenomic marine TA deficit is wrong.** That
   review contains no marine-vs-land TA statement. Makarova et al. (2009) and Pandey & Gerdes
   (2005) point the *other* way.
   - **Fix:** cite the binning-independent defensome-depletion work (Meaden et al. 2025;
     Beavogui et al. 2024), and frame the deficit as accessory/defence DNA, as the paper already
     partly does.
7. **The QC rule passes NA genomes** (`10_build_table.py:73`). Six genomes with fewer than 10 HE
   genes were retained against the pre-registered rule. **Rerun done** (`qc_fix_rerun.tsv`,
   n = 7,321): H1 +0.0064 [0.0016, 0.0113], interaction +0.024. The impact is negligible, and
   no decision changes.

### Minor

- **"All sensitivity CIs inside the SESOI" is false.** S_marine reaches 0.0155, and
  `decisions.md` already says "inconclusive".
- **"Nine of ten" sensitivity analyses** miscounts the pre-registered items. S_sum is on a scale
  where the SESOI is meaningless.
- **Timeline wording.** "Frozen before any confirmatory genome was annotated" contradicts the
  pre-registration, which says annotation was running at the freeze. Use the pre-registration's
  own wording.
- **Student report:** "λ = 0.94 means 94 % of the pattern follows the family tree" is incorrect.
  λ scales the phylogenetic covariance; it is not a variance share.
- **Negative control.** The ribosomal-protein negative control is not neutral in MAGs: Mise &
  Iwasaki (2022) report codon-usage-driven loss of ribosomal-protein genes from MAGs. Replace or
  supplement it with typical-codon-usage single-copy genes.
- **Smoke-test files.** `results/siteMCMC_cluster_robust.tsv` is a G = 2, R = 2 smoke test.
  Label it or delete it.
- **Sensitivity 10** omits the sporulation module.
- **Coordinate-precision parser** edge cases (2 genomes); site cells are centred, not
  grid-aligned.
- **Rounding.** 7 values differ from their source only through double rounding (e.g. MAG CI lower
  bound 0.0114 printed as 0.012; top-10-site share 42 % printed as 41 %).
- **Scope of "never reported".** Systematic MAG-vs-isolate functional divergence is documented
  (e.g. Rühlemann et al. 2025 preprint; Gupta & Almeida 2025). Claim novelty only for the
  *environment interaction*.

## Number audit summary

- 223 claims checked: 210 match, 6 differ only by rounding, 6 mismatch (items 2 and 5 above, plus
  three minor items), 1 unsourced ("about 900,000" GTDB genomes in the student report).
- The upstream fix of "three times" → "about twice" (commit a78a042) was verified.

## Scale (question raised by the user)

The confirmatory set, 7,327 species with a continuous predictor and phylogenetic control, is
**large for this question**:
- Prior work, with the published n verified:
  - TA and dormancy comparisons: Makarova 2009 (750 genomes), Fucich & Chen 2020 (71 strains),
    Galperin 2022 (180), Browne 2021 (1,358 Firmicutes). All of these are categorical, and most are
    non-phylogenetic.
  - The closest template is Hao et al. 2026 PNAS: 25,785 species, MCMCglmm, but with habitat
    clusters rather than sites.
  - Studies of 60,000–200,000 genomes exist for other traits, but none uses a site-level climate
    predictor.
- **The binding constraint is not genome count. It is the effective number of independent sites:**
  - only 368 sites for 4,492 MAGs;
  - the number of species that have both a MAG and an isolate genome (74 at HQ, 111 at medium
    quality).
