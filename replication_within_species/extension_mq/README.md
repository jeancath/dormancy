# Medium-quality extension of the within-species replication

**Status.** This is an *extension*, reported separately from the HQ result in `../README.md`.
- **Selection rules.** The same as the HQ analysis, but with the MIMAG medium-quality threshold
  (CheckM2 completeness ≥ 50 %, contamination < 10 %) instead of ≥ 90 % / < 5 %.
- **Scripts.** `05_select_mq.py` (selection), `run_annotation_list.sh` (annotation), then the same
  `07_run_downstream.sh`, `09_between_on_subset.R` and `10_completeness_strata.R` as the HQ run.
- **Sample.** 11,798 genomes annotated (7,206 HQ + 4,592 new), none lost. 601 genomes failed QC.
  **11,050 genomes, 2,515 species and 97 mixed species were analysed**, of which 995 species span
  more than one site.
- **Files:**
  - `results_mq.tsv`: every model, the between-species subset control and the completeness strata;
  - `fig_within_species_mq.png`;
  - `mixed_species_pairs_mq.tsv`;
  - `within_species_table_mq.tsv.gz`.

## Result

Species fixed effects, site-clustered 95 % CI:

| DORM_z ~ P0 | MQ, all (11,050 genomes) | MQ, completeness ≥ 90 % | MQ, ≥ 95 % | Between species (main) |
|---|---|---|---|---|
| Pooled within species | −0.002 [−0.009, +0.005] | — | — | +0.0065 |
| MAG slope | **−0.010 [−0.020, −0.000]** | −0.005 [−0.019, +0.010] | +0.009 [−0.002, +0.019] | +0.018 |
| Isolate slope | +0.001 [−0.007, +0.010] | +0.003 [−0.005, +0.012] | +0.004 [−0.004, +0.012] | −0.008 |
| MAG × P0 interaction | −0.011 [−0.024, +0.001] | −0.008 [−0.024, +0.009] | +0.005 [−0.008, +0.018] | +0.024 |
| Mixed species (97): interaction | +0.003 [−0.030, +0.036] | — | — | — |
| Mixed species (97): MAG − isolate gap | +0.032 [+0.003, +0.062] | — | — | — |

## What the extra genomes add

1. **The interaction stays absent.** With 71 % more species, the within-species MAG × P0
   interaction still excludes the between-species +0.024 at every completeness threshold. This
   agrees with the HQ result.

2. **The significantly negative MAG slope is a completeness artefact, not biology.**
   - Mechanism: it is driven by the stringent-response module (MAG slope −0.017, BH q = 0.02), and
     *relA/spoT* homologues are near single-copy. 12 % of MAGs have **zero** hits, against 0.2 %
     of isolates.
   - Evidence: restricted to ≥ 95 % completeness, the stringent-response MAG slope is +0.002 and
     the DORM_z MAG slope is +0.009.
   - **Lesson:** at medium quality, count-based modules measure completeness as much as gene
     content, even with completeness as a covariate. **The HQ analysis is the primary result.**

3. **Extra hibernation copies:**
   - MAG slope: positive again, +0.006 [+0.001, +0.012] for all genomes and +0.015
     [+0.004, +0.027] at ≥ 95 %;
   - isolate slope: flat (−0.002);
   - mixed-species interaction: +0.012 [−0.001, +0.024].

   This is the same partial MAG-side signal as in the HQ set. It is not an independent
   replication, because the HQ genomes are a subset of these.

4. **Mixed species.** The MAG − isolate gap at P0 = 0 (+0.032, p = 0.03) means MAGs of a species
   score slightly higher on DORM_z than its isolates, after adjusting for completeness and genome
   size. This is a level difference, not a slope difference. It is consistent with
   multi-copy-collapsing or contaminating bins adding paralogues, and should be read as a caution
   about MAG gene counts, not as dormancy biology.

## Opinion

The medium-quality data confirm the HQ conclusion that the composite reversal is not a
genome-source effect, and they show why the medium-quality tier cannot be the primary analysis:
completeness bias in single-copy modules is larger than the effects being studied. If the paper
cites this extension, it should cite it only as "same answer for the interaction with 2,515
species, and completeness artefacts below 90 %".
