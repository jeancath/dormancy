# Within-species replication of the MAG/isolate reversal 

**Status.** This is an exploratory replication of an exploratory finding. It is not
pre-registered and does not affect the confirmatory H1 decision.
- Files: `methods.md` (full methods), `results.tsv` (every model), `fig_within_species.png`.
- Fraction analysed: **100 % of the eligible high-quality pool.** 7,206 genomes in 1,503 species
  with ≥ 2 eligible genomes were selected and annotated. **7,098 genomes in 1,473 species** passed
  QC, including **71 species with both MAG and isolate genomes**.

## Question

Between species, DORM_z rises with climate variability (P0) in MAG-represented species (+0.018) and
falls in isolate-represented species (−0.008), an interaction of +0.024. Is this a property of
*genome source* (it should then reappear between MAG and isolate genomes of the same species), or
of *lineage*, i.e. never-cultured vs cultured species (it should then vanish within species)?

## Answer

**The reversal, read as an interaction, does not reappear within species.** Within species, the
MAG × P0 interaction excludes the between-species +0.024 at every completeness threshold. The
isolates' negative slope also vanishes within species.

**One component holds up on the MAG side:** extra hibernation-factor copies. Their within-species
slope in MAGs is positive and about as large as between species. This is the one place where a
MAG-specific effect cannot be ruled out.

All estimates use species fixed effects with site-clustered 95 % CIs. In the first column, "≥ 95 %"
is the subset with CheckM2 completeness ≥ 95 % (5,734 genomes, 1,113 species).

| DORM_z ~ P0 | Within species, all HQ | Within species, ≥ 95 % | Between species (main) |
|---|---|---|---|
| All genomes | +0.001 [−0.006, +0.009] | — | +0.0065 (H1) |
| MAG genomes | −0.008 [−0.023, +0.006] | +0.009 [−0.001, +0.019] | +0.018 |
| Isolate genomes | **+0.003 [−0.006, +0.011]** | **+0.004 [−0.005, +0.012]** | −0.008 |
| MAG × P0 interaction | **−0.011 [−0.028, +0.005]** | **+0.005 [−0.007, +0.018]** | +0.024 |
| Mixed species only: interaction (71 species, 413 genomes) | −0.011 [−0.036, +0.013] | — | — |
| Mixed species only: MAG − isolate gap | +0.014 [−0.015, +0.043] | — | — |
| Paired species means: ΔDORM_z on ΔP0 | +0.013 [−0.021, +0.048] | — | — |

Bold marks intervals that exclude the between-species value.

**The MAG slope for DORM_z depends on completeness.** Its sign changes with the threshold
(`10_completeness_strata.R`): −0.008 for all genomes, +0.009 at ≥ 95 %, +0.005 at ≥ 98 %.
- The negative value in the full set comes from missing single-copy stringent-response genes in
  less complete MAGs. 12 % of MAGs have zero *relA/spoT* hits, against 0.2 % of isolates (see the
  medium-quality extension, where this is significant).
- In the cleanest genomes, MAGs and isolates show a small, similar positive within-species slope
  (+0.004 to +0.009): no genome-type difference, but also no clean exclusion of +0.018 for MAGs at
  ≥ 95 %.

**Pooled slope.** Every model family agrees on the pooled within-species slope for DORM_z:
- within-P0 random slopes, with or without a site random effect;
- the replication plan's `(1 + P0 | species)` model;
- species fixed effects.

All lie between −0.003 and +0.001, inside the SESOI (±0.015), consistent with the practically-null
H1.

**Extra hibernation copies** (`log1p(hib_extra)`):

| Term | HQ, all | HQ ≥ 95 % | HQ ≥ 98 % | MQ extension, all | Between species |
|---|---|---|---|---|---|
| MAG slope | +0.009 [+0.001, +0.018] | +0.017 [+0.004, +0.030] | +0.006 [−0.000, +0.012] | +0.006 [+0.001, +0.012] | +0.013 |
| Isolate slope | −0.002 [−0.012, +0.008] | −0.002 [−0.011, +0.008] | −0.002 [−0.012, +0.008] | −0.002 [−0.012, +0.008] | −0.011 |

- The MAG half is stable in sign and similar in size to the between-species value.
- In the main HQ test it does not survive BH correction across ten responses (q = 0.26). The
  completeness strata are nested subsets, not independent replications.
- No other module has a MAG slope with p < 0.1 in the full HQ set.

## Control: does this species subset carry the reversal between species?

Species with ≥ 2 genomes are not a random subset: they are better-studied. Seven of the twelve
largest are pathogens or pathogen complexes (*Burkholderia ubonensis*, 250 genomes; *Pseudomonas
aeruginosa*; *Salmonella enterica*; *Escherichia coli*; *Listeria monocytogenes*; *B. mallei*;
*Klebsiella pneumoniae*). I therefore refit the between-species model on species means of the same single-type
species (1,402 species; `09_between_on_subset.R`). The DORM_z interaction is **+0.020 ± 0.029**,
which is the same size as the main analysis (+0.024) and opposite in sign to the within-species
value (−0.011). With 1,402 species means and no tree it is imprecise on its own (p = 0.49), and
with phylum fixed effects it shrinks to +0.005. This control does not prove that the subset
carries the reversal. It does show that the within-species null is not caused by a subset whose
between-species pattern already looked different.

## Sensitivity

`08_sensitivity_tables.py`. Species fixed effects, site-clustered 95 % CI:

| DORM_z ~ P0 | Main (7,098 genomes / 1,473 species) | ≤ 20 genomes per species (5,850 / 1,473) | SAGs removed (7,058 / 1,462) |
|---|---|---|---|
| MAG slope | −0.008 [−0.023, +0.006] | −0.009 [−0.025, +0.006] | −0.007 [−0.021, +0.008] |
| Isolate slope | +0.003 [−0.005, +0.011] | −0.000 [−0.009, +0.008] | +0.003 [−0.005, +0.011] |
| MAG × P0 interaction | −0.011 [−0.028, +0.005] | −0.009 [−0.026, +0.008] | −0.010 [−0.026, +0.007] |
| Interaction, mixed species only | −0.011 [−0.036, +0.013] | −0.015 [−0.044, +0.014] | −0.008 [−0.031, +0.016] |
| Extra hibernation copies, MAG slope | +0.009 [+0.001, +0.018] | +0.010 [+0.001, +0.019] | +0.010 [+0.001, +0.018] |

- Neither the large species nor the SAGs drive any estimate.
- Capping lost 3 mixed species (68 remain), because the random draw removed their only genome of
  one type.

## Interpretation (my opinion)

1. **The isolate half of the reversal is a lineage pattern.** Within species, isolate genomes from
   more variable climates have no fewer dormancy genes (DORM_z +0.003; extra hibernation copies
   −0.002). The between-species negative slope must therefore come from *which* cultured lineages
   live where. This matches the audit's marine test, where the between-species reversal was driven
   by isolates losing copies.
2. **The composite reversal is not a genome-source effect.** The within-species interaction
   excludes +0.024 at every completeness threshold.
3. **The MAG half is partly unresolved, and only for hibernation paralogues.** MAG genomes of one
   species carry more extra hibernation copies at more variable sites, by about as much as between
   species.
   - Two readings remain:
     - a real within-species association (strains from variable sites carry extra *hpf*/*raiA*
       copies);
     - a binning effect that collapses co-occurring strains into one bin and so adds paralogues.
   - The audit's strain-heterogeneity test argues weakly against the binning reading (the slope is
     equal in zero-heterogeneity MAGs), but it is not decisive.
   - A long-read or SAG comparison would settle it.
4. **Recommendation for the paper:**
   - Describe the DORM_z reversal as a between-lineage pattern (never-cultured vs cultured species).
   - Drop "population consensus" as the favoured reading of the composite.
   - Report the hibernation-paralogue result in MAGs as the one within-species signal, flagged as
     not multiplicity-robust.
5. **Caveat.** A within-species null also fits a lineage-level adaptation fixed per species. The
   design rules out assembly as the explanation of the *composite* reversal. It does not decide
   between "adaptation fixed per lineage" and "culturability correlates with habitat".

## Medium-quality extension

Reported separately in `extension_mq/`, once its annotation finishes: 11,788 genomes, 2,716
species and 111 mixed species before QC.
