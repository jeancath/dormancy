# Validation: phenotype calibration of the sporulation and Rpf calls

**This checks the dictionary against outside labels; it does not tune it.** Nothing in `dictionary/`, `scripts/` or the pre-registered
thresholds was changed. Every score comes from the unchanged pipeline (`scripts/03` → `05` → `08`),
run on genomes with independent phenotype labels.

Interpretation and recommendations for the main analysis: [`INTERPRETATION.md`](INTERPRETATION.md).

Details:
- [`sporulation_calibration.md`](sporulation_calibration.md): sporulation calibration.
- [`rpf_check.md`](rpf_check.md): Rpf and Sps check.

## Headline numbers

**Sporulation.** 543 Bacillota type-strain genomes (195 sporeformers, 348 not), with labels from
the Madin et al. 2020 phenotype compilation:

| | Value [95 % bootstrap CI] |
|---|---|
| AUC of `spo_strict` | **0.917** [0.891–0.940] (0.964 when ≥ 2 curated sources agree on the label) |
| AUC of `spo_relaxed` | 0.910 [0.882–0.933] |
| Pre-registered rule (Spo0A & ≥ 25): sensitivity | **0.903** [0.857–0.942] |
| specificity | **0.779** [0.735–0.822] |
| precision | 0.696 |
| MCC | **0.655** [0.594–0.713] |
| Youden-optimal threshold | ≥ 27 (sensitivity 0.867, specificity 0.871, MCC 0.724) |
| 71 non-Bacillota non-sporulators | all score 0: no false positive outside Bacillota |
| 12 non-endospore "sporulators" (actinomycete exospores) | all score 0, as intended |

**Results by lineage (pre-registered rule)**
- Negativicutes: MCC 0.95.
- Bacillales: sensitivity 0.95–1.00.
- **Clostridia: specificity only 0.31.** 64 of the 77 false positives are here. Most carry 25–33
  core genes and have "non-sporulating" labels that are doubtful (e.g. *Ruminococcus bromii*, shown
  to sporulate in 2018), or they are partial losers of sporulation.
- **Spore-forming Erysipelotrichia: 0 of 3 detected.** They score 18–19, a limitation the
  dictionary README already predicted.

**Rpf**
- Tier 1, literature counts: `n_rpf` matches the documented count **exactly in 23 of 24
  Actinomycetota**. This includes 5 in the *M. tuberculosis* complex, 4 in other pathogenic
  mycobacteria, 2–3 in corynebacteria and 5–6 in streptomycetes, plus the documented absences.
- Tier 2, Puspita 2013 BLASTp survey (95 Actinomycetota): **exact in 83 %, within ±1 in 98 %,
  Spearman ρ = 0.96, no missed Rpf-positive genome.**
- Sps matches Ravagnani 2005 exactly in 16 of 20 Bacillota; the other 4 are +1.
- Outside Actinomycetota, PF06737 also detects the "distantly related" Rpf-like protein of
  *S. aureus*, *S. epidermidis* and *O. oeni* (n_rpf = 1). Interpret `n_rpf` outside
  Actinomycetota with care.

## Is the pre-registered sporulation threshold supported?

**Yes, as a high-sensitivity rule, but it is not the best single cut.**
- Discrimination is strong. Spo0A plus 25/33 recovers 90 % of labelled sporeformers and never
  fires outside Bacillota.
- On these labels, ≥ 27 would trade 7 true positives for 32 fewer false positives.
- The threshold stays at 25, as the replication plan requires. ≥ 27 is a post-hoc optimum on noisy
  labels, and the pre-registration already includes a ≥ 27 (80 %) sensitivity analysis, which this
  result supports running.
- Clostridia and Erysipelotrichia are where the binary call is least reliable. Reporting the
  continuous score alongside the binary call, as planned, is justified.

## Limitations

- **Label sources: Galperin et al. 2022 Table S1 was not reachable.** PMC serves the supplement
  behind a browser challenge, and Europe PMC reports the article as not open access.
- **Label sources: Browne et al. 2021 contributed only 5 experimentally labelled genomes**, all
  ethanol-negative and all correctly called non-sporulators. Its 1,358-genome classification is a
  genome prediction and was not used as truth.
- **BacDive was not queried separately**, because its API needs registration.
- **The labels are compiled species-level phenotypes**, not strain-matched assays. "Spores not
  observed" is often recorded as "non-sporulating", so the specificity estimate is conservative.
  All label errors were kept in the metrics.
- **Sampling is stratified by order and capped**, so the numbers are not prevalence-weighted for
  Bacillota.
- **Bootstrap CIs ignore phylogenetic non-independence**, so they are somewhat too narrow.
- **Rpf tier 1 has 25 scored Actinomycetota** (target ≥ 30). Four Machowski 2014 genomes have
  suppressed NCBI assemblies that are no longer downloadable. Tier 2 adds 95, but from a
  bioinformatic survey.
- **The non-Actinomycetota Rpf absences rest on a 2005 general statement**, not per-genome
  documentation.
- **Not verified:**
  - the protein identity of the staphylococcal PF06737 hits (probably SceD/IsaA-type)
  - the cause of the missing Spo0A in *Thermanaeromonas toyohensis*
  - the true phenotype of the 33/33 Bacillaceae labelled non-sporulating

## Reproduce

Run from the repo root. Large files go in `/tmp/val/`.

```bash
pip install pyhmmer pyrodigal pandas requests scikit-learn matplotlib scipy
```

Download the label sources into `/tmp/val/src/`:
- Madin et al. `condensed_traits_NCBI.csv` from github.com/jmadin/bacteria-archaea-traits
  (`output/`), saved as `madin_traits.csv`.
- Browne et al. 2021 Table S1 (Europe PMC, PMC8340488), saved as `browne_s1.tsv`.
- Puspita et al. 2013 supplement (Europe PMC, PMC4070687), converted to `s1.txt` with
  `pdftotext -layout`.

Then run, in order:

```bash
python3 validation/01_build_labels.py /tmp/val/src /tmp/val
```

```bash
python3 validation/04_build_rpf_labels.py
```

```bash
python3 validation/04b_build_rpf_puspita.py /tmp/val/src/pusp/s1.txt
```

```bash
bash validation/02_run_pipeline.sh /tmp/val/spo /tmp/val/val_accessions.tsv spo 4
```

```bash
bash validation/02_run_pipeline.sh /tmp/val/rpf /tmp/val/rpf_accessions.tsv rpf 4
```

`04b_build_rpf_puspita.py` does not write an accessions file, so create the input for the tier-2
run first:

```bash
python3 -c "import pandas as pd; t=pd.read_csv('validation/rpf_labels_puspita.tsv',sep='\t'); pd.DataFrame({'accession':'RS_'+t.accession,'ncbi_assembly_name':'','ncbi_translation_table':11}).to_csv('/tmp/val/rpf2_accessions.tsv',sep='\t',index=False)"
```

```bash
bash validation/02_run_pipeline.sh /tmp/val/rpf2 /tmp/val/rpf2_accessions.tsv rpf2 4
```

```bash
python3 validation/03_calibrate_sporulation.py /tmp/val/spo spo
```

```bash
python3 validation/05_rpf_check.py /tmp/val/rpf rpf /tmp/val/rpf2 rpf2
```

The NCBI assembly mapping for the sporulation labels is cached in
`/tmp/val/assembly_map_cache.json`. Newer NCBI releases may pick different assemblies.
