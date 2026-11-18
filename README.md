# Dormancy-associated gene content and seasonal temperature amplitude across 7,327 bacterial species

[![DOI](https://zenodo.org/badge/1463415847.svg)](https://doi.org/10.5281/zenodo.24355693)

Code, derived data and results for the pre-registered study *"Dormancy-associated genomic repertoire
across bacterial taxa and the seasonal temperature amplitude of their sampling environments: a
pre-registered phylogenetic test across 7,327 species"*. The manuscript is in [`paper/`](paper/).

## Summary
**Question.** Seed-bank theory predicts that microbial dormancy is favoured where environments
fluctuate. Does the dormancy-associated gene content of individual bacterial genomes increase with
the temperature variability of their isolation sites?
- **Dormancy genes:** sporulation, resuscitation-promoting factors, ribosome hibernation and the
  stringent response.
- **Variability:** the 30-year monthly temperature SD, which at these sites is essentially the
  seasonal amplitude.
- **Models:** phylogenetic regression (Pagel's λ) on the GTDB r232 tree.

**Main results.**
- **Confirmatory.** The pooled slope (+0.0065 SD of the dormancy composite per SD of variability)
  is statistically equivalent to zero within the pre-registered margin of ±0.015.
  - The test is a cluster-robust jackknife over 0.5° cells.
  - p_TOST = 0.0005 for the originally specified procedure, and 0.012 for the averaged procedure
    adopted afterwards. Equivalence requires both.
- **Descriptive.** Gene content is strongly structured by relatedness: λ = 0.94, with 95 % of the
  residual variance phylogenetic.
- **Exploratory.**
  - The slope is positive in species represented by metagenome-assembled genomes (+0.017, 95 % CI
    0.004–0.031, not shown to exceed the margin) and weakly negative in species represented by
    cultured isolates.
  - The difference is robust to adjustment for phylum, habitat, latitude and mean temperature.
  - It was not detected within species, but within-species leverage is limited.
- **Marine genomes carry fewer toxin–antitoxin loci** than land genomes of the same size. The
  deficit is concentrated in metagenome-assembled genomes.

## Contents
| Path | What it holds |
|---|---|
| `paper/` | Manuscript (PDF and Markdown), figures, supplementary tables |
| `prereg/PREREGISTRATION.md` | The pre-registration, frozen before any confirmatory model was fitted |
| `provenance/prereg_freeze.txt` | Freeze commit, date and checksum showing the plan is unchanged |
| `scripts/` | Numbered analysis pipeline, models, calibration and sensitivity analyses |
| `dictionary/` | Dormancy gene dictionary: HMM library, cut-offs and evidence |
| `data_small/` | Derived tables: habitat classes, environmental variability, analysis set, pruned tree, final analysis table |
| `results/` | Model outputs, figures, results report (`REPORT.md`), deviation log (`deviations.md`), final genome-source table, reconciliation of superseded numbers |
| `validation/` | Phenotype calibration of the sporulation and Rpf calls |
| `replication_within_species/` | Within-species replication (7,098 high-quality genomes) |
| `audit/` | Computational audit: re-computation of every reported number, pre-registration audit, code review, refits |
| `docs/phase0/` | Literature search: queries, screened corpus and gap memo |

## Reproducing
The scripts download the large public inputs again into a directory outside the repository:
- GTDB R11-RS232 metadata and bac120 tree;
- NCBI assemblies;
- NOAA OISST v2.1, GHCN_CAMS and PREC/L;
- ESA OC-CCI v6.0.

Run the steps in this order:
1. `01`: classify isolation habitats.
2. `02`, `02b`, `02c`: compute environmental variability.
3. `03`: download genomes and call genes.
4. `05`: HMM search.
5. `06`, `06b`: growth rates.
6. `07`: prune the tree.
7. `08`: aggregate genes into modules.
8. `10`: build the analysis table.
9. `09`: fit the confirmatory models.
10. `09` with the exploratory configuration: fit the exploratory models.
11. `14`, `16`, `18`: extra tests, site robustness and the MAG/isolate dissection.
12. `13`, `17`: figures and decision table.
13. `22`–`34`: equivalence tests, their calibration and sensitivity analyses.

`data_small/final/` holds the final analysis table and the per-genome module counts, so the models
can be refitted without the raw inputs. `python3 audit/number_check_v2.py` re-computes every number
reported in the manuscript from the result files.

## Citation and licence
- **Citation:** see `CITATION.cff`. Please cite the article and this archive: version 1.0.0,
  doi:[10.5281/zenodo.24355694](https://doi.org/10.5281/zenodo.24355694) (all versions: doi:[10.5281/zenodo.24355693](https://doi.org/10.5281/zenodo.24355693)).
- **Licence:** code under MIT; data, results and text under CC BY 4.0 (`LICENSE`).
- **Use of AI tools:** the author found the idea, critiqued it, searched the literature, reviewed
  the work, and analysed and interpreted the results. A generative AI assistant also did some of
  this work, including writing code.
