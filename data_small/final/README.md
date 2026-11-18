# Final analysis inputs

| File | Rows | What it is |
|---|---|---|
| `analysis_table.tsv.gz` | 7,715 | Final analysis table from `scripts/10_build_table.py`, after QC. Filter `pilot == 0` to get the 7,327 confirmatory species. |
| `modules.tsv.gz` | 7,867 | Per-genome dormancy, TA and mobile-element module counts (`scripts/08_aggregate_modules.py`). |
| `grodon.tsv.gz` | 7,867 | gRodon2 predicted doubling time and number of highly expressed (HE) genes per genome. |

To reproduce every model without the large raw inputs, run from the repo root:

```bash
Rscript scripts/09_models.R data_small/final/analysis_table.tsv.gz data_small/analysis_tree.nwk confirmatory OUTPREFIX [CONFIG]
```
