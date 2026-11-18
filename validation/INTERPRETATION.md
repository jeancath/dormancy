# What the validation means for the main analysis (interpretation)

This is my reading of the results in `README.md`, `sporulation_calibration.md` and
`rpf_check.md`. It is interpretation, not new analysis. Where it goes beyond the numbers, it says
so. Nothing here changes the pre-registered plan.

## 1. The sporulation score measures sporulation *capacity* well, and the binary call less well

- **The continuous score separates the label groups strongly.** `spo_strict` gives an AUC of 0.92,
  or 0.96 on better-supported labels, and its distribution is clearly bimodal: almost every genome
  has either 0–2 or 25–33 core genes (figure, right panel). Very few genomes fall in between. So the
  score captures a real biological split: lineages that kept the sporulation programme versus
  lineages that lost it wholesale.
- **The errors sit on that boundary, not in the bulk.** Nearly all false positives keep 25–33
  genes, and most false negatives are either Erysipelotrichia, which have a genuinely reduced but
  functional programme, or genomes with 0–2 genes and doubtful "sporeformer" labels.
- **Implication for the seed-bank test.** The project asks whether dormancy *capacity* scales with
  environmental variability. For that question, the gene content is arguably the better measure, and
  culture labels are the noisier one. A strain that keeps 31/33 genes but was never seen to
  sporulate in a lab medium may still sporulate in nature; the reverse is much less likely.

## 2. Specific recommendations, all within the existing pre-registration

1. **Keep the pre-registered rule (Spo0A & ≥ 25) as the primary binary call.** It is the
   high-sensitivity choice and has no false positives outside Bacillota.
2. **Run the pre-registered ≥ 27 (80 %) sensitivity analysis and report it next to the primary
   result.** The calibration shows that it is the Youden optimum on these labels and that it
   removes the clearest genuine false positives (*Caldicellulosiruptor*, *Bhargavaea*).
   - If the confirmatory result holds at both 25 and 27, it is robust to the one place where
     the threshold matters.
   - If it holds only at 25, the gut and thermophile clostridia scoring 25–26 are driving the
     result, and that should be stated.
3. **Report the continuous `spo_strict` model as well.** Given the bimodality, a binary call and the
   continuous score should usually agree. Disagreement would point to the 17–27 band.
4. **Watch Clostridia and Erysipelotrichia.**
   - Clostridia are where labels and genome calls disagree most (specificity 0.31).
   - In Erysipelotrichia, the binary rule misses real sporeformers.
   - A lineage-stratified or phylogenetically controlled model handles both, and the plan already
     uses phylogenetic comparative methods. A per-class breakdown of the main result would be a
     cheap robustness check. It is exploratory, not pre-registered.
5. **Restrict `n_rpf` to Actinomycetota for interpretation.** It reproduces published Rpf counts
   very well there. Outside Actinomycetota, a value of 1 can be the distantly related Rpf-like
   protein of *Staphylococcus* or *Oenococcus* rather than a classical Rpf. The marine analysis
   should check whether any non-Actinomycetota genomes with `n_rpf` > 0 drive a result.
6. **The relaxed sporulation score adds little.** Its AUC is 0.910 against 0.917 for strict. It
   rescues known near-miss genes in true sporeformers (spoIIE, ylbJ, yabP) but raises the scores
   of label-negative genomes just as much. Keep it as the pre-registered sensitivity analysis, but
   do not expect it to change conclusions.

## 3. What this validation does not show

- **It does not validate the *environmental* half of the hypothesis.** Only the genome-side
  measurement was tested.
- **The labels are not ground truth.** They are compiled species phenotypes, and I argue that a
  large share of the "false positives" are label problems. I cannot prove that without assays, so
  the reported specificity (0.78) is best read as a lower bound on the rule's specificity for
  sporulation *capacity*. That is an interpretation, not a measurement.
- **Marine genomes are under-represented in the calibration set**, which is mostly type strains
  from soil, gut, thermal and saline sources. Performance on marine Bacillota is assumed to be
  similar, not shown. A modest check would be to look at the score distribution of marine
  Bacillota in the main dataset and confirm it is also bimodal.
- **Switching mode (responsive vs stochastic) is untouched.** As the dictionary README says, gene
  presence cannot measure it.

## 4. One-sentence summary

The genome sporulation and Rpf calls are reliable measures of dormancy *capacity* for the main
test. The pre-registered threshold is defensible, and the planned ≥ 27 sensitivity analysis is the
one that matters most. The weakest points are clostridial and erysipelotrichial lineages and
`n_rpf` outside Actinomycetota, and all of these can be handled with analyses the plan already
contains.
