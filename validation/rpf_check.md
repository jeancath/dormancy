# Rpf (and Sps) check against literature-documented gene counts

This compares `n_rpf`, `n_rpfB` and `n_sps` from the unchanged `scripts/08_aggregate_modules.py`
with gene counts published for named genomes. The dictionary rule is: a protein counts as Rpf when
its best transglycosylase hit is Pfam PF06737 and passes the GA cutoff (27 bits). Sps (Pfam 3D,
PF06725) is counted only in Bacillota.

Inputs:
- `rpf_labels.tsv` (tier 1) and `rpf_labels_puspita.tsv` (tier 2), built by
  `04_build_rpf_labels.py` and `04b_build_rpf_puspita.py`.

Outputs, from `05_rpf_check.py`:
- `rpf_check_table.tsv` and `rpf_check_table_puspita.tsv` (per-genome tables)
- `rpf_check_numbers.json` (summary numbers)
- `rpf_puspita_overlap_tier1.tsv` (genomes in both tiers)

## Sources

| Tier | Source | What it documents | Genomes used |
|---|---|---|---|
| 1 | Ravagnani, Finan & Young 2005, *BMC Genomics* 6:39, doi:10.1186/1471-2164-6-39 | Table 1A: Rpf-domain genes per genome (12 Actinobacteria). Table 1B: genomes with only a *distantly related* domain, which their Rpf HMM did not detect. Text: Rpf absent from *Bifidobacterium*, *Thermobifida* and *Tropheryma*, and "confined to the actinobacteria". Table 2A: Sps genes per genome in 18 Firmicutes. | 16 Actinomycetota, 20 Bacillota, 8 other-phylum controls |
| 1 | Machowski et al. 2014, *BMC Microbiol* 14:75, doi:10.1186/1471-2180-14-75 | Table 1: rpfA–E locus tags in 19 mycobacterial genomes | 13 further mycobacteria; 4 could not be downloaded (below) |
| 1 | Mukamolova et al. 1998, *PNAS*, doi:10.1073/pnas.95.15.8916 | *M. luteus* Rpf | 1 (also in Ravagnani) |
| 2 | Puspita et al. 2013, *Microbes Environ* 28:58, doi:10.1264/jsme2.me12122 | Table S1: rpf copies per genome, by BLASTp with the Rpf domain (Actinobacteria) | 95 finished genomes not in tier 1 |

**Tier 1 is literature documentation** (counts stated in papers, mostly from gene-by-gene analysis).
**Tier 2 is a bioinformatic survey** that uses a different method (BLASTp) from our Pfam HMM. Its
counts are not experimental truth, and agreement with it means the two methods agree, not that both
are right.

**Excluded rows**
- **Puspita table.** The PDF's accession column is misaligned for some rows; for example, a
  *Corynebacterium pseudotuberculosis* row points to the *Rhodococcus jostii* chromosome. Rows whose
  linked assembly did not carry the strain name in the table were dropped: 111 of 136 finished
  genomes were consistent, and 16 of those were already in tier 1.
- **Tier 1 download failures.** *Mycobacterium* sp. MCS, KMS and JLS and *M. gilvum* PYR-GCK
  (Machowski Table 1) have NCBI assemblies that are now suppressed and are no longer on the FTP
  server. They are listed in `rpf_labels.tsv` but **not scored**.

## Results, tier 1 (53 genomes scored)

**Actinomycetota (25 scored; 24 with a documented count).**

| | |
|---|---|
| exact count match | **23 / 24** |
| presence/absence correct | 24 / 24 |
| the one mismatch | *Mycolicibacterium vanbaalenii* PYR-1: documented 4, predicted 5 |
| documented absences (*Bifidobacterium longum* NCC2705, *Tropheryma whipplei* Twist, *Thermobifida fusca* YX) | all 0 |
| *Saccharopolyspora erythraea* (presence only in Ravagnani) | predicted 4; Puspita gives 4 |

Every multi-copy count reproduces exactly:
- *M. tuberculosis* H37Rv, H37Ra and CDC1551, *M. bovis* AF2122/97 and BCG Pasteur, and
  *M. africanum*: 5 (RpfA–E).
- *M. avium* 104 and K-10, *M. intracellulare*, *M. ulcerans*, *M. marinum*, *M. abscessus* and
  *M. smegmatis*: 4.
- *M. leprae*: 3. *C. diphtheriae*: 3. *C. glutamicum* and *C. efficiens*: 2. *M. luteus*: 1.
- *S. coelicolor*: 5. *S. avermitilis*: 6.

`n_rpfB` (Rpf + G5 + DUF348 in one protein) is 1 in every genome above except *M. luteus*, where
it is 0. That is consistent with RpfB being a single-copy gene in mycobacteria, corynebacteria and
streptomycetes. *M. luteus* has a single classical Rpf, and its RpfB-like protein (mentioned in the
replication plan) is not called RpfB by our architecture rule. That is a disagreement with the replication plan's
expectation, and **not verified** further.

**Bacillota (20 genomes, documented Rpf = 0).**
- 17 / 20 correctly give 0.
- **3 give n_rpf = 1:** *Staphylococcus aureus* N315 (score 89.8), *S. epidermidis* ATCC 12228
  (86.9) and *Oenococcus oeni* PSU-1 (33.1).
  - These are exactly the genomes Ravagnani lists in Table 1B as carrying one protein with a domain
    "distantly related" to Rpf, which their own Rpf HMM did not detect.
  - Current Pfam PF06737 does detect them above its gathering threshold of 27. In the staphylococci
    this is most likely the SceD/IsaA-type lytic transglycosylase family, but I **did not verify**
    the protein identity.
  - So these are not random false positives: the current PF06737 model is broader than the 2005
    Rpf definition. Outside Actinomycetota, `n_rpf` can therefore count these distant Rpf-like
    proteins. The main analysis should treat `n_rpf` outside Actinomycetota with care, or report it
    within Actinomycetota only.
- *Bifidobacterium longum* and *Tropheryma whipplei*, the Actinomycetota with "distantly related"
  proteins in Table 1B, both correctly score 0.

**Other phyla (8 controls).** *E. coli*, *P. putida*, *V. cholerae*, *Caulobacter*, *Deinococcus*,
*Synechococcus*, *Prochlorococcus* and *Pelagibacter* all give 0.

**Sps (YabE/YocH/YuiC-type, Bacillota only).**
- **16 / 20** Bacillota match the documented count exactly. That includes 0 for the two
  staphylococci and *O. oeni*.
- The other 4 over-count by exactly one:

  | Genome | Documented | Predicted |
  |---|---|---|
  | *B. cereus* ATCC 14579 | 5 | 6 |
  | *O. iheyensis* HTE831 | 4 | 5 |
  | *C. acetobutylicum* ATCC 824 | 2 | 3 |
  | *D. hafniense* DCB-2 | 1 | 2 |

  For *O. iheyensis*, Ravagnani separately lists one protein with a "distantly related" Sps domain,
  which would explain the extra copy. *D. hafniense* DCB-2 and *A. thermocellus* were draft genomes
  in 2005, so their counts may differ from the finished genomes.
- No genome outside Bacillota gets a non-zero `n_sps`, as expected from the rule. This only
  confirms the rule's design: 08 sets `n_sps` to 0 outside Bacillota.

## Results, tier 2 (Puspita 2013; 95 Actinomycetota)

| | |
|---|---|
| exact count match | **79 / 95 (83 %)** |
| within ±1 | 93 / 95 (98 %) |
| over-count / under-count | 15 / 1 |
| Spearman ρ (documented vs predicted) | **0.96** |
| presence/absence | TP 56, FN 0, FP 4, TN 35 |

The two methods agree on Rpf presence in every genome Puspita calls positive.

**Discordances**

- **Absent in Puspita but predicted present (4):** *Xylanimonas cellulosilytica* (2 predicted),
  *Conexibacter woesei*, *Sanguibacter keddieii* and *Thermomonospora curvata* (1 each).
- **Largest under-count:** *Beutenbergia cavernae*, documented 3, predicted 1.
- **Pattern:** most discordances are +1, i.e. our HMM finds one more Rpf protein than BLASTp did.
  That is the expected direction for a profile HMM, which is more sensitive to divergent members
  than a BLASTp search with a single domain sequence.
- **Agreement between the two published sources.** In the 15 genomes counted by both Puspita and
  tier 1 (`rpf_puspita_overlap_tier1.tsv`), 13 counts match. Puspita gives one fewer for
  *S. coelicolor* (4 vs 5) and *S. avermitilis* (5 vs 6). So the published counts themselves carry
  about ±1 uncertainty in some genomes.

## Conclusion

- **Within Actinomycetota, `n_rpf` reproduces documented Rpf counts almost perfectly.** Tier 1
  matches exactly in 23 of 24 genomes. Tier 2 matches exactly in 79 of 95 genomes and within ±1 in
  93, with no missed Rpf-positive genome.
- **Specificity outside Actinomycetota is limited by one known family.** The current PF06737
  model also detects the distantly related Rpf-like proteins of *Staphylococcus* and *Oenococcus*,
  so `n_rpf` > 0 outside Actinomycetota should not be read as classical Rpf without checking.
- **Limitations:**
  - Tier 1 has 25 scored Actinomycetota, short of the replication plan's target of ≥ 30 for literature
    counts. Tier 2 provides 95 more, but from a bioinformatic survey.
  - The non-Actinomycetota "absences" for Pseudomonadota, cyanobacteria and *Deinococcus* rest on
    Ravagnani's 2005 general statement, not a per-genome check. Lennon et al. 2025 (*mSystems*)
    note that Rpf may also occur outside Actinomycetota.
  - The dictionary and thresholds were not changed.
