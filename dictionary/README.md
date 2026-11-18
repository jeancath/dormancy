# Dormancy gene dictionary (Phase 1)

Profile-HMM dictionary for scoring bacterial dormancy gene content across prodigal-predicted
proteomes. It feeds the seed-bank test: does dormancy gene content scale with environmental
variability at the isolation site? (Lennon & Jones 2011, doi:10.1038/nrmicro2504)

## Files

| File | Content |
|---|---|
| `dormancy_dictionary.tsv` | 176 rows: one row per (gene, model). Columns: module, gene, model_source, model_acc, model_name, cutoff_type, cutoff_value, kofam_score_type, switching_class, evidence_ref (DOIs, `;`-separated), notes |
| `thresholds.tsv` | One row per HMM in the library (174): `model_name, threshold, score_type` (required by the pipeline), then `domain_threshold, noise_cutoff, cutoff_type, model_acc, model_source` |
| `sporulation_core.tsv` | The sporulation gene set: gene, set (galperin2022_core / extended), in_score (yes/no), model_names |
| `/home/user/data/hmm/dormancy_library.hmm` | 174 concatenated, unmodified HMMER3 models (72 Pfam, 58 NCBIfam, 44 KOfam). Each model keeps its own NAME/ACC/GA/TC/NC lines. KOfam models have NAME only (no ACC/GA/TC): their cutoffs are only in `thresholds.tsv`. `hmmstat` reads all 174 models without error. |
| `/home/user/data/hmm/{pfam,ncbifam,kofam}/` | The individual model files. `/home/user/data/hmm/raw/` holds `ko_list.gz` and `hmm_PGAP.tsv` for provenance. |
| `build_dictionary.py` | Regenerates the three TSVs and the library from the individual HMMs (all choices and notes live in this script). |
| `check_library.py` | Reference implementation of the cutoff and competitive-assignment rules below (pyhmmer). Usage: `python3 check_library.py proteome1.faa [...]` |

Source releases: Pfam 39.0 (individual HMMs from the InterPro API; versioned accessions such as
PF06737.21 are in `thresholds.tsv`). NCBIfam is hmm_PGAP release 20.0 (2026-06-23); it holds both
the former TIGRFAMs (`TIGR*`) and `NF*` models, and every file carries TC = GA. The KOfam
`ko_list` is dated 2026-10-26, with profiles taken from the matching `profiles.tar.gz`.

## Applying the cutoffs (all modules)

1. Run `pyhmmer.hmmsearch` with the whole library and no bit-cutoff option, at a permissive E-value
   (for example E=10). **Do not use `bit_cutoffs="gathering"/"trusted"` globally:** KOfam models
   have no GA or TC lines. Then filter each hit with `thresholds.tsv`:
   - `score_type = full`: keep the hit if its full-sequence score ≥ `threshold`. If
     `domain_threshold` is set, at least one domain must also be ≥ `domain_threshold`. This
     reproduces `--cut_ga` (Pfam) and `--cut_tc` (NCBIfam).
   - `score_type = domain` (KOfam only): keep the hit if its best single-domain score ≥ `threshold`.
     This is the KofamScan rule.
   - `noise_cutoff` (NC, Pfam/NCBIfam only) is used only for the relaxed sporulation score.
2. A gene is present if any of its models passes. Many genes have a primary model plus a
   KOfam (or second NCBIfam) alternative.
3. Copy number is the number of distinct proteins that pass, counted once per protein even when
   several models of the gene hit it.
4. **Competitive assignment, sporulation modules only.** Each protein counts only for the gene whose
   model gives it the highest bit score, among all models in `sporulation_core` and
   `sporulation_extended`. This is required for the paralog pairs SigF/SigG, SigE/SigK and
   SpoVAC/SpoVAEB. Example: C. difficile SigE scores 257 on the SigK model PRK05803 (TC 252), but
   374 on the SigE model, so it is assigned to SigE.
5. Rows marked ARCHITECTURE-ONLY (G5, DUF348, TGS, ACT, the RSH SYNTH and HD domains) are never
   counted alone. They only type proteins into classes (rules below).

## Module rules

### 1. Endospore sporulation (`sporulation_core`, `sporulation_extended`)

**Core set.** The core set comes from Galperin et al. 2022, *J Bacteriol* 204:e00079-22,
Table 1 (doi:10.1128/jb.00079-22). It uses the genes in bold, whose products built the
concatenated sporulation tree. The PMC rendering of Table 1 shows 37 bold entries. The text says
"40 genes" (41 with YmfB), and the difference comes from paralog groups such as sspA/B/C/D,
spoVB/ykvU, ylmC/ymxH, yitE/yqfU and ydfS/yetF. We count each bold entry once. YmfB (a ClpP
paralog) is left out because ClpP is universal.

Of the 37 entries, 33 are scored. Four are kept in the dictionary but not scored, because their
models cannot separate sporeformers from other bacteria:
- spoIIIJ: a YidC2 paralog (YidC is universal).
- ydfS/yetF: the YetF_C domain is widespread.
- yitE/yqfU: the YitT_membrane domain is widespread.
- spoIVFB: M50 site-2 proteases are ubiquitous. The KOfam K06402 model gave a false positive
  in *M. tuberculosis* (Rip3) and missed *C. difficile* SpoIVFB.

**Sporulation score.** The score is the number of distinct scored core genes present (0–33),
with spo0A among them. The score has two versions:
- *strict* (primary): hits must pass the GA/TC/KOfam cutoffs.
- *relaxed* (sensitivity analysis): also counts a missing gene when one of its Pfam/NCBIfam
  models scores between NC and TC and that model is the protein's best sporulation hit.
  TIGRFAM TCs were tuned on Bacillus-heavy data, and clostridial orthologues often fall just
  below them. For example, C. difficile YabP scores 74.7 against a TC of 75.4, and C. difficile
  YlbJ scores 288 against a TC of 304.7.

**Proposed rule.**
- **sporulation-capable:** Spo0A present AND strict core score ≥ 25/33 (≥ 75 %).
- **reduced / uncertain:** Spo0A present AND score 17–24 (50–74 %).
- **non-sporulator:** Spo0A absent, OR score < 17.

Why these numbers:
- Spo0A is necessary but not sufficient. Galperin 2022 found it in all 76 sporeformers among
  180 Firmicutes, and also in 40 non-sporeformers. Its absence is "an excellent predictor" of
  inability to sporulate.
- Abecasis et al. 2013 (doi:10.1128/jb.02110-12) used a 48-gene genomic signature. Known
  endospore formers kept more than 75 % of it and non-formers kept less than 10 %. Organisms they
  predicted to sporulate, such as *Eubacterium rectale* (62 %), fell in between. Hence 75 %, with
  a 50–74 % grey zone.
- Galperin 2022 found that clostridial non-sporeformers lost at least 29 of the 40 analysed
  genes. The spore-forming Erysipelotrichia (*C. innocuum*, *E. ramosum*) were each missing more
  than a dozen genes of the broader conserved set, mostly the SpoIIQ–SpoIIIA channel genes
  plus spoIIR/sspA or spoVT. Under our rule they would probably fall in the reduced class.
  That is a known limitation, and it is the reason the continuous score should also be reported.
- **Calibrate before the main analysis.** Compute an ROC on labelled genomes (Galperin 2022
  Table S1, 76 sporeformers and 104 others; Browne et al. 2021, doi:10.1186/s13059-021-02428-6,
  1,358 Firmicutes) and pre-register N. Run sensitivity analyses at 60 % and 80 % and with the
  relaxed score.
- The rule only makes sense within Bacillota (Firmicutes). Outside Bacillota every scored model
  should give 0, and a non-zero score there should be treated as a QC flag.

**Sanity check** (UniProt reference proteomes run through `check_library.py`):

| Proteome | Spo0A | strict | relaxed | missing |
|---|---|---|---|---|
| *B. subtilis* 168 | yes | 33/33 | 33 | – |
| *C. difficile* 630 | yes | 30/33 | 32 | spoIIQ (clostridial SpoIIQ is divergent), ylbJ and yabP (both between NC and TC) |
| *E. coli* K-12 | no | 0 | 0 | – |
| *M. tuberculosis* H37Rv | no | 0 | 0 | – |
| *S. pneumoniae* R6 | no | 0 | 0 | – |

**Extended set** (reported, not scored): spoIIM, spoIIIAA/AB/AD/AF, spoVD, spoVE, dpaA/B,
GerA-family A and C subunits, cotE. dpaAB is not a necessary marker (many clostridia lack it).
cotE is specific to Bacillus. Count Ger receptors as loci, not genes.

### 2. Resuscitation (`resuscitation`)

- **rpf:** proteins whose best transglycosylase hit is PF06737 (Rpf domain). This is mostly
  Actinobacteria. In the check, *M. tuberculosis* gave 5 hits (RpfA–E).
- **RpfB-type:** PF06737 plus G5 (PF07501) plus DUF348 (PF03990) in one protein. K21688
  (RpfB) also hits B. subtilis YabE through the shared DUF348–G5 region, so it is called RpfB
  only when PF06737 is in the same protein.
- **Sps (YabE/YocH/YuiC/YorM; Ravagnani et al. 2005):** the "Sps" domain is the Pfam 3D domain
  (PF06725). 3D is also part of MltA lytic transglycosylases in Proteobacteria (E. coli MltA
  hits), so Sps is counted only in Bacillota, or when MltA (PF03562) is absent. Experimental
  support is weak (an unpublished survival phenotype), so switching_class is `unclear`.
- PrkC-type muropeptide-sensing PASTA kinases (Shah et al. 2008, doi:10.1016/j.cell.2008.08.039)
  are left out on purpose. PknB/PrkC kinases are present in nearly all Bacillota and
  Actinobacteria, so their presence is uninformative.

### 3. Ribosome hibernation (`ribosome_hibernation`)

- **long HPF:** a PF02482 hit with PF16321 (C-terminal dimerisation domain) in the same protein.
- **short HPF / RaiA:** TIGR00741, K05808 or K05809, without PF16321.
- **RMF:** PF04957 or K03812 (Gammaproteobacteria, about 55 aa).
- **Balon** (Helena-Bueno et al. 2024, *Nature* 626:1125, doi:10.1038/s41586-024-07041-8) is
  present if any of four models hits:
  - NF052287 (AOC03_06830 family, the *Psychrobacter urativorans* protein; UniProt A0A0M3V8U3)
  - NF052289 (Rv2629 family)
  - NF052288 (MSMEG_1130 family)
  - PF18845 (Balon_mid)

  The three NCBIfam models cite the Balon paper (PMID 38355796). In the check, *M. tuberculosis*
  Rv2629 was detected. PF18845 alone covers only about 1,650 proteins, far fewer than the paper's
  "20 % of bacteria". If Balon prevalence looks low, build a dedicated profile: jackhmmer seeded
  with A0A0M3V8U3 plus the paper's homologue set against UniRef50 (3 iterations), excluding
  aeRF1/Pelota (PF03463–PF03465) and the other baeRF families, then set a trusted cutoff on
  genomes known to carry or lack Balon.
- **MRF** (NF047431): mycobacterial hibernation under zinc depletion.
- **YqjD/ElaB** (DUF883): weak evidence, `unclear`.
- **HflX:** near-universal, so record copy number only.
- **Score.** Count the hibernation systems present: long HPF, short HPF, RMF, Balon, MRF.

### 4. Stringent response (`stringent_response`)

These rules follow the RSH classification of Atkinson et al. 2011 (doi:10.1371/journal.pone.0023479).
- **long RSH:** a protein hit by TIGR00691, K00951 or K01139, or by PF04607 plus PF02824 (TGS).
  Analyse long-RSH copy number, not RelA vs SpoT, because that split is specific to
  Beta/Gammaproteobacteria.
- **SAS (small alarmone synthetase):** K07816, or a protein ≤ 350 aa with PF04607 but no
  TGS, HD_4 or ACT.
- **SAH (small alarmone hydrolase):** a protein ≤ 300 aa with PF13328 (HD_4) but no PF04607.
  No dedicated model exists; this rule is based on Ruwe et al. 2018 (doi:10.3389/fmicb.2018.00916).
- **DksA:** TIGR02420. The generic TraR/DksA zinc finger is not used.
- Watch for type-IV "aTA" SAS toxins (toxSAS), which also carry PF04607. A PF04607-only
  protein next to a small antitoxin should be flagged rather than counted as SAS.

### 5. Persistence / stress-dormancy candidates (`persistence_stress`)

- **HipA:** PF07804 plus PF13657, or K07154. **PF07805 does not exist in current Pfam:** the
  N-terminal domain is PF13657 (Couple_hipA). HipA_C also occurs in non-TA kinases (YjjJ,
  HipA_2), so prefer HipA adjacent to a HipB/HTH gene.
- **TisB:** TIGR04459 and PF13939 are in the library, but TisB is 29 aa. Prodigal's minimum gene
  length (about 90 nt) means it is usually not called, and the IstR-1/2 sRNA antitoxin cannot be
  found with protein HMMs (it would need Rfam/Infernal). It is listed but **not analysable**.
- **ObgE** (TIGR02729): Verstraeten et al. 2015 is the evidence, but the gene is essential and
  universal. Use it only as a QC check (expect exactly 1 copy); it is never scored.
- **SutA** (NF041490): a slow-growth/dormancy transcription factor in Pseudomonas
  (Gammaproteobacteria only).

### 6. Toxin–antitoxin systems (`TA_typeII_toxin`, `TA_typeII_antitoxin`): analyse SEPARATELY

TA loci are confounded with plasmids, prophages and mobile-element load. In addition, the TA →
persister link has been seriously weakened (Goormaghtigh et al. 2018, doi:10.1128/mBio.00640-18;
the (p)ppGpp–Lon–TA persister model of Maisonneuve et al. 2013 *Cell* was retracted). Analyse
TA as its own module, with mobile-element load as a covariate. Never merge it into a dormancy
score.

- **Pairing rule.** A TA locus is a toxin hit and an antitoxin hit within ≤ 3 genes of each
  other (|gene index difference| ≤ 3) on the same contig, usually the same strand, often
  overlapping or adjacent. Count paired loci per genome. Report orphan toxins and orphan
  antitoxins separately and do not analyse them as TA.
- **Promiscuous families** (count only when paired):
  - PIN (VapC; *M. tuberculosis* gives 46 PIN hits)
  - Fic/DOC (most Fic proteins are AMPylators or effectors)
  - RES
  - ParD/RHH
  - AbrB/MazE: also hits SpoVT and AbrB. In the B. subtilis check, all 3 "MazE" hits were
    regulators, not antitoxins (AbrB, Abh, SpoVT). The one PIN hit was YacL, a non-TA PIN-TRAM
    protein. The real MazF homologue EndoA pairs with EndoAI, a CopG-type antitoxin that
    PF04014 does not detect. Always exclude proteins whose best hit is TIGR02851 (SpoVT).
  - HipB (low-specificity HTH KO)
- **No usable model:** GNAT acetyltransferase toxins (TacT/AtaT/ItaT; Cheverton et al. 2016)
  share Acetyltransf_1/7 with dozens of non-toxin GNATs per genome, so they are excluded.
- CcdB/CcdA and many ParE/RelE loci are plasmid-borne. Record whether each contig is a
  chromosome or a plasmid (e.g. geNomad/MOB-suite) when the assembly allows it.

### 7. Covariate: mobile-element load (`covariate_mobile`)

- **Transposases:** the number of distinct proteins with ≥ 1 hit to any `transposase:*` model
  (union, deduplicated per protein), per 1000 CDS. The accessory models (IS3 OrfA HTH, rve_2,
  rve_3, IS5 N-terminus, Mu C-terminus) only contribute through this protein-level union.
- **Integrases/recombinases:** PF00589, PF00239 or PF07508 hits, minus proteins whose best hit is
  XerC (TIGR02224) or XerD (TIGR02225).
- **Pitfalls:**
  - Y1/Y2 IS200 TnpA also match domesticated RAYTs.
  - rve matches phage integrases.
  - Short-read assemblies collapse IS copies, so transposase counts are biased downward in
    fragmented genomes. Include N50/contig number as covariates, or run a sensitivity analysis on
    complete genomes only.

### 8. Covariate: general stress (`covariate_stress`)

- RpoS (K03087; TIGR02394 in Pseudomonadota) is an annotation covariate only, for VBNC and
  stationary-phase capacity (Lange & Hengge-Aronis 1991; Boaretti et al. 2003). RpoS is
  essentially restricted to Proteobacteria, so it is confounded with phylogeny.

## Responsive vs stochastic switching (`switching_class`)

The labels follow the theory of Kussell & Leibler 2005 (doi:10.1126/science.1114383): responsive
switching (sensing) pays off when environmental cues are reliable and sensing is cheap, and
stochastic switching (bet hedging) pays off under rare or unpredictable change.

| Class | Meaning | Families | Key evidence |
|---|---|---|---|
| responsive | entry is induced by a sensed signal | sporulation (starvation → phosphorelay → Spo0A~P), RSH/SAS/DksA, RpoS, TisB (SOS-induced), SutA | Burbulys 1991; Xiao 1991; Paul 2004; Dörr 2010; Lange 1991 |
| stochastic | entry is noise-driven in a subpopulation without a trigger | HipA (hipA7 type II persisters), ObgE (bet-hedging claim) | Balaban 2004, doi:10.1126/science.1099390; Rotem 2010, doi:10.1073/pnas.1004333107; Verstraeten 2015 |
| exit_trigger | resuscitation / germination | Rpf, RpfB, germinant receptors, GPR, HflX (100S disassembly) | Mukamolova 1998; Kana 2008; Paidhungat & Setlow 2000; Basu & Yap 2017 |
| hibernation | translational shutdown of the existing cell | HPF, RaiA, RMF, Balon, MRF | Ueta 2005/2013; Helena-Bueno 2024; Li 2018 |
| unclear | role or mode not established, or covariate rows | Sps, YqjD, most TA toxins/antitoxins, mobile elements | Goormaghtigh 2018 for TA |

**Honest caveat: gene presence cannot discriminate switching mode.** The same machinery can
produce both behaviours depending on network parameters and alleles, not gene content:
- The Spo0A phosphorelay responds to starvation, yet sporulation in clonal B. subtilis cultures
  is bistable and shows bet hedging (Chung et al. 1994, doi:10.1128/jb.176.7.1977-1984.1994;
  Veening et al. 2008, doi:10.1073/pnas.0700463105).
- High persistence comes from the hipA7 allele and expression thresholds (Rotem 2010), not from
  owning hipA.
- Germination and resuscitation can also be heterogeneous (Sturm & Dworkin 2015,
  doi:10.1016/j.cub.2015.07.018; scout hypothesis, Epstein 2009, doi:10.1038/4571083a).

`switching_class` is therefore a literature label for interpreting results. It is not a
measured trait. The defensible test is dormancy *capacity* (module presence or score) against
environmental variability. Claims about the responsive vs stochastic mode need phenotypic data
such as single-cell time-lapse or persister assays.

## General pitfalls

- **Short ORFs:** prodigal under-calls small proteins (TisB, RMF, SASP, SpoVM, HicA). Absence of
  a short gene is weak evidence.
- **Split genes:** sigK is interrupted by a skin element in B. subtilis and C. difficile. In
  genomic prodigal calls, each fragment can fall below TC. If a strict score is exactly one
  short, check for SigK fragments.
- **MAG/SAG incompleteness:** absences are uninformative in incomplete genomes. Restrict to
  ≥ 90 % completeness (CheckM2) or model the score with completeness as a covariate.
- **Phylogeny:** presence/absence is phylogenetically structured. Use phylogenetic comparative
  models (PGLS / phylogenetic logistic regression), and keep the sporulation analysis within
  Bacillota.
- **Dormancy forms not covered:** actinomycete exospores (whi/bld), myxospores, *Azotobacter*
  cysts, cyanobacterial akinetes and Bdellovibrio bdelloplasts. A genome scoring 0 here is not
  necessarily non-dormant.
- **Families without a usable model:**
  - GNAT TA toxins (TacT/AtaT)
  - IstR sRNA
  - TisB (model exists but the ORF is not called)
  - a dedicated SAH model (architecture rule used instead)
  - reliable SpoIVFB and clostridial SpoIIQ (models miss them)
  - SRA/protein D: a gammaproteobacterial hibernation factor, omitted because its NCBIfam file
    (NF019746) was not on the FTP server
  - Balon via Pfam alone has limited sensitivity (NCBIfam models added).
