# Adversarial literature second opinion: The Dormancy Rule

**Date:** 2026-11-17. **Reviewer role:** adversarial (find what undercuts novelty or interpretation).
**Search effort:** 55 Europe PMC REST queries (incl. SRC:PPR preprints and author runs), 8 Crossref
queries, 4 bioRxiv-API status checks (api.biorxiv.org has no keyword search, so bioRxiv was keyword-
searched through Europe PMC's PPR index). All queries and hit counts are in `literature_queries.tsv`;
study scales in `study_scale_comparison.tsv`. **Every DOI below was resolved against the Crossref
`/works/{DOI}` endpoint on 2026-11-17 and the title/year/journal matched** (see `doi_verified` column).
Abstracts, and full text where open, were read for each numeric claim; guesses are marked.

**Bottom line.** The core novelty claim survives: I found **no 2024–2026 paper or preprint that
regresses genome-level dormancy gene content (sporulation / Rpf / hibernation / stringent response /
TA) on a continuous measure of isolation-site temporal temperature variability under phylogenetic
control**. This agrees with the project's own Nov-2026 novelty memo. But several new items narrow the
novelty, and — more seriously — the three headline *interpretations* (MAG/isolate reversal, marine TA
deficit, Balon/cold) each have close competitors or outright contradictions that the paper should
engage before publication. The Lennon & Jones (2011) citation for the marine TA deficit is
**not supported** (detail at the end).

> **Lead auditor's note (added after review of this subagent report).**
>
> **DOI check.** I spot-checked the DOIs of items A1–A3 against Crossref. They resolve, and the
> titles match.
>
> **I rate A1 as *major*, not critical.**
> - What Mise & Iwasaki (2022) report: ribosomal-protein genes go *missing* from MAGs through a
>   codon-usage and binning effect.
> - Direction: a recovery loss that worsens at variable sites would push gene counts *down*. In
>   this dataset that is what MAG ribosomal HE genes do (−0.005 per SD, p = 1e-4). The hibernation
>   paralogues go *up* (+0.013).
> - Conclusion: a single generic recovery artefact cannot produce both signs, so this paper argues
>   *against* a recovery artefact behind the paralogue increase.
>
> **What it does undermine** is the use of ribosomal-protein genes as a *neutral* negative
> control, and the n_HE ≥ 20 QC rule, because both are known to be biased in MAGs. A better
> negative control would be single-copy genes with typical codon usage (e.g. tRNA synthetases).
>
> **The independent tests in this audit:**
> - `strain_heterogeneity_tests.tsv`: CheckM strain heterogeneity *falls* with P0 in MAGs, and the
>   MAG slope is as large in MAGs with zero heterogeneity;
> - `sag_test.tsv`: the marine paralogue reversal comes from isolates losing copies, not from MAGs
>   gaining them.
>
> Both of these point away from an explanation based on population consensus or binning, and
> towards differences between cultured and uncultured lineages. The within-species replication
> (`replication_within_species/`) is the decisive test.

---

## A. Threats to the INTERPRETATION of the MAG vs isolate reversal (most important)

### 1. [CRITICAL] The reversal is the headline, yet the "population-consensus / paralog" reading has a live, published mechanism that predicts exactly this sign and this module.
- **Mise & Iwasaki 2022, ISME Commun 2:118, doi:10.1038/s43705-022-00204-6 (verified).** Across
  >190,000 MAGs, ribosomal-protein genes are missing from 20–40% of near-complete MAGs because of
  codon-usage × tetranucleotide binning, **and the loss is worse in fast-growing lineages** — i.e.
  recovery bias is *correlated with life-history strategy*, which co-varies with dormancy. This is
  the single strongest argument that a MAG-specific positive slope concentrated in *accessory/
  paralogous* hibernation copies can be a recovery/compositional artifact rather than biology.
- **Meziti et al. 2021, AEM 87:e02593-20, doi:10.1128/AEM.02593-20 (verified).** ~95%-complete MAGs
  recover only 77% of core and 50% of population-variable genes. Paralogous HPF/RaiA/Balon copies are
  exactly the variable compartment.
- **Why it matters:** the paper's own §3.1/§3.6 negative controls (near-universal genes don't rise;
  clean/near-complete MAGs attenuate to +0.007 with CI still including +0.018) are honest but
  *underpowered*, and the literature gives a concrete, trait-correlated artifact that the controls
  cannot fully exclude. Severity is critical because this is the paper's flagship exploratory result.

### 2. [MAJOR] A 2025 preprint directly reports the MAG-vs-isolate functional gap in the same direction and attributes it to recovery, weakening any "never reported" framing.
- **Rühlemann et al. 2025, bioRxiv, doi:10.1101/2025.10.20.683395 (verified; preprint, not yet
  published — confirm before citing).** States plainly that *"compared to MAGs, paired isolate
  genomes recover more functional elements and signals of horizontal gene transfer."* The present
  study's "MAGs are gene-richer at variable sites" runs opposite in sign, so this is not a scoop, but
  it establishes that **systematic MAG-vs-isolate functional-content divergence is already an
  expected, documented phenomenon** — the novelty is the *environmental interaction*, not the
  divergence itself. The paper should claim only the interaction.
- **Gupta & Almeida 2025, Nat Commun 16:9959, doi:10.1038/s41467-025-64950-6 (verified).** MAGs and
  isolates of one species sample different diversity (214 MAG-only genes; >60% of MAGs new STs):
  supports "different clades/populations," a competing explanation for the reversal.

### 3. [MAJOR] A GTDB-scale study finds MAGs are only ~3.7% smaller than co-clustered isolates — cutting against a large systematic MAG/isolate gene-content gap.
- **Rodríguez-Gijón et al. 2022, Front Microbiol 12:761869, doi:10.3389/fmicb.2021.761869
  (verified).** 26,101 genomes; in ANI clusters containing both, MAGs averaged only 3.7% smaller than
  isolates. This is double-edged: it supports the paper (gross size bias is modest) but also means a
  *large* MAG/isolate slope difference needs a module-specific, not whole-genome, mechanism — which
  is precisely the paralog compartment that artifact arguments (threat 1) target.
- **Zhu et al. 2026, Adv Sci, doi:10.1002/advs.77591 (verified).** Shows method choice alone can flip
  average-genome-size–environment relationships and that a recent salt-stress "eco-evolutionary
  trade-off" may be a method artifact — a cautionary precedent that method-dependent trait–environment
  signs (exactly this study's MAG/isolate split) can be artifactual.

### 4. [MINOR] Culturing-selection / collider-bias reading is well supported but not unique.
- **Weissman, Hou & Fuhrman 2021, PNAS 118:e2016810118 (verified)** (culture collections biased to
  fast growers) and **Rodríguez del Río et al. 2026, BMC Genomics 27:655,
  doi:10.1186/s12864-026-13092-0 (verified)** (uncultivated soil taxa: smaller genomes, slower growth,
  biogeographically structured by latitude/climate) both support the collider story — but they equally
  support "different clades," so they cannot adjudicate between the paper's remaining candidates. The
  paper correctly flags that its is_isolate ~ P0×DORM model is just the mirror of the interaction.

---

## B. Threats to the MARINE TA DEFICIT and its interpretation

### 5. [MAJOR] The deficit's direction is contradicted by the one genome-size-corrected baseline (Makarova 2009), and confounded by a temperature effect (Danov 2024).
- **Makarova, Wolf & Koonin 2009, Biol Direct 4:19 (verified).** As a residual after genome size,
  *"aquatic" (marine + freshwater pooled) genomes carry MORE TA.* The paper does cite this as reversed,
  but an adversary notes the reversal rests on the paper splitting marine from freshwater — a split the
  paper asserts resolves it but does not fully demonstrate. State the freshwater result explicitly.
- **Danov et al. 2024, mBio 15:e01911-23 (Toxinome; verified).** Across 59,475 genomes, toxins are
  depleted **at both hot and cold temperature extremes.** Marine SST sites in this study are the cold
  tail (median SST SD 3 °C; cold-sea isolates feature in the Balon result), so part of a "marine" TA
  deficit could be a **temperature** effect, not a marine/accessory-DNA effect. Not previously flagged
  in the project memo; worth a covariate check.

### 6. [MAJOR, supportive-but-reframes] Marine defense-system depletion is real and binning-independent — but it is a *defensome/MGE* story, which argues the deficit is NOT about dormancy.
- **Meaden, Westra & Fineran 2025, Phil Trans R Soc B 380:20240069 (verified).** Contig-level,
  binning-independent: *"marine environments appear to be particularly depleted in defence systems."*
  This is the strongest evidence the marine deficit is real, not a binning artifact.
- **Beavogui et al. 2024, Nat Commun 15:2146 (verified).** ~65% of marine MAGs have a "limited
  defensome"; authors flag MAG under-recovery of MGE-borne defenses as a likely bias (supports the
  paper's "MAG estimate is inflated" claim).
- **Bograd et al. 2025, Genome Biol 26:163, doi:10.1186/s13059-025-03641-3 (verified).** Across 38,912
  isolate genomes + 6,073 metagenomes, plasmids, prophages and defense systems are depleted in *plant-
  associated* bacteria, taxon-independently. Combined with Meaden/Beavogui/Tesson, this shows a broad
  pattern: **accessory defense/MGE content tracks habitat generally.** The paper already concludes the
  TA deficit is "part of a general reduction of accessory defence and mobile DNA" — this literature
  strongly supports that framing and argues *against* reading it as reduced dormancy capacity. Keep
  the deficit as a defensome finding, not a dormancy finding.
- **Liu et al. 2025, Mar Drugs 23:436, doi:10.3390/md23110436 (verified)** is only descriptive (MAG-
  based, no terrestrial comparator, no size/MGE normalisation) — not a competitor.

**Net:** the *fact* of a marine TA/defense deficit is now well corroborated (minor novelty loss), but
the genome-size + MGE-corrected, isolate-vs-MAG decomposition remains unreplicated elsewhere. The
*interpretation* should lean fully into "accessory/mobile defensome," not dormancy.

---

## C. Threats to the BALON / cold association

### 7. [MAJOR] The Balon–cold link is an inference from a single organism in the discovery paper; no genome-scale test exists, so the paper's finding is novel BUT its mechanistic framing over-reads Helena-Bueno.
- **Helena-Bueno et al. 2024, Nature 626:1125, doi:10.1038/s41586-024-07041-8 (verified; full text
  read).** Balon is in 1,573 / 8,761 species (23/27 phyla) and multi-copy (2–4) in 603 species (38%).
  The paper discovered Balon in a cold *Psychrobacter* and inferred a cold/stress role, **but ran no
  temperature or habitat association across genomes** — the "cold-associated" label is organism-level,
  not comparative. So this study's cold-sea-isolate Balon result is genuinely new, but it should not
  be presented as "confirming" Helena-Bueno; Helena-Bueno made no such comparative claim.
- **[watch] Chan/Melnikov et al. 2025, bioRxiv doi:10.1101/2025.11.04.686544 (verified; also Research
  Square rs-8030510). I could NOT read the full text (bioRxiv 429 rate-limit; EuropePMC/RS full-text
  XML returned HTTP 500).** From the abstract only: across 46,015 genomes, hibernation factors
  undergo size reduction, extinction in >10% of bacteria, and replacement via HGT / de novo birth /
  fusion. This is the closest work to the paper's "reversal sits in hibernation paralogs" result and
  could partially explain paralog gain/loss patterns by HGT rather than environment. **Flagged as
  unverified-content; obtain full text before the paper relies on or distinguishes itself from it.**
- **Madru et al. 2026, Nat Commun 17:5751, doi:10.1038/s41467-026-72341-8 (verified)** and
  **Souza et al. 2026, EMBO J, doi:10.1038/s44318-026-00918-6 (verified)**: archaeal HPF/RaiA-family
  factors — context that hibernation-factor repertoires are modular and lineage-specific, reinforcing
  that paralog counts are evolutionarily labile (a caution for reading them as environmental tuning).

---

## D. Threats to NOVELTY of the core test and framing

### 8. [MINOR] The "predictable vs unpredictable variability" decomposition is not novel to ecology.
- Already conceded in the project memo (Colwell 1974; Marshall & Burgess 2015; Barneche 2018). The
  paper's §3.3 honestly reports it is non-identifiable in surface temperature (land r=1.00). Keep the
  claim narrowly to "first application to microbial genomes," and note it failed.

### 9. [MINOR] Several 2026 works show the same "phylogeny dominates environment" result for other traits — the substantive finding is part of a trend, not isolated.
- **Henderson et al. 2026, bioRxiv doi:10.64898/2026.06.03.729626 (verified; preprint).** Amino-acid
  auxotrophy: *"phylogeny is the dominant predictor… vastly outweighing environmental context."*
  Same phylogeny-first conclusion via phylogenetically-informed models. Good corroborating cite.
- **Favier et al. 2026, Ecol Lett 29:e70483, doi:10.1111/ele.70483 (verified).** Soil copiotroph
  thermal limits show "minimal influence of present-day bioclimatic variables" (251 isolates, 15
  sites, taxonomic random effects only). Corroborates a weak present-day climate signal but is itself
  weakly powered — do not lean on it as strong support.
- **Weissman et al. 2026, ISME J 20:wrag147, doi:10.1093/ismejo/wrag147 (verified).** Genomic
  copiotrophy signatures decouple from growth rate *along temperature gradients* — a confound for this
  study's temperature predictor and its growth-rate/rRNA covariates. Worth a sentence.

### 10. [MINOR] No prior genome-level site-random-effect model found — novelty of the methodological point holds, but the project should not overclaim.
- **Hao et al. 2026, PNAS 123:e2520645123 (verified; full text read):** 62,291 GTDB genomes obtained,
  25,785 species analysed; non-independence handled at the *habitat-cluster* level (MCMCglmm + PGLS),
  not per-site, because predictors are species-level. Confirms the project's claim that no comparable
  study faces the per-site pseudoreplication problem — but also shows the field's answer is to use
  species-level predictors, which is a legitimate alternative the paper should acknowledge rather than
  present site-random-effects as the only correct route.

---

## E. The Lennon & Jones (2011) check (task 4)

**Claim in paper (Discussion, "The marine TA deficit"):** the marine TA deficit "agrees with the
metagenomic observation of Lennon & Jones (2011)."

**Finding: NOT supported as a *marine* TA observation. The citation is misapplied.**
- Lennon & Jones 2011 (Nat Rev Microbiol 9:119–130, doi:10.1038/nrmicro2504; **verified**: title,
  year, pages 119–130 all match) is a review of microbial seed banks and dormancy. I read the full
  abstract/key-points, the complete reference list (139 refs, via the article page and the Crossref
  reference metadata), the figure titles, and all three open-access supplementary files.
- The review **does contain a metagenomic TA analysis** — but it is **not about marine depletion.**
  Supplementary S2 (Box 2 methods) states they tblastx'd *hipA/B, relB/E, mazE/F, yafQ, dinJ* and RpfC
  against metagenomes, computing "percent of genomes" by recA recovery and gene-length correction
  (method ref: Mou et al. 2008 Nature 451:708, coastal-ocean carbon processing — verified). Box 2
  therefore shows **how widespread TA/Rpf dormancy genes are across environments**, i.e. that dormancy
  genes are commonly detected in metagenomes. There is **no statement that marine bacteria are
  TA-depleted**; "marine" appears in the review only in the context of VBNC, Pelagibacter proteomics,
  and starvation-survival references, never a marine-vs-terrestrial TA contrast.
- **What was almost certainly intended:** the marine/aquatic TA-abundance claim belongs to
  **Pandey & Gerdes 2005, Nucleic Acids Res 33:966–976, doi:10.1093/nar/gki201 (verified)** — "TA loci
  are highly abundant in free-living but lost from host-associated prokaryotes" — and/or
  **Fucich & Chen 2020 (verified)** for the open-ocean picocyanobacterial depletion. Note that Pandey &
  Gerdes actually argue free-living (incl. marine) genomes are TA-*rich*, so if that is the intended
  reference it points *against* a marine deficit, matching Makarova 2009 (threat 5). Kearns & Shade
  2018 cite Lennon & Jones 2011 only as a general dormancy reference and cite Pandey & Gerdes 2005
  separately for the TA abundance claim — consistent with this being a citation slip.
- **Recommendation:** replace the Lennon & Jones (2011) citation at that sentence. Either (a) drop the
  "agrees with a prior metagenomic marine observation" clause, since no clean prior marine genome-level
  deficit exists (the paper's own isolate result and Meaden 2025 are the real support), or (b) cite
  Meaden 2025 (binning-independent marine defense depletion) and Fucich & Chen 2020, and state openly
  that Makarova 2009 / Pandey & Gerdes 2005 point the other way for free-living/aquatic genomes.

---

## F. Papers that SUPPORT the project's interpretations (cite these)
- **MAG recovery biases** (support "MAG estimate is an upper bound; paralog/accessory compartment is
  least robust"): Mise & Iwasaki 2022; Meziti 2021; Maguire 2020 (doi:10.1099/mgen.0.000436, verified);
  Nebauer et al. 2026 Microb Genom doi:10.1099/mgen.0.001736 (verified; GC-driven taxonomy bias in MAG
  recovery).
- **Culturing bias** (support collider framing): Weissman 2021; Rodríguez del Río 2026.
- **Marine defensome depletion is real** (support, binning-independent): Meaden 2025; Beavogui 2024;
  Bograd 2025; Tesson 2022 (doi:10.1038/s41467-022-30269-9, verified, isolate baseline).
- **Phylogeny ≫ environment for other traits** (support the null): Henderson 2026; Favier 2026.
- **Sporulation is energetically costly / last-resort** (supports "regulation not copy number"):
  Karakoç, Shoemaker & Lennon 2026, PNAS 123:e2524274123, doi:10.1073/pnas.2524274123 (verified).

## G. What I could not access / verify
- **Chan/Melnikov 2025 hibernation-factor-evolution preprint full text** (bioRxiv 429; EuropePMC/RS
  full-text XML HTTP 500). Abstract only — see threat 7. **Most important gap.**
- **Marshall & Burgess 2015; Barneche 2018** full texts: paywalled (consistent with the project memo).
- **Donhauser 2026 Commun Earth Environ** has no Crossref abstract; read via publisher page — 70 sites
  on 5 gradients; spore-formation and rapid resuscitation track precipitation seasonality at the
  *community* level (supports the seasonality angle, not a genome-level competitor).
- **bioRxiv keyword API:** api.biorxiv.org offers no text search; bioRxiv coverage relied on EuropePMC
  SRC:PPR, which indexes it but may lag a few days for the newest items.
