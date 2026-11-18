# Phase 0 gap memo: The Dormancy Rule

**Date:** 2026-11-17

**Search effort:**
- 111 Europe PMC query runs: 68 distinct query strings (title+abstract) plus 38 full-text and 5 preprint/recent-year variants.
- Forward citation chasing on 14 seed papers.
- 13 web searches, plus full-text checks of 10 key papers.
- 14,565 unique records were screened programmatically. About 1,000 titles and about 220 abstracts were read by hand, and 12 full texts were checked.

Files: `search_log.tsv`, `literature_corpus.tsv`, `citation_verification.md`.

## Headline: no one has done the full design. One set of results points the opposite way, and the design needs changes.

I found **no paper or preprint (to Nov 2026)** that tests whether the dormancy gene content of *individual genomes* scales with a *continuous, organism-independent measure of temporal environmental variability at the isolation site* under *phylogenetic regression*. Every close work fails on at least two of these four axes:

1. Unit of analysis is the genome.
2. The variability measure is continuous.
3. Phylogeny is controlled.
4. More than one dormancy mechanism is analysed.

**Verdict: PARTLY OPEN.** The question "do dormancy genes track environment?" has been asked many times, but only with categorical habitats or community-weighted means. The genome-level, continuous-variability, phylogenetically controlled, multi-mechanism test is **open**.

## The closest existing works

| Work | Unit / n | Environment measure | Phylo control | Mechanisms | What it did NOT do |
|---|---|---|---|---|---|
| **Makarova, Wolf & Koonin 2009** *Biol Direct* 4:19 | genome, 750 complete genomes | **categorical**: aquatic / terrestrial / host / specialized / multi-environment, plus temperature class | none (Archaea vs Bacteria factor only); genome-size allometry residuals | TA only | No continuous or variability index and no PGLS. **Found terrestrial and multi-environment genomes carry *fewer* TA than genome size predicts**, the opposite of a naive seed-bank expectation. |
| **Pandey & Gerdes 2005** *NAR* 33:966 | genome, 126 | categorical lifestyle (obligate intracellular vs free-living) | none; no genome-size correction | TA only | Makarova later suggested the absence pattern is largely genome-size allometry. |
| **Fucich & Chen 2020** *ISME J* 14:2843 | genome, 71 picocyanobacteria | categorical: open ocean vs coastal vs freshwater | none | type II TA | **This is a mini-version of our hypothesis in one clade**: coastal and freshwater "variable" habitats have more TA, but TA scales linearly with genome size and there are no statistics on variability. |
| **Browne et al. 2021** *Genome Biol* 22:204 | genome, 1,358 Firmicutes | categorical: host vs environment | phylogenetic mapping of losses, no regression | sporulation (66-gene signature) | No environmental variability; one phylum only. |
| **Donhauser et al. 2026** *Commun Earth Environ* 7:15 | **community** (70 soil samples, ASV-inferred trait CWMs) | **continuous BIO15 precipitation seasonality** | trait inference uses phylogeny; no PGLS | sporulation phenotype, plus RRN, genome size | **Closest on the predictor.** BIO15 is positively related to community sporulation. It is not genome-level, does not split predictable from unpredictable variability, and covers one mechanism. |
| **Piton et al. 2023** *Nat Microbiol* 8:2093 | community metagenomes, global soils | precipitation patterns, pH, C:N (random forest) | none | dormancy/sporulation genes inside trait axes | Community level. It reports dormancy genes co-vary with *large* genomes while seasonality and drought favour *small* genomes, a confound (see below). |
| **Kearns & Shade 2018** *ISME J* 12:2575 | community (succession) plus a genome-level side analysis | succession stage | none | sporulation, TA, Rpf | Genome level: genomes carrying dormancy genes have **more rRNA operons** (H = 1326.6; ρ = 0.78, TA 0.44). No environment test at genome level. |
| **Lennon & Jones 2011** *Nat Rev Microbiol* 9:119, Box 2 | metagenomes (ocean, sludge, soil, gut) | categorical biome | none | Spo0A/B, TA, Rpf | Descriptive percentages only. |
| **Parter, Kashtan & Alon 2007** *BMC Evol Biol* 7:169 | genome, 117 species | **categorical "habitat variability" classes** from NCBI lifestyle | none | metabolic network modularity | This is the methodological ancestor ("variable habitat → genome architecture"), not about dormancy. |
| **Danov et al. 2024** *mBio* 15:e01911-23 (Toxinome) | genome, 59,475 | categorical temperature class | none | all protein toxins incl. TA | TA/toxins are depleted in thermophiles and psychrophiles, which contradicts Makarova. No variability measure. |
| **Hao et al. 2026** *PNAS* 123:e2520645123 | species, 25,785 GTDB r207 | niche breadth (26 habitat clusters) | **yes (MCMCglmm, PGLS, λ)** | cooperation genes, not dormancy | This is the **design template** for our secondary predictor. It covers neither dormancy nor temporal variability. |
| **Favier et al. 2026** *Ecol Lett* | isolate, 251 Gammaproteobacteria | **isolation-site bioclimatic variables** | phylogeny-aware | thermal performance, not dormancy | Precedent for isolation-site climate. It found **minimal influence of present-day bioclimate** on traits, a power warning. |

Gene-family surveys that define our modules do not test environment at all:
- Galperin 2022 (sporulation core)
- Helena-Bueno 2024 (Balon in 1,573 of 8,761 species)
- Chan et al. 2025 bioRxiv (hibernation-factor loss in >10% of 46,015 genomes)
- Atkinson 2011 (RSH/SAS)
- Schroeckh 2006 (Rpf)

Community-metagenome climate studies relate dormancy or sporulation genes to **mean aridity**, not variability: Fierer 2012 PNAS, Tripathi 2017, Song 2019, Li 2022, Broderick 2025.

## The specific contribution (what is genuinely new)

1. **Organism-level test.** Genome dormancy modules are regressed on the physical variability at each genome's isolation site, not community-weighted means.
2. **Variability decomposed into predictable and unpredictable parts.** We separate seasonal amplitude (predictable, which favours *responsive* or anticipatory strategies; Bernhardt et al. 2020) from de-seasonalised residual variance and lag-1 autocorrelation (unpredictable, which bet-hedging theory says should favour seed banks). No existing microbial study makes this distinction. **This is the sharpest novelty claim, so lead with it.**
3. **Multi-mechanism modules in one phylogenetic framework.** These are sporulation, Rpf, hibernation, and stringent response. TA is analysed as a separate, confounded module that serves as a negative or positive control.
4. **Phylogenetic regression** (phylolm/phyloglm, Pagel's λ) for dormancy versus environment. This has only been done for other traits, such as genome size, cooperation and thermal limits.

## Threats that change the design

1. **Genome size is a mediator, not just a covariate.**
   - Seasonality and aridity select for *smaller* genomes (Piton 2023, Donhauser 2026, Liu 2023, Wang 2025).
   - TA (and probably other dormancy) gene counts scale super-linearly with genome size (exponent 1.64; Makarova).
   - Raw counts could therefore show a *negative* association with variability, and adjusting for genome size could absorb a real effect.
   - **Change:** add phylogenetic path analysis (`phylopath`), or report both the total effect (no genome-size term) and the direct effect (with it), and pre-register which one tests the hypothesis.
2. **rRNA operon count is also on the causal path.** Dormancy genes go with high RRN (Kearns & Shade), and RRN rises under high seasonality (Donhauser). Treat RRN the same way as genome size; do not use it as a plain nuisance covariate.
3. **TA is unreliable as a dormancy proxy.**
   - Published directions conflict: Makarova finds terrestrial and multi-environment genomes depleted and thermophiles enriched, while Danov finds thermophiles depleted.
   - TA systems act as phage defence and addiction modules (LeRoux & Laub 2022; Belot 2026), and TA deletions do not affect persistence (Goormaghtigh 2018).
   - Keep TA strictly separate. Pre-register it as the module expected to track mobile-element and defence load, not variability.
4. **Power and phylogenetic signal.**
   - Sporulation is nearly binary and restricted to Bacillota, Rpf is concentrated in Actinomycetota, and niche breadth has λ ≈ 0.88.
   - Isolation-site climate explained little in Favier 2026.
   - **Change:** fit mechanism-specific models within the clades where the trait actually varies (phyloglm within Bacillota for sporulation, within Actinomycetota for Rpf), and run a simulation-based power analysis before Phase 1.
5. **Host-associated and clinical genomes dominate GTDB.** Isolation-site SST or air temperature means nothing for a gut isolate. **Change:** restrict to free-living environmental isolates, or include host association as a stratum.
6. **The right physical variable may differ by habitat.** Precipitation seasonality (water), not temperature, drove sporulation in soils (Donhauser). **Change:** use a habitat-matched index: SST for marine and freshwater, and air temperature *plus* precipitation or soil moisture (e.g. ERA5-Land) for terrestrial sites.
7. **Scoop risk is moderate.**
   - The Lennon lab explicitly calls for genome and metagenome *rpf* surveys with environmental niche models (Lennon et al. 2025 mSystems), and is active on sporulation bioenergetics (Karakoç 2026 PNAS).
   - The Hao/West/Iwasaki pipeline (GTDB + PGLS + niche breadth) could add dormancy genes with little effort.
   - The Donhauser group already links BIO15 to sporulation.
   - Preprint early.
8. **Secondary predictor data exist but need mapping.**
   - Hao niche breadth is per species on **GTDB r207** (25,785 species), not r232.
   - von Meijenfeldt SNB is on **NCBI taxonomy** (16,632 species, 3,182 genera; mostly amplicon-based), so genus-level SNB is safer.
   - Both are downloaded to `/home/user/data/lit/`.

## Bottom line
Nobody has published the paper. The surviving, defensible contribution is: **the first genome-level, phylogenetically controlled test of the seed-bank prediction, using organism-independent physical variability that is decomposed into predictable and unpredictable components, across separate dormancy mechanisms, with genome size and rRNA operon count handled explicitly as mediators.**

Cite Makarova 2009, Fucich & Chen 2020, Browne 2021 and Donhauser 2026 prominently as partial precedents. Makarova's "terrestrial and multi-environment depleted" result should be stated up front as the existing (crude) evidence that runs against the hypothesis for TA.
