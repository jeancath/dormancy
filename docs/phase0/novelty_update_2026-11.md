# Novelty and bias update: The Dormancy Rule (Nov 2026)

**Date:** 2026-11-17
**Scope:** This is a targeted follow-up to `gap_memo.md`. It is driven by four preliminary findings (a–d) and covers six questions: a 2025–2026 scoop check, MAG gene-content biases, culturing bias, site non-independence, decomposing predictability, and marine TA depletion.

**Search effort:**
- About 45 Europe PMC REST queries, including preprints (SRC:PPR) and author-targeted runs on the Lennon, Shade, Fierer, Barberán, Westoby, Donhauser and Hao/West/Iwasaki groups.
- 4 web searches.
- Full-text checks of Hao et al. 2026 (PMC12974524), Beavogui et al. 2024 (PMC10924106), Meaden et al. 2025 (PMC12409356) and the Favier et al. 2026 bioRxiv preprint.
- The envPred package documentation, for the Marshall & Burgess method.

**Verification:** Every DOI below was resolved against Crossref on 2026-11-17, and the bibliographic details given here are the Crossref ones. Where the task brief named a paper with slightly different details, the correction is noted in place.

---

## 1. Scoop check, 2025–2026

**Verdict: nothing closer than `gap_memo.md` already lists.** I found no 2025–2026 paper or preprint that regresses genome-level sporulation, Rpf, hibernation, stringent-response or TA content on environmental variability, seasonality or climate under phylogenetic control. The nearest new items fall into three groups:

- community-level or experimental studies;
- reviews;
- comparative-trait papers on other traits.

Several should be cited, and two groups should be watched.

**Relevant new items (none is a scoop):**

- **Guigard L et al. 2026 Microbiol Mol Biol Rev 90:e00354-25 doi:10.1128/mmbr.00354-25.** A Shade-lab review that recasts dormancy as an activity spectrum shaped by "fluctuating conditions". It is conceptual and community-level, with no genome-by-environment test. It is a good framing citation.
- **Weissman JL et al. 2026 ISME J 20:wrag147 doi:10.1093/ismejo/wrag147.** Genomic growth-optimisation signatures (codon usage bias) are confounded by temperature, so genomic "copiotrophy" decouples from maximum growth rate along temperature gradients. This matters because our covariates (rRNA operon count, genome size) and our predictor (temperature) share this confound.
- **Hao C et al. 2026 bioRxiv doi:10.64898/2026.10.22.761225** ("Growth rate structures the dynamics of bacterial trait evolution"). Multivariate Ornstein–Uhlenbeck models over 283 species with measured growth rates and 11 traits, four of them environmental. It is not about dormancy, but it is the same group as Hao et al. 2026 PNAS and shows the group is moving into environment–trait coevolution. Watch it.
- **Santillan E et al. 2026 ISME J 20:wrag099 doi:10.1093/ismejo/wrag099.** CSR life-history classification of MAGs across a gradient of disturbance frequency, in bioreactors. It is experimental and community-level.
- **González-Sánchez A et al. 2026 Microb Ecol 89:114 doi:10.1007/s00248-026-02744-z.** In a desiccating desert wetland, the spore fraction (heat-resistant CFU) rose over 2007–2025. This is culture-based and community-level.
- **Donhauser J et al. 2025 Research Square preprint doi:10.21203/rs.3.rs-5145577/v1.** The preprint of the Donhauser 2026 *Commun Earth Environ* paper already in the gap memo. Nothing new.
- **Peng Z et al. 2025 Adv Sci 12:2405947 doi:10.1002/advs.202405947.** Life-history traits vs niche breadth. It is already in the corpus, has no temporal-variability predictor, and does not touch dormancy.

**Synthesis.** The genome-level gap stated in the gap memo still holds as of 17 Nov 2026. The groups with both the pipeline and the motive are unchanged: Hao/West/Iwasaki (GTDB + PGLS), Lennon (Rpf/sporulation) and Shade (dormancy spectrum). The Hao group's new OU-model preprint shows active work on environment–trait coevolution, so the scoop risk stays **moderate**, and a preprint should go up early.

---

## 2. Gene-content biases of MAGs vs isolate genomes

**Verified references.** Corrections to the task brief are marked.

- **Meziti A et al. 2021 Appl Environ Microbiol 87:e02593-20 doi:10.1128/AEM.02593-20.** *Correction:* this is 2021, not 2019. MAGs and isolates were taken from the same eight faecal samples (pathogenic *E. coli*). MAGs whose estimated completeness was about 95% recovered only **77% of population core genes and 50% of population-variable genes**. About 5% of MAG genes were foreign contamination, against a CheckM-estimated 1.5%.
- **Nelson WC, Tully BJ, Mobberley JM 2020 PeerJ 8:e10119 doi:10.7717/peerj.10119.** *Correction:* the title is "Biases in genome reconstruction from metagenomic data". No Nelson 2020 paper titled "Disentangling…" or "Mobile genetic elements in MAGs" was found. MAGs from an enrichment were compared with complete genomes of 10 isolates from the same enrichment, and the comparison was checked against Tara Oceans MAGs vs RefSeq. Repeats and regions of atypical composition went unbinned, and the unbinned genes were **strongly biased towards rRNA, tRNA, mobile-element functions and unknown genes**. Uneven plasmid distribution can mislead inference.
- **Maguire F et al. 2020 Microb Genom 6(10):e000436 doi:10.1099/mgen.0.000436.** First author is Finlay Maguire. A simulated metagenome of 30 genomes rich in genomic islands (GIs) and plasmids was run through 12 binning pipelines. **82–94% of chromosomal sequence was binned correctly, but only 38–44% of GIs and 1–29% of plasmid sequence.** No plasmid-borne AMR or virulence genes were recovered at all.
- **Chen LX et al. 2020 Genome Res 30:315–333 doi:10.1101/gr.258640.119.** Draft MAGs have gaps, local misassemblies, chimeras and contamination. Curation to complete (circular) MAGs is rare. GC-skew analysis of about 7,000 complete isolate genomes flags misassemblies at repeats even in references.
- **Mise K, Iwasaki W 2022 ISME Commun 2:118 doi:10.1038/s43705-022-00204-6.** Across more than 190,000 MAGs, ribosomal-protein genes are missing in **20–40% of near-complete MAGs**. The cause is codon-usage bias interacting with tetranucleotide binning. **The loss is worse in fast-growing lineages, so recovery bias is correlated with life-history strategy.**
- **Sczyrba A et al. 2017 Nat Methods 14:1063–1071 doi:10.1038/nmeth.4458** (CAMI). Assembly and binning perform well for single-strain species but are **"substantially affected by the presence of related strains"**.
- **Parks DH et al. 2017 Nat Microbiol 2:1533–1542 doi:10.1038/s41564-017-0012-7** and **Bowers RM et al. 2017 Nat Biotechnol 35:725–731 doi:10.1038/nbt.3893** (MIMAG). Multi-copy rRNA operons are often absent from MAGs. MIMAG therefore makes rRNA and tRNA presence a separate quality criterion from CheckM completeness.
- **Shaiber A, Eren AM 2019 mBio 10:e00725-19 doi:10.1128/mbio.00725-19.** Composite (multi-population) MAGs pass completeness and redundancy thresholds yet contaminate public repositories.
- **Evans JT, Denef VJ 2020 mSphere 5:e00971-19 doi:10.1128/msphere.00971-19.** Dereplication choices, and the strain heterogeneity they paper over, change what accessory content a MAG represents.
- **Beavogui A et al. 2024 Nat Commun 15:2146 doi:10.1038/s41467-024-46489-0.** Defence systems were surveyed in 7,759 high-quality MAGs (395 soil, 386 marine, 6,978 gut). The authors state that short-read binning misses MGE families, and that because defence genes sit on or near MGEs, their counts are "most likely a partial underestimated picture".
- **Gupta S, Almeida A 2025 Nat Commun 16:9959 doi:10.1038/s41467-025-64950-6.** In *K. pneumoniae*, MAGs added lineages missing from isolate collections and 214 MAG-only genes. Isolates and MAGs therefore sample different parts of a species' diversity, not just different quality tiers.

**Synthesis.** The evidence agrees that short-read MAGs systematically under-recover several kinds of DNA, even when CheckM completeness is above 90%:

- plasmids (the most severe loss, often less than 30% recovered);
- genomic islands and other horizontally acquired regions (about 40% recovered);
- multi-copy and repetitive loci (rRNA operons);
- population-variable accessory genes (about 50% recovered).

Defence systems and TA loci are concentrated in exactly these compartments, so their counts in MAGs are biased downwards. Two features make the bias a threat to a comparative analysis rather than just noise:

- **It is heterogeneous.** It depends on strain microdiversity, coverage and community complexity, all of which differ between habitats.
- **It can track the trait of interest.** In Mise & Iwasaki, faster growers lose more core genes. Rpf and sporulation co-vary with growth strategy.

Single-copy core genes, such as *relA/spoT* or the *hpf*/*rmf* hibernation factors, should be the most robust to the bias. Accessory and mobile modules (TA, some Rpf paralogs) should be the least robust.

---

## 3. Culturing bias: isolates vs the uncultured majority

**Verified references:**

- **Staley JT, Konopka A 1985 Annu Rev Microbiol 39:321–346 doi:10.1146/annurev.mi.39.100185.001541.** The origin of the "great plate count anomaly": plate counts are orders of magnitude below direct counts.
- **Lloyd KG et al. 2018 mSystems 3(5):e00055-18 doi:10.1128/msystems.00055-18.** Using primer-free metagenomic 16S, uncultured genera make up an estimated **81% of cells** and uncultured phyla 25%. Human-associated environments are the exception, dominated by cultured genera (45–97%).
- **Martiny AC 2019 ISME J 13:2125–2128 doi:10.1038/s41396-019-0410-3** argued that high proportions are culturable. **Steen AD et al. 2019 ISME J 13:3126–3130 doi:10.1038/s41396-019-0484-y** rebutted it: gap handling and amplicon primer bias inflate apparent culturability, so most taxa remain uncultured.
- **Weissman JL, Hou S, Fuhrman JA 2021 PNAS 118:e2016810118 doi:10.1073/pnas.2016810118.** Codon-based maximal growth rates over more than 200,000 genomes, MAGs and SAGs: **culture collections are "strongly biased toward organisms capable of rapid growth"**. Organisms fall into two growth classes.
- **Rodríguez del Río Á et al. 2026 BMC Genomics 27:655 doi:10.1186/s12864-026-13092-0.** 40,039 soil bins from the SMAG catalogue. Compared with genera that have cultured members, uncultivated soil prokaryotes have **smaller genomes, lower GC, slower growth** (the best predictor of cultivation status) and fewer biosynthesis and motility genes. Completely uncultivated genera are more abundant in tropical and arctic soils, so culturability is itself biogeographically structured.
- **Roller BR, Stoddard SF, Schmidt TM 2016 Nat Microbiol 1:16160 doi:10.1038/nmicrobiol.2016.160** and **Klappenbach JA et al. 2000 Appl Environ Microbiol 66:1328–1333 doi:10.1128/aem.66.4.1328-1333.2000.** rRNA operon copy number predicts maximum growth rate and copiotrophy. It is high in organisms that grow quickly on plates.
- **Browne HP et al. 2016 Nature 533:543–546 doi:10.1038/nature17645.** Targeted culturing of gut bacteria, including ethanol treatment that selects for spores, showed that **50–60% of gut genera form spores**. Isolation protocols can enrich spore-formers, and sequencing alone under-detects them.
- **Lewis WH et al. 2021 Nat Rev Microbiol 19:225–240 doi:10.1038/s41579-020-00458-8.** A review of why the majority stays uncultured: slow growth, dependencies and dormancy.

**Synthesis.** Isolates are not a random sample of species. Culturing preferentially captures fast-growing, high-rRNA-copy, larger-genome copiotrophs, and in some protocols spore-formers. The uncultured majority, represented only by MAGs, is slower-growing, more streamlined and more biosynthetically dependent. Each of these traits is causally linked to dormancy-gene content (gap memo, threats 1–2). Culturability is also geographically structured (Rodríguez del Río 2026). An isolate-only slope is therefore estimated within a selected subpopulation, and a MAG-only slope within a differently selected subpopulation that is also differently measured.

---

## 4. Spatial pseudoreplication and site non-independence

**How the reference studies handled it:**

- **Hao C et al. 2026 PNAS 123:e2520645123 doi:10.1073/pnas.2520645123.** Checked in the full text. The predictor is species-level niche breadth from ProkAtlas presence across 114 habitats, grouped into 26 habitat clusters. Non-independence is handled **at the habitat level**: compositionally similar habitats are clustered so they are not counted twice. Phylogeny enters through MCMCglmm and PGLS, and robustness through repeated 50% species subsampling. There is no site-level random effect, because genomes are not tied to individual sampling sites. **Their design therefore does not face our problem.**
- **Favier AI et al. 2026 Ecol Lett 29:e70483 doi:10.1111/ele.70483.** The bioRxiv preprint (doi:10.64898/2026.02.06.703311) was read through a summarising fetch tool. It reports 251 isolates from **15 sites** in California, with site-level bioclimatic PCs as fixed effects and **only taxonomic random effects** (`(1|Order/Genus)`). No site random effect or spatial term was reported. Their "minimal influence of bioclimate" therefore comes from about 15 effective climate observations. Check the published version before citing this detail.
- **Sriswasdi S, Yang CC, Iwasaki W 2017 Nat Commun 8:1162 doi:10.1038/s41467-017-01265-1** and **von Meijenfeldt FAB et al. 2023 Nat Ecol Evol 7:768–781 doi:10.1038/s41559-023-02027-7.** Habitat and niche breadth are summarised per taxon across many samples, which avoids per-genome site pseudoreplication by design.
- **Barberán A et al. 2012 Mol Ecol 21:1909–1917 doi:10.1111/j.1365-294x.2011.05383.x.** Community-trait (metagenome-level) analysis, where the sample is the unit. It does not address multiple taxa per site.
- **Barneche DR, Burgess SC, Marshall DJ 2018 Glob Ecol Biogeogr 27:890–898 doi:10.1111/geb.12748.** The closest structural analogue outside microbes: 1,078 observations of 288 marine fish species at 242 localities. Decomposed SST predictability enters a **Bayesian phylogenetic hierarchical model**. I did not check in the full text whether a locality effect was included.

**Methods references:**

- **Hurlbert SH 1984 Ecol Monogr 54:187–211 doi:10.2307/1942661.** Defines pseudoreplication.
- **Freckleton RP, Jetz W 2009 Proc R Soc B 276:21–30 doi:10.1098/rspb.2008.0905.** A comparative model with simultaneous phylogenetic and spatial covariance terms, estimating their relative contributions.
- **Hadfield JD, Nakagawa S 2010 J Evol Biol 23:494–508 doi:10.1111/j.1420-9101.2009.01915.x** and **Hadfield JD 2010 J Stat Softw 33(2) doi:10.18637/jss.v033.i02.** MCMCglmm with a phylogenetic random effect alongside further random effects, for example `random = ~ phylo + site`.
- **Li D et al. 2020 Methods Ecol Evol 11:1455–1463 doi:10.1111/2041-210X.13471** (`phyr`). Phylogenetic GLMMs with crossed species and site random effects, including phylogeny-by-site terms.
- **Dormann CF et al. 2007 Ecography 30:609–628 doi:10.1111/j.2007.0906-7590.05171.x.** A review of spatial-autocorrelation corrections: GLS with spatial covariance, spatial eigenvectors and autocovariates.
- **Roberts DR et al. 2017 Ecography 40:913–929 doi:10.1111/ecog.02881** and **Valavi R et al. 2019 Methods Ecol Evol 10:225–232 doi:10.1111/2041-210X.13107** (`blockCV`). Block cross-validation for spatial, hierarchical and phylogenetic structure.
- **Dewar AE, Belcher LJ, West SA 2025 Nat Rev Genet 26:395–405 doi:10.1038/s41576-024-00803-0.** A review of phylogenetic non-independence in comparative genomics. It is the West group's methods statement.

**Synthesis.** I found **no microbial comparative-genomics paper that models site as a random effect** alongside phylogeny for an isolation-site environmental predictor. The two templates sidestep the issue: Hao uses species-level predictors, and Favier uses taxonomic random effects only. Our design has the predictor measured per site, with up to dozens of species per site. The closest precedents are therefore Barneche et al. 2018 (phylogenetic hierarchical model over localities) and the general methods.

The recommended fit is a phylogenetic mixed model with **two random effects**, the phylogenetic covariance and an i.i.d. site intercept. Options:

- MCMCglmm with `~ phylo + site`;
- `phyr::pglmm` with `(1|sp__) + (1|site)`;
- `brms` with a phylogenetic covariance matrix plus `(1|site)`;
- a spatial covariance term in place of the site intercept, if sites cluster geographically (Freckleton & Jetz).

The model should be backed by **site-blocked or spatially blocked cross-validation and jackknife** (Roberts 2017; `blockCV`). The effective sample size for a site-level predictor is the number of sites, not the number of genomes. phylolm cannot fit the extra random effect, so the current phylolm fits should be treated as sensitivity analyses.

---

## 5. Collinearity of seasonal and residual variability, and decomposing predictability

**Verified references:**

- **Colwell RK 1974 Ecology 55:1148–1153 doi:10.2307/1940366.** Predictability (P) = constancy (C) + contingency (M). The ratio M/P is the fraction of predictability that comes from periodicity. These are information-theoretic indices on binned states, so they are scale-free.
- **Vasseur DA, Yodzis P 2004 Ecology 85:1146–1152 doi:10.1890/02-3122.** Environmental noise colour is the spectral exponent of a time series. This is the reference for the colour metric only; I did not check the paper's habitat comparisons.
- **Marshall DJ, Burgess SC 2015 Ecol Lett 18:174–181 doi:10.1111/ele.12402.** Global analysis of more than 800 marine invertebrates. Predictability is decomposed into **seasonality** and **colour**, and they relate differently to life-history traits.
  - The operational recipe, documented in the authors' `envPred` R package:
    1. Detrend the series.
    2. Compute a monthly-mean seasonal cycle.
    3. Seasonality = α/(α+β): the variance of the seasonal cycle over the total, where α is seasonal variance and β is residual variance.
    4. Colour θ = the negative slope of log spectral density against log frequency in the de-seasonalised residuals.
  - I verified the recipe from the package documentation, not the paper's full text, which is paywalled.
- **Barneche DR et al. 2018 Glob Ecol Biogeogr 27:890–898 doi:10.1111/geb.12748.** The same decomposition applied to SST and chlorophyll-a: mean, seasonality and colour, used as joint predictors in a phylogenetic model.
- **Tonkin JD et al. 2017 Ecology 98:1201–1216 doi:10.1002/ecy.1761.** Separate wavelet-based seasonality and Colwell-based predictability. Temporal β-diversity peaks where both are high, so the two axes are treated as distinct.
- **Vázquez DP et al. 2017 Biol Rev 92:22–42 doi:10.1111/brv.12216.** A review of climatic variability and extremes across levels of organisation. It is conceptual and gives no specific decomposition.
- **Bernhardt JR et al. 2020 Phil Trans R Soc B 375:20190454 doi:10.1098/rstb.2019.0454.** Life in fluctuating environments; the framing for predictable vs unpredictable fluctuation. Already verified in the gap memo.
- **Botero CA et al. 2015 PNAS 112:184–189 doi:10.1073/pnas.1408589111.** Theory: bet-hedging is favoured in a specific region of the space defined by **predictability × timescale**. This supports a predictability *fraction*, not an absolute residual variance, as the bet-hedging predictor.
- **Simons AM 2011 Proc R Soc B 278:1601–1609 doi:10.1098/rspb.2011.0176.** Empirical evidence for bet hedging is elusive, partly because the modes of response are hard to separate.
- **Gremer JR, Venable DL 2014 Ecol Lett 17:380–387 doi:10.1111/ele.12241.** Desert annuals: observed germination fractions match bet-hedging optima under rainfall variance. This is the clearest empirical seed-bank analogue.
- **Rubio de Casas R et al. 2017 New Phytol 214:1527–1536 doi:10.1111/nph.14498.** Phylogenetically informed global legume analysis: seed dormancy is determined by **seasonality** (predictable fluctuation) and seed size. This is a direct plant analogue in which dormancy tracks *predictable* rather than unpredictable variation.
- **Fusi M et al. 2024 Sci Rep 14:26174 doi:10.1038/s41598-024-77621-1.** Thermal tolerance tracks the predictability of diel fluctuations. Tropical sites had more predictable fluctuations than temperate ones, so predictability and latitude are confounded.
- **Dormann CF et al. 2013 Ecography 36:27–46 doi:10.1111/j.1600-0587.2012.07348.x.** Collinearity problems start around |r| > 0.7. Re-parameterisation is preferred over dropping a variable, and sequential residualisation is not recommended.
- **Freckleton RP 2002 J Anim Ecol 71:542–545 doi:10.1046/j.1365-2656.2002.00618.x.** Regressing on residuals of one predictor against another biases the estimates; use multiple regression instead.

**Synthesis.** A correlation of 0.82 between seasonal amplitude and residual variance is expected. Both scale with the overall magnitude of variability, which itself scales with latitude and continentality. The ecology literature gives a concrete way out that the gap-memo design did not use. Express predictability as **scale-free proportions**, not absolute variances:

- **seasonal fraction** α/(α+β), or Colwell's M/P;
- **residual colour** θ (the spectral exponent, or lag-1 autocorrelation of the residuals);
- a separate **total-variability magnitude** term (log total SD).

Marshall & Burgess 2015 and Barneche et al. 2018 use exactly this mean + seasonality + colour triplet in phylogenetic models. Do not residualise one variance on the other (Freckleton 2002; Dormann 2013).

The plant analogue (Rubio de Casas 2017) and the theory (Botero 2015) give the hypothesis two sides. Dormancy may track *predictable* seasonality (anticipatory dormancy) rather than unpredictable noise (bet-hedging). The pre-registered expectation should be stated in both directions.

---

## 6. Marine TA depletion: genome-level evidence 2015–2026

**Verified references:**

- **Makarova KS, Wolf YI, Koonin EV 2009 Biol Direct 4:19 doi:10.1186/1745-6150-4-19.** Baseline, with 750 complete genomes. As a residual after genome size, **"aquatic" genomes (freshwater and marine pooled) carry *more* TA**, and terrestrial genomes fewer. This conflicts with our finding (c) and should be addressed explicitly.
- **Fucich D, Chen F 2020 ISME J 14:2843–2850 doi:10.1038/s41396-020-00746-4.** Picocyanobacteria: open-ocean strains have few or no TA (0/38 *Prochlorococcus*), coastal and freshwater *Synechococcus* many. TA count scales with genome size.
- **Danov A et al. 2024 mBio 15:e01911-23 doi:10.1128/mbio.01911-23** (Toxinome). 59,475 genomes; reports temperature classes only. **I found no marine vs terrestrial comparison in it** (consistent with the gap memo).
- **Guan J et al. 2024 Nucleic Acids Res 52:D784–D790 doi:10.1093/nar/gkad962** (TADB 3.0). 211,697 predicted type I–VIII TA loci in 34,789 *complete* genomes, of which **69,019 (about a third) are linked to mobile genetic elements**. It gives no habitat comparison, but it quantifies the MGE share that MAG binning loses.
- **Akarsu H et al. 2019 PLoS Comput Biol 15:e1006946 doi:10.1371/journal.pcbi.1006946** (TASmania). More than 2 × 10^6 TA candidate loci over more than 41,000 EnsemblBacteria assemblies. It is a resource with no habitat test.
- **Liu C et al. 2025 Mar Drugs 23:436 doi:10.3390/md23110436.** TA diversity in marine genomes from the Global Ocean Microbiome Catalogue and the Mariana Trench (MEER), which are largely MAG-based. It is descriptive: TA types I–VII are prevalent, with novel combinations. **It has no terrestrial comparator and no genome-size or MGE normalisation.**
- **Beavogui A et al. 2024 Nat Commun 15:2146 doi:10.1038/s41467-024-46489-0.** In MAGs, the defence systems (which include abortive-infection and TA-like systems) of **marine genomes follow a geometric-like distribution, with about 65% of marine MAGs having a "limited defensome"**. Soil and gut MAGs follow a binomial-like distribution with a mode of 3–4 systems. The authors flag MAG under-recovery of MGE-borne defences as a likely bias.
- **Meaden S, Westra ER, Fineran PC 2025 Phil Trans R Soc B 380:20240069 doi:10.1098/rstb.2024.0069.** A contig-level survey of whole metagenomes, which **does not depend on binning**: "marine environments appear to be particularly depleted in defence systems". Defence abundance rises with viral density, and is highest in host-associated and hot environments.
- **Tesson F et al. 2022 Nat Commun 13:2561 doi:10.1038/s41467-022-30269-9.** In 21,000 complete genomes, defence-system content varies with genome size, viral threat and lifestyle.
- **Sheinman M et al. 2021 eLife 10:e62719 doi:10.7554/eLife.62719.** Identical sequences in distant genomes indicate HGT; marine bacteria show far lower effective HGT rates than gut bacteria. Beavogui cites this as roughly 10^3-fold.
- **Makarova KS et al. 2011 J Bacteriol 193:6039–6056 doi:10.1128/jb.05535-11** and **Jurėnas D et al. 2022 Nat Rev Microbiol 20:335–350 doi:10.1038/s41579-021-00661-1.** TA loci cluster in defence islands and on MGEs. Their evolution is dominated by HGT and addiction and defence roles.

**Synthesis.** I found **no 2015–2026 study that tests TA abundance in marine vs terrestrial bacteria at genome level with genome-size and MGE correction**, in either MAGs or isolates. What exists falls into four kinds:

1. One clade-level isolate study, showing open-ocean depletion (Fucich & Chen).
2. A descriptive marine MAG catalogue survey (Liu 2025).
3. Two defence-system surveys showing marine depletion: Beavogui (MAGs) and Meaden (contigs, binning-independent).
4. The older Makarova result, pointing the other way for "aquatic" as a whole.

The binning-independent Meaden result and the low marine HGT rates make a real marine deficit of mobile defence/TA loci plausible. The size of any MAG-based estimate, however, is confounded by MAG under-recovery of MGEs, so (c) is a partially novel contribution.

---

## Implications for interpreting findings a–d

**(a) Pooled near-null (inside the SESOI).**
- The pooled estimate averages two populations that are **sampled differently** (§3) and **measured differently** (§2). The equivalence test is valid as a statement about the GTDB-representative population as assembled. It should not be read as "no organism-level relationship" until the source-specific slopes are explained.
- The null agrees with Favier et al. 2026 (minimal present-day bioclimate signal). But Favier had about 15 sites and no site random effect, so it is weak corroboration.
- Report the pooled equivalence result **with** the source-stratified estimates. Avoid a headline of "practically null".

**(b) MAG + vs isolate − slopes.** Four non-exclusive explanations follow from the literature. Each can be tested.

1. **Site concentration and pseudoreplication.** About 370 MAG sites, with 41% of MAGs in the top 10, means the MAG slope rests on few effective climate observations. The site jackknife and one-species-per-site checks are good, but the main model should be a phylogenetic mixed model with a site random effect (§4).
2. **Differential recovery correlated with environment.** MAG gene recovery depends on strain microdiversity, coverage and community complexity (CAMI; Meziti; Maguire). If more variable sites yield less complex communities, or higher-coverage or more clonal populations, MAGs there will look gene-richer. That artefact gives a positive slope.
   - Add MAG completeness, contamination, N50 and contig count, and their interaction with variability, as covariates.
   - Restrict to high-contiguity MAGs (Beavogui used N50 ≥ 100 kb).
   - Contrast single-copy core modules (stringent response, hibernation), which should be robust, with accessory and mobile ones.
   - Mise & Iwasaki 2022 show recovery bias co-varies with growth strategy, which is causally linked to dormancy.
3. **Culturing selection (collider bias).** Isolation favours fast growers, high rRNA operon counts, large genomes and, under some protocols, spore-formers (Weissman 2021; Rodríguez del Río 2026; Browne 2016). If both dormancy capacity and site climate affect whether an organism gets isolated, conditioning on isolation can induce a spurious negative association. This is a hypothesis. A partial test is to model isolation status as a function of variability and dormancy, or to compare isolates with MAGs **within the same genera**.
4. **Different clades.** MAG and isolate species occupy different parts of the tree (Lloyd 2018; Gupta & Almeida 2025). Test the source × variability interaction inside one phylogenetic model, and within families that contain both sources.

**(c) Marine TA deficit strong in MAGs, weak in isolates.**
- The MAG-specific strength is the pattern expected if binning loses plasmid, GI and MGE-borne TA (Maguire: 1–29% plasmid and about 40% GI recovery; TADB 3.0: about a third of TA loci sit with MGEs). The loss should be worst in marine metagenomes, which are microdiverse and strain-rich (CAMI; Meziti: 50% of variable genes recovered). The MGE correction cannot fully remove this, because the MGEs themselves are under-recovered.
- Treat the **complete-isolate-genome estimate as primary** and the MAG estimate as an upper bound on the deficit.
- Binning-independent support for a real marine deficit exists (Meaden 2025, at contig level for defence systems; Sheinman 2021 on HGT), as does clade-level support (Fucich & Chen). The finding must still be reconciled with Makarova 2009 ("aquatic" enriched after genome-size correction). Splitting freshwater from marine is one likely resolution.

**(d) Seasonal vs residual variability, r = 0.82.**
- Expected (§5). It does not have to end the predictable-vs-unpredictable aim.
- Re-parameterise as in Marshall & Burgess 2015 and Barneche et al. 2018:
  - seasonal **fraction** α/(α+β), or Colwell's M/P;
  - residual **colour** θ, or lag-1 autocorrelation;
  - a separate total-magnitude term.
- Check the variance inflation factors of this triplet before giving up.
- Never residualise one variance on the other (Freckleton 2002).
- Record this as a documented deviation from the pre-registration.
- State the competing directional predictions, citing Rubio de Casas 2017 and Botero 2015: anticipatory dormancy tracks the seasonal fraction, while bet-hedging tracks low predictability or red-shifted residuals.

---

## Still-open claim

As of 17 Nov 2026, no published paper or preprint tests, **at the level of individual genomes across a GTDB-wide set of about 7,500 species and under phylogenetic regression**, whether sporulation, Rpf, ribosome-hibernation and stringent-response gene content scales with a **continuous, organism-independent measure of temporal temperature variability at the isolation site**. That claim remains novel.

Two parts of the gap memo's novelty statement must be narrowed:

1. **The predictable-vs-unpredictable decomposition** is not new to ecology (Colwell 1974; Marshall & Burgess 2015; Barneche et al. 2018; Tonkin et al. 2017). It is new only in its application to microbial genomes. In the current parameterisation it is **not identifiable** in this dataset (r = 0.82). It can be claimed only if a scale-free re-parameterisation (seasonal fraction + colour) proves separable.
2. **The substantive answer is currently "no detectable pooled relationship".** It sits alongside a sign reversal between MAGs and isolates. I found no prior report of such a reversal for a trait–environment slope, so it may itself be a publishable methodological finding.

The defensible claim is the first genome-level, phylogenetically controlled test of the microbial seed-bank prediction against isolation-site thermal variability. It shows that the pooled effect is practically null, and that the conclusion depends on genome provenance (MAG vs isolate) in ways predicted by known MAG-recovery and culturing biases. The marine TA deficit in complete isolate genomes, after genome-size and MGE correction, is a secondary finding that is partly novel. It refines Fucich & Chen 2020 and contradicts the pooled-"aquatic" result of Makarova 2009.
