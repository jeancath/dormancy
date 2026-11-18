#!/usr/bin/env python3
"""Build the dormancy gene dictionary, threshold table and concatenated HMM library.

Rows: (module, gene, source, acc, switching_class, evidence DOIs, notes, core_role)
core_role is only used for the sporulation_core.tsv side table:
  scored / unscored / extended / "" (not sporulation)
Model names and cutoffs are read from the downloaded HMM files (Pfam GA, NCBIfam TC)
and from KOfam ko_list (threshold, score_type).
"""
import csv, gzip, os, re, sys

HMM = '/home/user/data/hmm'
OUT = '/home/user/project-melad/dictionary'

# ---------------------------------------------------------------- references
GALP = '10.1128/jb.00079-22'
ABEC = '10.1128/jb.02110-12'
MEES = '10.1371/journal.pbio.1002341'
R = dict(
    burbulys='10.1016/0092-8674(91)90238-T', duncan93='10.1073/pnas.90.6.2325',
    duncan95='10.1126/science.270.5236.641', stragier88='10.1016/0092-8674(88)90407-2',
    stragier89='10.1126/science.2536191', karmazyn89='10.1101/gad.3.2.150',
    abanes02='10.1101/gad.1039902', chastanet07='10.1111/j.1365-2958.2007.05652.x',
    camp08='10.1111/j.1365-2958.2008.06289.x', meisner08='10.1073/pnas.0806301105',
    errington92='10.1099/00221287-138-12-2609', karow95='10.1073/pnas.92.6.2012',
    cutting90='10.1016/0092-8674(90)90362-I', rudner99='10.1073/pnas.96.26.14765',
    bagyan96='10.1128/jb.178.15.4500-4507.1996', tovar02='10.1128/JB.184.2.584-587.2002',
    mason86='10.1128/jb.167.1.174-178.1986', kunkel89='10.1101/gad.3.11.1735',
    popham91='10.1128/jb.173.24.7942-7949.1991', ramirez17='10.1371/journal.pgen.1007015',
    asai01='10.1099/00221287-147-4-919', roels92='10.1128/jb.174.2.575-585.1992',
    sanchez92='10.1128/jb.174.3.807-814.1992', daniel94='10.1016/S0022-2836(05)80027-0',
    daniel93='10.1006/jmbi.1993.1403', paidhungat00='10.1128/JB.182.9.2513-2519.2000',
    zheng88='10.1101/gad.2.8.1047', piggot86='10.1099/00221287-132-7-1875',
    mukamolova98='10.1073/pnas.95.15.8916', kana08='10.1111/j.1365-2958.2007.06078.x',
    ravagnani05='10.1186/1471-2164-6-39',
    wada90='10.1073/pnas.87.7.2657', ueta05='10.1111/j.1365-2443.2005.00903.x',
    ueta13='10.1111/gtc.12057', basu16='10.1093/nar/gkw180', zhang15='10.1038/nsmb.3103',
    basu17='10.1073/pnas.1709588114', yoshida12='10.1128/JB.00396-12',
    balon24='10.1038/s41586-024-07041-8', li18='10.1073/pnas.1804555115',
    atkinson11='10.1371/journal.pone.0023479', xiao91='10.1016/S0021-9258(19)67694-5',
    nanamiya08='10.1111/j.1365-2958.2007.06018.x', ruwe17='10.3389/fmicb.2017.01601',
    ruwe18='10.3389/fmicb.2018.00916', paul04='10.1016/j.cell.2004.07.009',
    korch03='10.1046/j.1365-2958.2003.03779.x', nguyen11='10.1126/science.1211037',
    moyed83='10.1128/jb.155.2.768-775.1983', germain13='10.1016/j.molcel.2013.08.045',
    balaban04='10.1126/science.1099390', rotem10='10.1073/pnas.1004333107',
    dorr10='10.1371/journal.pbio.1000317', verstraeten15='10.1016/j.molcel.2015.05.011',
    babin16='10.1073/pnas.1514412113',
    pedersen02='10.1046/j.1365-2958.2002.03027.x', zhang03='10.1016/S1097-2765(03)00402-7',
    fiebig10='10.1111/j.1365-2958.2010.07207.x', jorgensen09='10.1128/JB.01013-08',
    winther11='10.1073/pnas.1019587108', prysak09='10.1111/j.1365-2958.2008.06572.x',
    liu08='10.1073/pnas.0711949105', mutschler11='10.1371/journal.pbio.1001033',
    christensen06='10.1111/j.1365-2958.2006.05385.x', christensen04='10.1046/j.1365-2958.2003.03941.x',
    yamaguchi09='10.1074/jbc.M109.032904', bernard92='10.1016/0022-2836(92)90629-X',
    goormaghtigh18='10.1128/mBio.00640-18', gotfredsen98='10.1046/j.1365-2958.1998.00993.x',
    aizenman96='10.1073/pnas.93.12.6059', lehnherr93='10.1006/jmbi.1993.1521',
    roberts92='10.1128/jb.174.24.8119-8132.1992', brown09='10.1371/journal.ppat.1000706',
    skjerning19='10.1111/mmi.14150',
    siguier14='10.1111/1574-6976.12067', siguier06='10.1093/nar/gkj014',
    lange91='10.1111/j.1365-2958.1991.tb01825.x', boaretti03='10.1046/j.1462-2920.2003.00497.x',
)

def refs(*keys):
    return ';'.join(R.get(k, k) for k in keys)

# ---------------------------------------------------------------- dictionary rows
rows = []
def add(module, gene, source, acc, sw, ev, notes, role=''):
    rows.append(dict(module=module, gene=gene, source=source, acc=acc, sw=sw, ev=ev, notes=notes, role=role))

SIG = ('sigma-70 family paralogs (SigA/SigB/SigH/ECF); use the sporulation-specific equivalog '
       'and assign each protein to its best-scoring sporulation model (SigF vs SigG cross-hit)')
M = 'sporulation_core'
# gene, NCBIfam/Pfam primary, KOfam alt, refs, notes, role
core = [
 ('spo0A', 'NCBIfam:TIGR02875.1', 'K07699', refs('burbulys', GALP),
  'Master regulator; NECESSARY but not sufficient (Galperin 2022: present in all 76 sporeformers and in 40 non-sporeformers + 2 of unclear status among 180 Firmicutes; up to 3 paralogs in Clostridiales). '
  'K07699 is a response-regulator KO (OmpR-like receiver paralogs abound): TIGR02875 (TC seq 200 / dom 195) is the primary call', 'scored'),
 ('spoIIE', 'NCBIfam:TIGR02865.1', 'K06382', refs('duncan95'), 'PP2C phosphatase; paralogs RsbU/RsbX/PrpC share PP2C domain - equivalog TC excludes them', 'scored'),
 ('spoIIIJ', 'NCBIfam:NF002803.0', '', refs('errington92', GALP),
  'YidC/Oxa1 family insertase (YidC2); YidC paralogs are universal, PRK-derived model may hit non-sporeformer YidC2 (e.g. streptococci): NOT SCORED', 'unscored'),
 ('sigE', 'NCBIfam:TIGR02835.1', '', refs('stragier88'), SIG, 'scored'),
 ('sigF', 'NCBIfam:TIGR02885.1', '', refs('duncan93'), SIG, 'scored'),
 ('sigG', 'NCBIfam:TIGR02850.1', '', refs('karmazyn89'), SIG, 'scored'),
 ('spoIIAA', 'NCBIfam:TIGR02886.1', 'K06378', refs('duncan93'), 'STAS anti-anti-sigma; paralog RsbV/RsbS (TIGR00377 superfamily) - equivalog excludes', 'scored'),
 ('spoIIAB', 'NCBIfam:TIGR01925.2', 'K06379', refs('duncan93'), 'anti-sigma kinase; paralog RsbW', 'scored'),
 ('spoIIGA', 'NCBIfam:TIGR02854.1', 'K06383', refs('stragier88'), 'pro-SigE processing protease', 'scored'),
 ('spoIID', 'NCBIfam:TIGR02870.1', 'K06381', refs('abanes02', 'chastanet07'),
  'SpoIID domain (PF08486) also in non-sporulation amidases/LytB-like proteins; use equivalog', 'scored'),
 ('spoIIP', 'NCBIfam:TIGR02867.1', 'K06385', refs('abanes02', 'chastanet07'), 'most sporeformers carry two paralogs; count once', 'scored'),
 ('spoIIQ', 'NCBIfam:NF053738.1', 'K06386', refs('camp08', 'meisner08'),
  'M23-peptidase-like (degenerate); clostridial SpoIIQ is divergent and may be missed by Bacillus-trained models (C. difficile); M23 LytM paralogs abundant', 'scored'),
 ('spoIIIAC', 'NCBIfam:TIGR02848.1', 'K06392', refs('meisner08', 'camp08'), 'small membrane protein; lost in Erysipelotrichia sporeformers (Galperin 2022)', 'scored'),
 ('spoIIIAE', 'NCBIfam:TIGR02829.1', 'K06394', refs('meisner08', 'camp08'), '', 'scored'),
 ('spoIIIAG', 'NCBIfam:TIGR02830.1', 'K06396', refs('meisner08', 'camp08'), 'ring-forming channel protein', 'scored'),
 ('spoIIIAH', 'Pfam:PF12685', 'K06397', refs('meisner08', 'camp08'), 'YscJ-like ring; no TIGRFAM - Pfam SpoIIIAH GA + KOfam', 'scored'),
 ('spoIIR', 'NCBIfam:TIGR02837.1', 'K06387', refs('karow95'), 'lost in C. innocuum (Galperin 2022)', 'scored'),
 ('spoIVB', 'NCBIfam:TIGR02860.1', 'K06399', refs('cutting90'), 'S55 serine protease; PDZ-containing paralogs exist', 'scored'),
 ('spoVT', 'NCBIfam:TIGR02851.1', 'K04769', refs('bagyan96'),
  'N-terminal AbrB/MazE (PF04014) domain shared with AbrB and MazE antitoxins: never count a SpoVT hit as a MazE antitoxin', 'scored'),
 ('ydfS/yetF', 'Pfam:PF04239', '', refs(GALP),
  'YetF_C (UPF0702/DUF421) is widespread in non-sporeformers; only comparative (no verified experimental) evidence: NOT SCORED', 'unscored'),
 ('spoVAC', 'NCBIfam:TIGR02838.1', 'K06405', refs('tovar02'), 'paralog of SpoVAEB (same Pfam PF03862); best-hit assignment', 'scored'),
 ('spoVAD', 'NCBIfam:TIGR02845.1', 'K06406', refs('tovar02'), '', 'scored'),
 ('spoVAEB', 'NCBIfam:TIGR02839.1', 'K06407', refs('tovar02'), 'paralog of SpoVAC; best-hit assignment', 'scored'),
 ('sspA/B/C/D', 'Pfam:PF00269', 'K06418', refs('mason86', GALP),
  'alpha/beta-type SASP; multicopy; present in many Spo0A+ clostridial non-sporeformers (Galperin 2022) and absent in C. innocuum. Short (~65 aa): check prodigal calls', 'scored'),
 ('sigK', 'NCBIfam:TIGR02846.1', 'NCBIfam:NF004471.0', refs('stragier89'), SIG + '; sigK is interrupted by a skin prophage in B. subtilis and C. difficile (split ORFs in genomic prodigal calls). '
  'Full-sequence cutoff only (B. subtilis SigK splits into 2 HMMER domains, best domain < TC2). PRK05803 (NF004471, TC 252) rescues clostridial SigK (C. difficile 297 vs TIGR02846 TC 298.65) '
  'but SigE also scores >252 on PRK05803: competitive best-hit among the five sporulation-sigma models is mandatory', 'scored'),
 ('spoIIID', 'NCBIfam:TIGR02844.1', 'K06283', refs('kunkel89'), 'K06283 annotated DeoR-family: rely on TIGR02844 TC', 'scored'),
 ('spoIVFB', 'KOfam:K06402', '', refs('rudner99', 'cutting90'),
  'M50 site-2 protease family (PF02163) is ubiquitous and no TIGRFAM exists. Sanity test: K06402 gives a FALSE POSITIVE in M. tuberculosis (Rip3, 122.8 > 93.2) and misses C. difficile SpoIVFB: NOT SCORED', 'unscored'),
 ('spoVB/ykvU', 'NCBIfam:TIGR02900.1', 'K06409', refs('popham91', MEES),
  'MOP/MviN-like flippase superfamily (MurJ, polysaccharide exporters); equivalog TC required', 'scored'),
 ('yitE/yqfU', 'Pfam:PF02588', '', refs(GALP),
  'YitT_membrane (COG1284) is widespread in non-sporeformers; comparative evidence only: NOT SCORED', 'unscored'),
 ('ylmC/ymxH', 'NCBIfam:TIGR02888.1', '', refs(ABEC),
  'PRC-barrel (PF05239) superfamily; TIGR02888 restricted to endospore formers. Abecasis 2013 showed ylmC required for sporulation', 'scored'),
 ('ytaF', 'NCBIfam:TIGR02840.1', '', refs(ABEC), 'Abecasis 2013: new locus essential for sporulation', 'scored'),
 ('spoVV(ylbJ)', 'NCBIfam:TIGR02871.1', '', refs('ramirez17', MEES), 'DPA transporter; MEES 2016 ylbJ mutant 0% sporulation', 'scored'),
 ('yabP', 'NCBIfam:TIGR02892.1', '', refs('asai01', GALP),
  'YabP/YqfC family (paralog yqfC, TIGR02856 not included); evidence = sigma-E-dependent sporulation expression only (yabQ, not yabP, has the cortex phenotype)', 'scored'),
 ('yabQ', 'NCBIfam:TIGR02893.2', '', refs('asai01'), 'deletion blocks cortex formation', 'scored'),
 ('yqfD', 'NCBIfam:TIGR02876.1', 'K06438', refs(MEES), 'Meeske 2016: yqfD mutant ~0% sporulation, reduced sigma-G activity', 'scored'),
 ('spoIVA', 'NCBIfam:TIGR02836.1', 'K06398', refs('roels92'), 'coat basement-layer ATPase', 'scored'),
 ('gpr', 'NCBIfam:TIGR01441.2', 'K06012', refs('sanchez92'), 'germination protease degrading SASP; included in Galperin core', 'scored'),
]
for g, prim, ko, ev, notes, role in core:
    sw = 'exit_trigger' if g == 'gpr' else 'responsive'
    src, acc = prim.split(':')
    add(M, g, src, acc, sw, ev, notes, role)
    if ko:
        asrc, aacc = ko.split(':') if ':' in ko else ('KOfam', ko)
        add(M, g, asrc, aacc, sw, ev, 'alternative model for the same gene (gene present if ANY of its models passes cutoff). ' + notes, role)

M = 'sporulation_extended'
ext = [
 ('spoIIM', 'NCBIfam:TIGR02831.1', '', refs('chastanet07'), 'SpoIIM domain (PF01944) also in non-sporulating proteins; use equivalog. KOfam K06384 rejected: false positive in M. tuberculosis (74.4 > 41.43)'),
 ('spoIIIAA', 'NCBIfam:TIGR02858.1', 'K06390', refs('meisner08', 'camp08'), 'secretion ATPase-like; P-loop paralogs'),
 ('spoIIIAB', 'NCBIfam:TIGR02833.1', 'K06391', refs('meisner08', 'camp08'), ''),
 ('spoIIIAD', 'NCBIfam:TIGR02849.1', 'K06393', refs('meisner08', 'camp08'), ''),
 ('spoIIIAF', 'NCBIfam:TIGR02896.1', 'K06395', refs('meisner08', 'camp08'), ''),
 ('spoVD', 'NCBIfam:TIGR02214.1', 'K08384', refs('daniel94'), 'class B PBP; paralogs PBP2B/PBP3 - equivalog TC 688 required'),
 ('spoVE', 'NCBIfam:TIGR02615.1', '', refs('piggot86'), 'SEDS family (FtsW/RodA paralogs) - TIGR02615 is an "exception"-type model with high TC (448.8); do not use generic FtsW_RodA_SpoVE Pfam'),
 ('dpaA', 'NCBIfam:TIGR02853.1', 'K06410', refs('daniel93'), 'dipicolinate synthase; absent in many clostridial sporeformers which use EtfA-type DPA synthesis (Galperin 2022) - not a necessary marker'),
 ('dpaB', 'NCBIfam:TIGR02852.1', 'K06411', refs('daniel93'), 'flavoprotein paralogs (dfp/CoaBC) - equivalog only'),
 ('gerA_A', 'Pfam:PF03323', '', refs('paidhungat00'), 'GerA family A-subunit (germinant receptor); multiple paralogous operons (gerA/gerB/gerK) - count loci, not genes'),
 ('gerA_C', 'NCBIfam:TIGR02887.1', '', refs('paidhungat00'), 'Ger(x)C lipoprotein subunit; subfamily model'),
 ('cotE', 'Pfam:PF10628', 'K06328', refs('zheng88'), 'Bacillus-specific outer coat morphogenetic protein; absent in clostridia (not a general marker)'),
]
for g, prim, ko, ev, notes in ext:
    sw = 'exit_trigger' if g.startswith('ger') else 'responsive'
    src, acc = prim.split(':')
    add(M, g, src, acc, sw, ev, notes, 'extended')
    if ko:
        add(M, g, 'KOfam', ko, sw, ev, 'alternative model. ' + notes, 'extended')

M = 'resuscitation'
add(M, 'rpf', 'Pfam', 'PF06737', 'exit_trigger', refs('mukamolova98', 'kana08'),
    'Rpf (lysozyme-like) domain; Actinobacteria-restricted. Transglycosylas weakly matches other lytic transglycosylases (SLT PF01464, goose-type lysozyme): require best hit to PF06737 and report taxon')
add(M, 'rpfB', 'KOfam', 'K21688', 'exit_trigger', refs('kana08'),
    'RpfB-type = Rpf domain (PF06737) + G5 (PF07501) + DUF348 (PF03990) in one protein. K21688 ALSO HITS firmicute Sps protein YabE (B. subtilis 212 > 166.7) via the shared DUF348-G5 region: call RpfB only if PF06737 is in the same protein, otherwise reassign to sps')
add(M, 'rpfB_G5', 'Pfam', 'PF07501', 'exit_trigger', refs('kana08'),
    'ARCHITECTURE-ONLY domain (G5 occurs in many unrelated surface proteins incl. IgA proteases): never count alone')
add(M, 'rpfB_DUF348', 'Pfam', 'PF03990', 'exit_trigger', refs('kana08', 'ravagnani05'),
    'ARCHITECTURE-ONLY domain (G5-linked ubiquitin-like); present in RpfB and in firmicute YabE/Sps proteins: never count alone')
add(M, 'sps(yabE/yocH/yuiC)', 'Pfam', 'PF06725', 'unclear', refs('ravagnani05'),
    'Sps domain = Pfam 3D domain (Ravagnani 2005, COG3584). 3D also forms part of MltA lytic transglycosylases in Proteobacteria: count only in Bacillota or when MltA (PF03562) is absent. '
    'Experimental evidence weak (stationary-phase survival defect of yocH / quadruple sps mutant reported as unpublished data) - treat as candidate')

M = 'ribosome_hibernation'
add(M, 'hpf', 'NCBIfam', 'TIGR00741.1', 'hibernation', refs('ueta05', 'ueta13', 'basu16'),
    'HPF/YfiA(RaiA) family equivalog; hits both short HPF and RaiA paralogs in Gammaproteobacteria')
add(M, 'hpf', 'Pfam', 'PF02482', 'hibernation', refs('ueta05', 'ueta13', 'basu16'),
    'Ribosomal_S30AE domain (named sigma-54 modulation protein in older annotations); present in nearly all bacteria - presence weakly informative, use long/short typing')
add(M, 'hpf_long_Cterm', 'Pfam', 'PF16321', 'hibernation', refs('ueta13', 'basu16'),
    'long-HPF C-terminal dimerisation domain: long HPF = PF02482 + PF16321 in same protein (forms 100S without RMF; Firmicutes, Cyanobacteria etc.)')
add(M, 'hpf', 'KOfam', 'K05808', 'hibernation', refs('ueta05'), 'KOfam HPF (alternative)')
add(M, 'raiA', 'KOfam', 'K05809', 'hibernation', refs('ueta05'), 'RaiA/YfiA (E. coli) - antagonises 100S formation, still a 70S hibernation factor')
add(M, 'rmf', 'Pfam', 'PF04957', 'hibernation', refs('wada90', 'ueta05'), 'short (55 aa) gammaproteobacterial factor; check that prodigal calls short ORFs')
add(M, 'rmf', 'KOfam', 'K03812', 'hibernation', refs('wada90'), 'alternative model (KOfam profile_type trim)')
add(M, 'balon', 'NCBIfam', 'NF052287.1', 'hibernation', refs('balon24'),
    'AOC03_06830 family = Psychrobacter urativorans Balon (UniProt A0A0M3V8U3); NCBIfam built from Helena-Bueno 2024 (PMID 38355796); aeRF1-like fold')
add(M, 'balon', 'NCBIfam', 'NF052289.1', 'hibernation', refs('balon24'),
    'Rv2629 family = actinobacterial Balon homologue (paper shows mycobacterial homologues bind ribosomes)')
add(M, 'balon', 'NCBIfam', 'NF052288.1', 'hibernation', refs('balon24'), 'MSMEG_1130 family Balon homologue (Actinomycetota)')
add(M, 'balon', 'Pfam', 'PF18845', 'hibernation', refs('balon24'),
    'Balon_mid (formerly a baeRF family); InterPro lists only ~1,650 proteins, far fewer than the ~20% of bacteria claimed in the paper, so sensitivity is limited; '
    'if coverage looks low build a profile: jackhmmer (3 iters, E<=1e-10) seeded with A0A0M3V8U3 and the paper\'s homologue set against UniRef50, exclude aeRF1/Pelota (PF03463-PF03465) and baeRF families, calibrate a TC on known Balon+/- genomes')
add(M, 'mrf', 'NCBIfam', 'NF047431.1', 'hibernation', refs('li18'),
    'mycobacterial hibernation-factor recruiting GTPase (zinc-depletion hibernation); narrow taxonomic range')
add(M, 'yqjD/elaB', 'Pfam', 'PF05957', 'unclear', refs('yoshida12'),
    'DUF883 family (YqjD, ElaB, YgaM); only YqjD shown to associate with stationary-phase 100S ribosomes - family-level hit is weak evidence')
add(M, 'hflX', 'NCBIfam', 'TIGR03156.1', 'exit_trigger', refs('basu17', 'zhang15'),
    'splits hibernating 100S/70S ribosomes; near-universal in bacteria - presence uninformative, report copy number only')
add(M, 'hflX', 'KOfam', 'K03665', 'exit_trigger', refs('basu17'), 'alternative model')

M = 'stringent_response'
add(M, 'rsh_long', 'NCBIfam', 'TIGR00691.1', 'responsive', refs('xiao91', 'atkinson11'),
    'RelA/SpoT/Rel long RSH subfamily model (TC 499.2): hits RelA, SpoT and bifunctional Rel; does NOT hit SAS/SAH')
add(M, 'relA', 'KOfam', 'K00951', 'responsive', refs('xiao91', 'atkinson11'),
    'KOfam GTP pyrophosphokinase (RelA and firmicute Rel); RelA vs SpoT split is lineage-specific (Atkinson 2011): analyse as long-RSH count')
add(M, 'spoT', 'KOfam', 'K01139', 'responsive', refs('xiao91', 'atkinson11'), 'bifunctional SpoT; see relA note')
add(M, 'sas(relP/relQ)', 'KOfam', 'K07816', 'responsive', refs('nanamiya08', 'ruwe17', 'atkinson11'),
    'small alarmone synthetase (YjbM/RelQ, YwaC/RelP). Architecture rule as backup: PF04607 present, TGS/HD_4/ACT absent, length < 350 aa')
add(M, 'rsh_SYNTH', 'Pfam', 'PF04607', 'responsive', refs('atkinson11'),
    'synthetase domain present in long RSH AND SAS (and in aTA SAS toxins): use only for architecture typing')
add(M, 'rsh_HD', 'Pfam', 'PF13328', 'responsive', refs('atkinson11', 'ruwe18'),
    'RSH hydrolase HD domain; SAH (small alarmone hydrolase) = PF13328 present, PF04607 absent, length < 300 aa (Ruwe 2018). No dedicated SAH model exists')
add(M, 'rsh_TGS', 'Pfam', 'PF02824', 'responsive', refs('atkinson11'),
    'ARCHITECTURE-ONLY: TGS also in ThrRS and Obg/GTPases - defines long RSH only together with PF04607')
add(M, 'rsh_ACT', 'Pfam', 'PF13291', 'responsive', refs('atkinson11'), 'ARCHITECTURE-ONLY: ACT domains are ubiquitous in amino-acid metabolism enzymes')
add(M, 'dksA', 'NCBIfam', 'TIGR02420.1', 'responsive', refs('paul04'),
    'equivalog; excludes conjugative TraR and most DksA2 (zf-dskA_traR C4 zinc finger PF01258 is non-specific - not used)')

M = 'persistence_stress'
add(M, 'hipA', 'Pfam', 'PF07804', 'stochastic', refs('moyed83', 'germain13', 'korch03', 'balaban04'),
    'HipA_C kinase domain; also in non-TA HipA-like kinases (YjjJ, HipA_2 family) - prefer pairing with HipB/HTH; PF07805 (named in plan) does NOT exist in current Pfam, N-terminal domain is PF13657')
add(M, 'hipA', 'Pfam', 'PF13657', 'stochastic', refs('moyed83', 'germain13'), 'Couple_hipA (HipA N-terminal); also found alone in some proteins - use with PF07804')
add(M, 'hipA', 'KOfam', 'K07154', 'stochastic', refs('moyed83', 'germain13'), 'KOfam HipA (low threshold 28.17, domain score) - noisy')
add(M, 'tisB', 'NCBIfam', 'TIGR04459.1', 'responsive', refs('dorr10'),
    'SOS-induced type I toxin, 29 aa: prodigal minimum gene length (~90 nt) means it is usually NOT called; IstR-1/2 antisense sRNAs cannot be detected by protein HMMs (would need Rfam/Infernal). Enterobacterales only. Report but do not analyse')
add(M, 'tisB', 'Pfam', 'PF13939', 'responsive', refs('dorr10'), 'see TIGR04459 note')
add(M, 'obgE', 'NCBIfam', 'TIGR02729.1', 'stochastic', refs('verstraeten15'),
    'essential, universal GTPase - presence carries no information; kept only for QC (expected 1 copy per genome), EXCLUDED from scores')
add(M, 'sutA', 'NCBIfam', 'NF041490.1', 'responsive', refs('babin16'),
    'Pseudomonas slow-growth/dormancy transcription factor (Gammaproteobacteria only)')

M = 'TA_typeII_toxin'
tox = [
 ('relE/parE', 'PF05016', refs('pedersen02', 'fiebig10'), 'RelE/ParE superfamily (ParE_toxin) incl. many unpaired orphans'),
 ('mazF/pemK', 'PF02452', refs('zhang03', 'aizenman96'), 'PemK_toxin; also hits CcdB-like and Kid/PemK of plasmids'),
 ('hicA', 'PF07927', refs('jorgensen09'), 'dsRBD-like small toxin'),
 ('hipA', 'PF07804', refs('moyed83', 'germain13'), 'shared with persistence_stress module'),
 ('vapC', 'PF01850', refs('winther11'), 'PIN domain is HIGHLY PROMISCUOUS (non-TA nucleases, Nob1, rRNA processing): count only when paired with an antitoxin'),
 ('yafQ', 'PF15738', refs('prysak09'), 'RelE superfamily'),
 ('doc', 'PF02661', refs('liu08', 'lehnherr93'), 'Fic/DOC family is HIGHLY PROMISCUOUS (most Fic proteins are AMPylators/effectors, not TA): count only when paired with Phd/YefM'),
 ('zeta', 'PF06414', refs('mutschler11'), 'P-loop kinase; pair with epsilon antitoxin'),
 ('higB', 'PF05015', refs('christensen06'), 'HigB-like RelE fold'),
 ('higB', 'PF09907', refs('christensen06'), 'second HigB family'),
 ('yoeB', 'PF06769', refs('christensen04'), 'RelE superfamily'),
 ('mqsR', 'PF15723', refs('yamaguchi09', 'brown09'), ''),
 ('ccdB', 'PF01845', refs('bernard92'), 'mostly plasmid (F) encoded - plasmid-borne copies inflate counts'),
 ('RES', 'PF08808', refs('skjerning19'), 'RES domain also in non-TA proteins: pairing with Xre/HTH or Phd antitoxin required'),
]
for g, acc, ev, notes in tox:
    sw = 'stochastic' if g == 'hipA' else 'unclear'
    add(M, g, 'Pfam', acc, sw, ev + ';' + R['goormaghtigh18'],
        'TA module analysed SEPARATELY (confounded with mobile elements/plasmids). ' + notes)
add(M, 'tacT/ataT (GNAT)', 'none', '', 'unclear', refs('10.1016/j.molcel.2016.05.002'),
    'NO usable model: GNAT acetyltransferase toxins (TacT/AtaT/ItaT) share Acetyltransf_1/7 domains with hundreds of non-toxin GNATs per genome; excluded')

M = 'TA_typeII_antitoxin'
anti = [
 ('relB/dinJ', 'Pfam', 'PF04221', refs('gotfredsen98'), ''),
 ('mazE', 'Pfam', 'PF04014', refs('aizenman96'),
  'AbrB/MazE/SpoVT swapped-hairpin domain is PROMISCUOUS: also AbrB, SpoVT (sporulation core!), MraZ-like and many transcription regulators. Exclude proteins best-hit by TIGR02851 (SpoVT) and count only when paired with a toxin'),
 ('phd/yefM', 'Pfam', 'PF02604', refs('lehnherr93', 'christensen04'), ''),
 ('hicB', 'Pfam', 'PF15919', refs('jorgensen09'), ''),
 ('hicB', 'Pfam', 'PF05534', refs('jorgensen09'), 'HicB family incl. non-paired homologues'),
 ('parD', 'Pfam', 'PF03693', refs('roberts92', 'fiebig10'), 'ribbon-helix-helix; RHH regulators are common - require pairing'),
 ('ccdA', 'Pfam', 'PF07362', refs('bernard92'), 'plasmid-associated'),
 ('epsilon', 'Pfam', 'PF08998', refs('mutschler11'), ''),
 ('mqsA', 'Pfam', 'PF15731', refs('brown09'), ''),
 ('higA', 'NCBIfam', 'TIGR02607.1', refs('christensen06'), 'HigA HTH (XRE/HTH_3 is ubiquitous; this equivalog restricts to HigA)'),
 ('hipB', 'KOfam', 'K15773', refs('moyed83', 'germain13'), 'HTH XRE; low-specificity KO (trim profile) - only count adjacent to hipA'),
]
for g, src, acc, ev, notes in anti:
    add(M, g, src, acc, 'unclear', ev, 'TA module analysed SEPARATELY; locus counted only if a toxin is within 3 genes. ' + notes)

M = 'covariate_mobile'
tnp = [
 ('IS4/IS5_family_DDE', 'PF01609'), ('rve_integrase_core(IS3/IS21/IS481)', 'PF00665'),
 ('rve_2', 'PF13333'), ('rve_3', 'PF13683'), ('IS3_OrfA_HTH', 'PF01527'),
 ('IS110_C', 'PF02371'), ('IS110_N', 'PF01548'), ('IS66', 'PF03050'),
 ('DDE_Tnp_1_2', 'PF13586'), ('DDE_Tnp_1_3', 'PF13612'), ('DDE_Tnp_1_5', 'PF13737'),
 ('DDE_Tnp_1_6', 'PF13751'), ('ISL3', 'PF01610'), ('IS605_OrfB/TnpB', 'PF01385'),
 ('IS200_TnpA', 'PF01797'), ('IS200_Y2', 'PF04986'), ('IS256_mutator', 'PF00872'),
 ('IS1', 'PF03400'), ('IS5_N(DUF772)', 'PF05598'), ('DDE_Tnp_4', 'PF13359'),
 ('DDE_5', 'PF13546'), ('Tn3', 'PF01526'), ('Mu_DDE', 'PF02914'), ('Mu_C', 'PF09299'),
 ('IS240', 'PF13610'),
]
for g, acc in tnp:
    extra = ''
    if acc == 'PF00665':
        extra = ' rve also matches retroviral-type integrases of phages/Mu-like elements.'
    if acc in ('PF01797', 'PF04986'):
        extra = ' Y1/Y2 TnpA also match domesticated REP-associated tyrosine transposases (RAYTs), which are not mobile.'
    if acc == 'PF01527':
        extra = ' Accessory HTH of IS3 OrfA; counts only via protein-level union.'
    add(M, 'transposase:' + g, 'Pfam', acc, 'unclear', refs('siguier14', 'siguier06'),
        'COVARIATE (not a switching gene). Count distinct proteins with any transposase hit (union, dedup per protein), per 1000 CDS; fragmented short-read assemblies collapse IS copies - include assembly N50/contig count.' + extra)
for g, acc in [('integrase_tyrosine', 'PF00589'), ('resolvase_serine', 'PF00239'), ('recombinase_serine_large', 'PF07508')]:
    add(M, g, 'Pfam', acc, 'unclear', refs('siguier14'),
        'COVARIATE. Tyrosine-recombinase hits include chromosomal XerC/XerD: subtract proteins best-hit by TIGR02224/TIGR02225')
add(M, 'xerC_exclusion', 'NCBIfam', 'TIGR02224.1', 'unclear', refs('siguier14'), 'NEGATIVE FILTER: housekeeping XerC, removed from integrase counts')
add(M, 'xerD_exclusion', 'NCBIfam', 'TIGR02225.1', 'unclear', refs('siguier14'), 'NEGATIVE FILTER: housekeeping XerD, removed from integrase counts')

M = 'covariate_stress'
add(M, 'rpoS', 'KOfam', 'K03087', 'responsive', refs('lange91', 'boaretti03'),
    'ANNOTATION COVARIATE only. sigma-70 family: RpoD/SigB paralogs score below KOfam threshold 425.6 but verify; RpoS essentially restricted to Pseudomonadota')
add(M, 'rpoS', 'NCBIfam', 'TIGR02394.1', 'responsive', refs('lange91', 'boaretti03'), 'Pseudomonadota-restricted equivalog (alternative)')

# ---------------------------------------------------------------- resolve models
kol = {}
with gzip.open(f'{HMM}/raw/ko_list.gz', 'rt') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        kol[r['knum']] = r

def hmm_path(src, acc):
    if src == 'Pfam':
        return f'{HMM}/pfam/{acc}.hmm'
    if src == 'NCBIfam':
        return f'{HMM}/ncbifam/{acc}.HMM'
    if src == 'KOfam':
        return f'{HMM}/kofam/{acc}.hmm'
    return None

def parse_hmm(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('HMM '):
                break
            k = line[:6].strip()
            if k in ('NAME', 'ACC', 'DESC', 'LENG', 'GA', 'TC', 'NC'):
                d[k] = line[6:].strip()
    return d

FULL_ONLY = {'TIGR02846.1', 'NF004471.0'}  # sigK: domain split, see notes
def nc(h):
    return h['NC'].rstrip(';').split()[0] if 'NC' in h else ''
models = {}
for r in rows:
    if r['source'] == 'none':
        r.update(name='', cutoff_type='none', cutoff='', st='', thr_full='', thr_dom='')
        continue
    p = hmm_path(r['source'], r['acc'])
    if not os.path.exists(p):
        sys.exit(f'missing HMM {p}')
    h = parse_hmm(p)
    r['name'] = h['NAME']
    r['nc'] = nc(h)
    if r['source'] == 'Pfam':
        a, b = h['GA'].rstrip(';').split()
        r.update(cutoff_type='GA', cutoff=f'{a};{b}', st='', thr=a, score_type='full', thr_dom=b)
        r['acc'] = h['ACC']  # versioned accession
    elif r['source'] == 'NCBIfam':
        a, b = h['TC'].rstrip(';').split()
        r.update(cutoff_type='TC', cutoff=f'{a};{b}', st='', thr=a, score_type='full', thr_dom=b)
        if r['acc'] in FULL_ONLY:
            r['thr_dom'] = ''
    else:
        k = kol[r['acc']]
        if k['threshold'] == '-':
            sys.exit(f'{r["acc"]} has no KOfam threshold')
        r.update(cutoff_type='KOfam_threshold', cutoff=k['threshold'], st=k['score_type'],
                 thr=k['threshold'], score_type=k['score_type'], thr_dom='')
    models.setdefault(r['name'], (p, r))

# ---------------------------------------------------------------- write outputs
os.makedirs(OUT, exist_ok=True)
cols = ['module', 'gene', 'model_source', 'model_acc', 'model_name', 'cutoff_type', 'cutoff_value',
        'kofam_score_type', 'switching_class', 'evidence_ref', 'notes']
with open(f'{OUT}/dormancy_dictionary.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(cols)
    for r in rows:
        w.writerow([r['module'], r['gene'], r['source'], r['acc'], r['name'], r['cutoff_type'],
                    r['cutoff'], r['st'], r['sw'], r['ev'], r['notes']])

with open(f'{OUT}/thresholds.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['model_name', 'threshold', 'score_type', 'domain_threshold', 'noise_cutoff', 'cutoff_type', 'model_acc', 'model_source'])
    for name, (p, r) in models.items():
        w.writerow([name, r['thr'], r['score_type'], r['thr_dom'], r['nc'], r['cutoff_type'], r['acc'], r['source']])

# sporulation core side table (one row per gene)
seen = {}
for r in rows:
    if r['role'] in ('scored', 'unscored', 'extended'):
        seen.setdefault((r['gene'], r['role']), []).append(r['name'])
with open(f'{OUT}/sporulation_core.tsv', 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['gene', 'set', 'in_score', 'model_names'])
    for (g, role), names in seen.items():
        st = 'galperin2022_core' if role in ('scored', 'unscored') else 'extended'
        w.writerow([g, st, 'yes' if role == 'scored' else 'no', ','.join(names)])

# concatenated library
with open(f'{HMM}/dormancy_library.hmm', 'w') as out:
    for name, (p, r) in models.items():
        txt = open(p).read()
        if txt.count('\n//') != 1:
            sys.exit(f'{p} does not contain exactly one model')
        out.write(txt if txt.endswith('\n') else txt + '\n')

print('rows', len(rows), 'models', len(models),
      'core_scored', sum(1 for k in seen if k[1] == 'scored'),
      'core_unscored', sum(1 for k in seen if k[1] == 'unscored'),
      'extended', sum(1 for k in seen if k[1] == 'extended'))
