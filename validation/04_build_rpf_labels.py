#!/usr/bin/env python3
"""Literature-documented Rpf (and Sps) gene counts for the Rpf check.

Every row cites the paper that states the count, or presence/absence, for that genome. The main
source is Ravagnani, Finan & Young 2005, BMC Genomics 6:39, doi:10.1186/1471-2164-6-39:
  Table 1A  genomes with Rpf-domain genes (count per genome)
  Table 1B  genomes with only a domain *distantly related* to Rpf, which the Rpf HMM did not
            detect ("absent from ... Bifidobacterium, Thermobifida and Tropheryma") -> canonical
            Rpf = 0
  Table 2A  Firmicutes with Sps-domain genes (count per genome); text: "Rpf-like proteins are not
            found in firmicutes" -> Rpf = 0
Second source: Machowski et al. 2014, BMC Microbiol 14:75, doi:10.1186/1471-2180-14-75, Table 1
(rpf loci in 19 mycobacterial genomes). Other sources are given per row. `expected_rpf` = documented number of Rpf-domain genes (NA =
presence documented but number not stated); `expected_sps` likewise for Sps (Ravagnani Table 2A).
Accessions were matched to the strain named in the source and checked against NCBI Datasets
(organism and strain names).

usage: python3 validation/04_build_rpf_labels.py  -> validation/rpf_labels.tsv, /tmp/val/rpf_accessions.tsv
"""
import pandas as pd

RAV = 'Ravagnani2005 BMC Genomics 6:39 doi:10.1186/1471-2164-6-39'
MAC = 'Machowski2014 BMC Microbiol 14:75 doi:10.1186/1471-2180-14-75'
L = [
    # organism, strain, accession, phylum, expected_rpf, expected_sps, source, evidence
    ('Mycobacterium tuberculosis', 'H37Rv', 'GCF_000195955.2', 'Actinomycetota', 5, 0, RAV + '; Kana et al. 2008 Mol Microbiol doi:10.1111/j.1365-2958.2007.06078.x', 'Table 1A: 5 (RpfA-E)'),
    ('Mycobacterium bovis', 'AF2122/97', 'GCF_000195835.3', 'Actinomycetota', 5, 0, RAV, 'Table 1A: 5'),
    ('Mycobacterium leprae', 'TN', 'GCF_000195855.1', 'Actinomycetota', 3, 0, RAV, 'Table 1A: 3'),
    ('Mycobacterium avium subsp. paratuberculosis', 'K-10', 'GCF_000007865.1', 'Actinomycetota', 4, 0, RAV, 'Table 1A: 4 (NC_002944)'),
    ('Mycobacterium marinum', 'M', 'GCF_000018345.1', 'Actinomycetota', 4, 0, RAV, 'Table 1A: 4 (genome unfinished in 2005)'),
    ('Mycolicibacterium smegmatis', 'mc2 155', 'GCF_000015005.1', 'Actinomycetota', 4, 0, RAV, 'Table 1A: 4 (genome unfinished in 2005)'),
    ('Corynebacterium diphtheriae', 'NCTC 13129', 'GCF_000195815.1', 'Actinomycetota', 3, 0, RAV, 'Table 1A: 3'),
    ('Corynebacterium glutamicum', 'ATCC 13032', 'GCF_000011325.1', 'Actinomycetota', 2, 0, RAV + '; Hartmann et al. 2004', 'Table 1A: 2'),
    ('Corynebacterium efficiens', 'YS-314', 'GCF_000011305.1', 'Actinomycetota', 2, 0, RAV, 'Table 1A: 2'),
    ('Micrococcus luteus', 'NCTC 2665', 'GCF_000023205.1', 'Actinomycetota', 1, 0, RAV + '; Mukamolova et al. 1998 PNAS doi:10.1073/pnas.95.15.8916', 'Table 1A: 1 (genome not sequenced in 2005; count from Mukamolova 1998)'),
    ('Streptomyces coelicolor', 'A3(2)', 'GCF_000203835.1', 'Actinomycetota', 5, 0, RAV + '; Sexton et al. 2015 J Bacteriol', 'Table 1A: 5 (plus 2 distantly related, Table 1B)'),
    ('Streptomyces avermitilis', 'MA-4680', 'GCF_000009765.2', 'Actinomycetota', 6, 0, RAV, 'Table 1A: 6 (plus 3 distantly related, Table 1B)'),
    ('Saccharopolyspora erythraea', 'NRRL 2338', 'GCF_000062885.1', 'Actinomycetota', None, 0, RAV, 'text: Rpf-like proteins present in Saccharopolyspora (count not given)'),
    ('Bifidobacterium longum', 'NCC2705', 'GCF_000007525.1', 'Actinomycetota', 0, 0, RAV, 'Table 1B: 3 distantly related; text: Rpf absent from Bifidobacterium'),
    ('Tropheryma whipplei', 'Twist', 'GCF_000007485.1', 'Actinomycetota', 0, 0, RAV, 'Table 1B: 2 distantly related; text: Rpf absent from Tropheryma'),
    ('Thermobifida fusca', 'YX', 'GCF_000012405.1', 'Actinomycetota', 0, 0, RAV, 'text: Rpf absent from Thermobifida'),
    # Machowski et al. 2014 BMC Microbiol 14:75, doi:10.1186/1471-2180-14-75, Table 1: rpf loci per
    # genome (rpfA-E; a cell with two locus tags, e.g. MSMEG_4643/MSMEG_4640, is counted as two genes,
    # which matches Ravagnani's count of 4 for M. smegmatis)
    ('Mycobacterium tuberculosis', 'H37Ra', 'GCF_000016145.1', 'Actinomycetota', 5, 0, MAC, 'Table 1: rpfA-E'),
    ('Mycobacterium tuberculosis', 'CDC1551', 'GCF_000008585.1', 'Actinomycetota', 5, 0, MAC, 'Table 1: rpfA-E'),
    ('Mycobacterium africanum', 'GM041182', 'GCF_000253355.1', 'Actinomycetota', 5, 0, MAC, 'Table 1: rpfA-E'),
    ('Mycobacterium bovis BCG', 'Pasteur 1173P2', 'GCF_000009445.1', 'Actinomycetota', 5, 0, MAC, 'Table 1: rpfA-E'),
    ('Mycobacterium avium', '104', 'GCF_000014985.1', 'Actinomycetota', 4, 0, MAC, 'Table 1: rpfA,B,C,E (no rpfD)'),
    ('Mycobacterium intracellulare', 'ATCC 13950', 'GCF_000277125.1', 'Actinomycetota', 4, 0, MAC, 'Table 1: rpfA,B,C,E'),
    ('Mycobacterium ulcerans', 'Agy99', 'GCF_000013925.1', 'Actinomycetota', 4, 0, MAC, 'Table 1: rpfA,B,C,E'),
    ('Mycobacteroides abscessus', 'ATCC 19977', 'GCF_000069185.1', 'Actinomycetota', 4, 0, MAC, 'Table 1: rpfA,B,C,E'),
    ('Mycolicibacterium vanbaalenii', 'PYR-1', 'GCA_000015305.1', 'Actinomycetota', 4, 0, MAC, 'Table 1: rpfA,B, rpfE x2 (RefSeq copy suppressed)'),
    ('Mycobacterium sp.', 'MCS', 'GCA_000014165.1', 'Actinomycetota', 5, 0, MAC, 'Table 1: rpfA, rpfB x2, rpfE x2 (assembly suppressed)'),
    ('Mycobacterium sp.', 'KMS', 'GCA_000015405.1', 'Actinomycetota', 5, 0, MAC, 'Table 1: rpfA, rpfB x2, rpfE x2 (assembly suppressed)'),
    ('Mycobacterium sp.', 'JLS', 'GCA_000016005.1', 'Actinomycetota', 5, 0, MAC, 'Table 1: rpfA, rpfB x2, rpfE x2 (assembly suppressed)'),
    ('Mycolicibacterium gilvum', 'PYR-GCK', 'GCA_000016365.1', 'Actinomycetota', 4, 0, MAC, 'Table 1: rpfA,B, rpfE x2 (assembly suppressed)'),
    # Firmicutes: Rpf 0 (text), Sps counts Table 2A
    ('Bacillus anthracis', 'Ames', 'GCF_000007845.1', 'Bacillota', 0, 6, RAV, 'Table 2A: Sps 6'),
    ('Bacillus cereus', 'ATCC 10987', 'GCA_000008005.1', 'Bacillota', 0, 6, RAV, 'Table 2A: Sps 6 (RefSeq copy suppressed; GenBank used)'),
    ('Bacillus cereus', 'ATCC 14579', 'GCF_000007825.1', 'Bacillota', 0, 5, RAV, 'Table 2A: Sps 5'),
    ('Halalkalibacterium halodurans', 'C-125', 'GCF_000011145.1', 'Bacillota', 0, 3, RAV, 'Table 2A: Sps 3'),
    ('Bacillus subtilis', '168', 'GCF_000009045.1', 'Bacillota', 0, 4, RAV, 'Table 2A: Sps 4 (YabE, YocH, YuiC, YorM; YorM on SPbeta)'),
    ('Oceanobacillus iheyensis', 'HTE831', 'GCF_000011245.1', 'Bacillota', 0, 4, RAV, 'Table 2A: Sps 4 (+1 distantly related)'),
    ('Listeria innocua', 'Clip11262', 'GCF_000195795.1', 'Bacillota', 0, 2, RAV, 'Table 2A: Sps 2'),
    ('Listeria monocytogenes', 'EGD-e', 'GCF_000196035.1', 'Bacillota', 0, 2, RAV, 'Table 2A: Sps 2'),
    ('Enterococcus faecalis', 'V583', 'GCF_000007785.1', 'Bacillota', 0, 1, RAV, 'Table 2A: Sps 1'),
    ('Lactococcus lactis subsp. lactis', 'IL1403', 'GCF_000006865.1', 'Bacillota', 0, 1, RAV, 'Table 2A: Sps 1'),
    ('Clostridium acetobutylicum', 'ATCC 824', 'GCF_000008765.1', 'Bacillota', 0, 2, RAV, 'Table 2A: Sps 2'),
    ('Clostridium botulinum A', 'ATCC 3502', 'GCF_000063585.1', 'Bacillota', 0, 2, RAV, 'Table 2A: Sps 2 (genome unfinished in 2005)'),
    ('Clostridium perfringens', '13', 'GCF_000009685.1', 'Bacillota', 0, 3, RAV, 'Table 2A: Sps 3'),
    ('Clostridium tetani', 'E88', 'GCF_000007625.1', 'Bacillota', 0, 2, RAV, 'Table 2A: Sps 2'),
    ('Acetivibrio thermocellus', 'ATCC 27405', 'GCF_000015865.1', 'Bacillota', 0, 4, RAV, 'Table 2A: Sps 4 (draft AABG03 in 2005)'),
    ('Desulfitobacterium hafniense', 'DCB-2', 'GCF_000021925.1', 'Bacillota', 0, 1, RAV, 'Table 2A: Sps 1 (draft AAAW in 2005)'),
    ('Caldanaerobacter subterraneus subsp. tengcongensis', 'MB4', 'GCF_000007085.1', 'Bacillota', 0, 1, RAV, 'Table 2A: Sps 1'),
    ('Staphylococcus aureus', 'N315', 'GCF_000009645.1', 'Bacillota', 0, 0, RAV, 'Table 1B: 1 distantly related Rpf-like; no Rpf domain'),
    ('Staphylococcus epidermidis', 'ATCC 12228', 'GCF_000007645.1', 'Bacillota', 0, 0, RAV, 'Table 1B: 1 distantly related Rpf-like; no Rpf domain'),
    ('Oenococcus oeni', 'PSU-1', 'GCF_000014385.1', 'Bacillota', 0, 0, RAV, 'Table 1B: 1 distantly related (strain in 2005 was a draft; PSU-1 used)'),
    # non-Gram-positive controls: Rpf domain "confined to the actinobacteria" (Ravagnani 2005, text)
    ('Escherichia coli', 'K-12 MG1655', 'GCF_000005845.2', 'Pseudomonadota', 0, None, RAV, 'text: Rpf family confined to actinobacteria'),
    ('Pseudomonas putida', 'KT2440', 'GCF_000007565.2', 'Pseudomonadota', 0, None, RAV, 'text: confined to actinobacteria'),
    ('Vibrio cholerae', 'N16961', 'GCF_000006745.1', 'Pseudomonadota', 0, None, RAV, 'text: confined to actinobacteria'),
    ('Caulobacter vibrioides', 'CB15', 'GCF_000006905.1', 'Pseudomonadota', 0, None, RAV, 'text: confined to actinobacteria'),
    ('Deinococcus radiodurans', 'R1', 'GCF_000008565.1', 'Deinococcota', 0, None, RAV, 'text: confined to actinobacteria; Table 2B: 1 distantly related Sps-like'),
    ('Synechococcus sp.', 'CC9311', 'GCF_000014585.1', 'Cyanobacteriota', 0, None, RAV, 'text: confined to actinobacteria'),
    ('Prochlorococcus marinus', 'CCMP1986', 'GCF_000011465.1', 'Cyanobacteriota', 0, None, RAV, 'text: confined to actinobacteria'),
    ('Candidatus Pelagibacter ubique', 'HTCC1062', 'GCF_000012345.1', 'Pseudomonadota', 0, None, RAV, 'text: confined to actinobacteria'),
]
t = pd.DataFrame(L, columns=['organism', 'strain', 'accession', 'phylum', 'expected_rpf', 'expected_sps', 'source', 'evidence'])
t['expected_rpf'] = t.expected_rpf.astype('Int64'); t['expected_sps'] = t.expected_sps.astype('Int64')
t.to_csv('validation/rpf_labels.tsv', sep='\t', index=False)
a = pd.DataFrame({'accession': [('RS_' if x.startswith('GCF') else 'GB_') + x for x in t.accession],
                  'ncbi_assembly_name': '', 'ncbi_translation_table': 11})
a.to_csv('/tmp/val/rpf_accessions.tsv', sep='\t', index=False)
print(t.groupby('phylum').size())
