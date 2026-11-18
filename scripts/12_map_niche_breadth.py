#!/usr/bin/env python3
"""Secondary predictor: Hao et al. 2026 PNAS (Dataset S5) habitat niche breadth (number of habitat
clusters) per GTDB r207 species, mapped to r232 species: the r207 representative genome accession is
looked up in r232 metadata and assigned to the r232 species cluster it belongs to.
Also von Meijenfeldt 2023 social niche breadth at genus level (NCBI taxonomy names) matched to the
GTDB r232 genus name with suffixes stripped (approximate; flagged as such)."""
import pandas as pd, numpy as np
h = pd.read_csv('/home/user/data/lit/hao2026/hao2026_datasetS5_species_niche_breadth.tsv', sep='\t')
m = pd.read_csv('/home/user/data/gtdb/bac120_metadata.tsv.gz', sep='\t', usecols=['accession', 'gtdb_genome_representative'], low_memory=False)
h = h.merge(m, left_on='Species_gtdb_genome_representative', right_on='accession', how='left')
print('Hao species:', len(h), 'mapped to r232:', h.gtdb_genome_representative.notna().sum())
nb = h.dropna(subset=['gtdb_genome_representative']).groupby('gtdb_genome_representative').Niche_breadth_number_of_cluster.max()
a = pd.read_csv('data_small/analysis_set.tsv.gz', sep='\t', usecols=['accession', 'gtdb_genome_representative', 'gtdb_taxonomy'])
a['hao_niche_breadth'] = a.gtdb_genome_representative.map(nb)
# von Meijenfeldt SNB, genus level
s = pd.read_csv('/home/user/data/lit/vonmeijenfeldt2023/2023-02-01.Supplementary_Data_File_3.SNB.txt', sep='\t')
print(s.columns.tolist()[:12]); print(s.head(3).to_string()[:800])
a.to_csv('data_small/niche_breadth_map.tsv.gz', sep='\t', index=False, compression='gzip')
print('analysis species with Hao niche breadth:', a.hao_niche_breadth.notna().sum(), 'of', len(a))
# genus-level SNB: NCBI genus name vs GTDB genus with alphabetic suffix removed (e.g. Bacillus_A -> Bacillus)
g = s[(s['rank'] == 'genus') & s.SNB.notna()].copy()
g['genus'] = g['taxonomic lineage'].str.extract(r'genus\.([^;]+)$')[0].str.replace('Candidatus_', '', regex=False)
snb = g.groupby('genus').SNB.mean()
a['gtdb_genus'] = a.gtdb_taxonomy.str.extract(r'g__([^;]+)')[0].str.replace(r'_[A-Z]+$', '', regex=True)
a['snb_genus'] = a.gtdb_genus.map(snb)
a.drop(columns=['gtdb_taxonomy']).to_csv('data_small/niche_breadth_map.tsv.gz', sep='\t', index=False, compression='gzip')
print('analysis species with genus SNB:', a.snb_genus.notna().sum())
