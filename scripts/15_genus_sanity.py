#!/usr/bin/env python3
"""Genus-level sanity check of module calls against textbook phenotypes (not a calibration;
see validation/ for the phenotype-labelled calibration).
usage: 15_genus_sanity.py TABLE.tsv.gz OUT.tsv"""
import sys, pandas as pd
d = pd.read_csv(sys.argv[1], sep='\t', low_memory=False)
d['genus'] = d.gtdb_taxonomy.str.extract(r'g__([^;]+)')[0]
d['family'] = d.gtdb_taxonomy.str.extract(r'f__([^;]+)')[0]
# textbook expectation: 1 = endospore-forming genus, 0 = non-sporulating genus
exp = {'Bacillus': 1, 'Paenibacillus': 1, 'Clostridium': 1, 'Peribacillus': 1, 'Priestia': 1, 'Brevibacillus': 1,
       'Lysinibacillus': 1, 'Neobacillus': 1, 'Cytobacillus': 1, 'Alkalihalobacillus': 1, 'Virgibacillus': 1,
       'Sporosarcina': 1, 'Desulfotomaculum': 1, 'Thermoactinomyces': 1, 'Halobacillus': 1, 'Oceanobacillus': 1,
       'Staphylococcus': 0, 'Streptococcus': 0, 'Lactobacillus': 0, 'Enterococcus': 0, 'Listeria': 0,
       'Lactococcus': 0, 'Exiguobacterium': 0, 'Planococcus': 0, 'Carnobacterium': 0, 'Macrococcus': 0,
       'Weissella': 0, 'Leuconostoc': 0, 'Mammaliicoccus': 0}
rows = []
for g, e in exp.items():
    s = d[d.genus.str.fullmatch(g + r'(_[A-Z])?', na=False)]
    if len(s):
        rows.append(dict(genus=g, expected=e, n=len(s), pct_sporulator=round(100 * s.sporulator.mean(), 1),
                         median_spo=s.spo_strict.median(), pct_spo0A=round(100 * s.spo0A.mean(), 1)))
r = pd.DataFrame(rows)
print(r.to_string(index=False))
out = [r]
# phylum-level sporulator counts and Rpf by phylum
ph = d.groupby('phylum').agg(n=('accession', 'size'), sporulators=('sporulator', 'sum'), spo0A=('spo0A', 'sum'),
                             rpf_pct=('rpf_any', 'mean'), hib_mean=('n_hib', 'mean'), sr_mean=('n_sr', 'mean'))
ph = ph[ph.n >= 30].sort_values('n', ascending=False); ph['rpf_pct'] = (100 * ph.rpf_pct).round(1)
print(ph.round(2).to_string())
r.to_csv(sys.argv[2], sep='\t', index=False)
ph.round(3).to_csv(sys.argv[2].replace('.tsv', '_phylum.tsv'), sep='\t')
