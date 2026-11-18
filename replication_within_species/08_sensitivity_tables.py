#!/usr/bin/env python3
"""Sensitivity tables for the within-species models (run 04_models.R on each):
  cap20  - at most 20 genomes per species (random draw, seed 1), so a few very large species
           (e.g. Burkholderia ubonensis, n = 250) cannot dominate the within-species slope;
  nosag  - single-amplified genomes removed (a SAG is one cell, not a population consensus).
Species left with < 2 genomes are dropped and mixed_species is recomputed.
usage (repo root): python3 replication_within_species/08_sensitivity_tables.py TABLE.tsv.gz OUTDIR"""
import sys, pandas as pd
t = pd.read_csv(sys.argv[1], sep='\t', low_memory=False); o = sys.argv[2]
def fin(d, tag):
    d = d[d.groupby('species').accession.transform('size') >= 2].copy()
    d['mixed_species'] = d.groupby('species').is_mag.transform(lambda v: int(v.min() == 0 and v.max() == 1))
    d.to_csv(f'{o}/table_{tag}.tsv.gz', sep='\t', index=False, compression='gzip')
    print(tag, len(d), 'genomes', d.species.nunique(), 'species', d[d.mixed_species == 1].species.nunique(), 'mixed')
fin(t.sample(frac=1, random_state=1).groupby('species').head(20).sort_values('accession'), 'cap20')
fin(t[t.is_sag == 0], 'nosag')
