#!/bin/bash
# Combine per-chunk annotation outputs and run table -> models -> figure.
# usage (repo root): bash replication_within_species/07_run_downstream.sh "HITDIRS" SELECTED.tsv.gz OUTDIR LABEL
#   HITDIRS: space-separated list of hit directories (e.g. "/data/rep/hits /data/rep/hits_mq")
set -eu
D=$1; SEL=$2; O=$3; LAB=$4; R=replication_within_species
mkdir -p $O
python3 - "$D" "$O" <<'PY'
import sys, glob, pandas as pd
dirs, o = sys.argv[1].split(), sys.argv[2]
def cat(pat): return pd.concat([pd.read_csv(f, sep='\t') for d in dirs for f in sorted(glob.glob(f'{d}/{pat}'))])
h = cat('c*.tsv.gz'); n = cat('c*_nprot.tsv'); e = cat('c*.he.tsv')
h.to_csv(f'{o}/hits.tsv.gz', sep='\t', index=False)
n.drop_duplicates('accession').to_csv(f'{o}/nprot.tsv', sep='\t', index=False)
e.drop_duplicates('accession').to_csv(f'{o}/he.tsv', sep='\t', index=False)
print('genomes annotated:', n.accession.nunique(), 'with HE:', e.accession.nunique())
PY
python3 - "$SEL" "$O" <<'PY'
import sys, pandas as pd
d = pd.read_csv(sys.argv[1], sep='\t')
d['phylum'] = d.gtdb_taxonomy.str.extract(r'p__([^;]+)')[0]
d[['accession', 'phylum']].to_csv(f'{sys.argv[2]}/meta.tsv', sep='\t', index=False)
PY
python3 scripts/08_aggregate_modules.py $O/hits.tsv.gz $O/nprot.tsv $O/meta.tsv $O/modules.tsv
python3 $R/03_build_table.py $O/modules.tsv $O/he.tsv $SEL $O/table.tsv.gz
Rscript $R/04_models.R $O/table.tsv.gz $O/$LAB $LAB
python3 $R/06_figure.py $O/${LAB}_results.tsv $O/${LAB}_mixed_pairs.tsv $O/fig_$LAB.png
