#!/bin/bash
# Restart-safe chunked annotation driver (step 2). Calls the repo scripts unchanged.
# usage (repo root): bash replication_within_species/run_annotation.sh /data/rep 4
set -u
W=$1; NP=${2:-4}; R=replication_within_species
mkdir -p $W/chunks $W/hits
python3 - "$W" <<'PY'
import sys, pandas as pd
W = sys.argv[1]
d = pd.read_csv(f'{W}/selected_genomes.tsv.gz', sep='\t')
# translation table: GTDB metadata is not in the repo; 11 everywhere except Gracilibacteria (25),
# matching the main analysis, where the only non-11 tables were 5 Gracilibacteria genomes (25).
d['ncbi_translation_table'] = d.gtdb_taxonomy.str.contains('c__Gracilibacteria').map({True: 25, False: 11})
d['ncbi_assembly_name'] = ''
d = d.sort_values('accession')
for i in range(0, len(d), 500):
    d.iloc[i:i + 500][['accession', 'ncbi_assembly_name', 'ncbi_translation_table']].to_csv(
        f'{W}/chunks/c{i // 500:02d}.tsv', sep='\t', index=False)
PY
for c in $W/chunks/c*.tsv; do
  b=$(basename $c .tsv)
  [ -s $W/hits/$b.he.tsv ] && { echo "$b done"; continue; }
  echo "$(date +%T) $b fetch"; python3 scripts/03_fetch_and_call_genes.py $c $W/prot $NP > /dev/null
  echo "$(date +%T) $b hmm";   python3 scripts/05_hmm_search.py $W/hmm/dormancy_library.hmm dictionary/thresholds.tsv $W/prot $c $W/hits/$b.tsv.gz $NP > /dev/null
  echo "$(date +%T) $b he";    python3 $R/02_count_he.py $W/hmm/ribosomal.hmm $W/prot $c $W/hits/$b.he.tsv $NP
  tail -n +2 $c | cut -f1 | while read a; do rm -f $W/prot_fna/${a:0:10}/$a.fna.gz; done
done
echo "$(date +%T) ALL DONE"
