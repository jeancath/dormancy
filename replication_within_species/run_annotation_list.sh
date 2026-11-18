#!/bin/bash
# Same as run_annotation.sh but for an arbitrary accession list (used for the medium-quality
# extension). usage (repo root): bash replication_within_species/run_annotation_list.sh LIST.tsv WORKDIR TAG NPROC
set -u
L=$1; W=$2; T=$3; NP=${4:-4}; R=replication_within_species
mkdir -p $W/chunks_$T $W/hits_$T
python3 -c "
import pandas as pd; d=pd.read_csv('$L',sep='\t').sort_values('accession')
for i in range(0,len(d),500): d.iloc[i:i+500].to_csv('$W/chunks_$T/c%02d.tsv'%(i//500),sep='\t',index=False)"
for c in $W/chunks_$T/c*.tsv; do
  b=$(basename $c .tsv)
  [ -s $W/hits_$T/$b.he.tsv ] && { echo "$b done"; continue; }
  echo "$(date +%T) $b fetch"; python3 scripts/03_fetch_and_call_genes.py $c $W/prot $NP > /dev/null
  echo "$(date +%T) $b hmm";   python3 scripts/05_hmm_search.py $W/hmm/dormancy_library.hmm dictionary/thresholds.tsv $W/prot $c $W/hits_$T/$b.tsv.gz $NP > /dev/null
  echo "$(date +%T) $b he";    python3 $R/02_count_he.py $W/hmm/ribosomal.hmm $W/prot $c $W/hits_$T/$b.he.tsv $NP
  tail -n +2 $c | cut -f1 | while read a; do rm -f $W/prot_fna/${a:0:10}/$a.fna.gz; done
done
echo "$(date +%T) ALL DONE"
