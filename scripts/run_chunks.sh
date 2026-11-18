#!/bin/sh
# Annotate all remaining confirmatory genomes in chunks of 500 (restart-safe: finished chunks are kept).
cd /home/user/project-melad
i=2
while true; do
  N=batch$i
  if [ -s /home/user/data/hits/${N}_hits.tsv.gz ]; then i=$((i+1)); continue; fi
  python3 - "$N" <<'PY'
import pandas as pd, os, glob, sys
a=pd.read_csv('data_small/analysis_set.tsv.gz',sep='\t',usecols=['accession','pilot'])
a=a[a.pilot==0]
done=set()
for f in glob.glob('/home/user/data/hits/batch*_hits_nprot.tsv'):
    done|=set(pd.read_csv(f,sep='\t').accession)
ok=a.accession.map(lambda x: x not in done and os.path.exists(f'/home/user/data/genomes/proteins/{x[:10]}/{x}.faa.gz'))
a[ok][['accession']].head(500).to_csv(f'/home/user/data/hits/{sys.argv[1]}.tsv',sep='\t',index=False); print(sys.argv[1], min(ok.sum(),500), 'remaining', ok.sum())
PY
  n=$(($(wc -l < /home/user/data/hits/$N.tsv)-1))
  [ "$n" -le 0 ] && { rm /home/user/data/hits/$N.tsv; echo ALLDONE; break; }
  python3 scripts/05_hmm_search.py /home/user/data/hmm/dormancy_library.hmm dictionary/thresholds.tsv /home/user/data/genomes/proteins /home/user/data/hits/$N.tsv /home/user/data/hits/${N}_hits.tsv.gz 4 > /dev/null
  i=$((i+1))
done
