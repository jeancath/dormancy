#!/bin/sh
# usage: run_batch.sh NAME  -> annotates confirmatory+pilot genomes with proteins on disk not yet in any batch
cd /home/user/project-melad
N=$1
python3 - "$N" <<'PY'
import pandas as pd, os, glob, sys
a=pd.read_csv('data_small/analysis_set.tsv.gz',sep='\t',usecols=['accession','pilot'])
a=a[a.pilot==0]
done=set()
for f in glob.glob('/home/user/data/hits/batch*.tsv'):
    if not f.endswith('_nprot.tsv'): done|=set(pd.read_csv(f,sep='\t').accession)
ok=a.accession.map(lambda x: x not in done and os.path.exists(f'/home/user/data/genomes/proteins/{x[:10]}/{x}.faa.gz'))
a[ok][['accession']].to_csv(f'/home/user/data/hits/{sys.argv[1]}.tsv',sep='\t',index=False); print(sys.argv[1], ok.sum())
PY
python3 scripts/05_hmm_search.py /home/user/data/hmm/dormancy_library.hmm dictionary/thresholds.tsv /home/user/data/genomes/proteins /home/user/data/hits/$N.tsv /home/user/data/hits/${N}_hits.tsv.gz 2
echo ${N}DONE
