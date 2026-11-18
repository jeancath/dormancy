#!/bin/sh
# gRodon for all analysis genomes in chunks of 300; CDS files are deleted after each chunk (disk-light, restart-safe).
cd /home/user/project-melad
OUT=/home/user/data/hits/grodon.tsv
while true; do
  python3 - <<'PY'
import pandas as pd, os
a=pd.read_csv('data_small/analysis_set.tsv.gz',sep='\t',usecols=['accession'])
done=set(pd.read_csv('/home/user/data/hits/grodon.tsv',sep='\t').accession) if os.path.exists('/home/user/data/hits/grodon.tsv') else set()
todo=a[~a.accession.isin(done)].head(300)
todo.to_csv('/home/user/data/hits/growth_chunk.tsv',sep='\t',index=False); print('chunk',len(todo),'remaining',(~a.accession.isin(done)).sum(),flush=True)
PY
  n=$(($(wc -l < /home/user/data/hits/growth_chunk.tsv)-1))
  [ "$n" -le 0 ] && { echo GROWTHDONE; break; }
  rm -rf /home/user/data/genomes/cds; mkdir -p /home/user/data/genomes/cds
  python3 scripts/06_growth_rate.py /home/user/data/hits/growth_chunk.tsv /home/user/data/genomes/proteins /home/user/data/genomes/cds 4 > /dev/null
  Rscript scripts/06b_grodon.R /home/user/data/genomes/cds $OUT 4 > /dev/null 2>&1
  # genomes whose CDS could not be produced: record so the loop terminates
  python3 - <<'PY'
import pandas as pd, os
c=pd.read_csv('/home/user/data/hits/growth_chunk.tsv',sep='\t').accession
g=pd.read_csv('/home/user/data/hits/grodon.tsv',sep='\t')
miss=[a for a in c if a not in set(g.accession)]
if miss:
    pd.DataFrame({'accession':miss,'d_hours':'NA','lower':'NA','upper':'NA','CUBHE':'NA','ConsistencyHE':'NA','CUB':'NA','n_HE':'NA','err':'no_cds'}).to_csv('/home/user/data/hits/grodon.tsv',sep='\t',mode='a',header=False,index=False)
PY
done
