#!/usr/bin/env python3
# SUPERSEDED by scripts/23c_biosample_reconstruct.py, which rebuilds the log from scratch and overwrites its output.
"""Added 18 Nov 2026: write a committed provenance log for the BioSample collection-date retrieval from the
local caches written by 23_collection_dates.py (which are gitignored).
One row per BioSample record held in the cache:
  biosample            accession of the record as returned by NCBI (records are keyed by the returned accession,
                       so a record is linked to a genome only if this equals the genome's ncbi_biosample)
  requested_by_genome  True if some analysed genome lists this BioSample
  retrieval            'efetch by accession' or 'esearch UID then efetch (accession verified)'
  collection_date_raw  the raw collection_date attribute (empty if the record has none)
  status               'linked' or 'not requested: returned for a misresolved DDBJ/EBI request, discarded'
usage (repo root): python3 scripts/23b_biosample_provenance.py -> data_small/final/biosample_retrieval_log.tsv"""
import json, pandas as pd
how = json.load(open('data_small/final/.collection_dates_how.json')); raw = json.load(open('data_small/final/.collection_dates_raw.json'))
cd = pd.read_csv('data_small/final/collection_dates.tsv', sep='\t')
req = set(cd.ncbi_biosample.dropna())
rows = [dict(biosample=k, requested_by_genome=k in req, retrieval=how.get(k, ''), collection_date_raw=raw.get(k) or '',
             status='linked' if k in req else 'not requested: returned for a misresolved DDBJ/EBI request, discarded')
        for k in sorted(set(how) | set(raw))]
o = pd.DataFrame(rows); o.to_csv('data_small/final/biosample_retrieval_log.tsv', sep='\t', index=False)
print(len(o), o.status.value_counts().to_dict(), o[o.requested_by_genome].retrieval.value_counts().to_dict())
print('requested BioSamples missing from the log:', len(req - set(o.biosample)))
