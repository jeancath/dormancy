#!/usr/bin/env python3
"""Added 18 Nov 2026: rebuild the BioSample retrieval from the first request to the final link to a genome,
WITHOUT using the caches, so that every step can be audited and re-run.

Step 1 (as in the original run): request every BioSample accession by string with efetch, in batches of 200.
  For each returned record keep its accession and NCBI UID. A request is 'matched' if a returned record carries the
  requested accession. A returned record whose accession was not requested is a misresolved return: NCBI read the
  digits of a DDBJ/EBI accession as a UID (e.g. SAMD00000346 -> UID 346 -> SAMN00000346). It is attributed to the
  request whose digits equal its UID, and discarded.
Step 2 (the repair): every unmatched request is resolved to its UID with esearch [accn], fetched by UID, and accepted
  only if the returned record carries the requested accession.
Final link: a genome receives the collection_date of the record whose accession equals the genome's ncbi_biosample.
The raw dates are compared with those of the original run (cache), as a reproducibility check.
No e-mail address, API key or tool name is sent to NCBI; requests are paced at about 1 per second.
usage (repo root): python3 scripts/23c_biosample_reconstruct.py
  -> data_small/final/biosample_retrieval_log.tsv     one row per requested BioSample
     data_small/final/biosample_discarded_returns.tsv one row per misresolved record returned in step 1
     results/biosample_retrieval_flow.tsv             counts at every step"""
import json, os, re, time, urllib.request, xml.etree.ElementTree as ET
import pandas as pd

E = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
def get(url):
    for k in range(5):
        try:
            r = urllib.request.urlopen(url, timeout=90).read(); time.sleep(1.0); return r
        except Exception:
            time.sleep(5 * (k + 1))
    raise RuntimeError('failed: ' + url[:120])
def date_of(bs):
    for a in bs.iter('Attribute'):
        if a.get('harmonized_name') == 'collection_date' or a.get('attribute_name') == 'collection_date':
            return (a.text or '').strip()
    return None

t = pd.read_csv('data_small/final/analysis_table.tsv.gz', sep='\t', usecols=['accession', 'ncbi_biosample', 'pilot'])
req = sorted(set(t.ncbi_biosample.dropna()))
ng = t.groupby('ncbi_biosample').size(); nconf = t[t.pilot == 0].groupby('ncbi_biosample').size()
# ---- step 1
ret, batch_of = {}, {}
for i in range(0, len(req), 200):
    b = req[i:i + 200]
    for acc in b: batch_of[acc] = i // 200 + 1
    root = ET.fromstring(get(E + 'efetch.fcgi?db=biosample&retmode=xml&id=' + ','.join(b)))
    for bs in root.findall('BioSample'):
        ret[bs.get('accession')] = dict(uid=bs.get('id'), date=date_of(bs), batch=i // 200 + 1)
    print('step 1:', i + len(b), 'of', len(req), flush=True)
reqset = set(req)
matched = [a for a in req if a in ret]
wrong = {a: v for a, v in ret.items() if a not in reqset}
unmatched = [a for a in req if a not in ret]
digits = {a: int(re.sub(r'\D', '', a)) for a in unmatched}
by_uid = {}
for a, v in wrong.items(): by_uid.setdefault(int(v['uid']), []).append(a)
# ---- step 2
step2 = {}
for acc in unmatched:
    js = json.loads(get(E + 'esearch.fcgi?db=biosample&retmode=json&term=' + acc + '%5Baccn%5D'))
    uids = js['esearchresult']['idlist']
    r = dict(uids=','.join(uids), verified=False, date=None, returned='')
    if uids:
        root = ET.fromstring(get(E + 'efetch.fcgi?db=biosample&retmode=xml&id=' + ','.join(uids)))
        for bs in root.findall('BioSample'):
            r['returned'] = bs.get('accession')
            if acc in ({bs.get('accession')} | {x.text for x in bs.iter('Id') if x.text}):
                r.update(verified=True, date=date_of(bs))
    step2[acc] = r
    print('step 2:', acc, r['verified'], flush=True)
# ---- log
orig = json.load(open('data_small/final/.collection_dates_raw.json')) if os.path.exists('data_small/final/.collection_dates_raw.json') else {}
rows = []
for acc in req:
    s1 = ret.get(acc); s2 = step2.get(acc, {})
    mis = by_uid.get(digits.get(acc, -1), []) if acc in digits else []
    linked = bool(s1) or s2.get('verified', False)
    date = s1['date'] if s1 else s2.get('date')
    rows.append(dict(biosample=acc, n_genomes=int(ng.get(acc, 0)), n_confirmatory_genomes=int(nconf.get(acc, 0)),
                     step1_batch=batch_of[acc], step1_own_record_returned=bool(s1),
                     step1_misresolved_return=';'.join(mis), step2_esearch_uids=s2.get('uids', ''),
                     step2_returned_accession=s2.get('returned', ''), step2_accession_verified=s2.get('verified', ''),
                     final_status='linked' if linked else 'not retrieved',
                     final_method=('step 1: efetch by accession' if s1 else ('step 2: esearch UID then efetch, verified' if linked else '')),
                     collection_date_raw=date if date is not None else '',
                     same_as_original_run=(orig.get(acc) == date) if acc in orig else ''))
log = pd.DataFrame(rows); log.to_csv('data_small/final/biosample_retrieval_log.tsv', sep='\t', index=False)
att = {a: r for r in unmatched for a in by_uid.get(digits[r], [])}
pd.DataFrame([dict(returned_accession=a, returned_uid=v['uid'], step1_batch=v['batch'], attributed_request=att.get(a, ''),
                   status='discarded: not requested')
              for a, v in sorted(wrong.items())]).to_csv('data_small/final/biosample_discarded_returns.tsv', sep='\t', index=False)
conf = log[log.n_confirmatory_genomes > 0]
flow = [('BioSamples requested (all analysed genomes)', len(req)),
        ('  analysed genomes linked to them', int(ng.sum())),
        ('step 1: records returned', len(ret)),
        ('step 1: returned records matching a request', len(matched)),
        ('step 1: returned records NOT requested (misresolved, discarded)', len(wrong)),
        ('step 1: requests without their own record', len(unmatched)),
        ('  of which a misresolved record was returned for them', sum(1 for a in unmatched if by_uid.get(digits[a]))),
        ('  of which nothing was returned for them', sum(1 for a in unmatched if not by_uid.get(digits[a]))),
        ('step 2: requests resolved by esearch', len(step2)),
        ('step 2: records verified (accession matches)', sum(r['verified'] for r in step2.values())),
        ('final: BioSamples linked', int((log.final_status == 'linked').sum())),
        ('final: BioSamples not retrieved (missing)', int((log.final_status == 'not retrieved').sum())),
        ('final: confirmatory genomes', int(conf.n_confirmatory_genomes.sum())),
        ('final: confirmatory genomes with a linked BioSample', int(conf[conf.final_status == 'linked'].n_confirmatory_genomes.sum())),
        ('check: raw date identical to the original run (BioSamples)', int((log.same_as_original_run == True).sum())),
        ('check: raw date differs from the original run (BioSamples)', int((log.same_as_original_run == False).sum()))]
f = pd.DataFrame(flow, columns=['step', 'count']); f.to_csv('results/biosample_retrieval_flow.tsv', sep='\t', index=False)
print(f.to_string(index=False))
