#!/usr/bin/env python3
"""Fetch BioSample collection dates (cached in data_small/final/.collection_dates_raw.json after the first run) for the analysis genomes (NCBI E-utilities efetch, batches of 200,
<= 3 requests/s), parse the collection year, and write data_small/final/collection_dates.tsv.
Missing or unparseable dates (e.g. 'missing', 'not collected') are recorded as NA.
usage (repo root): python3 scripts/23_collection_dates.py"""
import re, time, urllib.request, pandas as pd, xml.etree.ElementTree as ET

t = pd.read_csv('data_small/final/analysis_table.tsv.gz', sep='\t', usecols=['accession', 'ncbi_biosample', 'pilot'])
ids = sorted(set(t.ncbi_biosample.dropna()))
import json, os
CACHE = 'data_small/final/.collection_dates_raw.json'
rec = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
for i in (range(0, len(ids), 200) if not rec else []):
    batch = ids[i:i + 200]
    url = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&retmode=xml&id=' + ','.join(batch))
    for attempt in range(4):
        try:
            xml = urllib.request.urlopen(url, timeout=60).read()
            break
        except Exception as e:
            time.sleep(3 * (attempt + 1))
    else:
        print('failed batch', i); continue
    root = ET.fromstring(xml)
    for bs in root.findall('BioSample'):
        acc = bs.get('accession')
        val = None
        for a in bs.iter('Attribute'):
            if a.get('harmonized_name') == 'collection_date' or a.get('attribute_name') == 'collection_date':
                val = (a.text or '').strip(); break
        rec[acc] = val
    time.sleep(0.4)
    print(i + len(batch), 'of', len(ids), flush=True)

import json
# Repair step (added 18 Nov 2026). efetch expects numeric BioSample UIDs; for some DDBJ (SAMD...) and EBI (SAMEA...)
# accessions, passing the accession string made NCBI interpret its digits as a UID and return a DIFFERENT record
# (e.g. SAMD00000346 -> SAMN00000346). Records are matched by their own accession, so no wrong date was ever
# assigned, but those genomes were left undated. Here each unmatched accession is resolved to its UID with esearch
# [accn], fetched by UID, and accepted only if the returned record carries the requested accession.
how = json.load(open('data_small/final/.collection_dates_how.json')) if os.path.exists('data_small/final/.collection_dates_how.json') else {}
for acc in rec:
    how.setdefault(acc, 'efetch by accession')
for acc in [x for x in ids if x not in rec]:
    try:
        js = json.loads(urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=biosample&retmode=json&term='
                                               + acc + '%5Baccn%5D', timeout=60).read())
        uids = js["esearchresult"]["idlist"]; time.sleep(1.0)
        if not uids:
            how[acc] = 'not found by esearch'; continue
        root = ET.fromstring(urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&retmode=xml&id='
                                                    + ','.join(uids), timeout=60).read()); time.sleep(1.0)
        for bs in root.findall('BioSample'):
            ids_in_record = {bs.get('accession')} | {i.text for i in bs.iter('Id') if i.text}
            if acc in ids_in_record:
                val = None
                for at in bs.iter('Attribute'):
                    if at.get('harmonized_name') == 'collection_date' or at.get('attribute_name') == 'collection_date':
                        val = (at.text or '').strip(); break
                rec[acc] = val; how[acc] = 'esearch UID then efetch (accession verified)'
        how.setdefault(acc, 'UID record did not match accession')
    except Exception as e:
        how[acc] = 'error: ' + str(e)[:60]
json.dump(rec, open('data_small/final/.collection_dates_raw.json', 'w'))
json.dump(how, open('data_small/final/.collection_dates_how.json', 'w'))
YR = r'(?<!\d)(1[89]\d\d|20\d\d)(?!\d)'
def parse(v):
    """Return (format, start_year, end_year). INSDC collection_date allows a full date, a year-month, a year,
    or a range written 'date1/date2'. Text such as 'missing' or 'not collected' is undated."""
    if not isinstance(v, str) or not v.strip():
        return 'absent', None, None
    ys = [int(y) for y in re.findall(YR, v)]
    if not ys:
        return 'non-date text', None, None
    if '/' in v and len(set(ys)) > 1:
        return 'range', min(ys), max(ys)
    if re.fullmatch(r'\d{4}', v.strip()):
        return 'year only', ys[0], ys[0]
    if re.fullmatch(r'\d{4}-\d{2}', v.strip()) or re.fullmatch(r'[A-Za-z]{3}-\d{4}', v.strip()):
        return 'year-month', ys[0], ys[0]
    return 'full date', ys[-1], ys[-1]
out = pd.DataFrame({'ncbi_biosample': list(rec), 'collection_date': list(rec.values())})
howmap = json.load(open('data_small/final/.collection_dates_how.json')) if os.path.exists('data_small/final/.collection_dates_how.json') else {}
out['retrieval'] = out.ncbi_biosample.map(howmap).fillna('efetch by accession')
parsed = out.collection_date.map(parse)
out['date_format'] = [x[0] for x in parsed]
out['year_start'] = [x[1] for x in parsed]
out['year_end'] = [x[2] for x in parsed]
out['collection_year'] = out.year_end
# in the 1991-2020 climate window only if the WHOLE stated period lies inside it
out['in_window'] = ((out.year_start >= 1991) & (out.year_end <= 2020)).astype(int)
res = t.merge(out, on='ncbi_biosample', how='left')
res['retrieval'] = res.retrieval.fillna('BioSample record not retrieved')
res['date_format'] = res.date_format.fillna('BioSample record not retrieved')
res['in_window'] = res.in_window.fillna(0).astype(int)
res.drop(columns='pilot').to_csv('data_small/final/collection_dates.tsv', sep='\t', index=False)
c = res[res.pilot == 0]
dated = c.collection_year.notna()
print('confirmatory genomes', len(c), '| dated', dated.sum(), f'({dated.mean():.1%})',
      '| in window', c.in_window.sum(), f'= {c.in_window.sum() / len(c):.1%} of all, {c.in_window.sum() / dated.sum():.1%} of dated')
print(c.retrieval.value_counts().to_string()); print(c.date_format.value_counts().to_string())
