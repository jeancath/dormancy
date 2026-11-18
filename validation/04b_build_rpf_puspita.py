#!/usr/bin/env python3
"""Second Rpf tier: Puspita et al. 2013, Microbes Environ 28:58-64, doi:10.1264/jsme2.me12122,
Table S1 ("putative number of rpf gene copies found in Actinobacteria finished and draft genome
sequences"). Counts there come from BLASTp with the Rpf domain, i.e. a method independent of the
Pfam PF06737 HMM used by our dictionary, but still bioinformatic, not experimental.

Only finished genomes (RefSeq NC_/CP accession given) are used. Each chromosome accession is
linked to its assembly through NCBI E-utilities (elink nuccore->assembly). Rows where the
linked assembly's organism does not carry the table's strain designation are dropped: the PDF
table has column misalignments (e.g. a C. pseudotuberculosis row pointing to the R. jostii
chromosome). Genomes already in validation/rpf_labels.tsv are left out here.

usage: python3 validation/04b_build_rpf_puspita.py S1_PDF_TEXT.txt  -> validation/rpf_labels_puspita.tsv
  (S1_PDF_TEXT.txt = `pdftotext -layout 28_58_s1.pdf`, supplement from Europe PMC PMC4070687)
"""
import sys, re, time
import pandas as pd
import requests

E = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
rows = []
for line in open(sys.argv[1], encoding='utf-8', newline='').read().split('\n'):
    t = re.sub(r'[\t\r ]', ' ', line).replace('­', '').replace('‐', '-')
    m = re.match(r'\s*(.+?)\s+(\w+aceae)\s+(\d{6,8})\s+(\d+)\s+(.*?)\s+(NC_\d+|CP\d+|NZ_\w+|Draft)\s', t)
    if m:
        rows.append(dict(name=re.sub(r'\s+', ' ', m.group(1)).strip(), family=m.group(2), genome_size=int(m.group(3)),
                         expected_rpf=int(m.group(4)), habitat=re.sub(r'\s+', ' ', m.group(5)), nuc=m.group(6)))
d = pd.DataFrame(rows)
f = d[d.nuc != 'Draft'].copy()
print(f'{len(d)} rows parsed, {len(f)} finished genomes')

asm, org = [], []
for nuc in f.nuc:
    a = o = None
    for k in range(3):
        try:
            r = requests.get(E + 'elink.fcgi', params=dict(dbfrom='nuccore', db='assembly', id=nuc, retmode='json'), timeout=60).json()
            ids = [l for ls in r['linksets'] for db in ls.get('linksetdbs', []) for l in db['links']]
            if ids:
                s = requests.get(E + 'esummary.fcgi', params=dict(db='assembly', id=','.join(ids[:5]), retmode='json'), timeout=60).json()['result']
                c = [(s[i].get('assemblyaccession'), s[i].get('organism')) for i in s.get('uids', [])]
                c = [x for x in c if x[0] and x[0].startswith('GCF')] or c
                if c:
                    a, o = c[0]
            break
        except Exception:
            time.sleep(2 ** k)
    asm.append(a); org.append(o); time.sleep(0.35)
f['accession'], f['asm_org'] = asm, org


def norm(s):
    return re.sub(r'[^a-z0-9]', '', str(s).lower())


def consistent(name, o):
    if not isinstance(o, str):
        return False
    toks = [t for t in re.split(r'[\s,]+', name.replace('--', '-')) if t]
    strain = [t for t in toks[2:] if re.search(r'\d', t)] or toks[2:]
    sp = toks[1] if len(toks) > 1 else ''
    on = norm(o)
    return any(norm(t) in on for t in strain if len(norm(t)) >= 2) and (norm(sp) in on or toks[0] in o)


f['consistent'] = [consistent(n, o) for n, o in zip(f.name, f.asm_org)]
f = f[f.consistent]
known = pd.read_csv('validation/rpf_labels.tsv', sep='\t')
overlap = f[f.accession.isin(known.accession)].merge(known[['accession', 'expected_rpf']], on='accession', suffixes=('_puspita', '_tier1'))
print('overlap with tier 1 (count agreement):'); print(overlap[['name', 'expected_rpf_puspita', 'expected_rpf_tier1']].to_string())
f = f[~f.accession.isin(known.accession)].drop_duplicates('accession')
f = f.assign(organism=f.name, phylum='Actinomycetota', source='Puspita2013 Microbes Environ 28:58 doi:10.1264/jsme2.me12122',
             evidence=[f'Table S1: {n} rpf copies (BLASTp); chromosome {c}' for n, c in zip(f.expected_rpf, f.nuc)])
f[['organism', 'accession', 'phylum', 'family', 'expected_rpf', 'habitat', 'source', 'evidence', 'asm_org']].to_csv(
    'validation/rpf_labels_puspita.tsv', sep='\t', index=False)
overlap.to_csv('validation/rpf_puspita_overlap_tier1.tsv', sep='\t', index=False)
print(len(f), 'new genomes;', f.expected_rpf.value_counts().sort_index().to_dict())
