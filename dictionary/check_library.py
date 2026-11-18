#!/usr/bin/env python3
"""Sanity check: run the dormancy library on reference proteomes using thresholds.tsv."""
import csv, sys, collections
import pyhmmer

LIB = '/home/user/data/hmm/dormancy_library.hmm'
D = '/home/user/project-melad/dictionary'
thr = {r['model_name']: r for r in csv.DictReader(open(f'{D}/thresholds.tsv'), delimiter='\t')}
dic = list(csv.DictReader(open(f'{D}/dormancy_dictionary.tsv'), delimiter='\t'))
core = list(csv.DictReader(open(f'{D}/sporulation_core.tsv'), delimiter='\t'))

with pyhmmer.plan7.HMMFile(LIB) as fh:
    hmms = list(fh)
assert len(hmms) == len(thr), (len(hmms), len(thr))
for h in hmms:
    n = h.name if isinstance(h.name, str) else h.name.decode()
    assert n in thr, n

alpha = pyhmmer.easel.Alphabet.amino()
for fa in sys.argv[1:]:
    with pyhmmer.easel.SequenceFile(fa, digital=True, alphabet=alpha) as sf:
        seqs = sf.read_block()
    spor_models = {r['model_name'] for r in dic if r['module'].startswith('sporulation')}
    best = {}  # protein -> (score, model) best sporulation model (competitive assignment)
    passing = collections.defaultdict(set)  # model -> proteins (strict)
    probable = collections.defaultdict(set)  # model -> proteins with NC <= score < TC
    for top in pyhmmer.hmmsearch(hmms, seqs, E=10, cpus=4):
        q = top.query.name if isinstance(top.query.name, str) else top.query.name.decode()
        t = thr[q]
        for hit in top:
            nm = hit.name if isinstance(hit.name, str) else hit.name.decode()
            if q in spor_models and (nm not in best or hit.score > best[nm][0]):
                best[nm] = (hit.score, q)
            if t['score_type'] == 'full':
                ok = hit.score >= float(t['threshold'])
                if ok and t['domain_threshold']:
                    ok = any(d.score >= float(t['domain_threshold']) for d in hit.domains)
            else:
                ok = max(d.score for d in hit.domains) >= float(t['threshold'])
            if ok:
                passing[q].add(nm)
            elif t['noise_cutoff'] and hit.score >= float(t['noise_cutoff']):
                probable[q].add(nm)
    # competitive assignment BY GENE: a protein counts only for the gene of its best-scoring sporulation model
    m2g = {r['model_name']: r['gene'] for r in dic if r['module'].startswith('sporulation')}
    for m in spor_models:
        passing[m] = {p for p in passing.get(m, set()) if m2g[best[p][1]] == m2g[m]}
        probable[m] = {p for p in probable.get(m, set()) if m2g[best[p][1]] == m2g[m]}
    gene_models = collections.defaultdict(set)
    for r in dic:
        gene_models[(r['module'], r['gene'])].add(r['model_name'])
    present = {k for k, ms in gene_models.items() if any(passing.get(m) for m in ms)}
    scored = [c['gene'] for c in core if c['in_score'] == 'yes']
    sc = sum(1 for g in scored if ('sporulation_core', g) in present)
    missing = [g for g in scored if ('sporulation_core', g) not in present]
    rel = [g for g in missing if any(probable.get(m) for m in gene_models[('sporulation_core', g)])]
    spo0a = ('sporulation_core', 'spo0A') in present
    print(f'\n== {fa.split("/")[-1]}  proteins={len(seqs)}  spo0A={spo0a} core_strict={sc}/{len(scored)} core_relaxed={sc+len(rel)}  missing={missing} probable(NC..TC)={rel}')
    bymod = collections.defaultdict(list)
    for (m, g) in sorted(present):
        if m != 'sporulation_core':
            n = len(set().union(*[passing.get(x, set()) for x in gene_models[(m, g)]]))
            bymod[m].append(f'{g}({n})')
    for m, gs in bymod.items():
        print(f'  {m}: {", ".join(gs)}')
