#!/usr/bin/env python3
"""Phase 5b: search the dormancy HMM library against stored proteomes with pyhmmer.

Cutoff logic per model (in order):
  1. thresholds.tsv row (model_name, threshold, score_type in {full, domain}) - KOfam style
  2. Pfam GA / NCBIfam TC lines in the HMM file (sequence and domain cutoffs)
  3. fallback: full-sequence E-value <= 1e-5 (flagged in output)

usage: 05_hmm_search.py LIB.hmm THRESHOLDS.tsv PROTEIN_DIR ACCESSION_LIST OUT.tsv.gz [nproc]
Output: one row per (protein, model) with E <= 1e-3; 'pass' = strict curated cutoff (GA/TC incl.
domain cutoff, or KOfam threshold), 'pass_nc' = score >= noise cutoff (relaxed sporulation score).
All hits are kept so that competitive best-hit assignment can be done downstream.
"""
import sys, os, gzip, glob
import multiprocessing as mp
import pandas as pd
import pyhmmer
from pyhmmer.easel import TextSequence, Alphabet

ALPHA = Alphabet.amino()
HMMS = None
THR = None


def init(lib, thr_path):
    global HMMS, THR
    with pyhmmer.plan7.HMMFile(lib) as fh:
        HMMS = list(fh)
    THR = {}
    if thr_path and os.path.exists(thr_path):
        t = pd.read_csv(thr_path, sep='\t')
        for r in t.itertuples():
            f = lambda v: float(v) if str(v) not in ('', 'nan', 'None') else None
            THR[str(r.model_name)] = (f(r.threshold), str(r.score_type), f(getattr(r, 'domain_threshold', None)),
                                      f(getattr(r, 'noise_cutoff', None)))


def name_of(h):
    n = h.name
    return n.decode() if isinstance(n, bytes) else n


def passes(hmm, hit):
    nm = name_of(hmm)
    best_dom = max((d.score for d in hit.domains), default=-1e9)
    if nm in THR and THR[nm][0] is not None:
        thr, st, dthr, nc = THR[nm]
        if st == 'full':
            ok = hit.score >= thr and (dthr is None or best_dom >= dthr)
        else:
            ok = best_dom >= thr
        return ok, ('NC' if (nc is not None and hit.score >= nc) else '')
    for cut, label in ((hmm.cutoffs.gathering, 'GA'), (hmm.cutoffs.trusted, 'TC')):
        if cut is not None and cut[0] is not None:
            return (hit.score >= cut[0] and best_dom >= cut[1]), label
    return hit.evalue <= 1e-5, 'evalue'


def read_seqs(path):
    seqs = []
    with gzip.open(path, 'rt') as fh:
        name, buf = None, []
        for line in fh:
            if line.startswith('>'):
                if name:
                    seqs.append(TextSequence(name=name.encode(), sequence=''.join(buf)).digitize(ALPHA))
                name = line[1:].strip(); buf = []
            else:
                buf.append(line.strip().replace('*', ''))
        if name:
            seqs.append(TextSequence(name=name.encode(), sequence=''.join(buf)).digitize(ALPHA))
    return seqs


def work(path):
    acc = os.path.basename(path).replace('.faa.gz', '')
    try:
        seqs = read_seqs(path)
    except Exception as e:
        return acc, 0, [], str(e)
    rows = []
    for hits in pyhmmer.hmmsearch(HMMS, seqs, cpus=1, E=1e-3, domE=1e-3):
        hmm = hits.query
        for hit in hits:
            ok, how = passes(hmm, hit)
            nm = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            p = nm.split('|')
            bd = max(hit.domains, key=lambda d: d.score)
            rows.append((acc, p[1], int(p[2]), int(p[3]), int(p[4]), int(p[5]), int(p[6]), name_of(hmm),
                         (hmm.accession.decode() if isinstance(hmm.accession, bytes) else hmm.accession) or '',
                         round(hit.score, 1), hit.evalue, round(bd.score, 1), bd.env_from, bd.env_to, int(ok), int(how == 'NC' or ok)))
    return acc, len(seqs), rows, ''


def main():
    lib, thr, pdir, acclist, out = sys.argv[1:6]
    nproc = int(sys.argv[6]) if len(sys.argv) > 6 else 4
    accs = pd.read_csv(acclist, sep='\t').accession.tolist()
    paths = [os.path.join(pdir, a[:10], f'{a}.faa.gz') for a in accs]
    paths = [p for p in paths if os.path.exists(p)]
    print(f'{len(paths)} proteomes', flush=True)
    allrows, nprot = [], []
    with mp.Pool(nproc, initializer=init, initargs=(lib, thr)) as pool:
        for i, (acc, n, rows, err) in enumerate(pool.imap_unordered(work, paths, chunksize=4), 1):
            allrows.extend(rows); nprot.append((acc, n, err))
            if i % 200 == 0:
                print(i, flush=True)
    cols = ['accession', 'contig', 'gene_idx', 'start', 'end', 'strand', 'plasmid', 'model', 'model_acc',
            'score', 'evalue', 'dom_score', 'env_from', 'env_to', 'pass', 'pass_nc']
    pd.DataFrame(allrows, columns=cols).to_csv(out, sep='\t', index=False, compression='gzip')
    pd.DataFrame(nprot, columns=['accession', 'n_proteins', 'error']).to_csv(
        out.replace('.tsv.gz', '_nprot.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
