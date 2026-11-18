#!/usr/bin/env python3
"""Ribosomal-protein HE gene count per genome, replicating scripts/06_growth_rate.py exactly
(Pfam GA cutoffs on the ribosomal library; a gene counts only if its CDS is >= 240 nt and a
multiple of 3, i.e. the CDS that 06 writes and gRodon's n_HE counts). Coordinates come from the
protein headers, so the nucleotide assembly is not needed. gRodon itself is not run: only n_HE
is needed (for the pre-registered QC rule n_HE < 20).
usage: 02_count_he.py RIBO.hmm PROT_DIR ACCESSION_LIST OUT.tsv [nproc]"""
import sys, os, gzip, multiprocessing as mp, pandas as pd, pyhmmer
from pyhmmer.easel import TextSequence, Alphabet
A = Alphabet.amino(); H = None
def init(lib):
    global H
    with pyhmmer.plan7.HMMFile(lib) as fh: H = list(fh)
def work(path):
    acc = os.path.basename(path)[:-7]; prots = {}; name = None; buf = []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name: prots[name] = ''.join(buf)
                name = line[1:].split()[0]; buf = []
            else: buf.append(line.strip())
    if name: prots[name] = ''.join(buf)
    seqs = [TextSequence(name=k.encode(), sequence=v.replace('*', '')).digitize(A) for k, v in prots.items()]
    he = set()
    for hits in pyhmmer.hmmsearch(H, seqs, cpus=1, bit_cutoffs='gathering'):
        for h in hits.reported: he.add(h.name.decode() if isinstance(h.name, bytes) else h.name)
    n = 0
    for k in he:
        s, e = int(k.split('|')[3]), int(k.split('|')[4]); L = e - s + 1
        n += int(L >= 240 and L % 3 == 0)
    return acc, n, len(he)
if __name__ == '__main__':
    lib, pdir, lst, out = sys.argv[1:5]; npr = int(sys.argv[5]) if len(sys.argv) > 5 else 4
    accs = pd.read_csv(lst, sep='\t').accession
    paths = [p for p in (os.path.join(pdir, a[:10], f'{a}.faa.gz') for a in accs) if os.path.exists(p)]
    with mp.Pool(npr, initializer=init, initargs=(lib,)) as pool:
        res = list(pool.imap_unordered(work, paths, chunksize=4))
    pd.DataFrame(res, columns=['accession', 'n_HE', 'n_ribo_hits_any_len']).to_csv(out, sep='\t', index=False)
