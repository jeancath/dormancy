#!/usr/bin/env python3
"""Growth-rate covariate: predicted minimal doubling time from codon usage bias (gRodon2,
Weissman et al. 2021 PNAS), using ribosomal proteins (Pfam GA cutoffs) as highly expressed genes.

Step 1 (this script): per genome, find ribosomal-protein genes with pyhmmer, extract all CDS
nucleotide sequences from the stored assembly using gene coordinates in the protein headers,
write CDS FASTA with ' HE' tag on highly-expressed genes.
Step 2: scripts/06b_grodon.R runs gRodon::predictGrowth in parallel.

usage: 06_growth_rate.py ACCESSION_LIST PROT_DIR OUT_DIR [nproc]
"""
import sys, os, gzip
import multiprocessing as mp
import pandas as pd
import pyhmmer
from pyhmmer.easel import TextSequence, Alphabet

LIB = '/home/user/data/hmm/ribo/ribosomal.hmm'
ALPHA = Alphabet.amino()
HMMS = None
COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def init():
    global HMMS
    with pyhmmer.plan7.HMMFile(LIB) as fh:
        HMMS = list(fh)


def fasta(path):
    out, name, buf = {}, None, []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    out[name] = ''.join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
    if name:
        out[name] = ''.join(buf)
    return out


def work(args):
    acc, pdir, odir = args
    out = os.path.join(odir, f'{acc}.cds.fna.gz')
    if os.path.exists(out):
        return acc, 'exists'
    faa = os.path.join(pdir, acc[:10], f'{acc}.faa.gz')
    fna = os.path.join(pdir + '_fna', acc[:10], f'{acc}.fna.gz')
    if not (os.path.exists(faa) and os.path.exists(fna)):
        return acc, 'missing'
    prots = fasta(faa)
    seqs = [TextSequence(name=k.encode(), sequence=v.replace('*', '')).digitize(ALPHA) for k, v in prots.items()]
    he = set()
    for hits in pyhmmer.hmmsearch(HMMS, seqs, cpus=1, bit_cutoffs='gathering'):
        for h in hits.reported:
            he.add(h.name.decode() if isinstance(h.name, bytes) else h.name)
    contigs = fasta(fna)
    with gzip.open(out + '.tmp', 'wt') as o:
        for k in prots:
            _, contig, idx, s, e, strand, _pl = k.split('|')
            s, e = int(s), int(e)
            nt = contigs.get(contig, '')[s - 1:e]
            if strand == '-1':
                nt = nt.translate(COMP)[::-1]
            if len(nt) < 240 or len(nt) % 3:
                continue
            o.write(f'>{contig}_{idx}{" HE" if k in he else ""}\n{nt}\n')
    os.replace(out + '.tmp', out)
    return acc, f'ok:{len(he)}'


def main():
    accs = pd.read_csv(sys.argv[1], sep='\t').accession.tolist()
    pdir, odir = sys.argv[2], sys.argv[3]
    n = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    os.makedirs(odir, exist_ok=True)
    with mp.Pool(n, initializer=init) as pool, open(os.path.join(odir, 'log.tsv'), 'a') as log:
        for i, (acc, st) in enumerate(pool.imap_unordered(work, [(a, pdir, odir) for a in accs], chunksize=4), 1):
            log.write(f'{acc}\t{st}\n')
            if i % 500 == 0:
                print(i, flush=True)


if __name__ == '__main__':
    main()
