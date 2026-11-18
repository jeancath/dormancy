#!/usr/bin/env python3
"""Phase 5a: download genome assemblies from NCBI and predict proteins with pyrodigal.

Uniform gene calling (same caller/version for isolates and MAGs) avoids annotation-pipeline
artefacts. Proteins are stored gzipped with headers:
  >{accession}|{contig}|{gene_idx}|{start}|{end}|{strand}|{plasmid_flag}
so that gene neighbourhoods (TA pairing) and plasmid location can be recovered later.

usage: 03_fetch_and_call_genes.py accessions.tsv OUTDIR [nproc]
accessions.tsv columns: accession (GTDB, e.g. RS_GCF_000005845.2), ncbi_assembly_name,
                        ncbi_translation_table
"""
import sys, os, re, gzip, io, time, random
import multiprocessing as mp
import pandas as pd
import requests
import pyrodigal

BASE = 'https://ftp.ncbi.nlm.nih.gov/genomes/all'
S = None


def session():
    global S
    if S is None:
        S = requests.Session()
    return S


def ftp_dir(acc):
    a = acc.split('_', 1)[1] if acc[:3] in ('RS_', 'GB_') else acc
    pre, num = a.split('_')
    digits = num.split('.')[0]
    return a, f'{BASE}/{pre}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}'


def get(url, tries=5, stream=False):
    for k in range(tries):
        try:
            r = session().get(url, timeout=120)
            if r.status_code == 200:
                return r
            if r.status_code == 404:
                return None
        except Exception:
            pass
        time.sleep(2 ** k + random.random())
    return None


def download(acc, asm_name):
    a, d = ftp_dir(acc)
    # list the parent directory to find the exact assembly folder (names can contain odd chars)
    r = get(d + '/')
    folder = None
    if r is not None:
        cands = re.findall(r'href="(' + re.escape(a) + r'_[^"/]+)/?"', r.text)
        if cands:
            folder = cands[0]
    if folder is None:
        # try the other prefix (GCA <-> GCF)
        alt = a.replace('GCF_', 'GCA_') if a.startswith('GCF_') else a.replace('GCA_', 'GCF_')
        _, d2 = ftp_dir(alt)
        r = get(d2 + '/')
        if r is not None:
            cands = re.findall(r'href="(' + re.escape(alt) + r'_[^"/]+)/?"', r.text)
            if cands:
                folder = cands[0]; d = d2
    if folder is None:
        return None
    r = get(f'{d}/{folder}/{folder}_genomic.fna.gz')
    if r is None:
        return None
    return gzip.decompress(r.content).decode()


def parse_fasta(txt):
    name, desc, seq = None, None, []
    for line in txt.splitlines():
        if line.startswith('>'):
            if name is not None:
                yield name, desc, ''.join(seq)
            parts = line[1:].split(None, 1)
            name = parts[0]; desc = parts[1] if len(parts) > 1 else ''; seq = []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, desc, ''.join(seq)


def work(args):
    acc, asm, ttable, outdir = args
    out = os.path.join(outdir, acc[:10], f'{acc}.faa.gz')
    fna = os.path.join(outdir + '_fna', acc[:10], f'{acc}.fna.gz')
    if os.path.exists(out) and os.path.getsize(out) > 100 and os.path.exists(fna):
        return acc, 'exists', 0
    os.makedirs(os.path.dirname(out), exist_ok=True)
    os.makedirs(os.path.dirname(fna), exist_ok=True)
    txt = download(acc, asm)
    if txt is None:
        return acc, 'download_failed', 0
    with gzip.open(fna + '.tmp', 'wt', compresslevel=6) as fh:
        fh.write(txt)
    os.replace(fna + '.tmp', fna)
    contigs = list(parse_fasta(txt))
    total = sum(len(s) for _, _, s in contigs)
    tt = int(ttable) if str(ttable).isdigit() else 11
    if total >= 100000:
        gf = pyrodigal.GeneFinder(meta=False)
        gf.train(*[s.encode() for _, _, s in contigs if len(s) >= 1], translation_table=tt)
    else:
        gf = pyrodigal.GeneFinder(meta=True)
    buf = io.StringIO()
    n = 0
    for cname, desc, seq in contigs:
        if len(seq) < 1:
            continue
        pl = 1 if re.search(r'plasmid', desc, re.I) else 0
        genes = gf.find_genes(seq.encode())
        for gi, g in enumerate(genes, 1):
            prot = g.translate(translation_table=tt if not gf.meta else None).rstrip('*')
            buf.write(f'>{acc}|{cname}|{gi}|{g.begin}|{g.end}|{g.strand}|{pl}\n{prot}\n')
            n += 1
    tmp = out + '.tmp'
    with gzip.open(tmp, 'wt', compresslevel=6) as fh:
        fh.write(buf.getvalue())
    os.replace(tmp, out)
    return acc, 'ok', n


def main():
    tab = pd.read_csv(sys.argv[1], sep='\t')
    outdir = sys.argv[2]
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    os.makedirs(outdir, exist_ok=True)
    jobs = [(r.accession, r.ncbi_assembly_name, r.ncbi_translation_table, outdir) for r in tab.itertuples()]
    log = open(os.path.join(outdir, 'fetch_log.tsv'), 'a')
    t0 = time.time()
    with mp.Pool(nproc) as pool:
        for i, (acc, status, n) in enumerate(pool.imap_unordered(work, jobs, chunksize=1), 1):
            log.write(f'{acc}\t{status}\t{n}\n'); log.flush()
            if i % 50 == 0:
                print(f'{i}/{len(jobs)} {time.time() - t0:.0f}s', flush=True)
    log.close()


if __name__ == '__main__':
    main()
