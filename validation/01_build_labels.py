#!/usr/bin/env python3
"""Build independent sporulation phenotype labels and map them to NCBI RefSeq assemblies.

Label sources (none was used to build dictionary/):
  * Madin et al. 2020, Sci Data 7:170, doi:10.1038/s41597-020-0497-4 ("bacteria-archaea-traits",
    condensed_traits_NCBI.csv). Its `sporulation` field is compiled from five phenotype sources:
    fierer (IJSEM species descriptions, Barberan et al. 2017), gold, patric, microbe-directory and
    protraits. ProTraits mixes text mining with genome-based prediction, so it is never allowed to
    create a label on its own; it can only veto one (a conflict drops the species).
  * Browne et al. 2021, Genome Biol 22:204, doi:10.1186/s13059-021-02428-6, Table S1: only the
    ethanol-resistance assay column (experimental), and only rows with an NCBI GCF accession.
    Their genome-signature column ("spore-former characterisation") is NOT used as truth.
  * Galperin et al. 2022 Table S1: NOT reachable (PMC serves the supplement behind a browser
    challenge; Europe PMC says not open access). Not used.

usage: python3 validation/01_build_labels.py SRC_DIR OUT_DIR
  SRC_DIR holds madin_traits.csv and browne_s1.tsv; writes
  validation/sporulation_labels.tsv and OUT_DIR/val_accessions.tsv (input of scripts/03).
"""
import sys, os, json, time, random
import pandas as pd
import requests

SRC, OUT = sys.argv[1], sys.argv[2]
SEED = 20261117
random.seed(SEED)
DS = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{}/dataset_report'
LEVEL = {'Complete Genome': 0, 'Chromosome': 1, 'Scaffold': 2, 'Contig': 3}
CURATED = ['fierer', 'gold', 'patric', 'microbe-directory']
BACILLOTA = {'Firmicutes', 'Tenericutes'}          # NCBI 2020 names; GTDB merges both into Bacillota
TT4_ORDERS = {'Mycoplasmatales', 'Entomoplasmatales'}  # genetic code 4 (UGA = Trp)
CAP_PER_ORDER_LABEL = 70
NON_BACILLOTA_NO = {'Proteobacteria': 25, 'Actinobacteria': 15, 'Bacteroidetes': 12,
                    'Cyanobacteria': 4, 'Deinococcus-Thermus': 3, 'Spirochaetes': 3,
                    'Chloroflexi': 2, 'Thermotogae': 2, 'Verrucomicrobia': 2, 'Euryarchaeota': 3}
NON_BACILLOTA_YES = 12   # QC only: non-endospore "sporulation" (exospores, myxospores)

x = pd.read_csv(os.path.join(SRC, 'madin_traits.csv'), low_memory=False,
                usecols=['species_tax_id', 'data_source', 'species', 'sporulation', 'phylum',
                         'class', 'order', 'family', 'ref_id'])
x = x[x.sporulation.notna() & x.species_tax_id.notna()]
x['species_tax_id'] = x.species_tax_id.astype(int)
per = x.groupby(['species_tax_id', 'data_source']).sporulation.agg(lambda v: '/'.join(sorted(set(v)))).unstack()
tax = x.drop_duplicates('species_tax_id').set_index('species_tax_id')[['species', 'phylum', 'class', 'order', 'family']]

rows = []
for tid, r in per.iterrows():
    cur = {s: r[s] for s in CURATED if s in r and isinstance(r[s], str)}
    if not cur:
        continue
    vals = set('/'.join(cur.values()).split('/'))
    pt = r.get('protraits')
    if isinstance(pt, str):
        vals |= set(pt.split('/'))
    if len(vals) != 1:
        continue                      # any disagreement -> drop
    lab = 1 if vals == {'yes'} else 0
    tier = 'A' if 'fierer' in cur else 'B'
    ev = '; '.join(f'{k}={v}' for k, v in cur.items()) + (f'; protraits={pt}' if isinstance(pt, str) else '')
    rows.append(dict(species_tax_id=tid, label=lab, tier=tier, n_sources=len(cur), evidence=ev))
lab = pd.DataFrame(rows).set_index('species_tax_id').join(tax).reset_index()
lab = lab[lab.species.notna() & ~lab.species.str.contains(r' sp\.|Candidatus|uncultured', regex=True)]
print('consistent labels:', len(lab), lab.groupby(lab.phylum.isin(BACILLOTA)).label.value_counts().to_dict())

# ---- sampling (fixed seed): Bacillota stratified by order x label; prefer tier A, more sources
def pick(df, n):
    df = df.assign(_r=[random.random() for _ in range(len(df))])
    df = df.sort_values(['tier', 'n_sources', '_r'], ascending=[True, False, True])
    return df.head(n)

bac = lab[lab.phylum.isin(BACILLOTA)]
sel = [pick(g, CAP_PER_ORDER_LABEL * 2) for _, g in bac.groupby(['order', 'label'])]  # 2x: mapping losses
nb = lab[~lab.phylum.isin(BACILLOTA)]
for ph, n in NON_BACILLOTA_NO.items():
    sel.append(pick(nb[(nb.phylum == ph) & (nb.label == 0)], n * 2))
sel.append(pick(nb[nb.label == 1], NON_BACILLOTA_YES * 2))
cand = pd.concat(sel)
print('candidates to map:', len(cand))

S = requests.Session()


def best_assembly(tid):
    for params in ({'filters.assembly_source': 'refseq', 'filters.is_type_material': 'true'},
                   {'filters.assembly_source': 'refseq', 'filters.reference_only': 'true'}):
        params = dict(params, page_size=200)
        for k in range(4):
            try:
                rr = S.get(DS.format(tid), params=params, timeout=60)
                if rr.status_code == 200:
                    break
            except Exception:
                pass
            time.sleep(2 ** k)
        else:
            continue
        reps = rr.json().get('reports', [])
        if not reps:
            continue
        reps.sort(key=lambda r: (LEVEL.get(r['assembly_info'].get('assembly_level'), 9),
                                 -float(r.get('assembly_stats', {}).get('contig_n50', 0) or 0)))
        r = reps[0]
        return dict(accession=r['accession'], ncbi_assembly_name=r['assembly_info'].get('assembly_name'),
                    assembly_level=r['assembly_info'].get('assembly_level'),
                    strain=(r['organism'].get('infraspecific_names') or {}).get('strain', ''),
                    organism=r['organism'].get('organism_name'),
                    type_material=int('type_material' in r))
    return None


cache_p = os.path.join(OUT, 'assembly_map_cache.json')
cache = json.load(open(cache_p)) if os.path.exists(cache_p) else {}
out = []
for i, r in enumerate(cand.to_dict('records'), 1):
    k = str(r['species_tax_id'])
    if k not in cache:
        cache[k] = best_assembly(r['species_tax_id'])
        time.sleep(0.25)
        if i % 50 == 0:
            json.dump(cache, open(cache_p, 'w')); print(i, flush=True)
    if cache[k]:
        out.append({**r, **cache[k]})
json.dump(cache, open(cache_p, 'w'))
m = pd.DataFrame(out)
m = m.drop_duplicates('accession')

# final caps after mapping
fin = []
for (o, l), g in m[m.phylum.isin(BACILLOTA)].groupby(['order', 'label']):
    fin.append(g.head(CAP_PER_ORDER_LABEL))
nbm = m[~m.phylum.isin(BACILLOTA)]
for ph, n in NON_BACILLOTA_NO.items():
    fin.append(nbm[(nbm.phylum == ph) & (nbm.label == 0)].head(n))
fin.append(nbm[nbm.label == 1].head(NON_BACILLOTA_YES))
m = pd.concat(fin)
m['label_source'] = 'Madin2020_tier' + m.tier
m['set'] = ['calibration' if p in BACILLOTA else ('specificity' if l == 0 else 'qc_nonendospore')
            for p, l in zip(m.phylum, m.label)]

# ---- Browne 2021 experimental (ethanol-resistance) rows with a GCF accession
b = pd.read_csv(os.path.join(SRC, 'browne_s1.tsv'), sep='\t')
b = b[b.etoh.notna() & b.acc.astype(str).str.startswith('GCF_')]
br = []
for r in b.itertuples():
    lab_b = 0 if r.etoh == 'no' else 1
    br.append(dict(accession=r.acc, label=lab_b, label_source='Browne2021_ethanol',
                   evidence=f'Table S1 ethanol resistance="{r.etoh}"; family {r.family}',
                   phylum='Firmicutes', order=r.taxorder, family=r.family, set='browne_experimental'))
m = pd.concat([m, pd.DataFrame(br)], ignore_index=True)

# fill names for Browne rows from NCBI
for i, r in m[m.label_source == 'Browne2021_ethanol'].iterrows():
    rr = S.get(f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{r.accession}/dataset_report', timeout=60).json()
    if not rr.get('reports'):
        rr = S.get(f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{r.accession}/dataset_report',
                   params={'filters.assembly_version': 'all_assemblies'}, timeout=60).json()
    if not rr.get('reports'):
        print('no NCBI report for', r.accession); m = m.drop(i); continue
    rep = rr['reports'][0]
    m.loc[i, 'organism'] = rep['organism']['organism_name']
    m.loc[i, 'strain'] = (rep['organism'].get('infraspecific_names') or {}).get('strain', '')
    m.loc[i, 'ncbi_assembly_name'] = rep['assembly_info']['assembly_name']
    m.loc[i, 'assembly_level'] = rep['assembly_info']['assembly_level']
    m.loc[i, 'type_material'] = int('type_material' in rep)

m['ncbi_translation_table'] = [4 if o in TT4_ORDERS else 11 for o in m.order]
m['gtdb_acc'] = 'RS_' + m.accession
m['organism'] = m.organism.fillna(m.species)
cols = ['organism', 'strain', 'accession', 'label', 'label_source', 'evidence', 'set', 'species',
        'species_tax_id', 'phylum', 'class', 'order', 'family', 'assembly_level', 'type_material']
m[cols].to_csv('validation/sporulation_labels.tsv', sep='\t', index=False)
m[['gtdb_acc', 'ncbi_assembly_name', 'ncbi_translation_table']].rename(columns={'gtdb_acc': 'accession'}).to_csv(
    os.path.join(OUT, 'val_accessions.tsv'), sep='\t', index=False)
print(m.groupby(['set', 'label']).size())
