#!/usr/bin/env bash
# Run the unchanged main-project pipeline (scripts/03, 05, 08) on the validation genomes.
# usage: bash validation/02_run_pipeline.sh WORKDIR ACCESSIONS.tsv NAME [NPROC]
#   ACCESSIONS.tsv: accession (RS_GCF_...), ncbi_assembly_name, ncbi_translation_table
#   writes WORKDIR/NAME_hits.tsv.gz, WORKDIR/NAME_hits_nprot.tsv, WORKDIR/NAME_modules.tsv
# Must be run from the repo root (08 reads dictionary/). Large files stay in WORKDIR.
set -euo pipefail
W=$1; ACC=$2; NAME=$3; NP=${4:-4}
mkdir -p "$W"
[ -s "$W/dormancy_library.hmm" ] || zcat dictionary/hmm/dormancy_library.hmm.gz > "$W/dormancy_library.hmm"

python3 scripts/03_fetch_and_call_genes.py "$ACC" "$W/prot" "$NP"

# accessions whose proteome exists
python3 - "$ACC" "$W" "$NAME" <<'EOF'
import sys, os, pandas as pd
acc, w, name = sys.argv[1:4]
t = pd.read_csv(acc, sep='\t')
ok = [a for a in t.accession if os.path.exists(os.path.join(w, 'prot', a[:10], f'{a}.faa.gz'))]
pd.DataFrame({'accession': ok}).to_csv(os.path.join(w, f'{name}_ok.tsv'), sep='\t', index=False)
print(f'{len(ok)}/{len(t)} proteomes available')
EOF

python3 scripts/05_hmm_search.py "$W/dormancy_library.hmm" dictionary/thresholds.tsv "$W/prot" \
    "$W/${NAME}_ok.tsv" "$W/${NAME}_hits.tsv.gz" "$NP"

# META for 08: phylum string must start with "Bacillota" for the Sps rule (GTDB naming)
python3 - "$W" "$NAME" <<'EOF'
import sys, pandas as pd
w, name = sys.argv[1:3]
lab = pd.concat([pd.read_csv(f, sep='\t') for f in ('validation/sporulation_labels.tsv', 'validation/rpf_labels.tsv')
                 if __import__('os').path.exists(f)])
lab['accession'] = [('RS_' if a.startswith('GCF') else 'GB_') + a for a in lab.accession.astype(str)]
bac = lab.phylum.isin(['Firmicutes', 'Tenericutes', 'Bacillota'])
lab['phylum'] = ['Bacillota' if b else p for b, p in zip(bac, lab.phylum)]
lab[['accession', 'phylum']].drop_duplicates('accession').to_csv(f'{w}/{name}_meta.tsv', sep='\t', index=False)
EOF

python3 scripts/08_aggregate_modules.py "$W/${NAME}_hits.tsv.gz" "$W/${NAME}_hits_nprot.tsv" \
    "$W/${NAME}_meta.tsv" "$W/${NAME}_modules.tsv"
