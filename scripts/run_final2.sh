#!/bin/sh
# Second-stage exploratory analyses, run after run_final.sh prints FINALDONE.
cd /home/user/project-melad
H=/home/user/data/hits
until grep -q FINALDONE $H/final.log; do sleep 60; done
Rscript scripts/18_mag_isolate.R $H/full_table.tsv.gz data_small/analysis_tree.nwk results/exploratory > results/mag_isolate_log.txt 2>&1
echo "FINAL2DONE $(date -u +%H:%M)"
