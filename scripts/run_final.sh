#!/bin/sh
# Final confirmatory + exploratory pipeline; waits for gRodon (GROWTHDONE), then runs everything.
cd /home/user/project-melad
H=/home/user/data/hits
until grep -q GROWTHDONE $H/growth.log; do sleep 60; done
echo "start $(date -u +%H:%M)"
python3 scripts/10_build_table.py $H/all_modules.tsv $H/grodon.tsv $H/full_table.tsv.gz
Rscript scripts/09_models.R $H/full_table.tsv.gz data_small/analysis_tree.nwk confirmatory results/confirmatory > results/confirmatory_log.txt 2>&1
Rscript scripts/09_models.R $H/full_table.tsv.gz data_small/analysis_tree.nwk confirmatory results/exploratory scripts/exploratory_config.tsv > results/exploratory_log.txt 2>&1
Rscript scripts/14_extra_tests.R $H/full_table.tsv.gz data_small/analysis_tree.nwk results/confirmatory > results/extra_tests_log.txt 2>&1
Rscript scripts/16_cluster_robust.R $H/full_table.tsv.gz data_small/analysis_tree.nwk results/siteK1 100 200 1 > results/siteK1_log.txt 2>&1
Rscript scripts/16_cluster_robust.R $H/full_table.tsv.gz data_small/analysis_tree.nwk results/siteK3 100 200 3 > results/siteK3_log.txt 2>&1
python3 scripts/15_genus_sanity.py $H/full_table.tsv.gz results/genus_sanity.tsv > /dev/null 2>&1
python3 scripts/13_figures.py $H/full_table.tsv.gz results/confirmatory_models.tsv results/figures results/exploratory_models.tsv results/siteK1_site_capped_betas.tsv > /dev/null 2>&1
python3 scripts/17_decisions.py results/confirmatory_models.tsv results/exploratory_models.tsv results/confirmatory_H6_wald.tsv results/siteK1_cluster_robust.tsv results/decisions.md > /dev/null 2>&1
echo "core done $(date -u +%H:%M)"
Rscript scripts/14_extra_tests.R $H/full_table.tsv.gz data_small/analysis_tree.nwk results/confirmatory mcmc > results/mcmc_modules_log.txt 2>&1
Rscript scripts/16_cluster_robust.R $H/full_table.tsv.gz data_small/analysis_tree.nwk results/siteMCMC 2 2 1 mcmc > results/siteMCMC_log.txt 2>&1
echo "FINALDONE $(date -u +%H:%M)"
