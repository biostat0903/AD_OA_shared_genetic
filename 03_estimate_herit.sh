#!/bin/bash
source activate ldsc

LDSC=/public/home/biostat03/biosoft/ldsc/ldsc.py
PROJ_PATH=/public/home/biostat03/project/oaProject
REF=/public/home/biostat03/gwas_ref_data/eur_w_ld_chr/

# Heritability
## OA
SUMM=${PROJ_PATH}/02_data/ldsc_summary_stat/ldsc_input/OA.sumstats.gz
HERIT=${PROJ_PATH}/03_result/01_h2/OA
python2 ${LDSC} --h2 ${SUMM}\
 --ref-ld-chr ${REF} --w-ld-chr ${REF} --out ${HERIT}
 ## AD
SUMM=${PROJ_PATH}/02_data/ldsc_summary_stat/ldsc_input/AD.sumstats.gz
HERIT=${PROJ_PATH}/03_result/01_h2/AD
python2 ${LDSC} --h2 ${SUMM}\
 --ref-ld-chr ${REF} --w-ld-chr ${REF} --out ${HERIT}

# Gene Correlation
SUMM_AD=${PROJ_PATH}/02_data/ldsc_summary_stat/ldsc_input/AD.sumstats.gz
SUMM_OA=${PROJ_PATH}/02_data/ldsc_summary_stat/ldsc_input/OA.sumstats.gz
CORR=${PROJ_PATH}/03_result/02_gene_correlation/OA_AD

cd ${PROJ_PATH}/02_data/ldsc_summary_stat/ldsc_input
python2 ${LDSC} \
 --rg ${SUMM_OA},${SUMM_AD} \
 --ref-ld-chr ${REF} \
 --w-ld-chr ${REF} \
 --out ${CORR} \
 --no-intercept

conda deactivate