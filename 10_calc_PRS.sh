#!/bin/bash
PACK_DIR=/public/home/biostat03/biosoft/
SPLITCHR=${PACK_DIR}DBSLMM/software/SPLITCHR.R
PLINK=${PACK_DIR}plink1.9
DBSLMM=${PACK_DIR}DBSLMM/DBSLMM_script.sh
PROJ_PATH=/public/home/biostat03/project/oaProject/
DATADIR=${PROJ_PATH}03_result/09_risk_model/PRS/

# # Split into chromosome
# summ=${DATADIR}summary/OA/OA.assoc.txt
# Rscript ${SPLITCHR} --summary ${summ}
# summ=${DATADIR}summary/AD/AD.assoc.txt
# Rscript ${SPLITCHR} --summary ${summ}

# Set parameters
BLOCK=${PACK_DIR}DBSLMM/block_data/EUR/chr
model=DBSLMM
type=t
col=1
index=r2
thread=1

# DBSLMM tuning version (without covariates)
## OA
herit=${PROJ_PATH}03_result/01_h2/OA.log
valp=${DATADIR}phenotype/OA_train.txt

# for chr in `seq 1 22`;do

# valg=${DATADIR}genotype/train/chr${chr}
# summ=${DATADIR}summary/OA/OA_chr${chr}
valg=${DATADIR}genotype/train/chr
summ=${DATADIR}summary/OA/OA
outpath=${DATADIR}DBSLMM_output/OA/
REF_GENOTYPE=/public/home/Datasets/1000GP/EUR/hm3/chr
sh ${DBSLMM} -D ${PACK_DIR} -p ${PLINK} -B ${BLOCK} -s ${summ} -H ${herit} -m ${model} -G ${valg} -P ${valp}\
             -R ${REF_GENOTYPE} -l ${col} -T ${type}  -i ${index} -t ${thread} -o ${outpath}
# done

## AD
herit=${PROJ_PATH}03_result/01_h2/AD.log
valp=${DATADIR}phenotype/AD_train.txt

for chr in `seq 1 22`;do

valg=${DATADIR}genotype/train/chr${chr}
summ=${DATADIR}summary/AD/AD_chr${chr}
outpath=${DATADIR}DBSLMM_output/AD/
sh ${DBSLMM} -D ${PACK_DIR} -p ${PLINK} -B ${BLOCK} -s ${summ} -H ${herit} -m ${model} -G ${valg} -P ${valp}\
             -l ${col} -T ${type}  -i ${index} -t ${thread} -o ${outpath}
done

# PLINK PRS
## OA
### Predict
bfilete=${DATADIR}genotype/train/test/chr
est=${DATADIR}DBSLMM_output/AD/
InterPred=${DATADIR}
## plink 1.9
plink=/your/path/plink1.9
${plink} --bfile ${bfilete}${chr} --score ${est}${chr}.assoc.dbslmm.txt 1 2 4 sum --out ${InterPred}${chr}
rm ${InterPred}${chr}.log


