
# 1.software path
projectpath=/public/home/biostat03/project/oaProject/
dbslmm_path=/public/home/biostat03/biosoft/DBSLMM/
TUNE=${dbslmm_path}software/TUNE.R
dbslmm=${dbslmm_path}scr/dbslmm

# 2.input/output path
summ_path=${projectpath}03_result/09_prs/PRS/summary/AD/
val_prefix=${projectpath}03_result/09_prs/PRS/genotype/valid/chr
valp=${projectpath}03_result/09_prs/id/val_pheno1.txt
out_path=${projectpath}03_result/09_prs/PRS/DBSLMM_output/AD/
val_cov=${projectpath}03_result/09_prs/id/val_covarites1.txt

# 3.model parameters
col=1
index=r2
thread=1

# 4. tuning
summchr=${summ_path}AD_chr22
summchr_prefix=`echo ${summchr##*/}`
summchr_prefix2=`echo ${summchr_prefix%_*}`

if [[ ! -n "$val_cov" ]]
then 
Rscript ${TUNE} --phenoPred ${out_path}${summchr_prefix2} --phenoVal ${valp},${col} \
--h2Range 0.8,1,1.2 --index ${index}
else 
Rscript ${TUNE} --phenoPred ${out_path}${summchr_prefix2} --phenoVal ${valp},${col} \
--h2Range 0.8,1,1.2 --index ${index} --cov ${val_cov}
fi

hbest=`cat ${out_path}${summchr_prefix2}_hbest.${index}`
for chr in `seq 1 22`
do
mv ${out_path}${summchr_prefix2}_chr${chr}_h2f${hbest}.dbslmm.txt ${out_path}${summchr_prefix2}_chr${chr}_best.dbslmm.txt
done

################################
# calculate PRS for AD Predict #
###############################
projectpath=/public/home/biostat03/project/oaProject/

type=test
bfile_prefix=${projectpath}03_result/09_prs/PRS/genotype/${type}/chr
out_path=${projectpath}03_result/09_prs/PRS/DBSLMM_output/
est_prefix=${out_path}AD/AD_chr
pred_prefix=${out_path}AD/AD_chr


# plink 1.9
for chr in `seq 1 3`
do
plink1.9 \
--bfile ${bfile_prefix}${chr} \
--score ${est_prefix}${chr}_best.dbslmm.txt 1 2 4 sum \
--out ${pred_prefix}${chr}_pred
done

rm ${pred_prefix}*_pred.log
rm ${pred_prefix}*_pred.nopred
##########################
for type in train valid2
do

projectpath=/public/home/biostat03/project/oaProject/
bfile_prefix=${projectpath}03_result/09_prs/PRS/genotype/${type}/chr
out_path=${projectpath}03_result/09_prs/PRS/DBSLMM_output/
est_prefix=${out_path}AD/AD_chr
pred_prefix=${out_path}AD_${type}/AD_chr

# plink 1.9
for chr in `seq 1 22`
do
plink1.9 \
--bfile ${bfile_prefix}${chr} \
--score ${est_prefix}${chr}_best.dbslmm.txt 1 2 4 sum \
--out ${pred_prefix}${chr}_pred
done

rm ${pred_prefix}*_pred.log
rm ${pred_prefix}*_pred.nopred

##########################
projectpath=/public/home/biostat03/project/oaProject/
bfile_prefix=/public/home/Datasets/ukb/geno/AFR/chr
out_path=${projectpath}03_result/09_prs/PRS/DBSLMM_output/
est_prefix=${out_path}AD/AD_chr

pred_prefix=${out_path}AD_AFR/AD_chr

# plink 1.9
for chr in `seq 1 22`
do
plink1.9 \
--bfile ${bfile_prefix}${chr} \
--score ${est_prefix}${chr}_best.dbslmm.txt 1 2 4 sum \
--out ${pred_prefix}${chr}_pred
done

rm ${pred_prefix}*_pred.log
rm ${pred_prefix}*_pred.nopred





