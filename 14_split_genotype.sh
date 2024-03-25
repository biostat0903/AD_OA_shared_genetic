
#
type=all
proj_path=/public/home/biostat03/project/oaProject/
PRS_path=${proj_path}03_result/09_prs/
snp_file=${PRS_path}PRS/genotype/AD_snp.txt

#
ukb_path=/public/home/Datasets/ukb/geno/EUR/
if [[ "$type" == "hm3" ]]
then
ukb_path=${ukb_path}/hm3_nosex/
out_train_path=${PRS_path}PRS/genotype/train_hm3/
out_valid_path=${PRS_path}PRS/genotype/valid_hm3/
out_valid2_path=${PRS_path}PRS/genotype/valid2_hm3/
out_test_path=${PRS_path}PRS/genotype/test_hm3/
else
out_train_path=${PRS_path}PRS/genotype/train/
out_valid_path=${PRS_path}PRS/genotype/valid/
out_valid2_path=${PRS_path}PRS/genotype/valid2/
out_test_path=${PRS_path}PRS/genotype/test/
fi

# train
eid_file=${PRS_path}id/summ_id.txt
for chr in `seq 1 22`;do
plink1.9 \
--bfile ${ukb_path}chr${chr} \
--keep ${eid_file} \
--make-bed \
--out ${out_train_path}chr${chr}
done

# valid1
eid_file=${PRS_path}id/val_id1.txt
for chr in `seq 1 22`;do
plink1.9 \
--bfile ${ukb_path}chr${chr} \
--keep ${eid_file} \
--make-bed \
--out ${out_valid_path}chr${chr}
done

# valid2
eid_file=${PRS_path}id/val_id2.txt
for chr in `seq 1 22`;do
plink1.9 \
--bfile ${ukb_path}chr${chr} \
--keep ${eid_file} \
--make-bed \
--out ${out_valid2_path}chr${chr}
done

# test
test_eid_file=${PRS_path}id/test_id.txt
for chr in `seq 1 22`;do
plink1.9 \
--bfile ${ukb_path}chr${chr} \
--keep ${test_eid_file} \
--make-bed \
--out ${out_test_path}chr${chr}
done

