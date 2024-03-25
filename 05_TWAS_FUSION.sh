#!/bin/bash

PROJ_PATH=/public/home/biostat03/project/oaProject
FUSION=/public/home/biostat03/biosoft/fusion_twas-master/
SUMM=${PROJ_PATH}/02_data/TWAS/FUSION_input
OUT=${PROJ_PATH}/03_result/04_twas


for var in OA AD
do

for i in `seq 1 22`
do

# Adipose_Subcutaneous
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Adipose_Subcutaneous.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Adipose_Subcutaneous_${i}.txt

# Adipose_Visceral
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Adipose_Visceral_Omentum.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Adipose_Visceral_${i}.txt

# Blood
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Whole_Blood.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Blood_${i}.txt

# Brain_Substantia_nigra.pos
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Brain_Substantia_nigra.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Brain_Substantia_nigra_${i}.txt

# Brain_Spinal_cord_cervical_c-1
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Brain_Spinal_cord_cervical_c-1.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Brain_Spinal_cord_cervical_${i}.txt

# Brain_Putamen_basal_ganglia
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Brain_Putamen_basal_ganglia.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Brain_Putamen_basal_ganglia_${i}.txt

# Brain_Nucleus_accumbens_basal_ganglia
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Brain_Nucleus_accumbens_basal_ganglia.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Brain_Nucleus_accumbens_basal_ganglia_${i}.txt

# Brain_Hypothalamus
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Brain_Hypothalamus.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Brain_Hypothalamus_${i}.txt

# Brain_Hippocampus
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Brain_Hippocampus.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Brain_Hippocampus_${i}.txt

# Brain_Frontal_Cortex_BA9
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Brain_Frontal_Cortex_BA9.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Brain_Frontal_Cortex_BA9_${i}.txt

# Brain_Cortex
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Brain_Cortex.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Brain_Cortex_${i}.txt

# Brain_Cerebellum
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Brain_Cerebellum.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Cerebellum_${i}.txt

# Brain_Cerebellar_Hemisphere
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Brain_Cerebellar_Hemisphere.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Brain_Cerebellar_Hemisphere_${i}.txt

# Brain_Caudate_basal_ganglia
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Brain_Caudate_basal_ganglia.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Brain_Caudate_basal_ganglia_${i}.txt

# Brain_Anterior_cingulate_cortex_BA24
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Brain_Anterior_cingulate_cortex_BA24.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Brain_Anterior_cingulate_cortex_BA24_${i}.txt

# Brain_Amygdala
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Brain_Amygdala.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Brain_Amygdala_${i}.txt

# Muscle
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Muscle_Skeletal.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Muscle_${i}.txt

# Colon_Sigmoid
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Colon_Sigmoid.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Colon_Sigmoid_${i}.txt

# Colon_Transverse
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Colon_Transverse.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Colon_Transverse_${i}.txt

# Liver
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Liver.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Liver_${i}.txt

# Lung
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Lung.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Lung_${i}.txt

# Nerve_Tibial
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Nerve_Tibial.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Nerve_Tibial_${i}.txt

# Pituitary
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Pituitary.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Pituitary_${i}.txt

# Small_Intestine_Terminal_Ileum
Rscript ${FUSION}FUSION.assoc_test.R \
--sumstats ${SUMM}/${var}.sumstats \
--weights ${FUSION}WEIGHTS/GTExv8.EUR.Small_Intestine_Terminal_Ileum.pos \
--weights_dir ${FUSION}WEIGHTS/ \
--ref_ld_chr ${FUSION}LDREF/1000G.EUR. \
--chr ${i} \
--out ${OUT}/${var}_Small_Intestine_Terminal_Ileum_${i}.txt

done

done