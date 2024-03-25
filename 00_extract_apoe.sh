
PLINK=/public/home/biostat03/biosoft/plink1.9
GENOTYPE=/public/home/Datasets/ukb/geno/plink/chr19
APOE=/public/home/biostat03/project/oaProject/03_result/00_population/apoe.txt
OUTAPOE=/public/home/biostat03/project/oaProject/03_result/00_population/apoe_genotype
${PLINK} --silent --bfile ${GENOTYPE} --extract ${APOE} --recode A --out ${OUTAPOE}

