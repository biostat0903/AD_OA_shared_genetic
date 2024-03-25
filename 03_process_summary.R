rm(list=ls())
gc()
# load packages
library(bigreadr)
library(data.table)
library(plyr)
library(dplyr)
library(stringr)

# set path
PROJ_PATH <- "/public/home/biostat03/project/oaProject/"
DATA_PROJ <- paste0(PROJ_PATH, "/02_data")
REF_PATH <- "/public/home/biostat03/gwas_ref_data/eur_w_ld_chr/"
KGP_PATH <- "/public/home/biostat03/gwas_ref_data/eur_1kgp_all/"
RAW_SUMM <- paste0(PROJ_PATH, "02_data/summary_stat/raw/")
SUMM_PATH <- paste0(PROJ_PATH, "02_data/summary_stat/cleaned/")

# load reference
hm3_snplist <- fread(paste0(REF_PATH, "/w_hm3.snplist"))
hm3_snp_info <- llply(c(1: 22), function(chr){
  
  ldsc_chr <- fread(paste0(REF_PATH, "/", chr, ".l2.ldscore.gz"))
}) %>% do.call("rbind", .) %>% select((CHR: BP))
hm3_snp_info_join <- left_join(hm3_snplist, hm3_snp_info, by = "SNP")


# munge function
munge.hm3 <- function(data_ldsc){
  
  data_ldsc_ref <- left_join(hm3_snplist, data_ldsc, by = c("SNP"))
  ## delete allele
  na_row1 <- which(data_ldsc_ref$A1.x != data_ldsc_ref$A2.y & data_ldsc_ref$A1.x != data_ldsc_ref$A1.y)
  na_row2 <- which(data_ldsc_ref$A2.x != data_ldsc_ref$A2.y & data_ldsc_ref$A2.x != data_ldsc_ref$A1.y)
  na_row <- c(na_row1, na_row2) %>% unique()
  if(length(na_row) != 0){
    
    data_ldsc_ref$N[na_row] <- NA
    data_ldsc_ref$Z[na_row] <- NA
  }
  ## flip allele
  fl_row <- which(data_ldsc_ref$A1.x == data_ldsc_ref$A2.y & data_ldsc_ref$A2.x == data_ldsc_ref$A1.y)
  if(length(fl_row) != 0){
    
    data_ldsc_ref$Z[fl_row] <- -1 * data_ldsc_ref$Z[fl_row]
  }
  data_ldsc_ref_f <- data.frame(SNP = data_ldsc_ref$SNP, 
                                N = data_ldsc_ref$N,
                                Z = data_ldsc_ref$Z,
                                A1 = data_ldsc_ref$A1.x,
                                A2 = data_ldsc_ref$A2.x)
  
  return(data_ldsc_ref_f)
}

# Process Summary

# LDSC input
# OA
## OA frontier
OA <- fread("/public/home/biostat03/project/oaProject/02_data/summary_stat/cleaned/OA_frontier.txt")
OA <- data.frame(SNP = OA$SNP,
                    N = 455221,
                    Z = OA$Beta/OA$SE,
                    A1 = OA$A1,
                    A2 = OA$A2)
OA <- munge.hm3(OA)
write.table(OA, row.names = F, quote = F, na = "NA", sep = "\t",
            file = paste0(DATA_PROJ, "/ldsc_summary_stat/ldsc_input/OA.sumstats"))
system(paste0("gzip -f ", DATA_PROJ, "/ldsc_summary_stat/ldsc_input/OA.sumstats"))

# AD
## nall ad
igap_ad <- fread(paste0(DATA_PROJ, "/summary_stat/raw/IGAP_stage_1.txt"))
igap_ad_ldsc <- data.frame(SNP = igap_ad$MarkerName, 
                           # N = 74046,
                           N = 63926,
                           Z = igap_ad$Beta/igap_ad$SE,
                           A1 = igap_ad$Effect_allele,
                           A2 = igap_ad$Non_Effect_allele)
igap_ad_hm3 <- munge.hm3(igap_ad_ldsc)
write.table(igap_ad_hm3, row.names = F, quote = F, na = "NA", sep = "\t",
            file = paste0(DATA_PROJ, "/ldsc_summary_stat/ldsc_input/AD.sumstats"))
system(paste0("gzip -f ", DATA_PROJ, "/ldsc_summary_stat/ldsc_input/AD.sumstats"))
