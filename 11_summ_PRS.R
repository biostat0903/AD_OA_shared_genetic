#######
library(bigreadr)
library(dplyr)

# set path
PROJ_PATH <- "/public/home/biostat03/project/oaProject/"
prs_path <- paste0(PROJ_PATH, "03_result/09_prs/PRS/")

### sum prs from DBSLMM
test_eid <- fread2(paste0(prs_path, "DBSLMM_output/AD/AD_chr22_pred.profile"))[,2]
prs_test_mat <- lapply(1:22, function(chrx){
  predx <- fread2(paste0(prs_path, "DBSLMM_output/AD/AD_chr", chrx, "_pred.profile"))
  scorex <- predx[match(test_eid, predx$IID), "SCORESUM"]
}) %>% Reduce("cbind", .)
prs_test <- rowSums(prs_test_mat)

train_eid <- fread2(paste0(prs_path, "DBSLMM_output/AD_train/AD_chr22_pred.profile"))[,2]
prs_train_mat <- lapply(1:22, function(chrx){
  predx <- fread2(paste0(prs_path, "DBSLMM_output/AD_train/AD_chr", chrx, "_pred.profile"))
  scorex <- predx[match(train_eid, predx$IID), "SCORESUM"]
}) %>% Reduce("cbind", .)
prs_train <- rowSums(prs_train_mat)

valid2_eid <- fread2(paste0(prs_path, "DBSLMM_output/AD_valid2/AD_chr22_pred.profile"))[,2]
prs_valid2_mat <- lapply(1:22, function(chrx){
  predx <- fread2(paste0(prs_path, "DBSLMM_output/AD_valid2/AD_chr", chrx, "_pred.profile"))
  scorex <- predx[match(valid2_eid, predx$IID), "SCORESUM"]
}) %>% Reduce("cbind", .)
prs_valid2 <- rowSums(prs_valid2_mat)

#
prs_df <- data.frame(eid = c(test_eid, train_eid, valid2_eid),
                     prs = c(prs_test, prs_train, prs_valid2),
                     group = c(rep("test", length(test_eid)),
                               rep("train", length(train_eid) + length(valid2_eid))))
saveRDS(prs_df, 
        file = paste0(PROJ_PATH, "03_result/09_prs/PRS/PRS_output/prs_df.rds"))

# ### sum target PRS from SCT
# test_eid <- fread2(paste0(PROJ_PATH, "03_result/09_prs/PRS/genotype/test/chr22.fam"))[,2]
# prs_test_path <- paste0(PROJ_PATH, "03_result/09_prs/PRS/CT_output/AD/")
# pred_test <- fread2(paste0(prs_test_path, "esteff_SCT_target_pred.profile"))
# prs_test <- pred_test[match(test_eid, pred_test$IID), "SCORESUM"]
# 
# # add prs 1
# train_eid <- fread2(paste0(PROJ_PATH, "03_result/09_prs/PRS/genotype/train/chr22.fam"))[,2]
# prs_train_path <- paste0(PROJ_PATH, "03_result/09_prs/PRS/CT_output/AD_train/")
# pred_train <- fread2(paste0(prs_train_path, "esteff_SCT_target_pred.profile"))
# prs_train <- pred_train[match(train_eid, pred_train$IID), "SCORESUM"]
# 
# #
# prs_df <- data.frame(eid = c(test_eid, train_eid),
#                      prs = c(prs_test, prs_train))
# 
# saveRDS(prs_df, 
#         file = paste0(PROJ_PATH, "03_result/09_prs/PRS/PRS_output/prs_target_df.rds"))

## combine with phenotype
library(bigreadr)
library(dplyr)

# load phenotype data
PROJ_PATH <- "/public/home/biostat03/project/oaProject/"
load(paste0(PROJ_PATH, "03_result/00_cohort/sample_info/tot_dff1227.RData"))

# load PRS data
PC_eur <- fread2(paste0(PROJ_PATH, "03_result/09_prs/id/PC_eur.txt"))
prs_df <- readRDS(paste0(PROJ_PATH, "03_result/09_prs/PRS/PRS_output/prs_df.rds"))
prs_target <- readRDS(paste0(PROJ_PATH, "03_result/09_prs/PRS/PRS_output/prs_target_df.rds"))
comb_df <- cbind(prs_df,
                 tot_dff[match(prs_df$eid, tot_dff$eid), 
                        c("Age", "Sex", "BMI_g", "oa", "nerv", "Time", "OS")],
                 PC_eur[match(prs_df$eid, PC_eur$eid), -1])
comb_df <- subset(comb_df, !is.na(comb_df$nerv))

# define prs group
quant10_prs <- quantile(comb_df$prs, probs = seq(0, 1, 0.1))
comb_df$prs_g10 <- ifelse(comb_df$prs <= quant10_prs[2], 1,
                          ifelse(comb_df$prs <= quant10_prs[3], 2,
                                 ifelse(comb_df$prs <= quant10_prs[4], 3,
                                        ifelse(comb_df$prs <= quant10_prs[5], 4,
                                               ifelse(comb_df$prs <= quant10_prs[6], 5,
                                                      ifelse(comb_df$prs <= quant10_prs[7], 6,
                                                             ifelse(comb_df$prs <= quant10_prs[8], 7,
                                                                    ifelse(comb_df$prs <= quant10_prs[9], 8,
                                                                           ifelse(comb_df$prs <= quant10_prs[10], 9, 10)))))))))
# 333
quant3_prs <- quantile(comb_df$prs, probs = seq(0, 1, length.out = 4))
comb_df$prs_g3 <- ifelse(comb_df$prs <= quant3_prs[2], 1,
                         ifelse(comb_df$prs <= quant3_prs[3], 2, 3))
# 343
quant3_prs2 <- quantile(comb_df$prs, probs = c(0, 0.3, 0.7, 1))
comb_df$prs_g32 <- ifelse(comb_df$prs <= quant3_prs2[2], 1,
                         ifelse(comb_df$prs <= quant3_prs2[3], 2, 3))
# 181
quant3_prs3 <- quantile(comb_df$prs, probs = c(0, 0.1, 0.9, 1))
comb_df$prs_g33 <- ifelse(comb_df$prs <= quant3_prs3[2], 1,
                          ifelse(comb_df$prs <= quant3_prs3[3], 2, 3))

# 262
quant3_prs4 <- quantile(comb_df$prs, probs = c(0, 0.2, 0.8, 1))
comb_df$prs_g34 <- ifelse(comb_df$prs <= quant3_prs4[2], 1,
                          ifelse(comb_df$prs <= quant3_prs4[3], 2, 3))
saveRDS(comb_df, 
        file = paste0(PROJ_PATH, "03_result/09_prs/PRS/PRS_output/comb_df.rds"))

### sum prs from DBSLMM AFR
train_eid <- fread2(paste0(prs_path, "DBSLMM_output/AD_AFR/AD_chr22_pred.profile"))[,2]
prs_afr_mat <- lapply(c(1:22), function(chrx){
  predx <- fread2(paste0(prs_path, "DBSLMM_output/AD_AFR/AD_chr", chrx, "_pred.profile"))
  scorex <- predx[match(train_eid, predx$IID), "SCORESUM"]
}) %>% Reduce("cbind", .)
prs_afr <- rowSums(prs_afr_mat)
prs_df <- data.frame(eid = train_eid,
                     prs = prs_afr)

saveRDS(prs_df, 
        file = paste0(PROJ_PATH, "03_result/09_prs/PRS/PRS_output/prs_afr_df.rds"))

# combine with phenotype
load(paste0(PROJ_PATH, "03_result/00_cohort/sample_info/tot_dat_AFR.RData"))
#
PC_afr_file <- paste0(PROJ_PATH, "03_result/09_prs/id/PC_afr.txt")
if (file.exists(PC_afr_file)) {
  PC_afr <- fread2(PC_afr_file)
} else {
  PC_afr <- c("eid", 
              paste0("22009-0.", c(1:20)))  %>%
    fread2("/public/home/Datasets/ukb/ukb47503.csv.gz", 
           select = .) %>% 
    filter(eid %in% prs_df$eid)            
  colnames(PC_afr) <- c("eid", paste0("PC", c(1:20)))
  fwrite2(PC_afr, file = paste0(PROJ_PATH, "03_result/09_prs/id/PC_afr.txt"), 
          sep = "\t")
}
#

comb_df <- cbind(prs_df,
                 tot_df[match(prs_df$eid, tot_df$eid), 
                        c("Age", "Sex", "BMI_g", "oa", "nerv")],
                 PC_afr[match(prs_df$eid, PC_afr$eid), -1])
comb_df <- subset(comb_df, !is.na(comb_df$nerv))

saveRDS(comb_df, 
        file = paste0(PROJ_PATH, "03_result/09_prs/PRS/PRS_output/comb_df_afr.rds"))

