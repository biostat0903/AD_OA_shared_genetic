library(bigreadr)
library(dplyr)

# Load data
PROJ_PATH <- "/public/home/biostat03/project/oaProject/"
load(paste0(PROJ_PATH, "03_result/00_cohort/sample_info/tot_dat_all0830.RData"))

# Split data
set.seed(20230528)
sample_id <- sample(1:10, size = nrow(tot_df), replace = T)

# Get validation set and test set
test_id <- tot_df$eid[sample_id == 1] %>% sort
val_id1 <- tot_df$eid[sample_id == 2] %>% sort
val_id2 <- tot_df$eid[sample_id == 3] %>% sort
summ_id <- tot_df$eid[sample_id %in% c(4:10)] %>% sort
# write.table(cbind(test_id, test_id), row.names = F, col.names = F,
#             file = paste0(PROJ_PATH, "03_result/09_prs/id/test_id.txt"))
# write.table(cbind(val_id1, val_id1), row.names = F, col.names = F,
#             file = paste0(PROJ_PATH, "03_result/09_prs/id/val_id1.txt"))
# write.table(cbind(val_id2, val_id2), row.names = F, col.names = F,
#             file = paste0(PROJ_PATH, "03_result/09_prs/id/val_id2.txt"))
# write.table(cbind(summ_id, summ_id), row.names = F, col.names = F,
#             file = paste0(PROJ_PATH, "03_result/09_prs/id/summ_id.txt"))

# Get genotype PC
PC_eur <- c("eid", 
             paste0("22009-0.", c(1:20)))  %>%
  fread2("/public/home/Datasets/ukb/ukb47503.csv.gz", 
         select = .) %>% 
  filter(eid %in% tot_df$eid)            
colnames(PC_eur) <- c("eid", paste0("PC", c(1:20)))
fwrite2(PC_eur, file = paste0(PROJ_PATH, "03_result/09_prs/id/PC_eur.txt"), 
        sep = "\t")

## Fit model in train set
use_y <- "nerv"
train_idx <- fread2(paste0(PROJ_PATH, "03_result/09_prs/PRS/genotype/train/chr22.fam"))[,2]
train_df <- tot_df[match(train_idx, tot_df$eid), ]
logit_prs_coef <- glm(train_df[[use_y]] ~ train_df$BMI_g + train_df$Sex + train_df$Age+ 
                        as.matrix(PC_eur[match(train_idx, PC_eur$eid), -1]), 
                      family = "binomial") %>% coef

# Valid two validation set
val_df1 <- cbind(1, cbind(tot_df[match(val_id1, tot_df$eid), c("Sex", "Age")], 
                          PC_eur[match(val_id1, PC_eur$eid), -1]))
val_cov1 <- as.matrix(val_df1) %*% as.matrix(logit_prs_coef) 
val_cov1 <- 1/(1+exp(-val_cov1))
auc_cov_oa1 <- pROC::auc(tot_df[match(val_id1, tot_df$eid), use_y], val_cov1[, 1])
val_df2 <- cbind(1, cbind(tot_df[match(val_id2, tot_df$eid), c("Sex", "Age")], 
                          PC_eur[match(val_id2, PC_eur$eid), -1]))
val_cov2 <- as.matrix(val_df2) %*% as.matrix(logit_prs_coef) 
val_cov2 <- 1/(1+exp(-val_cov2))
auc_cov_oa2 <- pROC::auc(tot_df[match(val_id2, tot_df$eid), use_y], val_cov2[, 1])
write.table(val_cov1, col.names = F, row.names = F, quote = F, 
            file = paste0(PROJ_PATH, "03_result/09_prs/id/val_cov1.txt"))
write.table(val_cov2, col.names = F, row.names = F, quote = F, 
            file = paste0(PROJ_PATH, "03_result/09_prs/id/val_cov2.txt"))

# output valid covariates for PRS
valid_idx1 <- fread2(paste0(PROJ_PATH, "03_result/09_prs/PRS/genotype/valid/chr22.fam"))[,2]
val_df1 <- cbind(1, cbind(tot_df[match(valid_idx1, tot_df$eid), c("Sex", "Age")], 
                          PC_eur[match(valid_idx1, PC_eur$eid), -1]))
write.table(val_df1, 
            col.names = F, row.names = F, quote = F, 
            file = paste0(PROJ_PATH, "03_result/09_prs/id/val_covarites1.txt"))

# output valid phenotype for PRS
write.table(tot_df[match(valid_idx1, tot_df$eid), "nerv"], 
            col.names = F, row.names = F, quote = F, 
            file = paste0(PROJ_PATH, "03_result/09_prs/id/val_pheno1.txt"))

# output valid covariates for PTRS
valid_idx2 <- fread2(paste0(PROJ_PATH, "03_result/09_prs/PRS/genotype/valid2_hm3/chr22.fam"))[,2]
val_df2 <- cbind(tot_df[match(valid_idx2, tot_df$eid), c("Sex", "Age")], 
                 PC_eur[match(valid_idx2, PC_eur$eid), -1])

write.table(val_df2, 
            col.names = F, row.names = F, quote = F, 
            file = paste0(PROJ_PATH, "03_result/09_prs/id/val_covarites2.txt"))

# output valid phenotype for PTRS
write.table(tot_df[match(valid_idx2, tot_df$eid), "nerv"], 
            col.names = F, row.names = F, quote = F, 
            file = paste0(PROJ_PATH, "03_result/09_prs/id/val_pheno2.txt"))


