#### use
library(bigreadr)
library(dplyr)
library(pROC)
library(trend)
# load data
PROJ_PATH <- "/public/home/biostat03/project/oaProject/"
comb_df <- readRDS(paste0(PROJ_PATH, "03_result/09_prs/PRS/PRS_output/comb_df.rds"))

####### PRS risk #######
fomu1 <- paste0("nerv ~ Age + Sex + BMI_g + ", 
                paste(paste0("PC", 1:20), collapse = " + "),
                " + prs")
logit_prs10 <- glm(fomu1, data = comb_df, family = "binomial")
summ_coef_prs <- summary(logit_prs10)$coefficients["prs", ]
RR_prs <- exp(summ_coef_prs[1]) %>% round(2)
lowCI_prs <- (summ_coef_prs[1] - 1.96*summ_coef_prs[2]) %>% exp %>% round(2)
highCI_prs <- (summ_coef_prs[1] + 1.96*summ_coef_prs[2]) %>% exp %>% round(2)
p_prs <- summ_coef_prs[4] %>% round(4)
print(c(paste0(RR_prs, "(",
               lowCI_prs, ",",
               highCI_prs, ")"),
        p_prs))

####### group PRS #######
comb_df$prs_g10 <- as.factor(comb_df$prs_g10)
fomu2 <- paste0("nerv ~ Age + Sex + BMI_g + ", 
                paste(paste0("PC", 1:20), collapse = " + "),
                " + prs_g10")
logit_prs10 <- glm(fomu2, data = comb_df, family = "binomial")
summ_coef_prs <- summary(logit_prs10)$coefficients[paste0("prs_g10", 2:10), ]
or_prs10_df <- data.frame(
  OR = c(1, summ_coef_prs[,1] %>% exp),
  lowCI = c(1, (summ_coef_prs[,1] - 1.96*summ_coef_prs[,2]) %>% exp),
  highCI = c(1, (summ_coef_prs[,1] + 1.96*summ_coef_prs[,2]) %>% exp))

print(paste0(or_prs10_df[10, 1] %>% round(2),
      "(",
      or_prs10_df[10, 2] %>% round(2),
      ",",
      or_prs10_df[10, 3] %>% round(2),
      ")"))
print(mk.test(or_prs10_df$OR, continuity = TRUE))
saveRDS(or_prs10_df, file = paste0(PROJ_PATH, "03_result/09_prs/PRS/PRS_output/or_prs10_df.rds"))

####
library(survival)
library(survminer)
nerv_surv <- Surv(comb_df$Time, comb_df$OS)
fomu_surv <- paste0("nerv_surv ~ Age + Sex + BMI_g + ", 
                    paste(paste0("PC", 1:20), collapse = " + "),
                    " + prs_g10")
cox_prs10 <- coxph(as.formula(fomu_surv), data = comb_df)
cox_coef_prs <- summary(cox_prs10)$conf.int[paste0("prs_g10", 2:10), c(1, 3, 4)]
surv_rr_prs10_df <- rbind(c(1, 1, 1), cox_coef_prs) %>% as.data.frame()
colnames(surv_rr_prs10_df) <- c("RR", "lowCI", "highCI")
saveRDS(surv_rr_prs10_df, file = paste0(PROJ_PATH, "03_result/09_prs/PRS/PRS_output/surv_rr_prs10_df.rds"))

####### OA by PRS #######
fomux <- paste0("nerv ~ Age + Sex + BMI_g + ", 
                paste(paste0("PC", 1:20), collapse = " + "),
                " + oa")
oa_sub_prs <- lapply(1:3, function(x){
  # use 181
  comb_dfx <- subset(comb_df, comb_df$prs_g34 == x)
  logit_oax <- glm(fomux, data = comb_dfx, family = "binomial")
  sum_oa <- summary(logit_oax)$coefficients["oa",]
  or <- sum_oa[1] %>% exp %>% round(2)
  or_ci_lower <- (sum_oa[1] - 1.96*sum_oa[2]) %>% exp %>% round(2)
  or_ci_upper <- (sum_oa[1] + 1.96*sum_oa[2]) %>% exp %>% round(2)
  p_val <- sum_oa[4] 
  #
  res <- data.frame(OR = or,
                    LowCI = or_ci_lower,
                    HighCI = or_ci_upper,
                    P_val = p_val,
                    n = sum(comb_dfx$nerv == 1))
}) %>% Reduce("rbind", .)

saveRDS(oa_sub_prs, 
        file = paste0(PROJ_PATH, "03_result/09_prs/PRS/PRS_output/oa_sub_prs262.rds"))

####### Compare model #######
## 0. split data
train_df <- comb_df[comb_df$group == "train",]
test_df <- comb_df[comb_df$group == "test",]

use_y <- "nerv"

## 1. base model
# 1.1 fit model
fomu_test1 <- paste0(use_y, " ~ BMI_g + Sex + Age")
logit_cov_coef <- glm(fomu_test1, 
                      data = train_df, family = "binomial") %>% coef

# 1.2 test model
test_cov_df <- cbind(1, test_df[, c("BMI_g", "Sex", "Age")])
test_cov <- as.matrix(test_cov_df) %*% as.matrix(logit_cov_coef) 
roc_cov <- pROC::roc(test_df[[use_y]], test_cov[, 1])

### 2. base + oa
# 2.1 fit model
fomu_test2 <- paste0(use_y, " ~ BMI_g + Sex + Age + oa")
logit_cov_oa_coef <- glm(fomu_test2, 
                         data = train_df, family = "binomial") %>% coef

# 2.2 test model
test_cov_oa_df <- cbind(1, test_df[, c("BMI_g", "Sex", "Age", "oa")])
test_cov_oa <- as.matrix(test_cov_oa_df) %*% as.matrix(logit_cov_oa_coef) 
roc_cov_oa <- pROC::roc(test_df[[use_y]], test_cov_oa[, 1])

### 3. base + oa + prs
# 3.1 fit model
fomu_test3 <- paste0(use_y, " ~ BMI_g + Sex + Age + ",
                     paste(paste0("PC", 1:20), collapse = " + "))
logit_oa_prs_coef <- glm(fomu_test3, 
                         data = train_df, family = "binomial") %>% coef

# 3.2 test model
test_cov_oa_pra_df <- cbind(1, test_df[, c("BMI_g", "Sex", "Age", 
                                           paste0("PC", 1:20))])
test_cov_oa_prs <- as.matrix(test_cov_oa_pra_df) %*% as.matrix(logit_oa_prs_coef) 

roc_cov_oa_prs <- pROC::roc(test_df[[use_y]], 
                            test_cov_oa_prs[, 1] + test_df$prs)

saveRDS(list("roc_cov" = roc_cov,
             "roc_cov_oa" = roc_cov_oa,
             "roc_cov_oa_prs" = roc_cov_oa_prs),
        file = paste0(PROJ_PATH, "03_result/09_prs/PRS/PRS_output/roc_list.rds"))

####### test in AFR #######
afr_test_df <- readRDS(paste0(PROJ_PATH, "03_result/09_prs/PRS/PRS_output/comb_df_afr.rds"))

# 1. test model 1
afr_test_cov_df <- cbind(1, afr_test_df[, c("BMI_g", "Sex", "Age")])
afr_test_cov <- as.matrix(afr_test_cov_df) %*% as.matrix(logit_cov_coef) 
afr_roc_cov <- pROC::roc(afr_test_df[[use_y]], afr_test_cov[, 1])

# 2. test model 2
afr_test_cov_oa_df <- cbind(1, afr_test_df[, c("BMI_g", "Sex", "Age", "oa")])
afr_test_cov_oa <- as.matrix(afr_test_cov_oa_df) %*% as.matrix(logit_cov_oa_coef) 
afr_roc_cov_oa <- pROC::roc(afr_test_df[[use_y]], afr_test_cov_oa[, 1])

# 3. test model 3
afr_test_cov_oa_pra_df <- cbind(1, afr_test_df[, c("BMI_g", "Sex", "Age", 
                                                   paste0("PC", 1:20))])
afr_test_cov_oa_prs <- as.matrix(afr_test_cov_oa_pra_df) %*% as.matrix(logit_oa_prs_coef) 

afr_roc_cov_oa_prs <- pROC::roc(afr_test_df[[use_y]], 
                                afr_test_cov_oa_prs[, 1] + afr_test_df$prs)

saveRDS(list("afr_roc_cov" = afr_roc_cov,
             "afr_roc_cov_oa" = afr_roc_cov_oa,
             "afr_roc_cov_oa_prs" = afr_roc_cov_oa_prs),
        file = paste0(PROJ_PATH, "03_result/09_prs/PRS/PRS_output/roc_list_afr.rds"))

