rm(list = ls())
# load packages
library(bigreadr)
library(dplyr)
library(plyr)
library(magrittr)

# Set parameters
SAMPLE_SIZE_EUR <- 310848
SAMPLE_SIZE_ASA <- 3538
SAMPLE_SIZE_AFR <- 6442

DATE_STAMP <- as.Date("1/1/2023", format = "%d/%m/%Y")
PHENO_FILE <- "/public/home/Datasets/ukb/ukb47503.csv.gz"
PROJ_PATH <- "/public/home/biostat03/project/oaProject/"

# Function 1: disease identification
disease.identify <- function(icd10_code, ancestry){
  
  #
  if (!all(nchar(icd10_code) == nchar(icd10_code)[1])) {
    stop("Input ICD10 codes not in the same level!")
  }
  
  pheno_str <- paste0("/public/home/Datasets/ukb/pheno/illness_mat_", ancestry, ".txt")
  ICD10_code <- c(paste0("41202-0.", 0:74),          # Diagnoses - main ICD10
                  paste0("41270-0.",0:222))          # Diagnoses - ICD10
  ICD10_diag_date <- c(paste0("41262-0.", 0:74),     # Date of first in-patient diagnosis - main ICD10
                       paste0("41280-0.",0:222))     # Date of first in-patient diagnosis - ICD10
  icd10_df <- fread2(pheno_str, select = c(ICD10_code, "eid"))
  icd10_date_df <- fread2(pheno_str, select = c(ICD10_diag_date, "eid"))
  
  disease_list <- plyr::alply(c(1: length(ICD10_code)), 1, function(x) {
    
    idx <- which(substr(icd10_df[, x], 1, nchar(icd10_code[1])) %in% c(icd10_code))
    diag <- data.frame(icd10_date_df$eid[idx], 
                       icd10_date_df[idx, x])
    return(diag)
  })
  disease_df <- plyr::ldply(disease_list, function(x) x)
  colnames(disease_df) <- c("X1", "eid", "date")
  disease_idx_uni <- unique(disease_df$eid)
  disease_df_uni <- plyr::adply(disease_idx_uni, 1, function(ii){
    
    tmp <- disease_df[disease_df$eid == ii, ]
    res <- tmp[which.min(tmp$date), c(2, 3)]
    return(res)
  })
  message(paste0("Number of cases: ", nrow(disease_df_uni)))
  message(paste0("Prevelence: ", round(nrow(disease_df_uni)/nrow(icd10_df), 4)))
  return(disease_df_uni)
}


# Function 2: disease identification with data stamp
disease.identify.stamp <- function(icd10_code, ancestry, STAMP){
  
  ## define sample size
  if (ancestry == "EUR"){
    
    assign("sample_size", SAMPLE_SIZE_EUR)
  }
  if (ancestry == "AFR"){
    
    assign("sample_size", SAMPLE_SIZE_AFR)
  }
  if (ancestry == "ASA"){
    
    assign("sample_size", SAMPLE_SIZE_ASA)
  }
  ## identify disease with date
  dis_dat <- disease.identify(icd10_code, ancestry) 
  eth_id <- fread2(paste0("/public/home/Datasets/ukb/pheno/eid_", ancestry, ".txt"))[,1]
  prev_label <- ifelse(eth_id %in% dis_dat$eid, 1, 0)
  incid_label <- ifelse(eth_id %in% dis_dat$eid[dis_dat$date > STAMP], 1, 0)
  # prev_label <- ifelse(c(1: sample_size) %in% dis_dat$idx, 1, 0)
  # incid_label <- ifelse(c(1: sample_size) %in% dis_dat$idx[dis_dat[, 3] > STAMP], 1, 0)
  if(all(incid_label == 0)){
    
    dis_label <- cbind.data.frame(eid = eth_id, prev_label)
    dis_label$prev_time <- dis_dat[match(dis_label$eid, dis_dat$eid), "date"]
    dis_label$prev_time[is.na(dis_label$prev_time)] <- STAMP
    # dis_label$prev_time[dis_label$prev_label == 1] <- dis_dat[dis_dat$idx %in% dis_label[dis_label$prev_label == 1, 1], 3]
  } else {
    
    dis_label <- cbind.data.frame(idx = c(1: sample_size), prev_label, incid_label)
    dis_label$prev_time <- STAMP
    dis_label$incid_time <- STAMP
    dis_label$prev_time[dis_label$prev_label == 1] <- dis_dat[dis_dat$idx %in% dis_label[dis_label$prev_label == 1, 1], 3]
    dis_label$incid_time[dis_label$incid_label == 1] <- dis_dat[dis_dat$idx %in% dis_label[dis_label$incid_label == 1, 1], 3]
    message(paste0("Number of case before DATE_STAMP: ", sum(is.na(incid_label))))
    message(paste0("Incidence: ", round(sum(incid_label == 1, na.rm = T)/(sample_size - sum(is.na(incid_label))), 4)))
  }
  
  return(dis_label)
}

# Load data
load(paste0(PROJ_PATH, "03_result/00_cohort/sample_info/cov_dat_all.RData"))
eur_id <- fread2(paste0("/public/home/Datasets/ukb/pheno/eid_EUR.txt"))
cov_dat <- cov_dat_all %>% filter(eid %in% eur_id[["V1"]])

############################
##### STEP 1: Define diseases and process variables
############################
## Osteoporosis
oa_eth_dat2 <- disease.identify.stamp(icd10_code = paste0("M8", 0:2),
                                     ancestry = "EUR",
                                     STAMP = DATE_STAMP)
## Alzheimer Disease 
nerv_dis <- data.frame(label = c("Alzheimer's disease"), 
                       icd10 = c("G30"))
nn = 1
nerv_eth_dat2 <- disease.identify.stamp(icd10_code = nerv_dis[nn, 2],
                                       ancestry = "EUR",
                                       STAMP = DATE_STAMP)
co_dis_df2 <- data.frame(eid = oa_eth_dat2[, 1],
                         oa = oa_eth_dat2[, 2],
                         oa_date = oa_eth_dat2[, 3],
                         nerv = nerv_eth_dat2[, 2],
                         nerv_date = nerv_eth_dat2[, 3])
## Delete samples suffered from AD first
del_cnd <- co_dis_df2$oa == 1 & co_dis_df2$nerv == 1 & co_dis_df2$oa_date > co_dis_df2$nerv_date
cat(sum(del_cnd), " should be deleted.\n")
if(sum(del_cnd) != 0){
  co_dis_df2 <- co_dis_df2[!del_cnd, ]
}
save(co_dis_df, file = paste0(PROJ_PATH, "03_result/00_cohort/sample_info/co_dis_df0830.RData"))

## Covariables
cov_dat1 <- cov_dat[match(co_dis_df2$eid, cov_dat$eid), ]
tot_df <- data.frame(nerv = co_dis_df2$nerv, 
                     BMI_g = cov_dat1$BMI_g,
                     Sex = cov_dat1$Sex,
                     Age = cov_dat1$Age,
                     oa = co_dis_df2$oa,
                     Diabetes = cov_dat1$Diabetes,
                     CVD = cov_dat1$CVD, 
                     Hearing = cov_dat1$Hearing, 
                     Fracture = cov_dat1$Fracture,
                     Activity = cov_dat1$Activity,
                     Dyslipidemia = cov_dat1$Dyslipidemia,
                     Income1 = as.factor(cov_dat1$Income1),
                     Education2 = as.factor(cov_dat1$Education2),
                     Employment = cov_dat1$Employment,
                     Nosmoke = cov_dat1$Nosmoke, 
                     Noalcohol = cov_dat1$Noalcohol,
                     sigma4 = as.factor(cov_dat1$sigma4), 
                     eid = cov_dat1$eid)

## Alzheimer Disease proxy
AD_proxy_df <- readRDS(paste0(PROJ_PATH, "03_result/00_cohort/sample_info/AD_proxy.rds"))

tot_df$AD_proxy <- AD_proxy_df$AD_proxy[match(tot_df$eid, AD_proxy_df$eid)]
tot_df$AD_proxy_AD <- tot_df$AD_proxy | tot_df$nerv
tot_df$AD_proxy_AD2 <- !tot_df$AD_proxy & tot_df$nerv
tot_df$AD_proxy_AD3 <- tot_df$AD_proxy & tot_df$nerv
tot_df$AD_proxy_AD4 <- tot_df$AD_proxy & !tot_df$nerv

save(tot_df, 
     file = paste0(PROJ_PATH, "03_result/00_cohort/sample_info/tot_dat_all.RData"))

############################
##### STEP 2: Main Analysis
############################
## univariate regression
load(paste0(PROJ_PATH, "03_result/00_cohort/sample_info/tot_dat_all0830.RData"))
tb_res <- alply(c(2:17)[-3], 1, function(ii){
  
  cnt <- with(tot_df, table(tot_df[, ii], nerv))
  cnt_prop <- with(tot_df, prop.table(table(tot_df[, ii], nerv), margin = 2))
  res_format <- data.frame(var = colnames(tot_df)[ii],
                           col1 = paste0(cnt[, 1], " (", round(cnt_prop[, 1]*100, 2), ")"),
                           col2 = paste0(cnt[, 2], " (", round(cnt_prop[, 2]*100, 2), ")"))
  return(res_format)
}) %>% do.call("rbind", .)
write.csv(tb_res, quote = F, row.names = F, 
          file = paste0(PROJ_PATH, "03_result/00_cohort/baseline_EUR.csv"))

aggregate(tot_df$Age, by = list(tot_df$nerv), mean)
aggregate(tot_df$Age, by = list(tot_df$nerv), sd)

sv_res <- alply(c(5:17), 1, function(nn){

  logit_mod <- glm(nerv ~ BMI_g + Sex + Age + tot_df[, nn],
                   data = tot_df, family = "binomial")
  
  or <- logit_mod$coefficients %>% exp %>% round(2)
  or_ci_lower <- (coefficients(logit_mod) - 1.96*summary(logit_mod)$coefficients[,2]) %>% 
    exp %>% round(2)
  or_ci_upper <- (coefficients(logit_mod) + 1.96*summary(logit_mod)$coefficients[,2]) %>% 
    exp %>% round(2)
  p_val <- coef(summary(logit_mod))[, 4] 
  res <- data.frame(OR_CI = paste0(or,
                                   "(",
                                   or_ci_lower,
                                   "–",
                                   or_ci_upper,
                                   ")"), 
                    P_val = p_val, 
                    Var = colnames(tot_df)[nn])
  res_format <- rbind(c("1.00", NA, colnames(tot_df)[nn]),
                      res[grep("tot_df", rownames(res)),])
  res_format <- rbind(rep(NA, 5),
                      res_format)
  return(res_format)
}) %>% do.call("rbind", .)
rownames(sv_res) <- NULL
sv_res$P_val <- as.numeric(sv_res$P_val)
write.table(sv_res, file = "sv_res.txt", sep = "\t", 
            row.names = F, col.names = T, quote = F)

## multivariable regression
sig_var <- sv_res$Var[which(sv_res$P_val < 0.1)] %>% unique()
if (length(sig_var) > 0) {
  mv_formula <- paste0("nerv ~ BMI_g + Sex + Age +",
                       paste(sig_var, collapse = " + "))
  logit_mod <-  glm(mv_formula,
                    data = tot_df, family = "binomial") 
  or <- logit_mod$coefficients %>% exp %>% round(2)
  or_ci_lower <- (coefficients(logit_mod) - 1.96*summary(logit_mod)$coefficients[,2]) %>% 
    exp %>% round(2)
  or_ci_upper <- (coefficients(logit_mod) + 1.96*summary(logit_mod)$coefficients[,2]) %>% 
    exp %>% round(2)
  p_val <- coef(summary(logit_mod))[, 4] 
  multi_res <- data.frame(OR_CI = paste0(or,
                                      "(",
                                      or_ci_lower,
                                      "~",
                                      or_ci_upper,
                                      ")"), 
                          P_val = p_val, 
                          Var = names(p_val))
  
  rownames(multi_res) <-NULL
} else {
  multi_res <- NA
}
nerv_result <- list("sv_res" = sv_res, 
                    "multi_res" = multi_res,
                    label = nerv_dis$label[nn])
save(nerv_result, file = paste0(PROJ_PATH, "03_result/nerv_result_main", eth, ".rda"))


############################
##### STEP 3: Sensitive Analysis
############################
###### sensitive analysis 1: + APOE ###### 
### univariate regression
sv_res <- alply(c(5: 16), 1, function(nn){
  
  #
  desp <- table(tot_df[, nn], tot_df$nerv) %>% as.data.frame()
  logit_mod <- glm(nerv ~ BMI_g + Sex + Age + sigma4 + tot_df[, nn],
                   data = tot_df, family = "binomial")
  
  or <- logit_mod$coefficients %>% exp %>% round(2)
  or_ci_lower <- (coefficients(logit_mod) - 1.96*summary(logit_mod)$coefficients[,2]) %>% 
    exp %>% round(2)
  or_ci_upper <- (coefficients(logit_mod) + 1.96*summary(logit_mod)$coefficients[,2]) %>% 
    exp %>% round(2)
  p_val <- coef(summary(logit_mod))[, 4] 
  res <- data.frame(OR_CI = paste0(or,
                                   "(",
                                   or_ci_lower,
                                   "~",
                                   or_ci_upper,
                                   ")"), 
                    P_val = p_val, 
                    Var = colnames(tot_df)[nn])
  res_format <- rbind(c("1.00", NA, colnames(tot_df)[nn]),
                      res[grep("tot_df", rownames(res)),])
  res_format <- cbind(matrix(desp[,3], ncol = 2),
                      res_format)
  res_format <- rbind(rep(NA, 5),
                      res_format)
  return(res_format)
}) %>% do.call("rbind", .)
rownames(sv_res) <- NULL
sv_res$P_val <- as.numeric(sv_res$P_val)
write.table(sv_res, file = "sv_res_sens1.txt", 
            sep = "\t", row.names = F, col.names = T, quote = F)

### multivariable regression
sig_var <- sv_res$Var[which(sv_res$P_val < 0.1)] 
if (length(sig_var) > 0) {
  mv_formula <- paste0("nerv ~ BMI_g + Sex + Age + sigma4 + ",
                       paste(unique(sig_var), collapse = " + "))
  logit_mod <-  glm(mv_formula,
                    data = tot_df, family = "binomial") 
  or <- logit_mod$coefficients %>% exp %>% round(2)
  or_ci_lower <- (coefficients(logit_mod) - 1.96*summary(logit_mod)$coefficients[,2]) %>% exp %>% round(2)
  or_ci_upper <- (coefficients(logit_mod) + 1.96*summary(logit_mod)$coefficients[,2]) %>% exp %>% round(2)
  p_val <- coef(summary(logit_mod))[, 4] 
  multi_res <- data.frame(OR_CI = paste0(or,
                                         "(",
                                         or_ci_lower,
                                         "~",
                                         or_ci_upper,
                                         ")"), 
                          P_val = p_val, 
                          Var = names(p_val))
  
  rownames(multi_res) <-NULL
} else {
  multi_res <- NA
}
nerv_result <- list("sv_res" = sv_res, 
                    "multi_res" = multi_res,
                    label = nerv_dis$label[nn])
save(nerv_result, file = paste0(PROJ_PATH, "03_result/nerv_result_APOE", eth, ".rda"))

###### sensitive analysis 2: Age > 60 ###### 
### univariate regression
PROJ_PATH <- "/public/home/biostat03/project/oaProject/"
load(paste0(PROJ_PATH, "03_result/00_cohort/sample_info/tot_dat_all0830.RData"))
tot_dff <- subset(tot_df, tot_df$Age >= 60)
sv_res <- alply(c(5:17), 1, function(nn){
  
  desp <- table(tot_dff[, nn], tot_dff$nerv) %>% as.data.frame()
  logit_mod <- glm(nerv ~ BMI_g + Sex + Age + tot_dff[, nn],
                   data = tot_dff, family = "binomial")
  
  or <- logit_mod$coefficients %>% exp %>% round(2)
  or_ci_lower <- (coefficients(logit_mod) - 1.96*summary(logit_mod)$coefficients[,2]) %>% 
    exp %>% round(2)
  or_ci_upper <- (coefficients(logit_mod) + 1.96*summary(logit_mod)$coefficients[,2]) %>% 
    exp %>% round(2)
  p_val <- coef(summary(logit_mod))[, 4] 
  res <- data.frame(OR_CI = paste0(or,
                                   "(",
                                   or_ci_lower,
                                   "~",
                                   or_ci_upper,
                                   ")"), 
                    P_val = p_val, 
                    Var = colnames(tot_dff)[nn])
  res_format <- rbind(c("1.00", NA, colnames(tot_dff)[nn]),
                      res[grep("tot_dff", rownames(res)),])
  res_format <- cbind(matrix(desp[,3], ncol = 2),
                      res_format)
  res_format <- rbind(rep(NA, 5),
                      res_format)
  return(res_format)
}) %>% do.call("rbind", .)
rownames(sv_res) <- NULL
sv_res$P_val <- as.numeric(sv_res$P_val)
write.table(sv_res, file = "sv_res_60.txt", sep = "\t", row.names = F, col.names = T, quote = F)

## multivariable regression
sig_var <- sv_res$Var[which(sv_res$P_val < 0.1)] %>% unique()
if (length(sig_var) > 0) {
  mv_formula <- paste0("nerv ~ BMI_g + Sex + Age +",
                       paste(sig_var, collapse = " + "))
  logit_mod <-  glm(mv_formula,
                    data = tot_dff, family = "binomial") 
  or <- logit_mod$coefficients %>% exp %>% round(2)
  or_ci_lower <- (coefficients(logit_mod) - 1.96*summary(logit_mod)$coefficients[,2]) %>% 
    exp %>% round(2)
  or_ci_upper <- (coefficients(logit_mod) + 1.96*summary(logit_mod)$coefficients[,2]) %>% 
    exp %>% round(2)
  p_val <- coef(summary(logit_mod))[, 4] 
  multi_res <- data.frame(OR_CI = paste0(or,
                                         "(",
                                         or_ci_lower,
                                         "~",
                                         or_ci_upper,
                                         ")"), 
                          P_val = p_val, 
                          Var = names(p_val))
  rownames(multi_res) <-NULL
} else {
  multi_res <- NA
}
nerv_result <- list("sv_res" = sv_res, 
                    "multi_res" = multi_res,
                    label = nerv_dis$label[nn])
save(nerv_result, file = paste0(PROJ_PATH, "03_result/nerv_result_60", eth, ".rda"))

###### sensitive 3: 10 folds CV ###### 
PROJ_PATH <- "/public/home/biostat03/project/oaProject/"
load(paste0(PROJ_PATH, "03_result/00_cohort/sample_info/tot_dat_all0830.RData"))

set.seed(20230528)
sample_id <- sample(1:10, size = nrow(tot_df), replace = T)
nerv_result_sens <- lapply(c(1: 10), function(nn){
  
  ### univariate regression
  tot_dff <- tot_df[nn != sample_id,]
  sv_res <- alply(c(5: 17), 1, function(nn){
    
    logit_mod <- glm(nerv ~ BMI_g + Sex + Age + tot_dff[, nn],
                     data = tot_dff, family = "binomial")
    or <- logit_mod$coefficients %>% exp %>% round(2)
    or_ci_lower <- (coefficients(logit_mod) - 1.96*summary(logit_mod)$coefficients[,2]) %>% 
      exp %>% round(2)
    or_ci_upper <- (coefficients(logit_mod) + 1.96*summary(logit_mod)$coefficients[,2]) %>%
      exp %>% round(2)
    p_val <- coef(summary(logit_mod))[, 4] 
    res <- data.frame(OR_CI = paste0(or,
                                     " (",
                                     or_ci_lower,
                                     "~",
                                     or_ci_upper,
                                     ")"), 
                      P_val = p_val, 
                      Var = colnames(tot_dff)[nn])
    return(res[grep("tot_dff", rownames(res)),])
  }) %>% do.call("rbind", .)
  rownames(sv_res) <- NULL
  
  ### multivariable regression
  sig_varx <- sv_res$Var[which(sv_res$P_val < 0.1)] %>% unique()
  if (length(sig_varx) > 0) {
    mv_formula <- paste0("nerv ~ BMI_g + Sex + Age + ",
                         paste(sig_varx, collapse = " + "))
    logit_mod <-  glm(mv_formula,
                      data = tot_df, family = "binomial") 
    or <- logit_mod$coefficients %>% exp %>% round(2)
    or_ci_lower <- (coefficients(logit_mod) - 1.96*summary(logit_mod)$coefficients[,2]) %>% exp %>% round(2)
    or_ci_upper <- (coefficients(logit_mod) + 1.96*summary(logit_mod)$coefficients[,2]) %>% exp %>% round(2)
    p_val <- coef(summary(logit_mod))[, 4]
    multi_res <- data.frame(OR_CI = paste0(or,
                                           " (",
                                           or_ci_lower,
                                           "~",
                                           or_ci_upper,
                                           ")"),
                            P_val = p_val,
                            Var = names(p_val))
    
    rownames(multi_res) <-NULL
  } else {
    multi_res <- NA
  }
  return(list("sv_res" = sv_res,
              "multi_res" = multi_res))
})
names(nerv_result_sens) <- paste0("CV", 1:10)
save(nerv_result_sens, file = paste0(PROJ_PATH, "03_result/nerv_result_CV", eth, ".rda"))

###### sensitive 4: ASA/AFR ###### 
load(paste0(PROJ_PATH, "03_result/00_cohort/sample_info/cov_dat_all.RData"))
eth_list <- c("AFR", "ASA")
lapply(eth_list, function(ethx){
  
  ethx_id <- fread2(paste0("/public/home/Datasets/ukb/pheno/eid_", ethx, ".txt"))
  cov_dat <- cov_dat_all %>%
  dplyr::filter(eid %in% ethx_id[["V1"]])
  
  ## Osteoporosis
  oa_ethx_dat <- disease.identify.stamp(icd10_code = paste0("M8", 0:2),
                                        ancestry = ethx,
                                        STAMP = DATE_STAMP)
  ## Alzheimer Disease
  nerv_dis <- data.frame(label = c("Alzheimer's disease"), 
                         icd10 = c("G30"))
  nn = 1
  nerv_ethx_dat <- disease.identify.stamp(icd10_code = nerv_dis[nn, 2],
                                          ancestry = ethx,
                                          STAMP = DATE_STAMP)
  
  co_dis_df <- data.frame(eid = oa_ethx_dat[, 1],
                          oa = oa_ethx_dat[, 2],
                          oa_date = oa_ethx_dat[, 3],
                          nerv = nerv_ethx_dat[, 2],
                          nerv_date = nerv_ethx_dat[, 3])
  del_cnd <- co_dis_df$oa == 1 & co_dis_df$nerv == 1 & co_dis_df$oa_date > co_dis_df$nerv_date
  cat(sum(del_cnd), " should be deleted.\n")
  if(sum(del_cnd) != 0){
    co_dis_df <- co_dis_df[!del_cnd, ]
  }
  
  cov_dat1 <- cov_dat[match(co_dis_df$eid, cov_dat$eid), ]
  tot_df <- data.frame(nerv = co_dis_df$nerv, 
                       BMI_g = cov_dat1$BMI_g,
                       Sex = cov_dat1$Sex,
                       Age = cov_dat1$Age,
                       oa = co_dis_df$oa,
                       Diabetes = cov_dat1$Diabetes,
                       CVD = cov_dat1$CVD, 
                       Hearing = cov_dat1$Hearing, 
                       Fracture = cov_dat1$Fracture,
                       Activity = cov_dat1$Activity,
                       Dyslipidemia = cov_dat1$Dyslipidemia,
                       Income1 = as.factor(cov_dat1$Income1),
                       Education2 = as.factor(cov_dat1$Education2),
                       Employment = cov_dat1$Employment,
                       Nosmoke = cov_dat1$Nosmoke, 
                       Noalcohol = cov_dat1$Noalcohol,
                       sigma4 = as.factor(cov_dat1$sigma4),
                       eid = cov_dat1$eid)
  
  # AD_proxy_df <- readRDS(paste0(PROJ_PATH, "03_result/00_cohort/sample_info/AD_proxy.rds"))
  # 
  # tot_df$AD_proxy <- AD_proxy_df$AD_proxy[match(tot_df$eid, AD_proxy_df$eid)]
  # tot_df$AD_proxy_AD <- tot_df$AD_proxy | tot_df$nerv
  # tot_df$AD_proxy_AD2 <- !tot_df$AD_proxy & tot_df$nerv
  # tot_df$AD_proxy_AD3 <- tot_df$AD_proxy & tot_df$nerv
  # tot_df$AD_proxy_AD4 <- tot_df$AD_proxy & !tot_df$nerv
  save(tot_df, 
       file = paste0(PROJ_PATH, "03_result/00_cohort/sample_info/tot_dat_", ethx, ".RData"))
  ethx <- "ASA"
  load(paste0(PROJ_PATH, "03_result/00_cohort/sample_info/tot_dat_", ethx, ".RData"))
  dim(tot_df)
  tb_res <- alply(c(2:17)[-3], 1, function(ii){
    
    cnt <- with(tot_df, table(tot_df[, ii], nerv))
    cnt_prop <- with(tot_df, prop.table(table(tot_df[, ii], nerv), margin = 2))
    res_format <- data.frame(var = colnames(tot_df)[ii],
                             col1 = paste0(cnt[, 1], " (", round(cnt_prop[, 1]*100, 2), ")"),
                             col2 = paste0(cnt[, 2], " (", round(cnt_prop[, 2]*100, 2), ")"))
    return(res_format)
  }) %>% do.call("rbind", .)
  write.csv(tb_res, quote = F, row.names = F, 
            file = paste0(PROJ_PATH, "03_result/00_cohort/baseline_", ethx, ".csv"))
  aggregate(tot_df$Age, by = list(tot_df$nerv), mean)
  aggregate(tot_df$Age, by = list(tot_df$nerv), sd)
  
  
  ### main analysis
  ## single variable
  sv_res <- alply(c(5:17), 1, function(nn){
    
    desp <- table(tot_df[, nn], tot_df$nerv) %>% as.data.frame()
    logit_mod <- glm(nerv ~ BMI_g + Sex + Age + tot_df[, nn],
                     data = tot_df, family = "binomial")
    
    or <- logit_mod$coefficients %>% exp %>% round(2)
    or_ci_lower <- (coefficients(logit_mod) - 1.96*summary(logit_mod)$coefficients[,2]) %>% 
      exp %>% round(2)
    or_ci_upper <- (coefficients(logit_mod) + 1.96*summary(logit_mod)$coefficients[,2]) %>% 
      exp %>% round(2)
    p_val <- coef(summary(logit_mod))[, 4] 
    res <- data.frame(OR_CI = paste0(or,
                                     "(",
                                     or_ci_lower,
                                     "~",
                                     or_ci_upper,
                                     ")"), 
                      P_val = p_val, 
                      Var = colnames(tot_df)[nn])
    res_format <- rbind(c("1.00", NA, colnames(tot_df)[nn]),
                        res[grep("tot_df", rownames(res)),])
    res_format <- cbind(matrix(desp[,3], ncol = 2),
                        res_format)
    res_format <- rbind(rep(NA, 5),
                        res_format)
    return(res_format)
  }) %>% do.call("rbind", .)
  rownames(sv_res) <- NULL
  sv_res$P_val <- as.numeric(sv_res$P_val)
  write.table(sv_res, file = "sv_res.txt", sep = "\t", row.names = F, col.names = T, quote = F)
  
  ## multivariable regression
  sig_var <- sv_res$Var[which(sv_res$P_val < 0.1)] %>% unique()
  if (length(sig_var) > 0) {
    mv_formula <- paste0("nerv ~ BMI_g + Sex + Age +",
                         paste(sig_var, collapse = " + "))
    logit_mod <-  glm(mv_formula,
                      data = tot_df, family = "binomial") 
    or <- logit_mod$coefficients %>% exp %>% round(2)
    or_ci_lower <- (coefficients(logit_mod) - 1.96*summary(logit_mod)$coefficients[,2]) %>% 
      exp %>% round(2)
    or_ci_upper <- (coefficients(logit_mod) + 1.96*summary(logit_mod)$coefficients[,2]) %>% 
      exp %>% round(2)
    p_val <- coef(summary(logit_mod))[, 4] 
    multi_res <- data.frame(OR_CI = paste0(or,
                                           "(",
                                           or_ci_lower,
                                           "~",
                                           or_ci_upper,
                                           ")"), 
                            P_val = p_val, 
                            Var = names(p_val))
    
    rownames(multi_res) <-NULL
  } else {
    multi_res <- NA
  }
  return(list("sv_res" = sv_res, 
              "multi_res" = multi_res))
})
save(nerv_result, file = paste0(PROJ_PATH, "03_result/nerv_result_main", ethx, ".rda"))
