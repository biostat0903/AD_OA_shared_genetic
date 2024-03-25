rm(list = ls())
# load packages
library(bigreadr)
library(plyr)
library(magrittr)

# Parameters
DATE_STAMP <- as.Date("1/1/2023", format = "%d/%m/%Y")
PHENO_FILE <- "/public/home/Datasets/ukb/ukb47503.csv.gz"
PROJ_PATH <- "/public/home/biostat03/project/oaProject/"

tab_var <- c("eid", 
             "738-0.0",               ## Average total household income before tax
             paste0("6138-0.", 0:5),  ## Education
             paste0("6142-0.", 0:6),  ## Employment status
             "20116-0.0",             ## Smoking status
             "884-0.0",               ## Number of days/week of moderate physical activity 10+ minutes
             "894-0.0",               ## Duration of moderate activity
             "904-0.0",               ## Number of days/week of vigorous physical activity 10+ minutes
             "914-0.0",               ## Duration of vigorous activity
             "20117-0.0",             ## Alcohol 
             "2443-0.0",              ## Diabetes diagnosed by doctor
             "20118-0.0",             ## Home area population density - urban or rural
             "30690-0.0",             ## cholesterol
             "30760-0.0",             ## HDL cholesterol (HDL)
             "30780-0.0",             ## LDL direct (LDL)
             "30870-0.0",             ## Triglycerides (TC)
             "2247-0.0",              ## Hearing difficulty/problems
             "2463-0.0",              ## Fractured/broken bones in last 5 years
             "100240-0.0",            ## Coffee consumed
             "6150-0.0",              ## Vascular/heart problems diagnosed by doctor
             "21001-0.0",             ## BMI
             "22001-0.0",             ## Genetic sex
             "31-0.0",                ## Sex
             "21003-0.0",             ## age when attended assessment centre
             "54-0.0"                 ## UK Biobank assessment centre
)  %>%
  fread2("/public/home/Datasets/ukb/ukb47503.csv.gz", 
         select = .)

################## data QC ##################
cnd_i <- tab_var[["22001-0.0"]] == tab_var[["31-0.0"]] & !is.na(tab_var[["22001-0.0"]])
cnd_ii <- tab_var$eid > 0 & !is.na(tab_var$eid)
qc_cnd <- cnd_i & cnd_ii

################## BMI group ##################
BMI_grp = ifelse(is.na(tab_var[["21001-0.0"]]), NA,
                 ifelse(tab_var[["21001-0.0"]] > 18.5 & 
                          tab_var[["21001-0.0"]]  < 30, 1, 
                        ifelse(tab_var[["21001-0.0"]] <= 18.5, 0, 
                               ifelse(tab_var[["21001-0.0"]] >= 30, 2, NA))))

################## Income group ##################
Income1 <- ifelse(is.na(tab_var[["738-0.0"]]), NA,
                  ifelse(tab_var[["738-0.0"]] > 0,
                         tab_var[["738-0.0"]], NA))
Income2 <- ifelse(is.na(Income1), NA,
                  ifelse(Income1 == 1, 1,
                         ifelse(Income1 %in% c(2, 3),2 ,3)))

################## education ##################
# values increase from 6138-0.0 to 6138-0.5
# take 6138-0.0 as the highest education status
# 1-2: College or above
# 3-6: High school or equivalent
# -7: Less than high school (33853828)
edu_cnd1 <- ifelse(is.na(tab_var[["6138-0.0"]]), NA,
                   ifelse(tab_var[["6138-0.0"]] == -3, NA,
                          ifelse(tab_var[["6138-0.0"]] == -7, 7, 
                                 tab_var[["6138-0.0"]])))
edu_cnd2 <- ifelse(is.na(tab_var[["6138-0.0"]]), NA,
                   ifelse(tab_var[["6138-0.0"]] %in% c(1:2), 3,
                          ifelse(tab_var[["6138-0.0"]] > 2, 2, 
                                 ifelse(tab_var[["6138-0.0"]] == -7, 1, NA))))

################## Employment ##################
# values increase from 6142-0.0 to 6142-0.6 
# take 6142-0.0 as the highest education status
# Employed Status including those in paid employment or self-employed, retired, 
# doing unpaid or voluntary work, or being full or part time students)
emp_cnd <- ifelse(is.na(tab_var[["6142-0.0"]]), NA,
                  ifelse(tab_var[["6142-0.0"]] %in% c(1, 2, 6, 7), 2,
                         ifelse(tab_var[["6142-0.0"]] %in% c(3:5, -7), 1, NA)))

################## no current smoke ##################
nosmoke_cnd <- ifelse(is.na(tab_var[["20116-0.0"]]), NA,
                      ifelse(tab_var[["20116-0.0"]] == 0, 1,
                             ifelse(tab_var[["20116-0.0"]] > 0, 0, NA)))

################## alcohol ##################
# never alcohol use
noalcohol_cnd <- ifelse(is.na(tab_var[["20117-0.0"]]), NA, 
                        ifelse(tab_var[["20117-0.0"]] == 0, 1,
                               ifelse(tab_var[["20117-0.0"]] > 0, 0, NA)))

################## activity ##################
num_moderate <- ifelse(tab_var[["884-0.0"]] < 0, NA, 
                       tab_var[["884-0.0"]])
num_vigorous <- ifelse(tab_var[["904-0.0"]] < 0, NA, 
                       tab_var[["904-0.0"]])
num_activity <- num_moderate > 5 & num_vigorous > 1
time_moderate <- ifelse(tab_var[["894-0.0"]] > 150, T, F)
time_vigorous <- ifelse(tab_var[["914-0.0"]] > 75, T, F)
act_mat <- cbind(num_activity, time_moderate, time_vigorous)
activity_cnd <- rep(0, nrow(act_mat))
activity_cnd[rowSums(act_mat, na.rm = T) > 0] <- 1
activity_cnd[rowSums(is.na(act_mat)) == 3] <- NA

################## Diabetes ##################
diabetes_cnd <- rep(0, nrow(tab_var))
diabetes_cnd[tab_var[["2443-0.0"]] == 1] <- 1
diabetes_cnd[tab_var[["2443-0.0"]] %in% c(-1,-3)] <- NA

################## CVD ##################
cvd_cnd <- rep(0, nrow(tab_var))
cvd_cnd[tab_var[["6150-0.0"]] %in% c(1, 3, 4)] <- 1
cvd_cnd[tab_var[["6150-0.0"]] %in% c(-7,-3)] <- NA

################## hearing ##################
hearing_cnd <- rep(0, nrow(tab_var))
hearing_cnd[tab_var[["2247-0.0"]] %in% c(1, 99)] <- 1
hearing_cnd[tab_var[["2247-0.0"]] %in% c(-1,-3)] <- NA

################## fracture ##################
fracture_cnd<- rep(0, nrow(tab_var))
fracture_cnd[tab_var[["2463-0.0"]] %in% 1] <- 1
fracture_cnd[tab_var[["2463-0.0"]] %in% c(-1,-3)] <- NA

################## dyslipidemia ##################
# "30690-0.0",             ## cholesterol  >= 6.2
# "30760-0.0",             ## HDL cholesterol (HDL) < 1.0
# "30780-0.0",             ## LDL direct (LDL)  >= 4.1
# "30870-0.0",             ## Triglycerides (TC)  >= 2.3
dyslipidemia_cnd <- ifelse(is.na(tab_var[["30690-0.0"]]) & is.na(tab_var[["30760-0.0"]]) & 
                             is.na(tab_var[["30780-0.0"]]) & is.na(tab_var[["30870-0.0"]]), NA,
                           ifelse(tab_var[["30690-0.0"]] >= 6.2, 1,
                                  ifelse(tab_var[["30760-0.0"]] < 1.0, 1, 
                                         ifelse(tab_var[["30780-0.0"]] >= 4.1, 1,
                                                ifelse(tab_var[["30870-0.0"]] >= 2.3, 1,
                                                       0)))))

################## APOE ##################
apoe_dat <- fread2(paste0(PROJ_PATH, "03_result/00_cohort/apoe_genotype.raw"))
apoe_dat$rs429358_C <- ifelse(apoe_dat$rs429358_C == 0,
                              "TT", ifelse(apoe_dat$rs429358_C == 1, "CT", "CC"))
apoe_dat$rs7412_T <- ifelse(apoe_dat$rs7412_T == 0,
                            "CC", ifelse(apoe_dat$rs7412_T == 1, "CT", "TT"))
# define APOE haplotype (PMID 29745836)
# rs429358_C / rs7412_C as sigma4
# rs429358_T / rs7412_C as sigma3
# rs429358_T / rs7412_T as sigma2
# rs429358_C / rs7412_T as sigma1 (rare, omitted)
# define APOE genotype (PMID 37062296)
# default Zero sigma4
apoe_dat$sigma4 <- 0
# Two sigma4: C/C+C/C -> CC+CC
apoe_dat$sigma4[apoe_dat$rs429358_C == "CC" &
                  apoe_dat$rs7412_T == "CC" ] <- 2
# One sigma4: sigma3/sigma4 -> T/C+C/C -> CT+CC
# One sigma4: sigma4/sigma3 -> C/C+T/C -> CT+CC
apoe_dat$sigma4[apoe_dat$rs429358_C == "CT" &
                  apoe_dat$rs7412_T == "CC" ] <- 1
# One sigma4: sigma2/sigma4 -> T/T+C/C -> CT+CT
# One sigma4: sigma4/sigma2 -> C/C+T/T -> CT+CT
# special case: sigma1/sigma3 -> C/T+T/C -> CT+CT (rare, omitted)
apoe_dat$sigma4[apoe_dat$rs429358_C == "CT" &
                  apoe_dat$rs7412_T == "CT" ] <- 1
# binary definition
apoe_dat$sigma4_b <- ifelse(apoe_dat$sigma4 == 0, 0, 1)
# old definition
apoe_dat$sigma4_old <- ifelse(apoe_dat$rs429358_C != "TT" & 
                                apoe_dat$rs7412_T == "CC", 
                              1, 0)

################## Output ##################
cov_dat_all <- data.frame(eid = tab_var[["eid"]], 
                          QC = qc_cnd, 
                          Age = tab_var[["21003-0.0"]], 
                          BMI = tab_var[["21001-0.0"]],
                          BMI_g = BMI_grp, 
                          Sex = tab_var[["22001-0.0"]],
                          Income1 = Income1,
                          Income2 = Income2,
                          Education1 = edu_cnd1,
                          Education2 = edu_cnd2,
                          Employment = emp_cnd,
                          Nosmoke = nosmoke_cnd, 
                          Activity = activity_cnd, 
                          Noalcohol = noalcohol_cnd,
                          Diabetes = diabetes_cnd,
                          CVD = cvd_cnd, 
                          Hearing = hearing_cnd,
                          Dyslipidemia = dyslipidemia_cnd,
                          Fracture = fracture_cnd)
cov_dat_all$sigma4 <- apoe_dat$sigma4[match(cov_dat_all$eid, apoe_dat$IID)]
cov_dat_all$sigma4_b <- apoe_dat$sigma4_b[match(cov_dat_all$eid, apoe_dat$IID)]
cov_dat_all$sigma4_old <- apoe_dat$sigma4_old[match(cov_dat_all$eid, apoe_dat$IID)]
save(cov_dat_all, 
     file = paste0(PROJ_PATH, "03_result/00_cohort/sample_info/cov_dat_all.RData"))

###############
#
library(bigreadr)
library(dplyr)
library(stringr)

#
PHENO_FILE <- "/public/home/Datasets/ukb/ukb47503.csv.gz"
PROJ_PATH <- "/public/home/biostat03/project/oaProject/"

#
parent_histry_code <- c(paste0("20107-", rep(0:3, each = 10), ".", rep(0:9, times = 4)),
                        paste0("20110-", rep(0:3, each = 11), ".", rep(0:10, times = 4)))
parent_histry_df <- fread2(PHENO_FILE)

parent_histry_df[is.na(parent_histry_df)] <- -999
AD_proxy_df <- lapply(2:ncol(parent_histry_df), function(x){
  ifelse(parent_histry_df[,x] == 10, 1, 0)
}) %>% Reduce("cbind", .)

AD_proxy <- data.frame(eid = parent_histry_df$eid,
                       AD_proxy = ifelse(rowSums(AD_proxy_df) > 0, 1, 0))
saveRDS(AD_proxy, file = paste0(PROJ_PATH, "03_result/00_cohort/sample_info/AD_proxy.rds"))


###############
## define last date with ICD recording
ICD10_diag_date <- c(paste0("41262-0.", 0:74),     # Date of first in-patient diagnosis - main ICD10
                     paste0("41280-0.",0:222))     # Date of first in-patient diagnosis - ICD10

all_eth <- c("EUR", "AFR", "ASA")
last_date_df <- lapply(all_eth, function(ancestryx){
  pheno_str <- paste0("/public/home/Datasets/ukb/pheno/illness_mat_", ancestryx, ".txt")
  icd10_date_df <- fread2(pheno_str, select = c(ICD10_diag_date, "eid"))
  
  last_date <- apply(icd10_date_df[,ICD10_diag_date], 1, function(x){
    x <- x[!is.na(x)]
    if (length(x) > 0) {
      return(max(x))
    } else {
      return(NA)
    }
  })
  return(data.frame(eid = icd10_date_df$eid,
                    last_date = last_date))
}) %>% Reduce("rbind", .)

saveRDS(last_date_df, file = paste0(PROJ_PATH, "03_result/00_cohort/sample_info/last_date_df.rds"))
