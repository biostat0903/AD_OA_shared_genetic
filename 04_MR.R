rm(list=ls())
gc()
# load packages
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ggplot2)
library(patchwork)

# Parameters
PROJ_PATH <- "/public/home/biostat03/project/oaProject/"
DATA_PROJ <- paste0(PROJ_PATH, "02_data/summary_stat/")
SUMM_PATH <- paste0(DATA_PROJ, "cleaned/")
setwd(paste0(PROJ_PATH, "03_result/03_mr/"))

# Process Summary Statistics


exp <- fread(paste0(SUMM_PATH, "OA_frontier.txt"))
exp <- exp[which(exp$P < 5E-8),]

exp_dat <- exp %>%
  format_data(
    type='exposure',
    snp_col = "SNP",
    beta_col = "Beta",
    se_col = "SE",
    eaf_col = "EAF",
    effect_allele_col ="A1",
    other_allele_col = "A2",
    pval_col = "P"
  )
exp_clump <- clump_data(exp_dat, clump_r2 = 0.01, clump_kb = 10000, pop = "EUR")

exp_clump$id.exposure <- "OA"

# out_dat <- read_outcome_data(
# snps = exp_clump$SNP,
# filename = paste0(SUMM_PATH, "OA_06.txt"),
# sep = "\t",
# snp_col = "SNP",
# beta_col = "Beta",
# se_col = "SE",
# effect_allele_col ="A1",
# other_allele_col = "A2",
# eaf_col = "EAF",
# pval_col = "P")

out_dat <- read_outcome_data(
  snps = exp_clump$SNP,
  filename = paste0(SUMM_PATH, "AD_06.txt"),
  sep = "\t",
  snp_col = "SNP",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col ="A1",
  other_allele_col = "A2",
  pval_col = "P"
)
out_dat <- get_eaf_from_1000G(out_dat, type = "outcome")
out_dat$id.outcome <- 'out'

# MR 
mydata <- harmonise_data(
  exposure_dat = exp_clump, 
  outcome_dat = out_dat,
  action= 2
)
mydata <- mydata[which(mydata$pval.outcome > 1E-5),]

res <- mr(mydata)

heterogeneity <- mr_heterogeneity(mydata)
pleiotropy_test <- mr_pleiotropy_test(mydata)

res$heterogeneity_MR <- heterogeneity$Q_pval[1]
res$heterogeneity_IVW <- heterogeneity$Q_pval[2]
res$pleiotropy <- pleiotropy_test$pval
res$outcome <- "AD"

# forest plot
res_single <- mr_singlesnp(mydata, all_method = c("mr_ivw"))
res_single$SNP[nrow(res_single)] <- "All-IVW"
p1 <- mr_forest_plot(res_single)
p1 <- p1[[1]] +
  theme_bw() +
  labs(x = "MR Effect Size for OA on AD", tag = "A") +
  guides(size = "none", colour = "none") +
  theme(legend.position = "none",
        plot.tag = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 15), 
        axis.text.x = element_text(size = 15), 
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 15), 
        axis.title.y = element_text(size = 18, face = "bold"))
p1
ggsave("p1.pdf", width = 7, height = 8)

# scatter plot
res <- res[1:3,]
res$method[3] <- "IVW"
p2 <- mr_scatter_plot(res, mydata)
p2 <- p2[[1]] +
  labs(x = "SNP Effect on OA",
       y = "SNP Effect on AD",
       color = "MR Method",
       tag = "B") +
  theme_bw() +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1)) +
  theme(legend.position = "right",
        legend.box = "vertical",
        plot.tag = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 15), 
        axis.text.x = element_text(size = 15), 
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 15), 
        axis.title.y = element_text(size = 18, face = "bold"))
p2
ggsave("p2.pdf", width = 6, height = 4)
  
# funnel plot
res_single <- mr_singlesnp(mydata, all_method = c("mr_ivw"))
p3 <- mr_funnel_plot(res_single)
p3 <- p3[[1]] +
  labs(color = "MR Method", tag = "C") +
  scale_color_discrete(labels = c("Inverse variance weighted" = "IVW")) +
  theme_bw() +
  theme(legend.position = "right",
        plot.tag = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 15), 
        axis.text.x = element_text(size = 15), 
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 15), 
        axis.title.y = element_text(size = 18, face = "bold"))
p3
ggsave("p3.pdf", width = 6, height = 4)

# Leaveoneout forestplot
res_loo <- mr_leaveoneout(mydata)
p4 <- mr_leaveoneout_plot(res_loo)
p4 <- p4[[1]] +
  theme_bw() +
  labs(x = "MR Leave-one-out Sensitivity Analysis for OA on AD",
       tag = "D") +
  guides(size = "none", colour = "none") +
  theme(legend.position = "none",
        plot.tag = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 15), 
        axis.text.x = element_text(size = 15), 
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 15), 
        axis.title.y = element_text(size = 18, face = "bold"))
p4
ggsave("p4.pdf", width = 8.5, height = 8)

# figs3
figs3 <- p1 | (p2 / p3)
figs3
ggsave("figs3.pdf", width = 13, height = 8)

# ## 反向MR
# exp <- fread(paste0(SUMM_PATH, "AD_06.txt"))
# exp <- exp[which(exp$P < 5E-8),]
# 
# exp_dat <- exp %>%
#   format_data(
#     type='exposure',
#     snp_col = "SNP",
#     beta_col = "Beta",
#     se_col = "SE",
#     eaf_col = "EAF",
#     effect_allele_col ="A1",
#     other_allele_col = "A2",
#     pval_col = "P"
#   )
# exp_dat <- get_eaf_from_1000G(exp_dat, type = "exposure")
# exp_clump <- clump_data(exp_dat, clump_r2 = 0.01, clump_kb = 10000, pop = "EUR")
# exp_clump$id.exposure <- "AD"
# 
# # out_dat <- read_outcome_data(
# # snps = exp_clump$SNP,
# # filename = paste0(SUMM_PATH, "OA_06.txt"),
# # sep = "\t",
# # snp_col = "SNP",
# # beta_col = "Beta",
# # se_col = "SE",
# # effect_allele_col ="A1",
# # other_allele_col = "A2",
# # eaf_col = "EAF",
# # pval_col = "P")
# 
# out_dat <- read_outcome_data(
#   snps = exp_clump$SNP,
#   filename = paste0(SUMM_PATH, "OA_frontier.txt"),
#   sep = "\t",
#   snp_col = "SNP",
#   beta_col = "Beta",
#   se_col = "SE",
#   effect_allele_col ="A1",
#   other_allele_col = "A2",
#   eaf_col = "EAF",
#   pval_col = "P"
# )
# out_dat$id.outcome <- 'out'
# 
# # MR 
# mydata <- harmonise_data(
#   exposure_dat = exp_clump, 
#   outcome_dat = out_dat,
#   action= 2
# )
# mydata <- mydata[which(mydata$pval.outcome > 1E-5),]
# 
# res <- mr(mydata)
# 
# heterogeneity <- mr_heterogeneity(mydata)
# pleiotropy_test <- mr_pleiotropy_test(mydata)
# 
# res$heterogeneity_MR <- heterogeneity$Q_pval[1]
# res$heterogeneity_IVW <- heterogeneity$Q_pval[2]
# res$pleiotropy <- pleiotropy_test$pval
# res$outcome <- "OA"
# 
# res_loo <- mr_leaveoneout(mydata)
# mr_leaveoneout_plot(res_loo)



