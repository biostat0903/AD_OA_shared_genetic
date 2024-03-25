
#
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(scales)
###
setwd("/public/home/biostat03/project/oaProject/03_result/09_prs/PRS/PRS_output/")

## Fig A
or_prs10_df <- readRDS("or_prs10_df.rds")
colnames(or_prs10_df) <- c("OR", "LowCI", "HighCI")
or_prs10_df$Group <- seq(10, 100, 10)
plt_A <- ggplot(or_prs10_df, aes(x = Group, y = OR))+
  geom_hline(yintercept = 1, 
             color = "grey50", 
             linetype = 2, 
             size = 1)+
  geom_point() + 
  geom_pointrange(aes(ymin = LowCI, ymax = HighCI),
                  stat='identity',
                  position = position_dodge(0.4),
                  size = 1)+
  scale_x_continuous(breaks = seq(0, 100, 20))+
  stat_smooth(method = "lm",
              se = F, linetype = 2, color = "red") +
  xlab("PRS(%)") + ylab("RR(95% CI)")+
  theme_bw() + 
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, face = "bold", color = "black"))


## Fig B
roc_list <- readRDS("roc_list.rds")
names(roc_list) <- c("Covariates",
                     "Covariates + OA",
                     "Covariates + OA + PRS")
auc_text <- lapply(roc_list, function(x){
  aucix <- pROC::ci(x) %>% round(4)
  paste0("AUC: ", 
         aucix[2] , "(",
         aucix[1], "~", 
         aucix[3], ")")}) %>% unlist
auc_text_df <- data.frame(x = rep(0.7, 3),
                          y = c(0.3, 0.25, 0.2),
                          auc_text)

roc_df <- lapply(names(roc_list), function(x){
  data.frame(Model = x,
             rev_specificity = 1-roc_list[[x]]$specificities,
             sensitivity = roc_list[[x]]$sensitivities)
}) %>% Reduce("rbind", .)

plt_B <- ggplot(roc_df) +
  geom_line(aes(x = rev_specificity, y = sensitivity, 
                color = Model, linetype = Model), 
            size = 1.3) +
  annotate("text", 
           x = auc_text_df[1, 1] , 
           y = auc_text_df[1, 2], 
           label = auc_text_df[1, 3],
           colour="#E41A1C") + 
  annotate("text", 
           x = auc_text_df[2, 1] , 
           y = auc_text_df[2, 2], 
           label = auc_text_df[2, 3],
           colour="#377EB8") + 
  annotate("text", 
           x = auc_text_df[3, 1] , 
           y = auc_text_df[3, 2], 
           label = auc_text_df[3, 3],
           colour="#4DAF4A") + 
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A")) + 
  xlab("1-Specificity") +
  ylab("Sensitivity") +
  theme_bw() + 
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"))

## Fig C
oa_sub_prs <- readRDS("oa_sub_prs181.rds")
colnames(oa_sub_prs) <- c("OR", "LowCI", "HighCI")
oa_sub_prs$Group <- c("Low", "Medium", "High") %>%
  factor(., levels = c("Low", "Medium", "High"))
plt_C <- ggplot(oa_sub_prs, aes(x = Group,y = OR))+
  geom_bar(aes(color = Group),
           stat = "identity", width = 0.5, size = 1,
           fill = NA)+
  scale_color_manual(values = brewer.pal(6, "Set3")[4:6]) + 
  geom_errorbar(aes(ymin = LowCI, ymax = HighCI),
                size = 0.8, width = 0.2, color = "grey10")+
  geom_hline(yintercept = 1, 
             color = "grey50", 
             linetype = 2, 
             size = 1)+
  xlab("PRS Group") + ylab("RR(95% CI)")+
  theme_bw() + 
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, face = "bold", color = "black"))


## Fig D
roc_list_afr <- readRDS("roc_list_afr.rds")
names(roc_list_afr) <- c("Covariates",
                     "Covariates + OA",
                     "Covariates + OA + PRS")
auc_text_afr <- lapply(roc_list_afr, function(x){
  aucix <- pROC::ci(x) %>% round(4)
  paste0("AUC: ", 
         aucix[2] , "(",
         aucix[1], "~", 
         aucix[3], ")")}) %>% unlist
auc_text_afr_df <- data.frame(x = rep(0.7, 3),
                          y = c(0.3, 0.25, 0.2),
                          auc_text_afr)
roc_df_afr <- lapply(names(roc_list_afr), function(x){
  data.frame(Model = x,
             rev_specificity = 1-roc_list_afr[[x]]$specificities,
             sensitivity = roc_list_afr[[x]]$sensitivities)
}) %>% Reduce("rbind", .)

plt_D <- ggplot(roc_df_afr) +
  geom_line(aes(x = rev_specificity, y = sensitivity, 
                color = Model, linetype = Model), 
            size = 1.3) +
  annotate("text", 
           x = auc_text_afr_df[1, 1] , 
           y = auc_text_afr_df[1, 2], 
           label = auc_text_afr_df[1, 3],
           colour="#E41A1C") + 
  annotate("text", 
           x = auc_text_afr_df[2, 1] , 
           y = auc_text_afr_df[2, 2], 
           label = auc_text_afr_df[2, 3],
           colour="#377EB8") + 
  annotate("text", 
           x = auc_text_afr_df[3, 1] , 
           y = auc_text_afr_df[3, 2], 
           label = auc_text_afr_df[3, 3],
           colour="#4DAF4A") + 
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A")) + 
  xlab("1-Specificity") +
  ylab("Sensitivity") +
  theme_bw() + 
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        # legend.position = "top",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"))

library(patchwork)
tiff("pred_plt.tiff",
     height = 8.5, width = 11, units = "in", compression = "lzw", res = 600)
# (((plt_A | plt_B) + 
#   plot_layout(widths = c(5.5, 6))) /
#   ((plt_C | plt_D) + 
#      plot_layout(widths = c(5.5, 6)))) + 
#   plot_annotation(tag_levels = "A") &
#   theme(plot.tag = element_text(size = 15, face = "bold", colour = "black"))

((plt_A | plt_C) /
    (plt_B | plt_D)) + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 15, face = "bold", colour = "black"))

dev.off()
