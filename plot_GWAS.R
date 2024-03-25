library(bigreadr)
library(dplyr)
library(stringr)
library(ggplot2)

setwd("/public/home/biostat03/project/oaProject/")
ref_EUR <- fread2("04_reference/SNP_BP_GRCh37.txt")
# calculate the length of each chromosome
chr_len <- ref_EUR %>% 
  group_by(CHR) %>% 
  summarise(chr_len = max(BP))

# set the start x of each chr
chr_pos <- chr_len  %>% 
  mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) 

# calculate the cumulative position of each SNP
SNP_info <- chr_pos %>%
  left_join(ref_EUR, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + total)

# set x axis
X_axis <-  SNP_info %>% 
  group_by(CHR) %>% 
  summarize(center=(max(BPcum) + min(BPcum))/2)

summ_AD <- fread2("02_data/summary_stat/raw/IGAP_stage_1.txt")
summ_OA <- fread2("02_data/summary_stat/cleaned/OA_frontier.txt")

SNP_info$P_AD <- summ_AD$Pvalue[match(SNP_info$SNP, summ_AD$MarkerName)] %>% as.numeric()
SNP_info$P_OA <- summ_OA$P[match(SNP_info$SNP, summ_OA$SNP)] %>% as.numeric()
SNP_info$P_AD[SNP_info$P_AD == 0] <- min(SNP_info$P_AD[SNP_info$P_AD != 0])
SNP_info$P_OA[SNP_info$P_OA == 0] <- min(SNP_info$P_OA[SNP_info$P_OA != 0])
SNP_info$CHR <- as.factor(SNP_info$CHR)

## Manhattan plot
# Manhattan plot for OA
Manh_plt_OA <- ggplot(SNP_info) +
  geom_point(aes(x = BPcum, y = -log10(P_OA), color = CHR),
             alpha = 0.8, size = 0.8)+
  scale_x_continuous(label = X_axis$CHR, 
                     breaks= X_axis$center, 
                     limits = range(SNP_info$BPcum)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = rep(c("#DE423D", "#3D5587"), 11)) +
  xlab("Chromosome") + ylab(bquote("-"~Log[10]~"(P value)"))+
  geom_hline(yintercept = -log10(5E-8), 
             color = 'red', linewidth = 0.8, linetype = 2) + 
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold", color = "black"),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# Manhattan plot for AD
Manh_plt_AD <- ggplot(SNP_info) +
  geom_point(aes(x = BPcum, y = -log10(P_AD), color = CHR),
             alpha = 0.8, size = 0.8)+
  scale_x_continuous(label = X_axis$CHR, 
                     breaks= X_axis$center, 
                     limits = range(SNP_info$BPcum)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = rep(c("#DE423D", "#3D5587"), 11)) +
  xlab("Chromosome") + ylab(bquote("-"~Log[10]~"(P value)"))+
  geom_hline(yintercept = -log10(5E-8), 
             color = 'red', linewidth = 0.8, linetype = 2) + 
  theme_bw() +
  theme(legend.position="none",
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold", color = "black"),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


## QQ plot
# QQ plot for OA
p_oa <- SNP_info$P_OA[!is.na(SNP_info$P_OA)]
qq_oa <- data.frame(obs = -log10(sort(p_oa, decreasing = FALSE)),
                    exp = -log10(ppoints(length(p_oa))))
qq_plt_oa <- ggplot(data = qq_oa, aes(exp, obs))+
  geom_point(alpha = 0.8, color = "grey60")+
  geom_abline(color = "red", linewidth = 0.8, linetype = 2)+
  xlab(bquote("Expected -"~Log[10]~"(P value)"))+
  ylab(bquote("Observed -"~Log[10]~"(P value)"))+
  theme_bw()+
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.background = element_blank())

# QQ plot for AD
p_ad <- SNP_info$P_AD[!is.na(SNP_info$P_AD)]
qq_ad <- data.frame(obs = -log10(sort(p_ad, decreasing = FALSE)),
                    exp = -log10(ppoints(length(p_ad))))
qq_plt_ad <- ggplot(data = qq_ad, aes(exp, obs))+
  geom_point(alpha = 0.8, color = "grey60")+
  geom_abline(color = "red", linewidth = 0.8, linetype = 2)+
  xlab(bquote("Expected -"~Log[10]~"(P value)"))+
  ylab(bquote("Observed -"~Log[10]~"(P value)"))+
  theme_bw()+
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.background = element_blank())

## output
library(patchwork)
tiff("05_plots/GWAS_plt.tiff", 
     height = 6, width = 13, units = "in", 
     res = 300, compression = "lzw")
((Manh_plt_OA / Manh_plt_AD) |
  (qq_plt_oa / qq_plt_ad)) + 
  plot_layout(widths = c(10, 3))+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 15, face = "bold"))
dev.off()

