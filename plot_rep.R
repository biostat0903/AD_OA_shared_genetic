library(ggplot2)
library(stringr)
setwd("/public/home/biostat03/project/oaProject/03_result")

eth <- "EUR"
load(paste0("nerv_result_CV", eth, ".rda"))
load(paste0("nerv_result_main", eth, ".rda"))

nerv_result_sens$main <- nerv_result
cv_sv_or <- lapply(nerv_result_sens, function(xx){
  sv_orx <- xx$sv_res[xx$sv_res$Var == "oa", 1]
  gsub("\\(", "_", sv_orx) %>%
    gsub("\\~", "_", .) %>%
    gsub("\\)", "", .) %>%
    str_split(., "_") %>% 
    unlist() %>% as.numeric()
}) %>% Reduce("rbind", .) %>% as.data.frame()

colnames(cv_sv_or) <- c("OR", "LowCI", "HighCI")
rownames(cv_sv_or) <- NULL
cv_sv_or$ID <- c(paste0("CV", 1:10), "Main") %>%
  factor(., 
         levels = c("Main", paste0("CV", 10:1)))
cv_sv_or$Group <- c(rep("CV", 10), "Main")

ggplot(cv_sv_or, aes(x = ID,y = OR, color = Group)) +
  geom_hline(yintercept = 1, color = "grey", lty = 5, lwd = 1)+
  geom_pointrange(aes(ymin = LowCI, ymax = HighCI),
                  stat='identity',position = position_dodge(0.4),
                  size = 1)+
  xlab("") + ylab("OR(95%CI)") +
  scale_color_manual(values = c("blue", "red"))+
  coord_flip()+
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 15, face = "bold.italic"),
        legend.position =  "none")

####
eth <- "EUR"
load(paste0("nerv_result_CV", eth, ".rda"))
load(paste0("nerv_result_main", eth, ".rda"))

nerv_result_sens$main <- nerv_result
cv_sv_or <- lapply(nerv_result_sens, function(xx){
  sv_orx <- xx$sv_res[xx$sv_res$Var == "oa", 1]
  gsub("\\(", "_", sv_orx) %>%
    gsub("\\~", "_", .) %>%
    gsub("\\)", "", .) %>%
    str_split(., "_") %>% 
    unlist() %>% as.numeric()
}) %>% Reduce("rbind", .) %>% as.data.frame()

colnames(cv_sv_or) <- c("OR", "LowCI", "HighCI")
rownames(cv_sv_or) <- NULL
cv_sv_or$ID <- c(paste0("CV", 1:10), "Main") %>%
  factor(., 
         levels = c("Main", paste0("CV", 10:1)))
cv_sv_or$Group <- c(rep("CV", 10), "Main")

ggplot(cv_sv_or, aes(x = ID,y = OR, color = Group)) +
  geom_hline(yintercept = 1, color = "grey", lty = 5, lwd = 1)+
  geom_pointrange(aes(ymin = LowCI, ymax = HighCI),
                  stat='identity',position = position_dodge(0.4),
                  size = 1)+
  xlab("") + ylab("OR(95%CI)") +
  scale_color_manual(values = c("blue", "red"))+
  coord_flip()+
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 15, face = "bold.italic"),
        legend.position =  "none")


###
cv_sv_or <- lapply(nerv_result_sens, function(xx){
  sv_orx <- xx$sv_res[xx$sv_res$Var == "oa", 1]
  gsub("\\(", "_", sv_orx) %>%
    gsub("\\~", "_", .) %>%
    gsub("\\)", "", .) %>%
    str_split(., "_") %>% 
    unlist() %>% as.numeric()
}) %>% Reduce("rbind", .) %>% as.data.frame()

colnames(cv_sv_or) <- c("OR", "LowCI", "HighCI")
rownames(cv_sv_or) <- NULL
cv_sv_or$ID <- c(paste0("Replication", 1:10), "Main") %>%
  factor(., 
         levels = c("Main", paste0("Replication", 10:1)))
cv_sv_or$Group <- c(rep("Replication", 10), "Main") %>%
  factor(., 
         levels = c("Replication", "Main"))

cv_sv_plt <- ggplot(cv_sv_or, aes(x = ID,y = OR, color = Group)) +
  geom_hline(yintercept = 1, color = "grey", lty = 5, lwd = 1)+
  geom_pointrange(aes(ymin = LowCI, ymax = HighCI),
                  stat='identity',position = position_dodge(0.4),
                  size = 1)+
  xlab("") + ylab("OR(95%CI)") +
  scale_color_manual(values = c("blue", "red"))+
  coord_flip()+
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 15, face = "bold.italic"),
        legend.position =  "none")

####
eth <- "EUR"
load(paste0("nerv_result_CV", eth, ".rda"))
load(paste0("nerv_result_main", eth, ".rda"))

nerv_result_sens$main <- nerv_result
cv_mv_or <- lapply(nerv_result_sens, function(xx){
  mv_orx <- xx$multi_res[xx$multi_res$Var == "oa", 1]
  gsub("\\(", "_", mv_orx) %>%
    gsub("\\~", "_", .) %>%
    gsub("\\)", "", .) %>%
    str_split(., "_") %>% 
    unlist() %>% as.numeric()
}) %>% Reduce("rbind", .) %>% as.data.frame()

colnames(cv_mv_or) <- c("OR", "LowCI", "HighCI")
rownames(cv_mv_or) <- NULL
cv_mv_or$ID <- c(paste0("Replication", 1:10), "Main") %>%
  factor(., 
         levels = c("Main", paste0("Replication", 10:1)))
cv_mv_or$Group <- c(rep("Replication", 10), "Main") %>%
  factor(., 
         levels = c("Replication", "Main"))

cv_mv_plt <- ggplot(cv_mv_or, aes(x = ID,y = OR, color = Group)) +
  geom_hline(yintercept = 1, color = "grey", lty = 5, lwd = 1)+
  geom_pointrange(aes(ymin = LowCI, ymax = HighCI),
                  stat='identity',position = position_dodge(0.4),
                  size = 1)+
  xlab("") + ylab("OR(95%CI)") +
  scale_color_manual(values = c("blue", "red"))+
  coord_flip()+
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 15, face = "bold.italic"),
        legend.position =  "none")

library(patchwork)

tiff("Rep_plot.tiff",
     height = 7, width = 9, units = "in", res = 600, compression = "lzw")
(cv_sv_plt | cv_mv_plt) + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 15, face = "bold"))
dev.off()

