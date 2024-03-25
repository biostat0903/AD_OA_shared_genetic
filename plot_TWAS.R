rm(list=ls())
gc()
# load packages
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(corrplot)

# Parameters
PROJ_PATH <- "/public/home/biostat03/project/oaProject/"
setwd(PROJ_PATH)
DATA_INPUT <- paste0(PROJ_PATH, "03_result/04_twas/")
WORK_PATH <- paste0(DATA_INPUT, "merge/")
SIGNIF_PATH <- paste0(DATA_INPUT, "signif/")
INTERSEC_PATH <- paste0(DATA_INPUT, "intersection/")

# set target 23 tissues
use_tissue <- c("Adipose_Subcutaneous", "Adipose_Visceral", "Blood", "Brain_Substantia_nigra",
                "Brain_Spinal_cord_cervical", "Brain_Putamen_basal_ganglia",
                "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Hypothalamus", "Brain_Hippocampus", 
                "Brain_Frontal_Cortex_BA9", "Brain_Cortex", "Cerebellum", "Brain_Cerebellar_Hemisphere",
                "Brain_Caudate_basal_ganglia", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Amygdala",
                "Muscle", "Colon_Sigmoid",
                "Colon_Transverse", "Liver", "Nerve_Tibial", "Pituitary",
                "Small_Intestine_Terminal_Ileum")
use_file_OA <- paste0("OA_", use_tissue, ".txt")
use_file_AD <- paste0("AD_", use_tissue, ".txt")

# extract common significant TWAS genes
TWAS_sig_OA <- lapply(use_file_OA, function(filex){
  fread2(paste0(SIGNIF_PATH, "/", filex))
}) %>% Reduce("rbind", .)

TWAS_sig_AD <- lapply(use_file_AD, function(filex){
  fread2(paste0(SIGNIF_PATH, "/", filex))
}) %>% Reduce("rbind", .)

common_gene <- intersect(unique(TWAS_sig_OA$genesymbol), 
                         unique(TWAS_sig_AD$genesymbol))
common_gene <- common_gene[!is.na(common_gene)]

# extract TWAS result on common genes
TWAS_common_OA <- lapply(use_file_OA, function(filex){
  dfx <- fread2(paste0(WORK_PATH, "/", filex))
  dfx <- subset(dfx, dfx$genesymbol %in% common_gene, 
                select = c("genesymbol", "LOC", "PANEL",
                           "TWAS.Z", "TWAS.P"))
}) %>% Reduce("rbind", .)

TWAS_common_AD <- lapply(use_file_AD, function(filex){
  dfx <- fread2(paste0(WORK_PATH, "/", filex))
  dfx <- subset(dfx, dfx$genesymbol %in% common_gene, 
                select = c("genesymbol", "LOC", "PANEL",
                           "TWAS.Z", "TWAS.P"))
}) %>% Reduce("rbind", .)

# add sig index
TWAS_common_OA$Tissue <- gsub("GTExv8.EUR.", "", TWAS_common_OA$PANEL)
TWAS_common_AD$Tissue <- gsub("GTExv8.EUR.", "", TWAS_common_AD$PANEL)
TWAS_sig_OA_index <- paste0(TWAS_sig_OA$genesymbol, TWAS_sig_OA$PANEL)
TWAS_sig_AD_index <- paste0(TWAS_sig_AD$genesymbol, TWAS_sig_AD$PANEL)
TWAS_common_OA_index <- paste0(TWAS_common_OA$genesymbol, TWAS_common_OA$PANEL)
TWAS_common_AD_index <- paste0(TWAS_common_AD$genesymbol, TWAS_common_AD$PANEL)
TWAS_common_OA$Sig <- ifelse(TWAS_common_OA_index %in% TWAS_sig_OA_index,
                             "*", "")
TWAS_common_AD$Sig <- ifelse(TWAS_common_AD_index %in% TWAS_sig_AD_index,
                             "*", "")
saveRDS(TWAS_common_OA, file = "04_twas/TWAS_common_OA.rds")
saveRDS(TWAS_common_AD, file = "04_twas/TWAS_common_AD.rds")

# plot heatmap
palette_1 <- COL2("RdYlBu", 200)[190:11]
# TWAS_common_OA_lim <- min(abs(range(TWAS_common_OA$TWAS.Z)))
# TWAS_common_AD_lim <- min(abs(range(TWAS_common_AD$TWAS.Z)))
# TWAS_common_OA$TWAS.Z[TWAS_common_OA$TWAS.Z > TWAS_common_OA_lim] <- TWAS_common_OA_lim
# TWAS_common_OA$TWAS.Z[TWAS_common_OA$TWAS.Z < -TWAS_common_OA_lim] <- -TWAS_common_OA_lim
# TWAS_common_AD$TWAS.Z[TWAS_common_AD$TWAS.Z > TWAS_common_AD_lim] <- TWAS_common_AD_lim
# TWAS_common_AD$TWAS.Z[TWAS_common_AD$TWAS.Z < TWAS_common_AD_lim] <- -TWAS_common_AD_lim

TWAS_common_plt_OA <- ggplot(TWAS_common_OA, aes(x = genesymbol, y = Tissue)) + 
  geom_raster(aes(fill = TWAS.Z)) + 
  geom_text(aes(label = Sig), color = "black", size = 4) +
  scale_fill_gradientn(colours = palette_1) + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.background = element_blank(),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

TWAS_common_plt_AD <- ggplot(TWAS_common_AD, aes(x = genesymbol, y = Tissue)) + 
  geom_raster(aes(fill = TWAS.Z)) + 
  geom_text(aes(label = Sig), color = "black", size = 4) +
  scale_fill_gradientn(colours = palette_1) + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.background = element_blank(),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

# output
library(patchwork)
tiff("05_plots/TWAS_common_plt.tiff",
     height = 7, width = 15, units = "in",
     compression = "lzw", res = 300)
(TWAS_common_plt_OA | TWAS_common_plt_AD)+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 15, face = "bold"))
dev.off()
