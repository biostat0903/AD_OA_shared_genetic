rm(list=ls())
gc()
# load packages
library(data.table)
library(dplyr)
library(purrr)

# Parameters
PROJ_PATH <- "/public/home/biostat03/project/oaProject/"
DATA_INPUT <- paste0(PROJ_PATH, "03_result/04_twas/")
WORK_PATH <- paste0(DATA_INPUT, "merge/")
SIGNIF_PATH <- paste0(DATA_INPUT, "signif/")
INTERSEC_PATH <- paste0(DATA_INPUT, "intersection/")

load(paste0(PROJ_PATH, "04_reference/Gencode_hg19.RData"))
region <- fread(paste0(PROJ_PATH, "04_reference/LAVA_s2500_m25_f1_w200.blocks.txt"))

gencode$LOC <- NA
for(i in 1:nrow(gencode)){
  
  LOC_in <- region$LOC[which(region$CHR == gencode$CHR[i] &
                               region$START < gencode$start[i] &
                               region$STOP > gencode$end[i])]
  LOC_between <- paste0(max(region$LOC[which(region$CHR == gencode$CHR[i] &
                                               region$START < gencode$start[i])]), "_",
                        min(region$LOC[which(region$CHR == gencode$CHR[i] &
                                               region$STOP > gencode$end[i])]))
  LOC <- ifelse(length(LOC_in) == 0,
                LOC_between,
                LOC_in)
  gencode$LOC[i] <- LOC
  
  if(i %% 1000 == 0){
    print(i)
  }
}
write.table(gencode, row.names = F, quote = F, sep = "\t",
            file = paste0(INTERSEC_PATH, "gencode.txt"))

gencode <- fread(paste0(INTERSEC_PATH, "gencode.txt"))
gencode <- gencode[, c(8:9, 11)]
names(gencode)[2] <- "ID"

for(dis in c("OA", "AD")){
  
  for(tissue in c("Adipose_Subcutaneous", "Adipose_Visceral")){
    
    data_merge <- data.frame()
    for(chr in 1:22){
      
      fusion_output <- fread(paste0(DATA_INPUT, dis, "_", tissue, "_", chr, ".txt"))
      if(chr == 6){
        
        fusion_mhc <- fread(paste0(DATA_INPUT, dis, "_", tissue, "_", chr, ".txt.MHC"))
        fusion_output <- rbind(fusion_output, fusion_mhc)
      }
      data_merge <- rbind(data_merge, fusion_output) 
      data_merge$ID <- substr(data_merge$ID, 1, 15)
      
    }
    data_merge <- merge(gencode, data_merge, by = "ID", all.y = T)
    write.table(data_merge, row.names = F, quote = F, sep = "\t",
                file = paste0(WORK_PATH, dis, "_", tissue, ".txt"))
    
  }
}

for(dis in c("OA", "AD")){
  
  for(tissue in c("Adipose_Subcutaneous", "Adipose_Visceral", "Blood", "Brain_Substantia_nigra",
                  "Brain_Spinal_cord_cervical", "Brain_Putamen_basal_ganglia",
                  "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Hypothalamus", "Brain_Hippocampus", 
                  "Brain_Frontal_Cortex_BA9", "Brain_Cortex", "Cerebellum", "Brain_Cerebellar_Hemisphere",
                  "Brain_Caudate_basal_ganglia", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Amygdala",
                  "Muscle", "Prostate", "Breast_Mammary", "EBV_transformed_lymphocytes", "Colon_Sigmoid",
                  "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis",
                  "Kidney_Cortex", "Liver", "Lung", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary",
                  "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid")){
    
    data_merge <- data.frame()
    for(chr in 1:22){
      
      fusion_output <- fread(paste0(DATA_INPUT, dis, "_", tissue, "_", chr, ".txt"))
      if(chr == 6){
        
        fusion_mhc <- fread(paste0(DATA_INPUT, dis, "_", tissue, "_", chr, ".txt.MHC"))
        fusion_output <- rbind(fusion_output, fusion_mhc)
      }
      data_merge <- rbind(data_merge, fusion_output) 
      data_merge$ID <- substr(data_merge$ID, 1, 15)
      
    }
    data_merge <- merge(gencode, data_merge, by = "ID", all.y = T)
    write.table(data_merge, row.names = F, quote = F, sep = "\t",
                file = paste0(WORK_PATH, dis, "_", tissue, ".txt"))
    
  }
}

LAVA <- fread("/public/home/biostat03/project/oaProject/03_result/02_gene_correlation/result_rg.csv")

# TWAS Significant results
results <- dir(WORK_PATH)
for(file in results){
  
  data_signif <- fread(paste0(WORK_PATH, file))
  data_signif$FDR_P <- p.adjust(data_signif$TWAS.P, method = "BH")
  data_signif <- data_signif[which(data_signif$FDR_P < 0.05),]
  write.table(data_signif, row.names = F, quote = F, sep = "\t",
              file = paste0(SIGNIF_PATH, file))
}

# Intersections across tissues
intersec <- combn(results, 2) %>% data.frame %>% t
full_intersec <- data.frame()
for(i in 1:nrow(intersec)){
  
  data1 <- fread(paste0(SIGNIF_PATH, intersec[i,1])) %>% data.frame
  data1$genesymbol[which(is.na(data1$genesymbol))] <- data1$ID[which(is.na(data1$genesymbol))]
  data2 <- fread(paste0(SIGNIF_PATH, intersec[i,2])) %>% data.frame
  data2$genesymbol[which(is.na(data2$genesymbol))] <- data2$ID[which(is.na(data2$genesymbol))]
  gene <- intersect(data1$genesymbol, data2$genesymbol)
  if(length(gene) != 0){
    
    data1[which(is.na(data1$genesymbol)),2] <- data1[which(is.na(data1$genesymbol)),1]
    rownames(data1) <- data1[,2] %>% as.character
    colnames(data1)[9:ncol(data1)] <- paste0(colnames(data1)[9:ncol(data1)], "_ver1")
    data1 <- data1[,-4:-5]
    data1 <- data1[gene,]
    
    data2[which(is.na(data2$genesymbol)),2] <- data2[which(is.na(data2$genesymbol)),1]
    rownames(data2) <- data2[,2] %>% as.character
    colnames(data2) <- paste0(colnames(data2), "_ver2")
    data2 <- data2[,-1:-8]
    data2 <- data2[gene,]
    
    data_intersec <- cbind(data1, data2)
    data_intersec$var1 <- gsub(".txt", "", intersec[i,1])
    data_intersec$var2 <- gsub(".txt", "", intersec[i,2])
    
    write.table(data_intersec, row.names = F, quote = F, sep = "\t",
                file = paste0(INTERSEC_PATH, gsub(".txt", "", intersec[i,1]), "_", intersec[i,2]))
    full_intersec <- rbind(full_intersec, data_intersec)
  }
}
write.table(full_intersec, row.names = F, quote = F, sep = "\t",
            file = paste0(INTERSEC_PATH, "full_result.txt"))

# 选tissue
tissue_select <- c(paste0("AD_", c("Adipose_Subcutaneous", "Adipose_Visceral", "Blood", "Brain_Substantia_nigra",
                                   "Brain_Spinal_cord_cervical", "Brain_Putamen_basal_ganglia",
                                   "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Hypothalamus", "Brain_Hippocampus", 
                                   "Brain_Frontal_Cortex_BA9", "Brain_Cortex", "Cerebellum", "Brain_Cerebellar_Hemisphere",
                                   "Brain_Caudate_basal_ganglia", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Amygdala",
                                   "Muscle", "Colon_Sigmoid",
                                   "Colon_Transverse", "Liver", "Nerve_Tibial", "Pituitary",
                                   "Small_Intestine_Terminal_Ileum")),
                   paste0("OA_", c("Adipose_Subcutaneous", "Adipose_Visceral", "Blood", "Brain_Substantia_nigra",
                                   "Brain_Spinal_cord_cervical", "Brain_Putamen_basal_ganglia",
                                   "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Hypothalamus", "Brain_Hippocampus", 
                                   "Brain_Frontal_Cortex_BA9", "Brain_Cortex", "Cerebellum", "Brain_Cerebellar_Hemisphere",
                                   "Brain_Caudate_basal_ganglia", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Amygdala",
                                   "Muscle", "Colon_Sigmoid",
                                   "Colon_Transverse", "Liver", "Nerve_Tibial", "Pituitary",
                                   "Small_Intestine_Terminal_Ileum")))
full_intersec <- full_intersec[which(full_intersec$var1 %in% tissue_select &
                                       full_intersec$var2 %in% tissue_select),]

# 与LAVA交集
full_intersec <- full_intersec %>%
  mutate(
    LOC2 = as.numeric(str_split(LOC, "_", simplify = TRUE)[, 1]),
    LOC3 = as.numeric(str_split(LOC, "_", simplify = TRUE)[, 2])
  )
full_intersec <- full_intersec[which(as.numeric(full_intersec$LOC) %in% LAVA$loc |
                                       as.numeric(full_intersec$LOC2) %in% LAVA$loc |
                                       as.numeric(full_intersec$LOC3) %in% LAVA$loc ),]
write.csv(full_intersec, paste0(INTERSEC_PATH, "full_intersec.csv"), row.names = F)

# OA与AD shared
full_intersec <- fread(paste0(INTERSEC_PATH, "full_result.txt"))
full_intersec <- full_intersec[which(full_intersec$var1 %in% tissue_select &
                                       full_intersec$var2 %in% tissue_select),]
axis_intersec <- full_intersec[which(substr(full_intersec$var1,1,2) != substr(full_intersec$var2,1,2)),]
write.csv(axis_intersec, paste0(INTERSEC_PATH, "LAVA_TWAS_intersec.csv"), row.names = F)

# LAVA整理
gencode <- fread(paste0(INTERSEC_PATH, "gencode.txt"))
LAVA <- LAVA[which(LAVA$rho > 0.8),]
gencode <- gencode %>%
  mutate(
    LOC2 = as.numeric(str_split(LOC, "_", simplify = TRUE)[, 1]),
    LOC3 = as.numeric(str_split(LOC, "_", simplify = TRUE)[, 2])
  )
gencode <- gencode[which(as.numeric(gencode$LOC) %in% LAVA$loc |
                           as.numeric(gencode$LOC2) %in% LAVA$loc |
                           as.numeric(gencode$LOC3) %in% LAVA$loc ),]
gencode <- gencode[,c(-2,-5,-12:-13)]
write.csv(gencode, paste0(INTERSEC_PATH, "LAVA_Gene.csv"), row.names = F)

#####################################################################################################################
# all TWAS Significant results combination
results <- dir(WORK_PATH)
all_results <- data.frame()
for(file in results){
  
  print(file)
  data_signif <- fread(paste0(WORK_PATH, file))
  data_signif$FDR_P <- p.adjust(data_signif$TWAS.P, method = "BH")
  data_signif <- data_signif[which(data_signif$FDR_P < 0.05),]
  all_results <- rbind(all_results, data_signif)
}
full_intersec <- full_intersec[which(full_intersec$var1 %in% tissue_select &
                                       full_intersec$var2 %in% tissue_select),]
axis_intersec <- full_intersec[which(substr(full_intersec$var1,1,2) != substr(full_intersec$var2,1,2)),]


#####################################################################################################################
# Table S6
result <- fread(paste0(INTERSEC_PATH, "LAVA_TWAS_intersec.csv"))

tables6 <- data.frame()
for(gene in unique(result$genesymbol)){
  
  a <- result[which(result$genesymbol == gene),]
  b <- a[order(a$FDR_P_ver2),]
  c <- a[order(a$FDR_P_ver1),]
  gene_min <- c(Gene = gene, tissue_oa = paste0(unique(b$var2), collapse = ", "), p_oa = b$FDR_P_ver2[1],
                tissue_ad = paste0(unique(c$var1), collapse = ", "), p_ad = c$FDR_P_ver1[1])
  tables6 <- rbind(tables6, gene_min)
}
names(tables6) <- c("Gene", "tissue_oa", "p_oa", "tissue_ad", "p_ad")
write.csv(tables6, paste0(INTERSEC_PATH, "tables6.csv"), row.names = F)

#####################################################################################################################
# Table SX4 / SX5
result <- fread(paste0(INTERSEC_PATH, "full_result.txt"))
tissue_select <- c(paste0("AD_", c("Adipose_Subcutaneous", "Adipose_Visceral", "Blood", "Brain_Substantia_nigra",
                                   "Brain_Spinal_cord_cervical", "Brain_Putamen_basal_ganglia",
                                   "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Hypothalamus", "Brain_Hippocampus", 
                                   "Brain_Frontal_Cortex_BA9", "Brain_Cortex", "Cerebellum", "Brain_Cerebellar_Hemisphere",
                                   "Brain_Caudate_basal_ganglia", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Amygdala",
                                   "Muscle", "Colon_Sigmoid",
                                   "Colon_Transverse", "Liver", "Nerve_Tibial", "Pituitary",
                                   "Small_Intestine_Terminal_Ileum")),
                   paste0("OA_", c("Adipose_Subcutaneous", "Adipose_Visceral", "Blood", "Brain_Substantia_nigra",
                                   "Brain_Spinal_cord_cervical", "Brain_Putamen_basal_ganglia",
                                   "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Hypothalamus", "Brain_Hippocampus", 
                                   "Brain_Frontal_Cortex_BA9", "Brain_Cortex", "Cerebellum", "Brain_Cerebellar_Hemisphere",
                                   "Brain_Caudate_basal_ganglia", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Amygdala",
                                   "Muscle", "Colon_Sigmoid",
                                   "Colon_Transverse", "Liver", "Nerve_Tibial", "Pituitary",
                                   "Small_Intestine_Terminal_Ileum")))
result <- result[which(result$var1 %in% tissue_select &
                         result$var2 %in% tissue_select),]
result <- result[,c(2, 7:38)]
a <- result[,c(1:16, 32)]
names(a) <- gsub("_ver1", "", names(a))
names(a)[17] <- "Tissue"
b <- result[,c(1, 17:31, 33)]
names(b) <- gsub("_ver2", "", names(b))
names(b)[17] <- "Tissue"
result <- rbind(a, b)

OA_result <- result[which(substr(result$Tissue, 1, 2) == "OA"),]
OA_result <- OA_result[!duplicated(OA_result),]
write.csv(OA_result, paste0(INTERSEC_PATH, "TWAS_OA.csv"), row.names = F)

AD_result <- result[which(substr(result$Tissue, 1, 2) == "AD"),]
AD_result <- AD_result[!duplicated(AD_result),]
write.csv(AD_result, paste0(INTERSEC_PATH, "TWAS_AD.csv"), row.names = F)


