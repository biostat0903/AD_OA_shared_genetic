rm(list=ls())
gc()
# load packages
library(coloc)
library(dplyr)
library(data.table)

# Parameters
PROJ_PATH <- "/public/home/biostat03/project/oaProject/"
EQTL_DATA <- "/public/home/Datasets/dbGaP/GTEx/V7_summary/GTEx_Analysis_v8_eQTL_sig/"
RAW_SUMM <- paste0(PROJ_PATH, "02_data/summary_stat/cleaned/")
SUMM_DATA <- paste0(PROJ_PATH, "02_data/Colocalization/COLOC_v8/")
EQTL_PATH <- paste0(SUMM_DATA, "GTEx_v8_cleaned/")
WORK_PATH <- paste0(PROJ_PATH, "03_result/05_colocalization/coloc_v8/")
LD_REF <- "/public/home/Datasets/1000GP/EUR/all/merge"
PLINK <- "/public/home/biostat03/biosoft/plink1.9"
KGP_PATH <- "/public/home/biostat03/gwas_ref_data/eur_1kgp_all/"

# Munge 1KGP function
munge.kgp <- function(data_ldsc, snplist, type){
  
  data_ldsc_ref <- left_join(snplist, data_ldsc, by = c("rs_id"))
  ## delete allele
  na_row1 <- which(data_ldsc_ref$A1.x != data_ldsc_ref$A2.y & data_ldsc_ref$A1.x != data_ldsc_ref$A1.y)
  na_row2 <- which(data_ldsc_ref$A2.x != data_ldsc_ref$A2.y & data_ldsc_ref$A2.x != data_ldsc_ref$A1.y)
  na_row <- c(na_row1, na_row2) %>% unique()
  if(type == "cc"){
    
    if(length(na_row) != 0){
      
      data_ldsc_ref$beta[na_row] <- NA
    }
    ## flip allele
    fl_row <- which(data_ldsc_ref$A1.x == data_ldsc_ref$A2.y & data_ldsc_ref$A2.x == data_ldsc_ref$A1.y)
    if(length(fl_row) != 0){
      
      data_ldsc_ref$beta[fl_row] <- -1 * data_ldsc_ref$beta[fl_row]
    }
    data_ldsc_ref_f <- data.frame(rs_id = data_ldsc_ref$rs_id,
                                  beta = data_ldsc_ref$beta,
                                  varbeta = data_ldsc_ref$varbeta,
                                  pvalues = data_ldsc_ref$pval,
                                  A1 = data_ldsc_ref$A1.x,
                                  A2 = data_ldsc_ref$A2.x)
  }
  if(type == "quant"){
    
    if(length(na_row) != 0){
      
      data_ldsc_ref$slope[na_row] <- NA
    }
    ## flip allele
    fl_row <- which(data_ldsc_ref$A1.x == data_ldsc_ref$A2.y & data_ldsc_ref$A2.x == data_ldsc_ref$A1.y)
    if(length(fl_row) != 0){
      
      data_ldsc_ref$slope[fl_row] <- -1 * data_ldsc_ref$slope[fl_row]
    }
    data_ldsc_ref_f <- data.frame(rs_id = data_ldsc_ref$rs_id,
                                  MAF = data_ldsc_ref$maf,
                                  ma_count = data_ldsc_ref$ma_count,
                                  pval_nominal = data_ldsc_ref$pval_nominal,
                                  beta = data_ldsc_ref$slope,
                                  varbeta = (data_ldsc_ref$slope_se)^2,
                                  A1 = data_ldsc_ref$A1.x,
                                  A2 = data_ldsc_ref$A2.x)
  }
  
  return(data_ldsc_ref_f)
}

## TWAS significant Genes
tissue <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Whole_Blood", "Brain_Substantia_nigra",
            "Brain_Spinal_cord_cervical_c-1", "Brain_Putamen_basal_ganglia",
            "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Hypothalamus", "Brain_Hippocampus", 
            "Brain_Frontal_Cortex_BA9", "Brain_Cortex", "Brain_Cerebellum", "Brain_Cerebellar_Hemisphere",
            "Brain_Caudate_basal_ganglia", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Amygdala",
            "Muscle_Skeletal", "Colon_Sigmoid",
            "Colon_Transverse", "Liver", "Nerve_Tibial", "Pituitary", "Small_Intestine_Terminal_Ileum")
tissue_n <- c(479, 393, 558, 100,
              115, 153,
              181, 156, 150,
              157, 183, 188, 157,
              172, 135, 119,
              588, 266,
              294, 178, 438, 219, 141)
load(paste0(PROJ_PATH, "04_reference/Gencode_hg19.RData"))
gene_pos <- gencode[,c(8,10,3,4)]
gene_pos <- gene_pos[!duplicated(gene_pos$genesymbol),]
gencode <- gencode[,8:9]
colnames(gencode) <- c("gene_symbol", "gene_id")
kgp_snplist <- fread(paste0(KGP_PATH, "kgp.frq"))
colnames(kgp_snplist)[2] <- "rs_id"

full_result <- data.frame()
# for(dis in c("OA", "AD")){
  dis = "AD"
  # for(i in tissue){
    gene_info <- read.table("/public/home/biostat03/project/oaProject/01_code/coloc_gene.txt", 
                            sep = "\t")
    gene_info
    snp_info <- matrix(NA, nrow(gene_info), 2)
    for (ii in 1: nrow(gene_info)){
      
      i = gene_info[ii, 2]
      target_gene = gene_info[ii, 1]
      eqtl <- fread(paste0(EQTL_PATH, i, "_cleaned.txt"))
      
      gwas <- fread(paste0(SUMM_DATA, dis, ".txt"))
      gwas <- gwas[!is.na(gwas$gwas_chr),]
      N <- ifelse(dis == "OA", 77052,
                  ifelse(dis == "AD", 21982, NA))
      
      # coloc_result <- data.frame()
      # Gene <- unique(eqtl$gene_symbol)
      # for(target_gene in Gene){
      
      win_size <- 1000000
      eqtl_mat <- eqtl[which(eqtl$gene_symbol == target_gene),]
      gwas_mat <- gwas[which(gwas$gwas_chr == gene_pos$CHR[gene_pos$genesymbol == target_gene] &
                               gwas$gwas_bp > (gene_pos$start[gene_pos$genesymbol == target_gene] - win_size) &
                               gwas$gwas_bp < (gene_pos$start[gene_pos$genesymbol == target_gene] + win_size)),]
      
      snp_intersec <- intersect(eqtl_mat$rs_id, gwas_mat$rs_id) %>% intersect(kgp_snplist$rs_id)
      
      # if(length(snp_intersec) > 10){
        
        intersec_snplist <- kgp_snplist %>% filter(rs_id %in% snp_intersec)
        
        gwas_mat <- gwas_mat[!duplicated(gwas_mat$rs_id),]
        gwas_mat <- munge.kgp(gwas_mat, intersec_snplist, type = "cc")
        
        eqtl_mat <- eqtl_mat[!duplicated(eqtl_mat$rs_id),]
        eqtl_mat <- munge.kgp(eqtl_mat, intersec_snplist, type = "quant")
        
        gwas_na <- gwas_mat$rs_id[is.na(gwas_mat$beta)]
        eqtl_na <- eqtl_mat$rs_id[is.na(eqtl_mat$beta)]
        na_snp <- union(gwas_na, eqtl_na)
        
        if(length(na_snp) > 0){
          
          gwas_mat <- gwas_mat[-which(gwas_mat$rs_id %in% na_snp),]
          eqtl_mat <- eqtl_mat[-which(eqtl_mat$rs_id %in% na_snp),]
        }
        
        dataset1 = list(snp = gwas_mat$rs_id,
                        beta = gwas_mat$beta,
                        varbeta = gwas_mat$varbeta,
                        pvalues = gwas_mat$pval, 
                        type = "cc", 
                        N = N)
        dataset2 = list(snp = eqtl_mat$rs_id,
                        pvalues = eqtl_mat$pval_nominal, 
                        varbeta = eqtl_mat$varbeta,
                        beta = eqtl_mat$beta,
                        type = "quant",
                        N = tissue_n[which(i == tissue)],
                        MAF = eqtl_mat$MAF)
        
        result <- coloc.abf(dataset1, dataset2,
                            p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)
        
        lead_index <- which.max(result[[2]]$SNP.PP.H4)
        snp_info[ii, 1] <- as.vector(result[[2]][lead_index, 1])
        snp_info[ii, 2] <- as.vector(result[[2]][lead_index, 11])
    }
    
   write.csv(snp_info, file = "/public/home/biostat03/project/oaProject/01_code/coloc_snp_AD.csv", 
   )
    
    
    # }
#     eqtl <- fread(paste0(EQTL_PATH, i, "_cleaned.txt"))
#     
#     gwas <- fread(paste0(SUMM_DATA, dis, ".txt"))
#     gwas <- gwas[!is.na(gwas$gwas_chr),]
#     N <- ifelse(dis == "OA", 77052,
#                 ifelse(dis == "AD", 21982, NA))
#     
#     # coloc_result <- data.frame()
#     Gene <- unique(eqtl$gene_symbol)
#     # for(target_gene in Gene){
#       
#       win_size <- 1000000
#       eqtl_mat <- eqtl[which(eqtl$gene_symbol == target_gene),]
#       gwas_mat <- gwas[which(gwas$gwas_chr == gene_pos$CHR[gene_pos$genesymbol == target_gene] &
#                                gwas$gwas_bp > (gene_pos$start[gene_pos$genesymbol == target_gene] - win_size) &
#                                gwas$gwas_bp < (gene_pos$start[gene_pos$genesymbol == target_gene] + win_size)),]
#       
#       snp_intersec <- intersect(eqtl_mat$rs_id, gwas_mat$rs_id) %>% intersect(kgp_snplist$rs_id)
#       
#       if(length(snp_intersec) > 10){
#         
#         intersec_snplist <- kgp_snplist %>% filter(rs_id %in% snp_intersec)
#         
#         gwas_mat <- gwas_mat[!duplicated(gwas_mat$rs_id),]
#         gwas_mat <- munge.kgp(gwas_mat, intersec_snplist, type = "cc")
#         
#         eqtl_mat <- eqtl_mat[!duplicated(eqtl_mat$rs_id),]
#         eqtl_mat <- munge.kgp(eqtl_mat, intersec_snplist, type = "quant")
#         
#         gwas_na <- gwas_mat$rs_id[is.na(gwas_mat$beta)]
#         eqtl_na <- eqtl_mat$rs_id[is.na(eqtl_mat$beta)]
#         na_snp <- union(gwas_na, eqtl_na)
#         
#         if(length(na_snp) > 0){
#           
#           gwas_mat <- gwas_mat[-which(gwas_mat$rs_id %in% na_snp),]
#           eqtl_mat <- eqtl_mat[-which(eqtl_mat$rs_id %in% na_snp),]
#         }
#         
#         dataset1 = list(snp = gwas_mat$rs_id,
#                         beta = gwas_mat$beta,
#                         varbeta = gwas_mat$varbeta,
#                         pvalues = gwas_mat$pval, 
#                         type = "cc", 
#                         N = N)
#         dataset2 = list(snp = eqtl_mat$rs_id,
#                         pvalues = eqtl_mat$pval_nominal, 
#                         varbeta = eqtl_mat$varbeta,
#                         beta = eqtl_mat$beta,
#                         type = "quant",
#                         N = tissue_n[which(i == tissue)],
#                         MAF = eqtl_mat$MAF)
#         
#         result <- coloc.abf(dataset1, dataset2,
#                             p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)
#         result <- result$summary %>%
#           data.frame %>% t
#         result <- cbind(result, disease = dis,
#                         tissue = i,
#                         Gene = target_gene)
#         # coloc_result <- rbind(coloc_result, result)
#       }
#     }
#     write.table(coloc_result, row.names = F, quote = F, sep = "\t",
#                 file = paste0(WORK_PATH, "TWAS/",  dis, "_", i,  "_coloc_result.txt"))
#     full_result <- rbind(full_result, coloc_result)
#     print(paste0(dis, "_", i))
#   }
# }
write.table(full_result, row.names = F, quote = F, sep = "\t",
            file = paste0(WORK_PATH, "TWAS/", "full_result.txt"))

# Table S7
result <- fread(paste0(WORK_PATH, "TWAS/", "full_result.txt"))
result <- result[result$PP.H4.abf > 0.7,]
result_oa <- result[result$disease == "OA",]
result_ad <- result[result$disease == "AD",]
result <- inner_join(result_oa, result_ad, by = "Gene", relationship = "many-to-many")

# Fig2 D E
# Blood and Adipose
tissue <- "Adipose_Subcutaneous"
tissue <- "Adipose_Visceral_Omentum"
eqtl <- fread(paste0(EQTL_PATH, tissue, "_cleaned.txt"))
eqtl <- eqtl[which(eqtl$gene_symbol == "PSMC3"),]
eqtl <- eqtl[order(eqtl$BP),]
write.table(eqtl, row.names = F, quote = F, sep = "\t",
            file = paste0("/public/home/biostat03/project/oaProject/05_plots", "eqtl_", tissue, ".txt"))
