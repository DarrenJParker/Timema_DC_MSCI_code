### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
### libraries
# 
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("edgeR")

library("edgeR")
library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(cowplot)
library(stringr)
library(gtable)
library(pheatmap)
library(RColorBrewer)
require(vegan)
library(pvclust)
library(raster)
library("SuperExactTest")
library(Rtsne)
library(dplyr)
library(hash)


print (sessionInfo())
# # 
# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS 26.0.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Europe/London
# tzcode source: internal
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] hash_2.2.6.3         dplyr_1.1.4          Rtsne_0.17           SuperExactTest_1.1.0
# [5] raster_3.6-32        sp_2.2-0             pvclust_2.2-0        vegan_2.7-1         
# [9] permute_0.9-8        RColorBrewer_1.1-3   pheatmap_1.0.13      gtable_0.3.6        
# [13] stringr_1.5.1        cowplot_1.1.3        lattice_0.22-6       ggplot2_3.5.2       
# [17] gridExtra_2.3        VennDiagram_1.7.3    futile.logger_1.4.3  edgeR_4.4.2         
# [21] limma_3.62.2        
# 
# loaded via a namespace (and not attached):
#   [1] generics_0.1.4       futile.options_1.0.1 stringi_1.8.7        magrittr_2.0.3      
# [5] Matrix_1.7-1         formatR_1.14         BiocManager_1.30.27  mgcv_1.9-1          
# [9] scales_1.4.0         codetools_0.2-20     cli_3.6.5            rlang_1.1.6         
# [13] splines_4.4.2        withr_3.0.2          tools_4.4.2          parallel_4.4.2      
# [17] locfit_1.5-9.12      lambda.r_1.2.4       vctrs_0.6.5          R6_2.6.1            
# [21] lifecycle_1.0.4      MASS_7.3-61          cluster_2.1.6        pkgconfig_2.0.3     
# [25] terra_1.8-60         pillar_1.10.2        glue_1.8.0           Rcpp_1.0.14         
# [29] statmod_1.5.1        tibble_3.2.1         tidyselect_1.2.1     rstudioapi_0.17.1   
# [33] farver_2.1.2         nlme_3.1-166         compiler_4.4.2       
# 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### filter and normalisation (code)

### takes: DGE structure, cpm cutoff (min cpm), number of samples to apply (i.e. at least this many samples must have expression greater than the cpm cutoff)
filt_and_norm <- function(y,cpm_cut,cut_in_Nsams){
  
  cat("\nNumber number of genes / samples in orig data\n")
  print(dim(y)) ### number of genes / samples
  head(cpm(y)) 
  keep <- rowSums(cpm(y)>cpm_cut) >= cut_in_Nsams 
  y <- y[keep,]
  
  y$samples$lib.size <- colSums(y$counts) # if filter recalc lib size
  y <- calcNormFactors(y)
  cat("\nLib norm factors\n")
  print(y$samples)	
  cat("\nNumber number of genes / samples after filtering\n")
  print(dim(y))
  
  return(y)
  
}


## returns full table and sig DE genes as a vector
get_DE_genes <- function(fita,FDRa){
  TT1 = topTags(fita, n =3000000000)
  TT2 = TT1$table
  temp_sig <- subset(TT2, TT2$FDR <= FDRa)	
  sig_genes = temp_sig$genes
  sig_logFC = temp_sig$logFC
  
  N_sig_genes <- length(sig_genes)
  cat("Number of sig genes: ", N_sig_genes )
  
  r_list <- list("table" = TT2, "S_gene_list" = sig_genes, "S_logFC_list" = sig_logFC )
  return(r_list)
}



setwd("../data/read_counts")


#### 
dat_to_Tce_raw <- read.csv("to_Tce_H2E.counts_wGLSLLG.csv", check.names=FALSE, stringsAsFactors=FALSE)
dat_to_Tpa_raw <- read.csv("to_Tpa_H2E.counts_wGLSL.csv", check.names=FALSE, stringsAsFactors=FALSE)
dat_to_Tps_raw <- read.csv("to_Tps_H2E.counts_wGLSL.csv", check.names=FALSE, stringsAsFactors=FALSE)

### tidy up colnames (delete genome name)
colnames(dat_to_Tce_raw) <- gsub("_to_Tce_LRv5a_mtDNAv350_HTseq", "", colnames(dat_to_Tce_raw))
colnames(dat_to_Tpa_raw) <- gsub("_to_Tpa_LRv5a_mtDNAv350_HTseq", "", colnames(dat_to_Tpa_raw))
colnames(dat_to_Tps_raw) <- gsub("_to_Tps_LRv5b_mtDNAv350_HTseq", "", colnames(dat_to_Tps_raw))

### remove '-' from sample names
colnames(dat_to_Tce_raw)     <- gsub("-", "", colnames(dat_to_Tce_raw))
colnames(dat_to_Tpa_raw)     <- gsub("-", "", colnames(dat_to_Tpa_raw))
colnames(dat_to_Tps_raw)     <- gsub("-", "", colnames(dat_to_Tps_raw))



#########################################################################################################$$$$$########
### gene count by scaffold

dat_to_Tpa_raw$scaf_keep <- ifelse(dat_to_Tpa_raw$scaf == "Tpa_LRv5a_scf1", dat_to_Tpa_raw$scaf, 
                                   ifelse(dat_to_Tpa_raw$scaf == "Tpa_LRv5a_scf2", dat_to_Tpa_raw$scaf, 
                                          ifelse(dat_to_Tpa_raw$scaf == "Tpa_LRv5a_scf3", dat_to_Tpa_raw$scaf, 
                                                 ifelse(dat_to_Tpa_raw$scaf == "Tpa_LRv5a_scf4", dat_to_Tpa_raw$scaf, 
                                                        ifelse(dat_to_Tpa_raw$scaf == "Tpa_LRv5a_scf5", dat_to_Tpa_raw$scaf, 
                                                               ifelse(dat_to_Tpa_raw$scaf == "Tpa_LRv5a_scf6", dat_to_Tpa_raw$scaf, 
                                                                      ifelse(dat_to_Tpa_raw$scaf == "Tpa_LRv5a_scf7", dat_to_Tpa_raw$scaf, 
                                                                             ifelse(dat_to_Tpa_raw$scaf == "Tpa_LRv5a_scf8", dat_to_Tpa_raw$scaf, 
                                                                                    ifelse(dat_to_Tpa_raw$scaf == "Tpa_LRv5a_scf9", dat_to_Tpa_raw$scaf, 
                                                                                           ifelse(dat_to_Tpa_raw$scaf == "Tpa_LRv5a_scf10", dat_to_Tpa_raw$scaf, 
                                                                                                  ifelse(dat_to_Tpa_raw$scaf == "Tpa_LRv5a_scf11", dat_to_Tpa_raw$scaf,
                                                                                                         ifelse(dat_to_Tpa_raw$scaf == "Tpa_LRv5a_scf12", dat_to_Tpa_raw$scaf, 
                                                                                                                ifelse(dat_to_Tpa_raw$scaf == "Tpa_LRv5a_scf13", dat_to_Tpa_raw$scaf, 
                                                                                                                       ifelse(dat_to_Tpa_raw$scaf == "Tpa_LRv5a_scf14", dat_to_Tpa_raw$scaf, NA))))))))))))))

dat_to_Tps_raw$scaf_keep <- ifelse(dat_to_Tps_raw$scaf == "Tps_LRv5b_scf1", dat_to_Tps_raw$scaf, 
                                   ifelse(dat_to_Tps_raw$scaf == "Tps_LRv5b_scf2", dat_to_Tps_raw$scaf, 
                                          ifelse(dat_to_Tps_raw$scaf == "Tps_LRv5b_scf3", dat_to_Tps_raw$scaf, 
                                                 ifelse(dat_to_Tps_raw$scaf == "Tps_LRv5b_scf4", dat_to_Tps_raw$scaf, 
                                                        ifelse(dat_to_Tps_raw$scaf == "Tps_LRv5b_scf5", dat_to_Tps_raw$scaf, 
                                                               ifelse(dat_to_Tps_raw$scaf == "Tps_LRv5b_scf6", dat_to_Tps_raw$scaf, 
                                                                      ifelse(dat_to_Tps_raw$scaf == "Tps_LRv5b_scf7", dat_to_Tps_raw$scaf, 
                                                                             ifelse(dat_to_Tps_raw$scaf == "Tps_LRv5b_scf8", dat_to_Tps_raw$scaf, 
                                                                                    ifelse(dat_to_Tps_raw$scaf == "Tps_LRv5b_scf9", dat_to_Tps_raw$scaf, 
                                                                                           ifelse(dat_to_Tps_raw$scaf == "Tps_LRv5b_scf10", dat_to_Tps_raw$scaf, 
                                                                                                  ifelse(dat_to_Tps_raw$scaf == "Tps_LRv5b_scf11", dat_to_Tps_raw$scaf,
                                                                                                         ifelse(dat_to_Tps_raw$scaf == "Tps_LRv5b_scf12", dat_to_Tps_raw$scaf, NA))))))))))))


write.csv(table(dat_to_Tce_raw$LG, useNA = "always"), "Tce_N_genes_all.csv")
write.csv(table(dat_to_Tpa_raw$scaf_keep, useNA = "always"), "Tpa_N_genes_all.csv")
write.csv(table(dat_to_Tps_raw$scaf_keep, useNA = "always"), "Tps_N_genes_all.csv")



##################################################################################################################################################
#### filter rRNA and tRNA out  

dat_to_Tce_f1 <- dat_to_Tce_raw[!grepl("^.RNA_", dat_to_Tce_raw$Gene_name), ]
dat_to_Tpa_f1 <- dat_to_Tpa_raw[!grepl("^.RNA_", dat_to_Tpa_raw$Gene_name), ]
dat_to_Tps_f1 <- dat_to_Tps_raw[!grepl("^.RNA_", dat_to_Tps_raw$Gene_name), ]

write.csv(table(dat_to_Tce_f1$LG, useNA = "always"), "Tce_N_genes_norRNAtRNA.csv")
write.csv(table(dat_to_Tpa_f1$scaf_keep, useNA = "always"), "Tpa_N_genes_norRNAtRNA.csv")
write.csv(table(dat_to_Tps_f1$scaf_keep, useNA = "always"), "Tps_N_genes_norRNAtRNA.csv")


#### filter mtDNA genes out  
dat_to_Tce <- dat_to_Tce_f1[!grepl("mtDNA_v350", dat_to_Tce_f1$scaf), ]
dat_to_Tpa <- dat_to_Tpa_f1[!grepl("mtDNA_v350", dat_to_Tpa_f1$scaf), ]
dat_to_Tps <- dat_to_Tps_f1[!grepl("mtDNA_v350", dat_to_Tps_f1$scaf), ]


write.csv(table(dat_to_Tce$LG, useNA = "always"), "Tce_N_genes_norRNAtRNAmtDNA.csv")
write.csv(table(dat_to_Tpa$scaf_keep, useNA = "always"), "Tpa_N_genes_norRNAtRNAmtDNA.csv")
write.csv(table(dat_to_Tps$scaf_keep, useNA = "always"), "Tps_N_genes_norRNAtRNAmtDNA.csv")

#### filter non-main scaf genes out 

dat_to_Tce <- subset(dat_to_Tce_f1, ! is.na(dat_to_Tce_f1$LG))
dat_to_Tpa <- subset(dat_to_Tpa_f1, ! is.na(dat_to_Tpa_f1$scaf_keep))
dat_to_Tps <- subset(dat_to_Tps_f1, ! is.na(dat_to_Tps_f1$scaf_keep))

table(dat_to_Tce$LG, useNA = "always")
table(dat_to_Tpa$scaf_keep, useNA = "always")
table(dat_to_Tps$scaf_keep, useNA = "always")

############################################################################################################################################################
### Get the samples for this project

count_samps <- function(df){
  N_samp        <- as.data.frame(colnames(df)[2:length(colnames(df))])
  N_samp$sp     <- str_split_fixed(as.character(N_samp[,1]), "_", 5)[,1]
  N_samp$sex    <- str_split_fixed(as.character(N_samp[,1]), "_", 5)[,2]
  N_samp$tiss   <- str_split_fixed(as.character(N_samp[,1]), "_", 5)[,3]
  N_samp$stage  <- str_split_fixed(as.character(N_samp[,1]), "_", 5)[,4]
  N_samp$group  <- paste(N_samp$sp, N_samp$sex,   N_samp$tiss,   N_samp$stage, sep = "_" )
  out_tab <- table(N_samp$group)
  return(out_tab)
}


### Tce
dat_to_Tce_ex <- dat_to_Tce[,  unique(c(grep("_Ha_", colnames(dat_to_Tce)),  grep("_WB_", colnames(dat_to_Tce)), grep("_rep", colnames(dat_to_Tce))))] ## excluded
dat_to_Tce    <- dat_to_Tce[, -unique(c(grep("_Ha_", colnames(dat_to_Tce)),  grep("_WB_", colnames(dat_to_Tce)), grep("_rep", colnames(dat_to_Tce))))] ## kept
colnames(dat_to_Tce_ex)
colnames(dat_to_Tce)


### Tpa
dat_to_Tpa_ex <- dat_to_Tpa[,  unique(c(grep("_Ha_", colnames(dat_to_Tpa)),  grep("_WB_", colnames(dat_to_Tpa)), grep("_rep", colnames(dat_to_Tpa))))] ## excluded
dat_to_Tpa    <- dat_to_Tpa[, -unique(c(grep("_Ha_", colnames(dat_to_Tpa)),  grep("_WB_", colnames(dat_to_Tpa)), grep("_rep", colnames(dat_to_Tpa))))] ## kept
colnames(dat_to_Tpa_ex)
colnames(dat_to_Tpa)

### Tps
dat_to_Tps_ex <- dat_to_Tps[,  unique(c(grep("_Ha_", colnames(dat_to_Tps)),  grep("_WB_", colnames(dat_to_Tps)), grep("_rep", colnames(dat_to_Tps)), grep("_ind", colnames(dat_to_Tps))))] ## excluded
dat_to_Tps    <- dat_to_Tps[, -unique(c(grep("_Ha_", colnames(dat_to_Tps)),  grep("_WB_", colnames(dat_to_Tps)), grep("_rep", colnames(dat_to_Tps)), grep("_ind", colnames(dat_to_Tps))))] ## kept
colnames(dat_to_Tps_ex)
colnames(dat_to_Tps)

count_samps(dat_to_Tce) 
count_samps(dat_to_Tpa)
count_samps(dat_to_Tps)


### filter samples as discussed in EdgeR_tSNE_sample_plots.R


## Tce
dat_to_Tce <- dat_to_Tce[,!colnames(dat_to_Tce) %in% c("Tce_F_DG_Ad_181009_lib181009DG", "Tce_M_Fe_Ad_181015_lib181015Fe")]
dat_to_Tce <- dat_to_Tce[,  -unique(c(grep("Tms_F_A_Ad_15", colnames(dat_to_Tce)),   
                                      grep("Tms_M_A_Ad_15", colnames(dat_to_Tce)),   
                                      grep("Tms_F_Fe_Ad_15", colnames(dat_to_Tce)),  
                                      grep("Tms_M_Fe_Ad_15", colnames(dat_to_Tce)),
                                      grep("Tms_F_FB_Ad_15", colnames(dat_to_Tce)),
                                      grep("Tms_M_FB_Ad_15", colnames(dat_to_Tce)),
                                      grep("Tms_F_Ta_Ad_15", colnames(dat_to_Tce)),
                                      grep("Tms_M_Ta_Ad_15", colnames(dat_to_Tce)),
                                      grep("Tms_F_B_Ad_150015_libATPQ415",  colnames(dat_to_Tce)),
                                      grep("Tms_F_Go_Ad_GN118_libATPQ684",  colnames(dat_to_Tce))))]
count_samps(dat_to_Tce)


## Tpa
dat_to_Tpa <- dat_to_Tpa[,!colnames(dat_to_Tpa) %in% c("Tpa_F_B_Ad_191009_libATPQ241")]
count_samps(dat_to_Tpa)

## Tps - nothing
count_samps(dat_to_Tps)



######################################################################################################################################################################
### add XA chr info

# Tce - genome fragmented - but can use HiC linkage groups - X == LG8
# Tpa - good genome - X = scf1, all males look correct
# Tps - Notes - good genome - X = scf3, all males look correct

dat_to_Tce$XA <- ifelse(dat_to_Tce$LG   == "HiCLG3",         "X", "A") 
dat_to_Tpa$XA <- ifelse(dat_to_Tpa$scaf == "Tpa_LRv5a_scf1", "X", "A")
dat_to_Tps$XA <- ifelse(dat_to_Tps$scaf == "Tps_LRv5b_scf3", "X", "A")

head(dat_to_Tce)

#############################################################################################################################
### 

### filter for genes expressed in male AND females IN BOTH sex and asex species (to keep it comparable)

filt_and_norm_male_female_FPKM_SA <- function(y,cpm_cut,cut_in_Nsams,gene_lens,N_fem_samps_sex, N_male_samps_sex, N_fem_samps_asex, N_male_samps_asex){
  cat("\nNumber number of genes / samples in orig data\n")
  print(dim(y)) ### number of genes / samples
  #print(head(rpkm(y, gene.length=gene_lens)))
  sex_fem_lib_Ns   <- seq(1, N_fem_samps_sex)
  sex_male_lib_Ns  <- seq(N_fem_samps_sex + 1, N_fem_samps_sex + N_male_samps_sex)
  asex_fem_lib_Ns  <- seq(N_fem_samps_sex + N_male_samps_sex + 1, N_fem_samps_sex + N_male_samps_sex +N_fem_samps_asex )
  asex_male_lib_Ns <- seq(N_fem_samps_sex + N_male_samps_sex + N_fem_samps_asex + 1, N_fem_samps_sex + N_male_samps_sex +  N_fem_samps_asex + N_male_samps_asex)  
  
  
  print(sex_fem_lib_Ns)
  print(sex_male_lib_Ns)
  print(asex_fem_lib_Ns)
  print(asex_male_lib_Ns)
  
  keep <- 
    rowSums(rpkm(y[,sex_fem_lib_Ns ], gene.length=gene_lens,log=FALSE)> cpm_cut) >= cut_in_Nsams & 
    rowSums(rpkm(y[,sex_male_lib_Ns], gene.length=gene_lens,log=FALSE)> cpm_cut) >= cut_in_Nsams &
    rowSums(rpkm(y[,asex_fem_lib_Ns ], gene.length=gene_lens,log=FALSE)> cpm_cut) >= cut_in_Nsams & 
    rowSums(rpkm(y[,asex_male_lib_Ns], gene.length=gene_lens,log=FALSE)> cpm_cut) >= cut_in_Nsams 
  
  y <- y[keep,]
  
  y$samples$lib.size <- colSums(y$counts) 
  y <- calcNormFactors(y)
  cat("\nLib norm factors\n")
  print(y$samples)	
  cat("\nNumber number of genes / samples after filtering\n")
  print(dim(y))
  
  return(y)
}







#########################################################################################################################################################################################
### to_Tps 

y_TpsTdi_A_Ad_to_Tps_UF  <- DGEList(counts=dat_to_Tps[,c(grep("_F_A_Ad_",  colnames(dat_to_Tps)),  grep("_M_A_Ad_",  colnames(dat_to_Tps))) ], genes=dat_to_Tps$Gene_name)
y_TpsTdi_B_Ad_to_Tps_UF  <- DGEList(counts=dat_to_Tps[,c(grep("_F_B_Ad_",  colnames(dat_to_Tps)),  grep("_M_B_Ad_",  colnames(dat_to_Tps))) ], genes=dat_to_Tps$Gene_name)
y_TpsTdi_DG_Ad_to_Tps_UF <- DGEList(counts=dat_to_Tps[,c(grep("_F_DG_Ad_", colnames(dat_to_Tps)),  grep("_M_DG_Ad_", colnames(dat_to_Tps))) ], genes=dat_to_Tps$Gene_name)
y_TpsTdi_FB_Ad_to_Tps_UF <- DGEList(counts=dat_to_Tps[,c(grep("_F_FB_Ad_", colnames(dat_to_Tps)),  grep("_M_FB_Ad_", colnames(dat_to_Tps))) ], genes=dat_to_Tps$Gene_name)
y_TpsTdi_Fe_Ad_to_Tps_UF <- DGEList(counts=dat_to_Tps[,c(grep("_F_Fe_Ad_", colnames(dat_to_Tps)),  grep("_M_Fe_Ad_", colnames(dat_to_Tps))) ], genes=dat_to_Tps$Gene_name)
y_TpsTdi_Gu_Ad_to_Tps_UF <- DGEList(counts=dat_to_Tps[,c(grep("_F_Gu_Ad_", colnames(dat_to_Tps)),  grep("_M_Gu_Ad_", colnames(dat_to_Tps))) ], genes=dat_to_Tps$Gene_name)
y_TpsTdi_Ta_Ad_to_Tps_UF <- DGEList(counts=dat_to_Tps[,c(grep("_F_Ta_Ad_", colnames(dat_to_Tps)),  grep("_M_Ta_Ad_", colnames(dat_to_Tps))) ], genes=dat_to_Tps$Gene_name)
y_TpsTdi_GoTe_Ad_to_Tps_UF <- DGEList(counts=dat_to_Tps[,c(grep("_F_Go_Ad_", colnames(dat_to_Tps)),  grep("_M_Te_Ad_", colnames(dat_to_Tps)))], genes=dat_to_Tps$Gene_name) 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### filter low expressed gene by FPKM or TPM, then TMM normalisation 

# FPKM 

FPKM_filt_value = 0 ### don't want to filter genes for tiss spec stuff

y_TpsTdi_A_Ad_to_Tps_F_FPKM    <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_A_Ad_to_Tps_UF, FPKM_filt_value,    0, dat_to_Tps$total_exon_length, length(grep("Tps_F_A_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_A_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_A_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_A_Ad_", colnames(dat_to_Tps))))
y_TpsTdi_B_Ad_to_Tps_F_FPKM    <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_B_Ad_to_Tps_UF, FPKM_filt_value,    0, dat_to_Tps$total_exon_length, length(grep("Tps_F_B_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_B_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_B_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_B_Ad_", colnames(dat_to_Tps))))
y_TpsTdi_DG_Ad_to_Tps_F_FPKM   <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_DG_Ad_to_Tps_UF, FPKM_filt_value,   0, dat_to_Tps$total_exon_length, length(grep("Tps_F_DG_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_DG_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_DG_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_DG_Ad_", colnames(dat_to_Tps))))
y_TpsTdi_FB_Ad_to_Tps_F_FPKM   <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_FB_Ad_to_Tps_UF, FPKM_filt_value,   0, dat_to_Tps$total_exon_length, length(grep("Tps_F_FB_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_FB_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_FB_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_FB_Ad_", colnames(dat_to_Tps))))
y_TpsTdi_Fe_Ad_to_Tps_F_FPKM   <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_Fe_Ad_to_Tps_UF, FPKM_filt_value,   0, dat_to_Tps$total_exon_length, length(grep("Tps_F_Fe_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_Fe_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_Fe_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_Fe_Ad_", colnames(dat_to_Tps))))
y_TpsTdi_Gu_Ad_to_Tps_F_FPKM   <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_Gu_Ad_to_Tps_UF, FPKM_filt_value,   0, dat_to_Tps$total_exon_length, length(grep("Tps_F_Gu_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_Gu_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_Gu_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_Gu_Ad_", colnames(dat_to_Tps))))
y_TpsTdi_Ta_Ad_to_Tps_F_FPKM   <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_Ta_Ad_to_Tps_UF, FPKM_filt_value,   0, dat_to_Tps$total_exon_length, length(grep("Tps_F_Ta_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_Ta_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_Ta_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_Ta_Ad_", colnames(dat_to_Tps))))
y_TpsTdi_GoTe_Ad_to_Tps_F_FPKM <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_GoTe_Ad_to_Tps_UF, FPKM_filt_value, 0, dat_to_Tps$total_exon_length, length(grep("Tps_F_Go_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_Te_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_Go_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_Te_Ad_", colnames(dat_to_Tps)))) ### Te,Go together


#########################################################################################################################################################################################
### to_Tpa 

y_TpaTge_A_Ad_to_Tpa_UF  <- DGEList(counts=dat_to_Tpa[,c(grep("_F_A_Ad_",  colnames(dat_to_Tpa)),  grep("_M_A_Ad_",  colnames(dat_to_Tpa))) ], genes=dat_to_Tpa$Gene_name)
y_TpaTge_B_Ad_to_Tpa_UF  <- DGEList(counts=dat_to_Tpa[,c(grep("_F_B_Ad_",  colnames(dat_to_Tpa)),  grep("_M_B_Ad_",  colnames(dat_to_Tpa))) ], genes=dat_to_Tpa$Gene_name)
y_TpaTge_DG_Ad_to_Tpa_UF <- DGEList(counts=dat_to_Tpa[,c(grep("_F_DG_Ad_", colnames(dat_to_Tpa)),  grep("_M_DG_Ad_", colnames(dat_to_Tpa))) ], genes=dat_to_Tpa$Gene_name)
y_TpaTge_FB_Ad_to_Tpa_UF <- DGEList(counts=dat_to_Tpa[,c(grep("_F_FB_Ad_", colnames(dat_to_Tpa)),  grep("_M_FB_Ad_", colnames(dat_to_Tpa))) ], genes=dat_to_Tpa$Gene_name)
y_TpaTge_Fe_Ad_to_Tpa_UF <- DGEList(counts=dat_to_Tpa[,c(grep("_F_Fe_Ad_", colnames(dat_to_Tpa)),  grep("_M_Fe_Ad_", colnames(dat_to_Tpa))) ], genes=dat_to_Tpa$Gene_name)
y_TpaTge_Gu_Ad_to_Tpa_UF <- DGEList(counts=dat_to_Tpa[,c(grep("_F_Gu_Ad_", colnames(dat_to_Tpa)),  grep("_M_Gu_Ad_", colnames(dat_to_Tpa))) ], genes=dat_to_Tpa$Gene_name)
y_TpaTge_Ta_Ad_to_Tpa_UF <- DGEList(counts=dat_to_Tpa[,c(grep("_F_Ta_Ad_", colnames(dat_to_Tpa)),  grep("_M_Ta_Ad_", colnames(dat_to_Tpa))) ], genes=dat_to_Tpa$Gene_name)
y_TpaTge_GoTe_Ad_to_Tpa_UF <- DGEList(counts=dat_to_Tpa[,c(grep("_F_Go_Ad_", colnames(dat_to_Tpa)),  grep("_M_Te_Ad_", colnames(dat_to_Tpa)))], genes=dat_to_Tpa$Gene_name) ### Te Go together

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### filter low expressed gene by FPKM or TPM, then TMM normalisation 

# FPKM 

y_TpaTge_A_Ad_to_Tpa_F_FPKM  <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_A_Ad_to_Tpa_UF, FPKM_filt_value, 0, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_A_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_A_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_A_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_A_Ad_", colnames(dat_to_Tpa))))
y_TpaTge_B_Ad_to_Tpa_F_FPKM  <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_B_Ad_to_Tpa_UF, FPKM_filt_value, 0, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_B_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_B_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_B_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_B_Ad_", colnames(dat_to_Tpa))))
y_TpaTge_DG_Ad_to_Tpa_F_FPKM <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_DG_Ad_to_Tpa_UF, FPKM_filt_value, 0, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_DG_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_DG_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_DG_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_DG_Ad_", colnames(dat_to_Tpa))))
y_TpaTge_FB_Ad_to_Tpa_F_FPKM <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_FB_Ad_to_Tpa_UF, FPKM_filt_value, 0, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_FB_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_FB_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_FB_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_FB_Ad_", colnames(dat_to_Tpa))))
y_TpaTge_Fe_Ad_to_Tpa_F_FPKM <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_Fe_Ad_to_Tpa_UF, FPKM_filt_value, 0, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_Fe_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_Fe_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_Fe_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_Fe_Ad_", colnames(dat_to_Tpa))))
y_TpaTge_Gu_Ad_to_Tpa_F_FPKM <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_Gu_Ad_to_Tpa_UF, FPKM_filt_value, 0, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_Gu_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_Gu_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_Gu_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_Gu_Ad_", colnames(dat_to_Tpa))))
y_TpaTge_Ta_Ad_to_Tpa_F_FPKM <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_Ta_Ad_to_Tpa_UF, FPKM_filt_value, 0, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_Ta_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_Ta_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_Ta_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_Ta_Ad_", colnames(dat_to_Tpa))))
y_TpaTge_GoTe_Ad_to_Tpa_F_FPKM <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_GoTe_Ad_to_Tpa_UF, FPKM_filt_value, 0, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_Go_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_Te_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_Go_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_Te_Ad_", colnames(dat_to_Tpa)))) ### Te Go together


#########################################################################################################################################################################################
### to_Tce 

y_TceTms_A_Ad_to_Tce_UF  <- DGEList(counts=dat_to_Tce[,c(grep("_F_A_Ad_",  colnames(dat_to_Tce)),  grep("_M_A_Ad_",  colnames(dat_to_Tce))) ], genes=dat_to_Tce$Gene_name)
y_TceTms_B_Ad_to_Tce_UF  <- DGEList(counts=dat_to_Tce[,c(grep("_F_B_Ad_",  colnames(dat_to_Tce)),  grep("_M_B_Ad_",  colnames(dat_to_Tce))) ], genes=dat_to_Tce$Gene_name)
y_TceTms_DG_Ad_to_Tce_UF <- DGEList(counts=dat_to_Tce[,c(grep("_F_DG_Ad_", colnames(dat_to_Tce)),  grep("_M_DG_Ad_", colnames(dat_to_Tce))) ], genes=dat_to_Tce$Gene_name)
y_TceTms_FB_Ad_to_Tce_UF <- DGEList(counts=dat_to_Tce[,c(grep("_F_FB_Ad_", colnames(dat_to_Tce)),  grep("_M_FB_Ad_", colnames(dat_to_Tce))) ], genes=dat_to_Tce$Gene_name)
y_TceTms_Fe_Ad_to_Tce_UF <- DGEList(counts=dat_to_Tce[,c(grep("_F_Fe_Ad_", colnames(dat_to_Tce)),  grep("_M_Fe_Ad_", colnames(dat_to_Tce))) ], genes=dat_to_Tce$Gene_name)
y_TceTms_Gu_Ad_to_Tce_UF <- DGEList(counts=dat_to_Tce[,c(grep("_F_Gu_Ad_", colnames(dat_to_Tce)),  grep("_M_Gu_Ad_", colnames(dat_to_Tce))) ], genes=dat_to_Tce$Gene_name)
y_TceTms_Ta_Ad_to_Tce_UF <- DGEList(counts=dat_to_Tce[,c(grep("_F_Ta_Ad_", colnames(dat_to_Tce)),  grep("_M_Ta_Ad_", colnames(dat_to_Tce))) ], genes=dat_to_Tce$Gene_name)
y_TceTms_GoTe_Ad_to_Tce_UF <- DGEList(counts=dat_to_Tce[,c(grep("_F_Go_Ad_", colnames(dat_to_Tce)),  grep("_M_Te_Ad_", colnames(dat_to_Tce)))], genes=dat_to_Tce$Gene_name) ### Te, Go together

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### filter low expressed gene by FPKM or TPM, then TMM normalisation 

# FPKM 

y_TceTms_A_Ad_to_Tce_F_FPKM      <- filt_and_norm_male_female_FPKM_SA(y_TceTms_A_Ad_to_Tce_UF, FPKM_filt_value, 0, dat_to_Tce$total_exon_length, length(grep("Tce_F_A_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_A_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_A_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_A_Ad_", colnames(dat_to_Tce))))
y_TceTms_B_Ad_to_Tce_F_FPKM      <- filt_and_norm_male_female_FPKM_SA(y_TceTms_B_Ad_to_Tce_UF, FPKM_filt_value, 0, dat_to_Tce$total_exon_length, length(grep("Tce_F_B_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_B_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_B_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_B_Ad_", colnames(dat_to_Tce))))
y_TceTms_DG_Ad_to_Tce_F_FPKM     <- filt_and_norm_male_female_FPKM_SA(y_TceTms_DG_Ad_to_Tce_UF, FPKM_filt_value, 0, dat_to_Tce$total_exon_length, length(grep("Tce_F_DG_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_DG_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_DG_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_DG_Ad_", colnames(dat_to_Tce))))
y_TceTms_FB_Ad_to_Tce_F_FPKM     <- filt_and_norm_male_female_FPKM_SA(y_TceTms_FB_Ad_to_Tce_UF, FPKM_filt_value, 0, dat_to_Tce$total_exon_length, length(grep("Tce_F_FB_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_FB_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_FB_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_FB_Ad_", colnames(dat_to_Tce))))
y_TceTms_Fe_Ad_to_Tce_F_FPKM     <- filt_and_norm_male_female_FPKM_SA(y_TceTms_Fe_Ad_to_Tce_UF, FPKM_filt_value, 0, dat_to_Tce$total_exon_length, length(grep("Tce_F_Fe_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_Fe_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_Fe_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_Fe_Ad_", colnames(dat_to_Tce))))
y_TceTms_Gu_Ad_to_Tce_F_FPKM     <- filt_and_norm_male_female_FPKM_SA(y_TceTms_Gu_Ad_to_Tce_UF, FPKM_filt_value, 0, dat_to_Tce$total_exon_length, length(grep("Tce_F_Gu_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_Gu_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_Gu_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_Gu_Ad_", colnames(dat_to_Tce))))
y_TceTms_Ta_Ad_to_Tce_F_FPKM     <- filt_and_norm_male_female_FPKM_SA(y_TceTms_Ta_Ad_to_Tce_UF, FPKM_filt_value, 0, dat_to_Tce$total_exon_length, length(grep("Tce_F_Ta_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_Ta_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_Ta_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_Ta_Ad_", colnames(dat_to_Tce))))
y_TceTms_GoTe_Ad_to_Tce_F_FPKM   <- filt_and_norm_male_female_FPKM_SA(y_TceTms_GoTe_Ad_to_Tce_UF, FPKM_filt_value, 0, dat_to_Tce$total_exon_length, length(grep("Tce_F_Go_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_Te_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_Go_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_Te_Ad_", colnames(dat_to_Tce)))) ### Te and Go together





#################################################################################################################################################################
### get FPKM and TPMs


get_FPKM <- function(y_F, full_df, female_prefix, male_prefix){
  y_F_genenames_df = as.data.frame(y_F$genes)
  y_F_genenames_df$genes <- as.character(y_F_genenames_df$genes)
  subset_df <- subset(full_df, full_df$Gene_name %in% y_F_genenames_df$genes)
  FPKM_df   <- rpkm(y_F, gene.length=subset_df$total_exon_length, normalized.lib.sizes=TRUE, log=FALSE)
  
  col_out = c()
  for(i in colnames(FPKM_df)){
    col_out <- c(col_out, paste(i, "FPKM", sep = "_"))
  }
  colnames(FPKM_df) <- col_out
  
  out_df   <- cbind(subset_df[,1], subset_df[,(length(subset_df) - 4): length(subset_df)], FPKM_df)
  col_out_2 <- colnames(out_df )
  
  #### calc mean FPKM for males and females
  
  subData_F <- out_df[,grepl(female_prefix,colnames(out_df))]
  print("cols for females")
  print(colnames(subData_F))
  
  subData_M <- out_df[,grepl(male_prefix,colnames(out_df))]
  print("cols for males")
  print(colnames(subData_M))
  
  out_df$female_FPKM_mean <-rowMeans(subData_F)
  out_df$male_FPKM_mean   <-rowMeans(subData_M)	
  out_df$FPKM_mean <- (out_df$female_FPKM_mean + out_df$male_FPKM_mean) / 2	
  colnames(out_df) <- c(col_out_2, paste(female_prefix, "meanFPKM", sep = "_"), paste(male_prefix, "meanFPKM", sep = "_"),  paste(str_replace(male_prefix, "M_", ""), "meanFPKM", sep = "_"))
  return(out_df )		
  
}

#########################################################################################################################################################################################
### extracting sep - but same genes run in same model
### to_Tps 
### _A_
Tps_A_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_A_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tps_F_A_Ad_", "Tps_M_A_Ad_")
Tps_A_Ad_to_Tps_F_FPKM$log2MF <- log2(Tps_A_Ad_to_Tps_F_FPKM$Tps_M_A_Ad__meanFPKM / Tps_A_Ad_to_Tps_F_FPKM$Tps_F_A_Ad__meanFPKM)
Tdi_A_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_A_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tdi_F_A_Ad_", "Tdi_M_A_Ad_")
Tdi_A_Ad_to_Tps_F_FPKM$log2MF <- log2(Tdi_A_Ad_to_Tps_F_FPKM$Tdi_M_A_Ad__meanFPKM / Tdi_A_Ad_to_Tps_F_FPKM$Tdi_F_A_Ad__meanFPKM)

#write.csv(Tdi_A_Ad_to_Tps_F_FPKM, "temp.csv" )
length(Tps_A_Ad_to_Tps_F_FPKM$log2MF)
length(Tdi_A_Ad_to_Tps_F_FPKM$log2MF)

### _B_
Tps_B_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_B_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tps_F_B_Ad_", "Tps_M_B_Ad_")
Tps_B_Ad_to_Tps_F_FPKM$log2MF <- log2(Tps_B_Ad_to_Tps_F_FPKM$Tps_M_B_Ad__meanFPKM / Tps_B_Ad_to_Tps_F_FPKM$Tps_F_B_Ad__meanFPKM)
Tdi_B_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_B_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tdi_F_B_Ad_", "Tdi_M_B_Ad_")
Tdi_B_Ad_to_Tps_F_FPKM$log2MF <- log2(Tdi_B_Ad_to_Tps_F_FPKM$Tdi_M_B_Ad__meanFPKM / Tdi_B_Ad_to_Tps_F_FPKM$Tdi_F_B_Ad__meanFPKM)

head(Tdi_B_Ad_to_Tps_F_FPKM)
head(Tps_B_Ad_to_Tps_F_FPKM)

### _DG_
Tps_DG_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_DG_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tps_F_DG_Ad_", "Tps_M_DG_Ad_")
Tps_DG_Ad_to_Tps_F_FPKM$log2MF <- log2(Tps_DG_Ad_to_Tps_F_FPKM$Tps_M_DG_Ad__meanFPKM / Tps_DG_Ad_to_Tps_F_FPKM$Tps_F_DG_Ad__meanFPKM)
Tdi_DG_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_DG_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tdi_F_DG_Ad_", "Tdi_M_DG_Ad_")
Tdi_DG_Ad_to_Tps_F_FPKM$log2MF <- log2(Tdi_DG_Ad_to_Tps_F_FPKM$Tdi_M_DG_Ad__meanFPKM / Tdi_DG_Ad_to_Tps_F_FPKM$Tdi_F_DG_Ad__meanFPKM)

head(Tdi_DG_Ad_to_Tps_F_FPKM)
head(Tps_DG_Ad_to_Tps_F_FPKM)

### _FB_
Tps_FB_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_FB_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tps_F_FB_Ad_", "Tps_M_FB_Ad_")
Tps_FB_Ad_to_Tps_F_FPKM$log2MF <- log2(Tps_FB_Ad_to_Tps_F_FPKM$Tps_M_FB_Ad__meanFPKM / Tps_FB_Ad_to_Tps_F_FPKM$Tps_F_FB_Ad__meanFPKM)
Tdi_FB_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_FB_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tdi_F_FB_Ad_", "Tdi_M_FB_Ad_")
Tdi_FB_Ad_to_Tps_F_FPKM$log2MF <- log2(Tdi_FB_Ad_to_Tps_F_FPKM$Tdi_M_FB_Ad__meanFPKM / Tdi_FB_Ad_to_Tps_F_FPKM$Tdi_F_FB_Ad__meanFPKM)

head(Tdi_FB_Ad_to_Tps_F_FPKM)
head(Tps_FB_Ad_to_Tps_F_FPKM)

### _Fe_
Tps_Fe_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_Fe_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tps_F_Fe_Ad_", "Tps_M_Fe_Ad_")
Tps_Fe_Ad_to_Tps_F_FPKM$log2MF <- log2(Tps_Fe_Ad_to_Tps_F_FPKM$Tps_M_Fe_Ad__meanFPKM / Tps_Fe_Ad_to_Tps_F_FPKM$Tps_F_Fe_Ad__meanFPKM)
Tdi_Fe_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_Fe_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tdi_F_Fe_Ad_", "Tdi_M_Fe_Ad_")
Tdi_Fe_Ad_to_Tps_F_FPKM$log2MF <- log2(Tdi_Fe_Ad_to_Tps_F_FPKM$Tdi_M_Fe_Ad__meanFPKM / Tdi_Fe_Ad_to_Tps_F_FPKM$Tdi_F_Fe_Ad__meanFPKM)

head(Tdi_Fe_Ad_to_Tps_F_FPKM)
head(Tps_Fe_Ad_to_Tps_F_FPKM)

### _Gu_
Tps_Gu_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_Gu_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tps_F_Gu_Ad_", "Tps_M_Gu_Ad_")
Tps_Gu_Ad_to_Tps_F_FPKM$log2MF <- log2(Tps_Gu_Ad_to_Tps_F_FPKM$Tps_M_Gu_Ad__meanFPKM / Tps_Gu_Ad_to_Tps_F_FPKM$Tps_F_Gu_Ad__meanFPKM)
Tdi_Gu_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_Gu_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tdi_F_Gu_Ad_", "Tdi_M_Gu_Ad_")
Tdi_Gu_Ad_to_Tps_F_FPKM$log2MF <- log2(Tdi_Gu_Ad_to_Tps_F_FPKM$Tdi_M_Gu_Ad__meanFPKM / Tdi_Gu_Ad_to_Tps_F_FPKM$Tdi_F_Gu_Ad__meanFPKM)

head(Tdi_Gu_Ad_to_Tps_F_FPKM)
head(Tps_Gu_Ad_to_Tps_F_FPKM)

### _Ta_
Tps_Ta_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_Ta_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tps_F_Ta_Ad_", "Tps_M_Ta_Ad_")
Tps_Ta_Ad_to_Tps_F_FPKM$log2MF <- log2(Tps_Ta_Ad_to_Tps_F_FPKM$Tps_M_Ta_Ad__meanFPKM / Tps_Ta_Ad_to_Tps_F_FPKM$Tps_F_Ta_Ad__meanFPKM)
Tdi_Ta_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_Ta_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tdi_F_Ta_Ad_", "Tdi_M_Ta_Ad_")
Tdi_Ta_Ad_to_Tps_F_FPKM$log2MF <- log2(Tdi_Ta_Ad_to_Tps_F_FPKM$Tdi_M_Ta_Ad__meanFPKM / Tdi_Ta_Ad_to_Tps_F_FPKM$Tdi_F_Ta_Ad__meanFPKM)

head(Tdi_Ta_Ad_to_Tps_F_FPKM)
head(Tps_Ta_Ad_to_Tps_F_FPKM)

### _GoTe_
Tps_GoTe_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_GoTe_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tps_F_Go_Ad_", "Tps_M_Te_Ad_")
Tps_GoTe_Ad_to_Tps_F_FPKM$log2MF <- log2(Tps_GoTe_Ad_to_Tps_F_FPKM$Tps_M_Te_Ad__meanFPKM / Tps_GoTe_Ad_to_Tps_F_FPKM$Tps_F_Go_Ad__meanFPKM)
Tdi_GoTe_Ad_to_Tps_F_FPKM        <- get_FPKM(y_TpsTdi_GoTe_Ad_to_Tps_F_FPKM, dat_to_Tps, "Tdi_F_Go_Ad_", "Tdi_M_Te_Ad_")
Tdi_GoTe_Ad_to_Tps_F_FPKM$log2MF <- log2(Tdi_GoTe_Ad_to_Tps_F_FPKM$Tdi_M_Te_Ad__meanFPKM / Tdi_GoTe_Ad_to_Tps_F_FPKM$Tdi_F_Go_Ad__meanFPKM)

head(Tdi_GoTe_Ad_to_Tps_F_FPKM)
head(Tps_GoTe_Ad_to_Tps_F_FPKM)


#########################################################################################################################################################################################
### to_Tpa 
### _A_
Tpa_A_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_A_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tpa_F_A_Ad_", "Tpa_M_A_Ad_")
Tpa_A_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tpa_A_Ad_to_Tpa_F_FPKM$Tpa_M_A_Ad__meanFPKM / Tpa_A_Ad_to_Tpa_F_FPKM$Tpa_F_A_Ad__meanFPKM)
Tge_A_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_A_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tge_F_A_Ad_", "Tge_M_A_Ad_")
Tge_A_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tge_A_Ad_to_Tpa_F_FPKM$Tge_M_A_Ad__meanFPKM / Tge_A_Ad_to_Tpa_F_FPKM$Tge_F_A_Ad__meanFPKM)

head(Tge_A_Ad_to_Tpa_F_FPKM)
head(Tpa_A_Ad_to_Tpa_F_FPKM)

### _B_
Tpa_B_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_B_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tpa_F_B_Ad_", "Tpa_M_B_Ad_")
Tpa_B_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tpa_B_Ad_to_Tpa_F_FPKM$Tpa_M_B_Ad__meanFPKM / Tpa_B_Ad_to_Tpa_F_FPKM$Tpa_F_B_Ad__meanFPKM)
Tge_B_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_B_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tge_F_B_Ad_", "Tge_M_B_Ad_")
Tge_B_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tge_B_Ad_to_Tpa_F_FPKM$Tge_M_B_Ad__meanFPKM / Tge_B_Ad_to_Tpa_F_FPKM$Tge_F_B_Ad__meanFPKM)

head(Tge_B_Ad_to_Tpa_F_FPKM)
head(Tpa_B_Ad_to_Tpa_F_FPKM)

### _DG_
Tpa_DG_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_DG_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tpa_F_DG_Ad_", "Tpa_M_DG_Ad_")
Tpa_DG_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tpa_DG_Ad_to_Tpa_F_FPKM$Tpa_M_DG_Ad__meanFPKM / Tpa_DG_Ad_to_Tpa_F_FPKM$Tpa_F_DG_Ad__meanFPKM)
Tge_DG_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_DG_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tge_F_DG_Ad_", "Tge_M_DG_Ad_")
Tge_DG_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tge_DG_Ad_to_Tpa_F_FPKM$Tge_M_DG_Ad__meanFPKM / Tge_DG_Ad_to_Tpa_F_FPKM$Tge_F_DG_Ad__meanFPKM)

head(Tge_DG_Ad_to_Tpa_F_FPKM)
head(Tpa_DG_Ad_to_Tpa_F_FPKM)

### _FB_
Tpa_FB_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_FB_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tpa_F_FB_Ad_", "Tpa_M_FB_Ad_")
Tpa_FB_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tpa_FB_Ad_to_Tpa_F_FPKM$Tpa_M_FB_Ad__meanFPKM / Tpa_FB_Ad_to_Tpa_F_FPKM$Tpa_F_FB_Ad__meanFPKM)
Tge_FB_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_FB_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tge_F_FB_Ad_", "Tge_M_FB_Ad_")
Tge_FB_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tge_FB_Ad_to_Tpa_F_FPKM$Tge_M_FB_Ad__meanFPKM / Tge_FB_Ad_to_Tpa_F_FPKM$Tge_F_FB_Ad__meanFPKM)

head(Tge_FB_Ad_to_Tpa_F_FPKM)
head(Tpa_FB_Ad_to_Tpa_F_FPKM)

### _Fe_
Tpa_Fe_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_Fe_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tpa_F_Fe_Ad_", "Tpa_M_Fe_Ad_")
Tpa_Fe_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tpa_Fe_Ad_to_Tpa_F_FPKM$Tpa_M_Fe_Ad__meanFPKM / Tpa_Fe_Ad_to_Tpa_F_FPKM$Tpa_F_Fe_Ad__meanFPKM)
Tge_Fe_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_Fe_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tge_F_Fe_Ad_", "Tge_M_Fe_Ad_")
Tge_Fe_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tge_Fe_Ad_to_Tpa_F_FPKM$Tge_M_Fe_Ad__meanFPKM / Tge_Fe_Ad_to_Tpa_F_FPKM$Tge_F_Fe_Ad__meanFPKM)

head(Tge_Fe_Ad_to_Tpa_F_FPKM)
head(Tpa_Fe_Ad_to_Tpa_F_FPKM)

### _Gu_
Tpa_Gu_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_Gu_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tpa_F_Gu_Ad_", "Tpa_M_Gu_Ad_")
Tpa_Gu_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tpa_Gu_Ad_to_Tpa_F_FPKM$Tpa_M_Gu_Ad__meanFPKM / Tpa_Gu_Ad_to_Tpa_F_FPKM$Tpa_F_Gu_Ad__meanFPKM)
Tge_Gu_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_Gu_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tge_F_Gu_Ad_", "Tge_M_Gu_Ad_")
Tge_Gu_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tge_Gu_Ad_to_Tpa_F_FPKM$Tge_M_Gu_Ad__meanFPKM / Tge_Gu_Ad_to_Tpa_F_FPKM$Tge_F_Gu_Ad__meanFPKM)

head(Tge_Gu_Ad_to_Tpa_F_FPKM)
head(Tpa_Gu_Ad_to_Tpa_F_FPKM)

### _Ta_
Tpa_Ta_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_Ta_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tpa_F_Ta_Ad_", "Tpa_M_Ta_Ad_")
Tpa_Ta_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tpa_Ta_Ad_to_Tpa_F_FPKM$Tpa_M_Ta_Ad__meanFPKM / Tpa_Ta_Ad_to_Tpa_F_FPKM$Tpa_F_Ta_Ad__meanFPKM)
Tge_Ta_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_Ta_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tge_F_Ta_Ad_", "Tge_M_Ta_Ad_")
Tge_Ta_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tge_Ta_Ad_to_Tpa_F_FPKM$Tge_M_Ta_Ad__meanFPKM / Tge_Ta_Ad_to_Tpa_F_FPKM$Tge_F_Ta_Ad__meanFPKM)

head(Tge_Ta_Ad_to_Tpa_F_FPKM)
head(Tpa_Ta_Ad_to_Tpa_F_FPKM)


### _GoTe_
Tpa_GoTe_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_GoTe_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tpa_F_Go_Ad_", "Tpa_M_Te_Ad_")
Tpa_GoTe_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tpa_GoTe_Ad_to_Tpa_F_FPKM$Tpa_M_Te_Ad__meanFPKM / Tpa_GoTe_Ad_to_Tpa_F_FPKM$Tpa_F_Go_Ad__meanFPKM)
Tge_GoTe_Ad_to_Tpa_F_FPKM        <- get_FPKM(y_TpaTge_GoTe_Ad_to_Tpa_F_FPKM, dat_to_Tpa, "Tge_F_Go_Ad_", "Tge_M_Te_Ad_")
Tge_GoTe_Ad_to_Tpa_F_FPKM$log2MF <- log2(Tge_GoTe_Ad_to_Tpa_F_FPKM$Tge_M_Te_Ad__meanFPKM / Tge_GoTe_Ad_to_Tpa_F_FPKM$Tge_F_Go_Ad__meanFPKM)

head(Tge_GoTe_Ad_to_Tpa_F_FPKM)
head(Tpa_GoTe_Ad_to_Tpa_F_FPKM)

#########################################################################################################################################################################################
### to_Tce 
### _A_
Tce_A_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_A_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tce_F_A_Ad_", "Tce_M_A_Ad_")
Tce_A_Ad_to_Tce_F_FPKM$log2MF <- log2(Tce_A_Ad_to_Tce_F_FPKM$Tce_M_A_Ad__meanFPKM / Tce_A_Ad_to_Tce_F_FPKM$Tce_F_A_Ad__meanFPKM)
Tms_A_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_A_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tms_F_A_Ad_", "Tms_M_A_Ad_")
Tms_A_Ad_to_Tce_F_FPKM$log2MF <- log2(Tms_A_Ad_to_Tce_F_FPKM$Tms_M_A_Ad__meanFPKM / Tms_A_Ad_to_Tce_F_FPKM$Tms_F_A_Ad__meanFPKM)

head(Tms_A_Ad_to_Tce_F_FPKM)
head(Tce_A_Ad_to_Tce_F_FPKM)

### _B_
Tce_B_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_B_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tce_F_B_Ad_", "Tce_M_B_Ad_")
Tce_B_Ad_to_Tce_F_FPKM$log2MF <- log2(Tce_B_Ad_to_Tce_F_FPKM$Tce_M_B_Ad__meanFPKM / Tce_B_Ad_to_Tce_F_FPKM$Tce_F_B_Ad__meanFPKM)
Tms_B_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_B_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tms_F_B_Ad_", "Tms_M_B_Ad_")
Tms_B_Ad_to_Tce_F_FPKM$log2MF <- log2(Tms_B_Ad_to_Tce_F_FPKM$Tms_M_B_Ad__meanFPKM / Tms_B_Ad_to_Tce_F_FPKM$Tms_F_B_Ad__meanFPKM)

head(Tms_B_Ad_to_Tce_F_FPKM)
head(Tce_B_Ad_to_Tce_F_FPKM)

### _DG_
Tce_DG_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_DG_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tce_F_DG_Ad_", "Tce_M_DG_Ad_")
Tce_DG_Ad_to_Tce_F_FPKM$log2MF <- log2(Tce_DG_Ad_to_Tce_F_FPKM$Tce_M_DG_Ad__meanFPKM / Tce_DG_Ad_to_Tce_F_FPKM$Tce_F_DG_Ad__meanFPKM)
Tms_DG_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_DG_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tms_F_DG_Ad_", "Tms_M_DG_Ad_")
Tms_DG_Ad_to_Tce_F_FPKM$log2MF <- log2(Tms_DG_Ad_to_Tce_F_FPKM$Tms_M_DG_Ad__meanFPKM / Tms_DG_Ad_to_Tce_F_FPKM$Tms_F_DG_Ad__meanFPKM)

head(Tms_DG_Ad_to_Tce_F_FPKM)
head(Tce_DG_Ad_to_Tce_F_FPKM)

### _FB_
Tce_FB_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_FB_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tce_F_FB_Ad_", "Tce_M_FB_Ad_")
Tce_FB_Ad_to_Tce_F_FPKM$log2MF <- log2(Tce_FB_Ad_to_Tce_F_FPKM$Tce_M_FB_Ad__meanFPKM / Tce_FB_Ad_to_Tce_F_FPKM$Tce_F_FB_Ad__meanFPKM)
Tms_FB_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_FB_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tms_F_FB_Ad_", "Tms_M_FB_Ad_")
Tms_FB_Ad_to_Tce_F_FPKM$log2MF <- log2(Tms_FB_Ad_to_Tce_F_FPKM$Tms_M_FB_Ad__meanFPKM / Tms_FB_Ad_to_Tce_F_FPKM$Tms_F_FB_Ad__meanFPKM)

head(Tms_FB_Ad_to_Tce_F_FPKM)
head(Tce_FB_Ad_to_Tce_F_FPKM)

### _Fe_
Tce_Fe_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_Fe_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tce_F_Fe_Ad_", "Tce_M_Fe_Ad_")
Tce_Fe_Ad_to_Tce_F_FPKM$log2MF <- log2(Tce_Fe_Ad_to_Tce_F_FPKM$Tce_M_Fe_Ad__meanFPKM / Tce_Fe_Ad_to_Tce_F_FPKM$Tce_F_Fe_Ad__meanFPKM)
Tms_Fe_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_Fe_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tms_F_Fe_Ad_", "Tms_M_Fe_Ad_")
Tms_Fe_Ad_to_Tce_F_FPKM$log2MF <- log2(Tms_Fe_Ad_to_Tce_F_FPKM$Tms_M_Fe_Ad__meanFPKM / Tms_Fe_Ad_to_Tce_F_FPKM$Tms_F_Fe_Ad__meanFPKM)

head(Tms_Fe_Ad_to_Tce_F_FPKM)
head(Tce_Fe_Ad_to_Tce_F_FPKM)

### _Gu_
Tce_Gu_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_Gu_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tce_F_Gu_Ad_", "Tce_M_Gu_Ad_")
Tce_Gu_Ad_to_Tce_F_FPKM$log2MF <- log2(Tce_Gu_Ad_to_Tce_F_FPKM$Tce_M_Gu_Ad__meanFPKM / Tce_Gu_Ad_to_Tce_F_FPKM$Tce_F_Gu_Ad__meanFPKM)
Tms_Gu_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_Gu_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tms_F_Gu_Ad_", "Tms_M_Gu_Ad_")
Tms_Gu_Ad_to_Tce_F_FPKM$log2MF <- log2(Tms_Gu_Ad_to_Tce_F_FPKM$Tms_M_Gu_Ad__meanFPKM / Tms_Gu_Ad_to_Tce_F_FPKM$Tms_F_Gu_Ad__meanFPKM)

head(Tms_Gu_Ad_to_Tce_F_FPKM)
head(Tce_Gu_Ad_to_Tce_F_FPKM)

### _Ta_
Tce_Ta_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_Ta_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tce_F_Ta_Ad_", "Tce_M_Ta_Ad_")
Tce_Ta_Ad_to_Tce_F_FPKM$log2MF <- log2(Tce_Ta_Ad_to_Tce_F_FPKM$Tce_M_Ta_Ad__meanFPKM / Tce_Ta_Ad_to_Tce_F_FPKM$Tce_F_Ta_Ad__meanFPKM)
Tms_Ta_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_Ta_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tms_F_Ta_Ad_", "Tms_M_Ta_Ad_")
Tms_Ta_Ad_to_Tce_F_FPKM$log2MF <- log2(Tms_Ta_Ad_to_Tce_F_FPKM$Tms_M_Ta_Ad__meanFPKM / Tms_Ta_Ad_to_Tce_F_FPKM$Tms_F_Ta_Ad__meanFPKM)

head(Tms_Ta_Ad_to_Tce_F_FPKM)
head(Tce_Ta_Ad_to_Tce_F_FPKM)

### _GoTe_
Tce_GoTe_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_GoTe_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tce_F_Go_Ad_", "Tce_M_Te_Ad_")
Tce_GoTe_Ad_to_Tce_F_FPKM$log2MF <- log2(Tce_GoTe_Ad_to_Tce_F_FPKM$Tce_M_Te_Ad__meanFPKM / Tce_GoTe_Ad_to_Tce_F_FPKM$Tce_F_Go_Ad__meanFPKM)
Tms_GoTe_Ad_to_Tce_F_FPKM        <- get_FPKM(y_TceTms_GoTe_Ad_to_Tce_F_FPKM, dat_to_Tce, "Tms_F_Go_Ad_", "Tms_M_Te_Ad_")
Tms_GoTe_Ad_to_Tce_F_FPKM$log2MF <- log2(Tms_GoTe_Ad_to_Tce_F_FPKM$Tms_M_Te_Ad__meanFPKM / Tms_GoTe_Ad_to_Tce_F_FPKM$Tms_F_Go_Ad__meanFPKM)

head(Tms_GoTe_Ad_to_Tce_F_FPKM)
head(Tce_GoTe_Ad_to_Tce_F_FPKM)


###################################################################################################################################################
###### get average FPKM

extract_avEXP = function(dat, sex){
  
  full_df_name <- deparse(substitute(dat))
  sp   <- strsplit(full_df_name, "_")[[1]][1]
  tiss <- strsplit(full_df_name, "_")[[1]][2]
  
  print(full_df_name)
  print(sp)
  print(tiss)
  
  df_out = as.data.frame(cbind(
    dat[,1],
    dat$scaf,
    dat$XA,
    dat$mid_point,
    as.character(eval(parse(text=paste('dat','$', sp, '_', sex, '_', tiss, '_Ad__meanFPKM', sep='')))),
    rep(sp, length(dat[,1])),
    rep(tiss, length(dat[,1])),
    rep(sex, length(dat[,1]))
  ))
  
  colnames(df_out) <- c("gene_id", "scaf", "XA", "mid_point", "meanFPKM", "sp", "tiss", "sex")
  
  df_out$rep_m <- ifelse(df_out$sp == "Tps", "sex", 
                         ifelse(df_out$sp == "Tdi", "asex", 
                                ifelse(df_out$sp == "Tpa", "sex", 
                                       ifelse(df_out$sp == "Tge", "asex", 
                                              ifelse(df_out$sp == "Tce", "sex", 
                                                     ifelse(df_out$sp == "Tms", "asex", 
                                                            ifelse(df_out$sp == "Tcm", "sex", 
                                                                   ifelse(df_out$sp == "Tsi", "asex",                     
                                                                          ifelse(df_out$sp == "Tbi", "sex", 
                                                                                 ifelse(df_out$sp == "Tte", "asex",  "ERROR")))))))))) 
  df_out$meanFPKM <- as.numeric(df_out$meanFPKM)
  df_out$sex_XA   <- paste(df_out$sex, df_out$XA, sep = "_")
  df_out$LG <- dat$LG
  return(df_out)
  
}


extract_avEXP_GoTe = function(dat, sex, tiss){
  
  full_df_name <- deparse(substitute(dat))
  sp   <- strsplit(full_df_name, "_")[[1]][1]
  
  print(full_df_name)
  print(sp)
  print(tiss)
  
  df_out = as.data.frame(cbind(
    dat[,1],
    dat$scaf,
    dat$XA,
    dat$mid_point,
    as.character(eval(parse(text=paste('dat','$', sp, '_', sex, '_', tiss, '_Ad__meanFPKM', sep='')))),
    rep(sp, length(dat[,1])),
    rep(tiss, length(dat[,1])),
    rep(sex, length(dat[,1]))
  ))
  
  colnames(df_out) <- c("gene_id", "scaf", "XA", "mid_point", "meanFPKM", "sp", "tiss", "sex")
  
  df_out$rep_m <- ifelse(df_out$sp == "Tps", "sex", 
                         ifelse(df_out$sp == "Tdi", "asex", 
                                ifelse(df_out$sp == "Tpa", "sex", 
                                       ifelse(df_out$sp == "Tge", "asex", 
                                              ifelse(df_out$sp == "Tce", "sex", 
                                                     ifelse(df_out$sp == "Tms", "asex", 
                                                            ifelse(df_out$sp == "Tcm", "sex", 
                                                                   ifelse(df_out$sp == "Tsi", "asex",                     
                                                                          ifelse(df_out$sp == "Tbi", "sex", 
                                                                                 ifelse(df_out$sp == "Tte", "asex",  "ERROR")))))))))) 
  df_out$meanFPKM <- as.numeric(df_out$meanFPKM)
  df_out$sex_XA   <- paste(df_out$sex, df_out$XA, sep = "_")
  df_out$LG <- dat$LG
  return(df_out)
  
}


#### calc Tau
## following Kryuchkova-Mostacci, Nadezda, and Marc Robinson-Rechavi. 2017. “A Benchmark of Gene Expression Tissue-Specificity Metrics.” Briefings in Bioinformatics 18 (2): 205–14.


fTau <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(max(x)!=0)
      {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    } 
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  } 
  return(res)
}



Tau_calc <- function(df){
  
  exp_df <- df[2:length(colnames(df))] ### don't want first col
  print("Using for Tau cals:")
  print(colnames(exp_df))
  
  Tau_vals <- c()
  for(i in 1:length(exp_df[,1])){
    exp_v <- exp_df[i,]
    #print(exp_v)
    exp_v <- as.numeric(exp_v)
    T_val <- fTau(exp_v)
    #print(T_val)
    Tau_vals <- c(Tau_vals, T_val)
  }
  
  df$Tau <-  Tau_vals
  
  return(df)
}



### TpsTdi
TpsTdi_Ad_to_Tps_F_avFPKM_all <- as.data.frame(cbind(
  extract_avEXP(Tps_A_Ad_to_Tps_F_FPKM, "M"),
  extract_avEXP(Tps_A_Ad_to_Tps_F_FPKM, "F"),
  extract_avEXP(Tps_B_Ad_to_Tps_F_FPKM, "M"),
  extract_avEXP(Tps_B_Ad_to_Tps_F_FPKM, "F"),  
  extract_avEXP(Tps_DG_Ad_to_Tps_F_FPKM, "M"),
  extract_avEXP(Tps_DG_Ad_to_Tps_F_FPKM, "F"),  
  extract_avEXP(Tps_FB_Ad_to_Tps_F_FPKM, "M"),
  extract_avEXP(Tps_FB_Ad_to_Tps_F_FPKM, "F"),  
  extract_avEXP(Tps_Fe_Ad_to_Tps_F_FPKM, "M"),
  extract_avEXP(Tps_Fe_Ad_to_Tps_F_FPKM, "F"),  
  extract_avEXP(Tps_Gu_Ad_to_Tps_F_FPKM, "M"),
  extract_avEXP(Tps_Gu_Ad_to_Tps_F_FPKM, "F"),  
  extract_avEXP(Tps_Ta_Ad_to_Tps_F_FPKM, "M"),
  extract_avEXP(Tps_Ta_Ad_to_Tps_F_FPKM, "F"), 
  extract_avEXP_GoTe(Tps_GoTe_Ad_to_Tps_F_FPKM, "F", "Go"),  
  extract_avEXP_GoTe(Tps_GoTe_Ad_to_Tps_F_FPKM, "M", "Te"),
  extract_avEXP(Tdi_A_Ad_to_Tps_F_FPKM, "M"),
  extract_avEXP(Tdi_A_Ad_to_Tps_F_FPKM, "F"),
  extract_avEXP(Tdi_B_Ad_to_Tps_F_FPKM, "M"),
  extract_avEXP(Tdi_B_Ad_to_Tps_F_FPKM, "F"),  
  extract_avEXP(Tdi_DG_Ad_to_Tps_F_FPKM, "M"),
  extract_avEXP(Tdi_DG_Ad_to_Tps_F_FPKM, "F"),  
  extract_avEXP(Tdi_FB_Ad_to_Tps_F_FPKM, "M"),
  extract_avEXP(Tdi_FB_Ad_to_Tps_F_FPKM, "F"),  
  extract_avEXP(Tdi_Fe_Ad_to_Tps_F_FPKM, "M"),
  extract_avEXP(Tdi_Fe_Ad_to_Tps_F_FPKM, "F"),  
  extract_avEXP(Tdi_Gu_Ad_to_Tps_F_FPKM, "M"),
  extract_avEXP(Tdi_Gu_Ad_to_Tps_F_FPKM, "F"),  
  extract_avEXP(Tdi_Ta_Ad_to_Tps_F_FPKM, "M"),
  extract_avEXP(Tdi_Ta_Ad_to_Tps_F_FPKM, "F"), 
  extract_avEXP_GoTe(Tdi_GoTe_Ad_to_Tps_F_FPKM, "F", "Go"),  
  extract_avEXP_GoTe(Tdi_GoTe_Ad_to_Tps_F_FPKM, "M", "Te")
))


##########################################################################################################################
## Tpa
TpaTge_Ad_to_Tpa_F_avFPKM_all <- as.data.frame(cbind(
  extract_avEXP(Tpa_A_Ad_to_Tpa_F_FPKM, "M"),
  extract_avEXP(Tpa_A_Ad_to_Tpa_F_FPKM, "F"),
  extract_avEXP(Tpa_B_Ad_to_Tpa_F_FPKM, "M"),
  extract_avEXP(Tpa_B_Ad_to_Tpa_F_FPKM, "F"),  
  extract_avEXP(Tpa_DG_Ad_to_Tpa_F_FPKM, "M"),
  extract_avEXP(Tpa_DG_Ad_to_Tpa_F_FPKM, "F"),  
  extract_avEXP(Tpa_FB_Ad_to_Tpa_F_FPKM, "M"),
  extract_avEXP(Tpa_FB_Ad_to_Tpa_F_FPKM, "F"),  
  extract_avEXP(Tpa_Fe_Ad_to_Tpa_F_FPKM, "M"),
  extract_avEXP(Tpa_Fe_Ad_to_Tpa_F_FPKM, "F"),  
  extract_avEXP(Tpa_Gu_Ad_to_Tpa_F_FPKM, "M"),
  extract_avEXP(Tpa_Gu_Ad_to_Tpa_F_FPKM, "F"),  
  extract_avEXP(Tpa_Ta_Ad_to_Tpa_F_FPKM, "M"),
  extract_avEXP(Tpa_Ta_Ad_to_Tpa_F_FPKM, "F"), 
  extract_avEXP_GoTe(Tpa_GoTe_Ad_to_Tpa_F_FPKM, "F", "Go"),  
  extract_avEXP_GoTe(Tpa_GoTe_Ad_to_Tpa_F_FPKM, "M", "Te"),
  extract_avEXP(Tge_A_Ad_to_Tpa_F_FPKM, "M"),
  extract_avEXP(Tge_A_Ad_to_Tpa_F_FPKM, "F"),
  extract_avEXP(Tge_B_Ad_to_Tpa_F_FPKM, "M"),
  extract_avEXP(Tge_B_Ad_to_Tpa_F_FPKM, "F"),  
  extract_avEXP(Tge_DG_Ad_to_Tpa_F_FPKM, "M"),
  extract_avEXP(Tge_DG_Ad_to_Tpa_F_FPKM, "F"),  
  extract_avEXP(Tge_FB_Ad_to_Tpa_F_FPKM, "M"),
  extract_avEXP(Tge_FB_Ad_to_Tpa_F_FPKM, "F"),  
  extract_avEXP(Tge_Fe_Ad_to_Tpa_F_FPKM, "M"),
  extract_avEXP(Tge_Fe_Ad_to_Tpa_F_FPKM, "F"),  
  extract_avEXP(Tge_Gu_Ad_to_Tpa_F_FPKM, "M"),
  extract_avEXP(Tge_Gu_Ad_to_Tpa_F_FPKM, "F"),  
  extract_avEXP(Tge_Ta_Ad_to_Tpa_F_FPKM, "M"),
  extract_avEXP(Tge_Ta_Ad_to_Tpa_F_FPKM, "F"), 
  extract_avEXP_GoTe(Tge_GoTe_Ad_to_Tpa_F_FPKM, "F", "Go"),  
  extract_avEXP_GoTe(Tge_GoTe_Ad_to_Tpa_F_FPKM, "M", "Te")
))



########################################################
## Tce

TceTms_Ad_to_Tce_F_avFPKM_all <- as.data.frame(cbind(
  extract_avEXP(Tce_A_Ad_to_Tce_F_FPKM, "M"),
  extract_avEXP(Tce_A_Ad_to_Tce_F_FPKM, "F"),
  extract_avEXP(Tce_B_Ad_to_Tce_F_FPKM, "M"),
  extract_avEXP(Tce_B_Ad_to_Tce_F_FPKM, "F"),  
  extract_avEXP(Tce_DG_Ad_to_Tce_F_FPKM, "M"),
  extract_avEXP(Tce_DG_Ad_to_Tce_F_FPKM, "F"),  
  extract_avEXP(Tce_FB_Ad_to_Tce_F_FPKM, "M"),
  extract_avEXP(Tce_FB_Ad_to_Tce_F_FPKM, "F"),  
  extract_avEXP(Tce_Fe_Ad_to_Tce_F_FPKM, "M"),
  extract_avEXP(Tce_Fe_Ad_to_Tce_F_FPKM, "F"),  
  extract_avEXP(Tce_Gu_Ad_to_Tce_F_FPKM, "M"),
  extract_avEXP(Tce_Gu_Ad_to_Tce_F_FPKM, "F"),  
  extract_avEXP(Tce_Ta_Ad_to_Tce_F_FPKM, "M"),
  extract_avEXP(Tce_Ta_Ad_to_Tce_F_FPKM, "F"), 
  extract_avEXP_GoTe(Tce_GoTe_Ad_to_Tce_F_FPKM, "F", "Go"),  
  extract_avEXP_GoTe(Tce_GoTe_Ad_to_Tce_F_FPKM, "M", "Te"),
  extract_avEXP(Tms_A_Ad_to_Tce_F_FPKM, "M"),
  extract_avEXP(Tms_A_Ad_to_Tce_F_FPKM, "F"),
  extract_avEXP(Tms_B_Ad_to_Tce_F_FPKM, "M"),
  extract_avEXP(Tms_B_Ad_to_Tce_F_FPKM, "F"),  
  extract_avEXP(Tms_DG_Ad_to_Tce_F_FPKM, "M"),
  extract_avEXP(Tms_DG_Ad_to_Tce_F_FPKM, "F"),  
  extract_avEXP(Tms_FB_Ad_to_Tce_F_FPKM, "M"),
  extract_avEXP(Tms_FB_Ad_to_Tce_F_FPKM, "F"),  
  extract_avEXP(Tms_Fe_Ad_to_Tce_F_FPKM, "M"),
  extract_avEXP(Tms_Fe_Ad_to_Tce_F_FPKM, "F"),  
  extract_avEXP(Tms_Gu_Ad_to_Tce_F_FPKM, "M"),
  extract_avEXP(Tms_Gu_Ad_to_Tce_F_FPKM, "F"),  
  extract_avEXP(Tms_Ta_Ad_to_Tce_F_FPKM, "M"),
  extract_avEXP(Tms_Ta_Ad_to_Tce_F_FPKM, "F"), 
  extract_avEXP_GoTe(Tms_GoTe_Ad_to_Tce_F_FPKM, "F", "Go"),  
  extract_avEXP_GoTe(Tms_GoTe_Ad_to_Tce_F_FPKM, "M", "Te")
))


#### 

### TpsTdi
TpsTdi_Ad_to_Tps_F_avFPKM <- as.data.frame(cbind(
  extract_avEXP(Tps_A_Ad_to_Tps_F_FPKM, "M")$gene_id,
  extract_avEXP(Tps_A_Ad_to_Tps_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tps_A_Ad_to_Tps_F_FPKM, "F")$meanFPKM,
  extract_avEXP(Tps_B_Ad_to_Tps_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tps_B_Ad_to_Tps_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tps_DG_Ad_to_Tps_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tps_DG_Ad_to_Tps_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tps_FB_Ad_to_Tps_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tps_FB_Ad_to_Tps_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tps_Fe_Ad_to_Tps_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tps_Fe_Ad_to_Tps_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tps_Gu_Ad_to_Tps_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tps_Gu_Ad_to_Tps_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tps_Ta_Ad_to_Tps_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tps_Ta_Ad_to_Tps_F_FPKM, "F")$meanFPKM, 
  extract_avEXP_GoTe(Tps_GoTe_Ad_to_Tps_F_FPKM, "F", "Go")$meanFPKM,  
  extract_avEXP_GoTe(Tps_GoTe_Ad_to_Tps_F_FPKM, "M", "Te")$meanFPKM,
  extract_avEXP(Tdi_A_Ad_to_Tps_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tdi_A_Ad_to_Tps_F_FPKM, "F")$meanFPKM,
  extract_avEXP(Tdi_B_Ad_to_Tps_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tdi_B_Ad_to_Tps_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tdi_DG_Ad_to_Tps_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tdi_DG_Ad_to_Tps_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tdi_FB_Ad_to_Tps_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tdi_FB_Ad_to_Tps_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tdi_Fe_Ad_to_Tps_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tdi_Fe_Ad_to_Tps_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tdi_Gu_Ad_to_Tps_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tdi_Gu_Ad_to_Tps_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tdi_Ta_Ad_to_Tps_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tdi_Ta_Ad_to_Tps_F_FPKM, "F")$meanFPKM, 
  extract_avEXP_GoTe(Tdi_GoTe_Ad_to_Tps_F_FPKM, "F", "Go")$meanFPKM,  
  extract_avEXP_GoTe(Tdi_GoTe_Ad_to_Tps_F_FPKM, "M", "Te")$meanFPKM
))



colnames(TpsTdi_Ad_to_Tps_F_avFPKM) <- c(
  "gene_id",
  "Tps_A_Ad_M",
  "Tps_A_Ad_F",
  "Tps_B_Ad_M",
  "Tps_B_Ad_F",  
  "Tps_DG_Ad_M",
  "Tps_DG_Ad_F",  
  "Tps_FB_Ad_M",
  "Tps_FB_Ad_F",  
  "Tps_Fe_Ad_M",
  "Tps_Fe_Ad_F",  
  "Tps_Gu_Ad_M",
  "Tps_Gu_Ad_F",  
  "Tps_Ta_Ad_M",
  "Tps_Ta_Ad_F", 
  "Tps_GoTe_Ad_F",  
  "Tps_GoTe_Ad_M",
  "Tdi_A_Ad_M",
  "Tdi_A_Ad_F",
  "Tdi_B_Ad_M",
  "Tdi_B_Ad_F",  
  "Tdi_DG_Ad_M",
  "Tdi_DG_Ad_F",  
  "Tdi_FB_Ad_M",
  "Tdi_FB_Ad_F",  
  "Tdi_Fe_Ad_M",
  "Tdi_Fe_Ad_F",  
  "Tdi_Gu_Ad_M",
  "Tdi_Gu_Ad_F",  
  "Tdi_Ta_Ad_M",
  "Tdi_Ta_Ad_F", 
  "Tdi_GoTe_Ad_F",  
  "Tdi_GoTe_Ad_M")

colnames(TpsTdi_Ad_to_Tps_F_avFPKM)
TpsTdi_Ad_to_Tps_F_avFPKM_wTau <- Tau_calc(TpsTdi_Ad_to_Tps_F_avFPKM)

### TpaTge
TpaTge_Ad_to_Tpa_F_avFPKM <- as.data.frame(cbind(
  extract_avEXP(Tpa_A_Ad_to_Tpa_F_FPKM, "M")$gene_id,
  extract_avEXP(Tpa_A_Ad_to_Tpa_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tpa_A_Ad_to_Tpa_F_FPKM, "F")$meanFPKM,
  extract_avEXP(Tpa_B_Ad_to_Tpa_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tpa_B_Ad_to_Tpa_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tpa_DG_Ad_to_Tpa_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tpa_DG_Ad_to_Tpa_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tpa_FB_Ad_to_Tpa_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tpa_FB_Ad_to_Tpa_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tpa_Fe_Ad_to_Tpa_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tpa_Fe_Ad_to_Tpa_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tpa_Gu_Ad_to_Tpa_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tpa_Gu_Ad_to_Tpa_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tpa_Ta_Ad_to_Tpa_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tpa_Ta_Ad_to_Tpa_F_FPKM, "F")$meanFPKM, 
  extract_avEXP_GoTe(Tpa_GoTe_Ad_to_Tpa_F_FPKM, "F", "Go")$meanFPKM,  
  extract_avEXP_GoTe(Tpa_GoTe_Ad_to_Tpa_F_FPKM, "M", "Te")$meanFPKM,
  extract_avEXP(Tge_A_Ad_to_Tpa_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tge_A_Ad_to_Tpa_F_FPKM, "F")$meanFPKM,
  extract_avEXP(Tge_B_Ad_to_Tpa_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tge_B_Ad_to_Tpa_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tge_DG_Ad_to_Tpa_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tge_DG_Ad_to_Tpa_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tge_FB_Ad_to_Tpa_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tge_FB_Ad_to_Tpa_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tge_Fe_Ad_to_Tpa_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tge_Fe_Ad_to_Tpa_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tge_Gu_Ad_to_Tpa_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tge_Gu_Ad_to_Tpa_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tge_Ta_Ad_to_Tpa_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tge_Ta_Ad_to_Tpa_F_FPKM, "F")$meanFPKM, 
  extract_avEXP_GoTe(Tge_GoTe_Ad_to_Tpa_F_FPKM, "F", "Go")$meanFPKM,  
  extract_avEXP_GoTe(Tge_GoTe_Ad_to_Tpa_F_FPKM, "M", "Te")$meanFPKM
))



colnames(TpaTge_Ad_to_Tpa_F_avFPKM) <- c(
  "gene_id",
  "Tpa_A_Ad_M",
  "Tpa_A_Ad_F",
  "Tpa_B_Ad_M",
  "Tpa_B_Ad_F",  
  "Tpa_DG_Ad_M",
  "Tpa_DG_Ad_F",  
  "Tpa_FB_Ad_M",
  "Tpa_FB_Ad_F",  
  "Tpa_Fe_Ad_M",
  "Tpa_Fe_Ad_F",  
  "Tpa_Gu_Ad_M",
  "Tpa_Gu_Ad_F",  
  "Tpa_Ta_Ad_M",
  "Tpa_Ta_Ad_F", 
  "Tpa_GoTe_Ad_F",  
  "Tpa_GoTe_Ad_M",
  "Tge_A_Ad_M",
  "Tge_A_Ad_F",
  "Tge_B_Ad_M",
  "Tge_B_Ad_F",  
  "Tge_DG_Ad_M",
  "Tge_DG_Ad_F",  
  "Tge_FB_Ad_M",
  "Tge_FB_Ad_F",  
  "Tge_Fe_Ad_M",
  "Tge_Fe_Ad_F",  
  "Tge_Gu_Ad_M",
  "Tge_Gu_Ad_F",  
  "Tge_Ta_Ad_M",
  "Tge_Ta_Ad_F", 
  "Tge_GoTe_Ad_F",  
  "Tge_GoTe_Ad_M")

colnames(TpaTge_Ad_to_Tpa_F_avFPKM)
TpaTge_Ad_to_Tpa_F_avFPKM_wTau <- Tau_calc(TpaTge_Ad_to_Tpa_F_avFPKM)



### TceTms
TceTms_Ad_to_Tce_F_avFPKM <- as.data.frame(cbind(
  extract_avEXP(Tce_A_Ad_to_Tce_F_FPKM, "M")$gene_id,
  extract_avEXP(Tce_A_Ad_to_Tce_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tce_A_Ad_to_Tce_F_FPKM, "F")$meanFPKM,
  extract_avEXP(Tce_B_Ad_to_Tce_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tce_B_Ad_to_Tce_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tce_DG_Ad_to_Tce_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tce_DG_Ad_to_Tce_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tce_FB_Ad_to_Tce_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tce_FB_Ad_to_Tce_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tce_Fe_Ad_to_Tce_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tce_Fe_Ad_to_Tce_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tce_Gu_Ad_to_Tce_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tce_Gu_Ad_to_Tce_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tce_Ta_Ad_to_Tce_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tce_Ta_Ad_to_Tce_F_FPKM, "F")$meanFPKM, 
  extract_avEXP_GoTe(Tce_GoTe_Ad_to_Tce_F_FPKM, "F", "Go")$meanFPKM,  
  extract_avEXP_GoTe(Tce_GoTe_Ad_to_Tce_F_FPKM, "M", "Te")$meanFPKM,
  extract_avEXP(Tms_A_Ad_to_Tce_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tms_A_Ad_to_Tce_F_FPKM, "F")$meanFPKM,
  extract_avEXP(Tms_B_Ad_to_Tce_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tms_B_Ad_to_Tce_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tms_DG_Ad_to_Tce_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tms_DG_Ad_to_Tce_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tms_FB_Ad_to_Tce_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tms_FB_Ad_to_Tce_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tms_Fe_Ad_to_Tce_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tms_Fe_Ad_to_Tce_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tms_Gu_Ad_to_Tce_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tms_Gu_Ad_to_Tce_F_FPKM, "F")$meanFPKM,  
  extract_avEXP(Tms_Ta_Ad_to_Tce_F_FPKM, "M")$meanFPKM,
  extract_avEXP(Tms_Ta_Ad_to_Tce_F_FPKM, "F")$meanFPKM, 
  extract_avEXP_GoTe(Tms_GoTe_Ad_to_Tce_F_FPKM, "F", "Go")$meanFPKM,  
  extract_avEXP_GoTe(Tms_GoTe_Ad_to_Tce_F_FPKM, "M", "Te")$meanFPKM
))

colnames(TceTms_Ad_to_Tce_F_avFPKM) <- c(
  "gene_id",
  "Tce_A_Ad_M",
  "Tce_A_Ad_F",
  "Tce_B_Ad_M",
  "Tce_B_Ad_F",  
  "Tce_DG_Ad_M",
  "Tce_DG_Ad_F",  
  "Tce_FB_Ad_M",
  "Tce_FB_Ad_F",  
  "Tce_Fe_Ad_M",
  "Tce_Fe_Ad_F",  
  "Tce_Gu_Ad_M",
  "Tce_Gu_Ad_F",  
  "Tce_Ta_Ad_M",
  "Tce_Ta_Ad_F", 
  "Tce_GoTe_Ad_F",  
  "Tce_GoTe_Ad_M",
  "Tms_A_Ad_M",
  "Tms_A_Ad_F",
  "Tms_B_Ad_M",
  "Tms_B_Ad_F",  
  "Tms_DG_Ad_M",
  "Tms_DG_Ad_F",  
  "Tms_FB_Ad_M",
  "Tms_FB_Ad_F",  
  "Tms_Fe_Ad_M",
  "Tms_Fe_Ad_F",  
  "Tms_Gu_Ad_M",
  "Tms_Gu_Ad_F",  
  "Tms_Ta_Ad_M",
  "Tms_Ta_Ad_F", 
  "Tms_GoTe_Ad_F",  
  "Tms_GoTe_Ad_M")

colnames(TceTms_Ad_to_Tce_F_avFPKM)
TceTms_Ad_to_Tce_F_avFPKM_wTau <- Tau_calc(TceTms_Ad_to_Tce_F_avFPKM)


#### export


TpsTdi_Ad_to_Tps_F_Tau <- as.data.frame(cbind(TpsTdi_Ad_to_Tps_F_avFPKM_wTau$gene_id, TpsTdi_Ad_to_Tps_F_avFPKM_wTau$Tau))
colnames(TpsTdi_Ad_to_Tps_F_Tau ) <- c("gene_id", "Tau")
write.csv(TpsTdi_Ad_to_Tps_F_Tau, "TpsTdi_Ad_to_Tps_wTau.csv", row.names = F)


TpaTge_Ad_to_Tpa_F_Tau <- as.data.frame(cbind(TpaTge_Ad_to_Tpa_F_avFPKM_wTau$gene_id, TpaTge_Ad_to_Tpa_F_avFPKM_wTau$Tau))
colnames(TpaTge_Ad_to_Tpa_F_Tau ) <- c("gene_id", "Tau")
write.csv(TpaTge_Ad_to_Tpa_F_Tau, "TpaTge_Ad_to_Tpa_wTau.csv", row.names = F)


TceTms_Ad_to_Tce_F_Tau <- as.data.frame(cbind(TceTms_Ad_to_Tce_F_avFPKM_wTau$gene_id, TceTms_Ad_to_Tce_F_avFPKM_wTau$Tau))
colnames(TceTms_Ad_to_Tce_F_Tau ) <- c("gene_id", "Tau")
write.csv(TceTms_Ad_to_Tce_F_Tau, "TceTms_Ad_to_Tce_wTau.csv", row.names = F)






