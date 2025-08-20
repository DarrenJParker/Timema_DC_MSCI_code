/*§#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
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
# 
# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Monterey 12.6.8
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
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
#   [1] hash_2.2.6.3         dplyr_1.1.3          Rtsne_0.16           SuperExactTest_1.1.0 raster_3.6-23       
# [6] sp_2.1-0             pvclust_2.2-0        vegan_2.6-4          permute_0.9-7        RColorBrewer_1.1-3  
# [11] pheatmap_1.0.12      gtable_0.3.4         stringr_1.5.0        cowplot_1.1.1        lattice_0.21-8      
# [16] ggplot2_3.4.3        gridExtra_2.3        VennDiagram_1.7.3    futile.logger_1.4.3  edgeR_3.42.4        
# [21] limma_3.56.2         BiocManager_1.30.22 
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.3           generics_0.1.3       futile.options_1.0.1 stringi_1.7.12       magrittr_2.0.3      
# [6] Matrix_1.5-4.1       formatR_1.14         mgcv_1.8-42          fansi_1.0.4          scales_1.2.1        
# [11] codetools_0.2-19     cli_3.6.1            rlang_1.1.1          munsell_0.5.0        splines_4.3.1       
# [16] withr_2.5.1          tools_4.3.1          parallel_4.3.1       colorspace_2.1-0     locfit_1.5-9.8      
# [21] lambda.r_1.2.4       vctrs_0.6.3          R6_2.5.1             lifecycle_1.0.3      MASS_7.3-60         
# [26] cluster_2.1.4        pkgconfig_2.0.3      terra_1.7-46         pillar_1.9.0         glue_1.6.2          
# [31] Rcpp_1.0.11          tibble_3.2.1         tidyselect_1.2.0     nlme_3.1-162         compiler_4.3.1      
# > 
### DJP code
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



setwd("data/read_counts")


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

FPKM_filt_value = 2

y_TpsTdi_A_Ad_to_Tps_F_FPKM    <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_A_Ad_to_Tps_UF, FPKM_filt_value,    2, dat_to_Tps$total_exon_length, length(grep("Tps_F_A_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_A_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_A_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_A_Ad_", colnames(dat_to_Tps))))
y_TpsTdi_B_Ad_to_Tps_F_FPKM    <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_B_Ad_to_Tps_UF, FPKM_filt_value,    2, dat_to_Tps$total_exon_length, length(grep("Tps_F_B_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_B_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_B_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_B_Ad_", colnames(dat_to_Tps))))
y_TpsTdi_DG_Ad_to_Tps_F_FPKM   <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_DG_Ad_to_Tps_UF, FPKM_filt_value,   2, dat_to_Tps$total_exon_length, length(grep("Tps_F_DG_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_DG_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_DG_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_DG_Ad_", colnames(dat_to_Tps))))
y_TpsTdi_FB_Ad_to_Tps_F_FPKM   <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_FB_Ad_to_Tps_UF, FPKM_filt_value,   2, dat_to_Tps$total_exon_length, length(grep("Tps_F_FB_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_FB_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_FB_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_FB_Ad_", colnames(dat_to_Tps))))
y_TpsTdi_Fe_Ad_to_Tps_F_FPKM.  <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_Fe_Ad_to_Tps_UF, FPKM_filt_value,   2, dat_to_Tps$total_exon_length, length(grep("Tps_F_Fe_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_Fe_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_Fe_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_Fe_Ad_", colnames(dat_to_Tps))))
y_TpsTdi_Gu_Ad_to_Tps_F_FPKM.  <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_Gu_Ad_to_Tps_UF, FPKM_filt_value,   2, dat_to_Tps$total_exon_length, length(grep("Tps_F_Gu_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_Gu_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_Gu_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_Gu_Ad_", colnames(dat_to_Tps))))
y_TpsTdi_Ta_Ad_to_Tps_F_FPKM   <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_Ta_Ad_to_Tps_UF, FPKM_filt_value,   2, dat_to_Tps$total_exon_length, length(grep("Tps_F_Ta_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_Ta_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_Ta_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_Ta_Ad_", colnames(dat_to_Tps))))
y_TpsTdi_GoTe_Ad_to_Tps_F_FPKM <- filt_and_norm_male_female_FPKM_SA(y_TpsTdi_GoTe_Ad_to_Tps_UF, FPKM_filt_value, 2, dat_to_Tps$total_exon_length, length(grep("Tps_F_Go_Ad_", colnames(dat_to_Tps))), length(grep("Tps_M_Te_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_F_Go_Ad_", colnames(dat_to_Tps))), length(grep("Tdi_M_Te_Ad_", colnames(dat_to_Tps)))) ### Te,Go together


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

FPKM_filt_value = 2

y_TpaTge_A_Ad_to_Tpa_F_FPKM  <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_A_Ad_to_Tpa_UF, FPKM_filt_value, 2, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_A_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_A_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_A_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_A_Ad_", colnames(dat_to_Tpa))))
y_TpaTge_B_Ad_to_Tpa_F_FPKM  <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_B_Ad_to_Tpa_UF, FPKM_filt_value, 2, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_B_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_B_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_B_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_B_Ad_", colnames(dat_to_Tpa))))
y_TpaTge_DG_Ad_to_Tpa_F_FPKM <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_DG_Ad_to_Tpa_UF, FPKM_filt_value, 2, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_DG_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_DG_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_DG_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_DG_Ad_", colnames(dat_to_Tpa))))
y_TpaTge_FB_Ad_to_Tpa_F_FPKM <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_FB_Ad_to_Tpa_UF, FPKM_filt_value, 2, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_FB_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_FB_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_FB_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_FB_Ad_", colnames(dat_to_Tpa))))
y_TpaTge_Fe_Ad_to_Tpa_F_FPKM <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_Fe_Ad_to_Tpa_UF, FPKM_filt_value, 2, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_Fe_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_Fe_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_Fe_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_Fe_Ad_", colnames(dat_to_Tpa))))
y_TpaTge_Gu_Ad_to_Tpa_F_FPKM <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_Gu_Ad_to_Tpa_UF, FPKM_filt_value, 2, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_Gu_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_Gu_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_Gu_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_Gu_Ad_", colnames(dat_to_Tpa))))
y_TpaTge_Ta_Ad_to_Tpa_F_FPKM <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_Ta_Ad_to_Tpa_UF, FPKM_filt_value, 2, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_Ta_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_Ta_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_Ta_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_Ta_Ad_", colnames(dat_to_Tpa))))
y_TpaTge_GoTe_Ad_to_Tpa_F_FPKM <- filt_and_norm_male_female_FPKM_SA(y_TpaTge_GoTe_Ad_to_Tpa_UF, FPKM_filt_value, 2, dat_to_Tpa$total_exon_length, length(grep("Tpa_F_Go_Ad_", colnames(dat_to_Tpa))), length(grep("Tpa_M_Te_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_F_Go_Ad_", colnames(dat_to_Tpa))), length(grep("Tge_M_Te_Ad_", colnames(dat_to_Tpa)))) ### Te Go together


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

FPKM_filt_value = 2

y_TceTms_A_Ad_to_Tce_F_FPKM      <- filt_and_norm_male_female_FPKM_SA(y_TceTms_A_Ad_to_Tce_UF, FPKM_filt_value, 2, dat_to_Tce$total_exon_length, length(grep("Tce_F_A_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_A_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_A_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_A_Ad_", colnames(dat_to_Tce))))
y_TceTms_B_Ad_to_Tce_F_FPKM      <- filt_and_norm_male_female_FPKM_SA(y_TceTms_B_Ad_to_Tce_UF, FPKM_filt_value, 2, dat_to_Tce$total_exon_length, length(grep("Tce_F_B_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_B_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_B_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_B_Ad_", colnames(dat_to_Tce))))
y_TceTms_DG_Ad_to_Tce_F_FPKM     <- filt_and_norm_male_female_FPKM_SA(y_TceTms_DG_Ad_to_Tce_UF, FPKM_filt_value, 2, dat_to_Tce$total_exon_length, length(grep("Tce_F_DG_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_DG_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_DG_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_DG_Ad_", colnames(dat_to_Tce))))
y_TceTms_FB_Ad_to_Tce_F_FPKM     <- filt_and_norm_male_female_FPKM_SA(y_TceTms_FB_Ad_to_Tce_UF, FPKM_filt_value, 2, dat_to_Tce$total_exon_length, length(grep("Tce_F_FB_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_FB_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_FB_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_FB_Ad_", colnames(dat_to_Tce))))
y_TceTms_Fe_Ad_to_Tce_F_FPKM     <- filt_and_norm_male_female_FPKM_SA(y_TceTms_Fe_Ad_to_Tce_UF, FPKM_filt_value, 2, dat_to_Tce$total_exon_length, length(grep("Tce_F_Fe_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_Fe_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_Fe_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_Fe_Ad_", colnames(dat_to_Tce))))
y_TceTms_Gu_Ad_to_Tce_F_FPKM     <- filt_and_norm_male_female_FPKM_SA(y_TceTms_Gu_Ad_to_Tce_UF, FPKM_filt_value, 2, dat_to_Tce$total_exon_length, length(grep("Tce_F_Gu_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_Gu_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_Gu_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_Gu_Ad_", colnames(dat_to_Tce))))
y_TceTms_Ta_Ad_to_Tce_F_FPKM     <- filt_and_norm_male_female_FPKM_SA(y_TceTms_Ta_Ad_to_Tce_UF, FPKM_filt_value, 2, dat_to_Tce$total_exon_length, length(grep("Tce_F_Ta_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_Ta_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_Ta_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_Ta_Ad_", colnames(dat_to_Tce))))
y_TceTms_GoTe_Ad_to_Tce_F_FPKM   <- filt_and_norm_male_female_FPKM_SA(y_TceTms_GoTe_Ad_to_Tce_UF, FPKM_filt_value, 2, dat_to_Tce$total_exon_length, length(grep("Tce_F_Go_Ad_", colnames(dat_to_Tce))), length(grep("Tce_M_Te_Ad_", colnames(dat_to_Tce))), length(grep("Tms_F_Go_Ad_", colnames(dat_to_Tce))), length(grep("Tms_M_Te_Ad_", colnames(dat_to_Tce)))) ### Te and Go together





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



##################################################################################################################################################################
##################################################################################################################################################################
plot_MF_box <- function(df_sex, df_asex){
  sex_df_name  <- deparse(substitute(df_sex))
  asex_df_name <- deparse(substitute(df_asex))
  sex_sp  <- str_split(sex_df_name, "_")[[1]][1]
  asex_sp <- str_split(asex_df_name, "_")[[1]][1]  
  print(sex_df_name)
  print(asex_df_name) 
  print(sex_sp)
  print(asex_sp)
  
  df <- as.data.frame(cbind(
    c(df_sex$log2MF, df_asex$log2MF),
    c(df_sex$XA,     df_asex$XA),
    c(rep(sex_sp, length(df_sex$log2MF)),  rep(asex_sp, length(df_asex$log2MF)))))
  colnames(df) <- c("log2MF", "XA", "sp")
  df$rep_m <- ifelse(df$sp == "Tps", "sex", 
                     ifelse(df$sp == "Tdi", "asex", 
                            ifelse(df$sp == "Tpa", "sex", 
                                   ifelse(df$sp == "Tge", "asex", 
                                          ifelse(df$sp == "Tce", "sex", 
                                                 ifelse(df$sp == "Tms", "asex", 
                                                        ifelse(df$sp == "Tcm", "sex", 
                                                               ifelse(df$sp == "Tsi", "asex",                     
                                                                      ifelse(df$sp == "Tbi", "sex", 
                                                                             ifelse(df$sp == "Tte", "asex",  "ERROR")))))))))) 
  
  df$group <- paste(df$rep_m, df$XA, sep = "_")
  df$group_o <- ordered(df$group, c("sex_A", "sex_X", "asex_A", "asex_X"))  
  df$log2MF <- as.numeric(as.character(df$log2MF))
  
  print(head(df)) 
  
  P_out <- ggplot(df, aes(group_o, log2MF)) + 
    theme_classic() +
    geom_boxplot(aes(fill = factor(group_o)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
    coord_cartesian(ylim=c(-4,4)) +
    ylab ("log2(Male FPKM / Female FPKM)")  +
    xlab ("Group") + 
    scale_fill_manual(values=c("white", "darkorange", "white", "darkorange"))  + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
    ggtitle(paste(strsplit(sex_df_name, "_")[[1]][1] , strsplit(sex_df_name, "_")[[1]][2]))
  P_out <- P_out + theme(legend.position = "none")
  
  return(P_out)
  
}

#### 
plot_MF_box_df <- function(df_sex, df_asex){
  sex_df_name  <- deparse(substitute(df_sex))
  asex_df_name <- deparse(substitute(df_asex))
  sex_sp  <- str_split(sex_df_name, "_")[[1]][1]
  asex_sp <- str_split(asex_df_name, "_")[[1]][1] 
  sex_tiss  <- str_split(sex_df_name, "_")[[1]][2] 
  asex_tiss <- str_split(asex_df_name, "_")[[1]][2] 
  print(sex_df_name)
  print(asex_df_name) 
  print(sex_sp)
  print(asex_sp)
  print(sex_tiss)
  print(asex_tiss)
  
  df <- as.data.frame(cbind(
    c(df_sex$log2MF, df_asex$log2MF),
    c(df_sex$XA,     df_asex$XA),
    c(rep(sex_sp, length(df_sex$log2MF)),  rep(asex_sp, length(df_asex$log2MF))),
    c(rep(sex_tiss, length(df_sex$log2MF)),  rep(asex_tiss, length(df_asex$log2MF))),
    c(rep(paste(sex_sp, "_", asex_sp, sep = ""), length(df_sex$log2MF)),  rep(paste(sex_sp, "_", asex_sp, sep = ""), length(df_asex$log2MF)))
  ))
  
  colnames(df) <- c("log2MF", "XA", "sp", "tiss", "sp_pair")
  df$rep_m <- ifelse(df$sp == "Tps", "sex", 
                     ifelse(df$sp == "Tdi", "asex", 
                            ifelse(df$sp == "Tpa", "sex", 
                                   ifelse(df$sp == "Tge", "asex", 
                                          ifelse(df$sp == "Tce", "sex", 
                                                 ifelse(df$sp == "Tms", "asex", 
                                                        ifelse(df$sp == "Tcm", "sex", 
                                                               ifelse(df$sp == "Tsi", "asex",                     
                                                                      ifelse(df$sp == "Tbi", "sex", 
                                                                             ifelse(df$sp == "Tte", "asex",  "ERROR")))))))))) 
  
  df$group <- paste(df$rep_m, df$XA, sep = "_")
  df$group_o <- ordered(df$group, c("sex_A", "sex_X", "asex_A", "asex_X"))  
  df$log2MF <- as.numeric(as.character(df$log2MF))
  
  print(head(df)) 
  return(df)
  
}

#### 
MF_box_wilcox <- function(df_sex, df_asex){
  sex_df_name  <- deparse(substitute(df_sex))
  asex_df_name <- deparse(substitute(df_asex))
  sex_sp  <- str_split(sex_df_name, "_")[[1]][1]
  asex_sp <- str_split(asex_df_name, "_")[[1]][1] 
  sex_tiss  <- str_split(sex_df_name, "_")[[1]][2] 
  asex_tiss <- str_split(asex_df_name, "_")[[1]][2] 
  print(sex_df_name)
  print(asex_df_name) 
  print(sex_sp)
  print(asex_sp)
  print(sex_tiss)
  print(asex_tiss)
  
  df <- as.data.frame(cbind(
    c(df_sex$log2MF, df_asex$log2MF),
    c(df_sex$XA,     df_asex$XA),
    c(rep(sex_sp, length(df_sex$log2MF)),  rep(asex_sp, length(df_asex$log2MF))),
    c(rep(sex_tiss, length(df_sex$log2MF)),  rep(asex_tiss, length(df_asex$log2MF))),
    c(rep(paste(sex_sp, "_", asex_sp, sep = ""), length(df_sex$log2MF)),  rep(paste(sex_sp, "_", asex_sp, sep = ""), length(df_asex$log2MF)))    
  ))
  
  colnames(df) <- c("log2MF", "XA", "sp", "tiss", "sp_pair")
  df$rep_m <- ifelse(df$sp == "Tps", "sex", 
                     ifelse(df$sp == "Tdi", "asex", 
                            ifelse(df$sp == "Tpa", "sex", 
                                   ifelse(df$sp == "Tge", "asex", 
                                          ifelse(df$sp == "Tce", "sex", 
                                                 ifelse(df$sp == "Tms", "asex", 
                                                        ifelse(df$sp == "Tcm", "sex", 
                                                               ifelse(df$sp == "Tsi", "asex",                     
                                                                      ifelse(df$sp == "Tbi", "sex", 
                                                                             ifelse(df$sp == "Tte", "asex",  "ERROR")))))))))) 
  
  df$group <- paste(df$rep_m, df$XA, sep = "_")
  df$group_o <- ordered(df$group, c("sex_A", "sex_X", "asex_A", "asex_X"))  
  df$log2MF <- as.numeric(as.character(df$log2MF))
  
  print(head(df)) 
  
  
  ###sp tiss comp sex_med asex_med pval
  
  Wilcox_sex_X_sex_A   <- wilcox.test(subset(df, df$group == "sex_X")$log2MF,  subset(df, df$group == "sex_A")$log2MF)$p.value
  Wilcox_asex_X_asex_A <- wilcox.test(subset(df, df$group == "asex_X")$log2MF,  subset(df, df$group == "asex_A")$log2MF)$p.value
  Wilcox_sex_A_asex_A  <- wilcox.test(subset(df, df$group == "sex_A")$log2MF,  subset(df, df$group == "asex_A")$log2MF)$p.value  
  Wilcox_sex_X_asex_X  <- wilcox.test(subset(df, df$group == "sex_X")$log2MF,  subset(df, df$group == "asex_X")$log2MF)$p.value
  
  
  sex_X_medMF  <- median(subset(df, df$group == "sex_X")$log2MF)
  asex_X_medMF <- median(subset(df, df$group == "asex_X")$log2MF)  
  sex_A_medMF  <- median(subset(df, df$group == "sex_A")$log2MF)
  asex_A_medMF <- median(subset(df, df$group == "asex_A")$log2MF)    
  
  
  wilcox_df <- as.data.frame(rbind(
    c(sex_sp, sex_tiss, "sex_X_sex_A",   sex_X_medMF, sex_A_medMF, Wilcox_sex_X_sex_A ),
    c(sex_sp, sex_tiss, "asex_X_asex_A", asex_X_medMF, asex_A_medMF, Wilcox_asex_X_asex_A),
    c(sex_sp, sex_tiss, "sex_A_asex_A",  sex_A_medMF, asex_A_medMF, Wilcox_sex_A_asex_A),
    c(sex_sp, sex_tiss, "sex_X_asex_X",  sex_X_medMF, asex_X_medMF, Wilcox_sex_X_asex_X))
  )
  
  colnames(wilcox_df) <- c("sp", "tiss", "comp", "med_comp1", "med_comp2", "wilcox_p")
  
  return(wilcox_df)
  
}

###############################################


png(filename = "MF_box_to_Tps_F_FPKM.png", width = 8, height = 12, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_box(Tps_A_Ad_to_Tps_F_FPKM , Tdi_A_Ad_to_Tps_F_FPKM  ),
  plot_MF_box(Tps_B_Ad_to_Tps_F_FPKM , Tdi_B_Ad_to_Tps_F_FPKM  ),
  plot_MF_box(Tps_DG_Ad_to_Tps_F_FPKM , Tdi_DG_Ad_to_Tps_F_FPKM  ),
  plot_MF_box(Tps_FB_Ad_to_Tps_F_FPKM , Tdi_FB_Ad_to_Tps_F_FPKM  ),
  plot_MF_box(Tps_Fe_Ad_to_Tps_F_FPKM , Tdi_Fe_Ad_to_Tps_F_FPKM),
  plot_MF_box(Tps_Gu_Ad_to_Tps_F_FPKM , Tdi_Gu_Ad_to_Tps_F_FPKM),
  plot_MF_box(Tps_Ta_Ad_to_Tps_F_FPKM , Tdi_Ta_Ad_to_Tps_F_FPKM),
  plot_MF_box(Tps_GoTe_Ad_to_Tps_F_FPKM , Tdi_GoTe_Ad_to_Tps_F_FPKM),
  ncol = 3)
dev.off()
getwd() ## where has my plot gone...



png(filename = "MF_box_to_Tpa_F_FPKM.png", width = 8, height = 12, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_box(Tpa_A_Ad_to_Tpa_F_FPKM , Tge_A_Ad_to_Tpa_F_FPKM  ),
  plot_MF_box(Tpa_B_Ad_to_Tpa_F_FPKM , Tge_B_Ad_to_Tpa_F_FPKM  ),
  plot_MF_box(Tpa_DG_Ad_to_Tpa_F_FPKM , Tge_DG_Ad_to_Tpa_F_FPKM  ),
  plot_MF_box(Tpa_FB_Ad_to_Tpa_F_FPKM , Tge_FB_Ad_to_Tpa_F_FPKM  ),
  plot_MF_box(Tpa_Fe_Ad_to_Tpa_F_FPKM , Tge_Fe_Ad_to_Tpa_F_FPKM  ),
  plot_MF_box(Tpa_Gu_Ad_to_Tpa_F_FPKM , Tge_Gu_Ad_to_Tpa_F_FPKM  ),
  plot_MF_box(Tpa_Ta_Ad_to_Tpa_F_FPKM , Tge_Ta_Ad_to_Tpa_F_FPKM  ), 
  plot_MF_box(Tpa_GoTe_Ad_to_Tpa_F_FPKM , Tge_GoTe_Ad_to_Tpa_F_FPKM), 
  ncol = 3)
dev.off()
getwd() ## where has my plot gone...


png(filename = "MF_box_to_Tce_F_FPKM.png", width = 8, height = 12, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_box(Tce_A_Ad_to_Tce_F_FPKM , Tms_A_Ad_to_Tce_F_FPKM  ),
  plot_MF_box(Tce_B_Ad_to_Tce_F_FPKM , Tms_B_Ad_to_Tce_F_FPKM  ),
  plot_MF_box(Tce_DG_Ad_to_Tce_F_FPKM , Tms_DG_Ad_to_Tce_F_FPKM  ),
  plot_MF_box(Tce_FB_Ad_to_Tce_F_FPKM , Tms_FB_Ad_to_Tce_F_FPKM  ),
  plot_MF_box(Tce_Fe_Ad_to_Tce_F_FPKM , Tms_Fe_Ad_to_Tce_F_FPKM  ),
  plot_MF_box(Tce_Gu_Ad_to_Tce_F_FPKM , Tms_Gu_Ad_to_Tce_F_FPKM  ),
  plot_MF_box(Tce_Ta_Ad_to_Tce_F_FPKM , Tms_Ta_Ad_to_Tce_F_FPKM  ), 
  plot_MF_box(Tce_GoTe_Ad_to_Tce_F_FPKM , Tms_GoTe_Ad_to_Tce_F_FPKM),  
  ncol = 3)
dev.off()
getwd() ## where has my plot gone...




#############################################################################################################################################
######## plot all non-rep tissues together


### Tce

Tce_somatic_MF_df <- rbind(
  plot_MF_box_df(Tce_A_Ad_to_Tce_F_FPKM , Tms_A_Ad_to_Tce_F_FPKM  ),
  plot_MF_box_df(Tce_B_Ad_to_Tce_F_FPKM , Tms_B_Ad_to_Tce_F_FPKM  ),
  plot_MF_box_df(Tce_DG_Ad_to_Tce_F_FPKM , Tms_DG_Ad_to_Tce_F_FPKM  ),
  plot_MF_box_df(Tce_FB_Ad_to_Tce_F_FPKM , Tms_FB_Ad_to_Tce_F_FPKM  ),
  plot_MF_box_df(Tce_Fe_Ad_to_Tce_F_FPKM , Tms_Fe_Ad_to_Tce_F_FPKM),
  plot_MF_box_df(Tce_Gu_Ad_to_Tce_F_FPKM , Tms_Gu_Ad_to_Tce_F_FPKM),
  plot_MF_box_df(Tce_Ta_Ad_to_Tce_F_FPKM , Tms_Ta_Ad_to_Tce_F_FPKM))

head(Tce_somatic_MF_df )


Tce_somatic_MF_plot <- ggplot(Tce_somatic_MF_df, aes(tiss, log2MF)) + 
  theme_classic() +
  geom_boxplot(aes(fill = factor(group_o)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
  coord_cartesian(ylim=c(-3,3)) +
  ylab ("log2(Male FPKM / Female FPKM)")  +
  xlab ("Tissue") + 
  scale_fill_manual(values=c("white", "darkorange","grey", "#2297E6"))   + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
  ggtitle("Tce") + theme(legend.position = "none")

Tce_somatic_MF_plot_legend <- ggplot(Tce_somatic_MF_df, aes(tiss, log2MF)) + 
  theme_classic() +
  geom_boxplot(aes(fill = factor(group_o)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
  coord_cartesian(ylim=c(-3,3)) +
  ylab ("log2(Male FPKM / Female FPKM)")  +
  xlab ("Tissue") + 
  scale_fill_manual(values=c("white", "darkorange","grey", "#2297E6"))   + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
  ggtitle("Tce")



## Tpa

Tpa_somatic_MF_df <- rbind(
  plot_MF_box_df(Tpa_A_Ad_to_Tpa_F_FPKM , Tge_A_Ad_to_Tpa_F_FPKM  ),
  plot_MF_box_df(Tpa_B_Ad_to_Tpa_F_FPKM , Tge_B_Ad_to_Tpa_F_FPKM  ),
  plot_MF_box_df(Tpa_DG_Ad_to_Tpa_F_FPKM , Tge_DG_Ad_to_Tpa_F_FPKM  ),
  plot_MF_box_df(Tpa_FB_Ad_to_Tpa_F_FPKM , Tge_FB_Ad_to_Tpa_F_FPKM  ),
  plot_MF_box_df(Tpa_Fe_Ad_to_Tpa_F_FPKM , Tge_Fe_Ad_to_Tpa_F_FPKM),
  plot_MF_box_df(Tpa_Gu_Ad_to_Tpa_F_FPKM , Tge_Gu_Ad_to_Tpa_F_FPKM),
  plot_MF_box_df(Tpa_Ta_Ad_to_Tpa_F_FPKM , Tge_Ta_Ad_to_Tpa_F_FPKM))

head(Tpa_somatic_MF_df )


Tpa_somatic_MF_plot <- ggplot(Tpa_somatic_MF_df, aes(tiss, log2MF)) + 
  theme_classic() +
  geom_boxplot(aes(fill = factor(group_o)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
  coord_cartesian(ylim=c(-3,3)) +
  ylab ("log2(Male FPKM / Female FPKM)")  +
  xlab ("Tissue") + 
  scale_fill_manual(values=c("white", "darkorange","grey", "#2297E6"))  + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
  ggtitle("Tpa") + theme(legend.position = "none")

Tpa_somatic_MF_plot_legend <- ggplot(Tpa_somatic_MF_df, aes(tiss, log2MF)) + 
  theme_classic() +
  geom_boxplot(aes(fill = factor(group_o)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
  coord_cartesian(ylim=c(-3,3)) +
  ylab ("log2(Male FPKM / Female FPKM)")  +
  xlab ("Tissue") + 
  scale_fill_manual(values=c("white", "darkorange","grey", "#2297E6"))   + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
  ggtitle("Tpa") 

## Tps 

Tps_somatic_MF_df <- rbind(
  plot_MF_box_df(Tps_A_Ad_to_Tps_F_FPKM , Tdi_A_Ad_to_Tps_F_FPKM  ),
  plot_MF_box_df(Tps_B_Ad_to_Tps_F_FPKM , Tdi_B_Ad_to_Tps_F_FPKM  ),
  plot_MF_box_df(Tps_DG_Ad_to_Tps_F_FPKM , Tdi_DG_Ad_to_Tps_F_FPKM  ),
  plot_MF_box_df(Tps_FB_Ad_to_Tps_F_FPKM , Tdi_FB_Ad_to_Tps_F_FPKM  ),
  plot_MF_box_df(Tps_Fe_Ad_to_Tps_F_FPKM , Tdi_Fe_Ad_to_Tps_F_FPKM),
  plot_MF_box_df(Tps_Gu_Ad_to_Tps_F_FPKM , Tdi_Gu_Ad_to_Tps_F_FPKM),
  plot_MF_box_df(Tps_Ta_Ad_to_Tps_F_FPKM , Tdi_Ta_Ad_to_Tps_F_FPKM))

head(Tps_somatic_MF_df )


Tps_somatic_MF_plot <- ggplot(Tps_somatic_MF_df, aes(tiss, log2MF)) + 
  theme_classic() +
  geom_boxplot(aes(fill = factor(group_o)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
  coord_cartesian(ylim=c(-3,3)) +
  ylab ("log2(Male FPKM / Female FPKM)")  +
  xlab ("Tissue") + 
  scale_fill_manual(values=c("white", "darkorange","grey", "#2297E6"))   + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
  ggtitle("Tps") + theme(legend.position = "none")


Tps_somatic_MF_plot_legend <- ggplot(Tps_somatic_MF_df, aes(tiss, log2MF)) + 
  theme_classic() +
  geom_boxplot(aes(fill = factor(group_o)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
  coord_cartesian(ylim=c(-3,3)) +
  ylab ("log2(Male FPKM / Female FPKM)")  +
  xlab ("Tissue") + 
  scale_fill_manual(values=c("white", "darkorange","grey", "#2297E6"))   + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
  ggtitle("Tps") 


pdf("TceTpaTps_MF_non-rep.pdf", width = 7, height = 10)
plot_grid(
  Tce_somatic_MF_plot,
  Tpa_somatic_MF_plot,
  Tps_somatic_MF_plot, ncol = 1
)
dev.off()
getwd() ## where has my plot gone....?

pdf("TceTpaTps_MF_non-rep_legend.pdf", width = 8, height = 10)
plot_grid(
  Tce_somatic_MF_plot_legend,
  Tpa_somatic_MF_plot_legend,
  Tps_somatic_MF_plot_legend, ncol = 1
)
dev.off()
getwd() ## where has my plot gone....?




#### wilcox (all rep and non-rep)

Tce_MF_wilcox_df <- rbind(
  MF_box_wilcox(Tce_A_Ad_to_Tce_F_FPKM , Tms_A_Ad_to_Tce_F_FPKM  ),
  MF_box_wilcox(Tce_B_Ad_to_Tce_F_FPKM , Tms_B_Ad_to_Tce_F_FPKM  ),
  MF_box_wilcox(Tce_DG_Ad_to_Tce_F_FPKM , Tms_DG_Ad_to_Tce_F_FPKM  ),
  MF_box_wilcox(Tce_FB_Ad_to_Tce_F_FPKM , Tms_FB_Ad_to_Tce_F_FPKM  ),
  MF_box_wilcox(Tce_Fe_Ad_to_Tce_F_FPKM , Tms_Fe_Ad_to_Tce_F_FPKM),
  MF_box_wilcox(Tce_Gu_Ad_to_Tce_F_FPKM , Tms_Gu_Ad_to_Tce_F_FPKM),
  MF_box_wilcox(Tce_Ta_Ad_to_Tce_F_FPKM , Tms_Ta_Ad_to_Tce_F_FPKM),
  MF_box_wilcox(Tce_GoTe_Ad_to_Tce_F_FPKM , Tms_GoTe_Ad_to_Tce_F_FPKM))


Tpa_MF_wilcox_df <- rbind(
  MF_box_wilcox(Tpa_A_Ad_to_Tpa_F_FPKM , Tge_A_Ad_to_Tpa_F_FPKM  ),
  MF_box_wilcox(Tpa_B_Ad_to_Tpa_F_FPKM , Tge_B_Ad_to_Tpa_F_FPKM  ),
  MF_box_wilcox(Tpa_DG_Ad_to_Tpa_F_FPKM , Tge_DG_Ad_to_Tpa_F_FPKM  ),
  MF_box_wilcox(Tpa_FB_Ad_to_Tpa_F_FPKM , Tge_FB_Ad_to_Tpa_F_FPKM  ),
  MF_box_wilcox(Tpa_Fe_Ad_to_Tpa_F_FPKM , Tge_Fe_Ad_to_Tpa_F_FPKM),
  MF_box_wilcox(Tpa_Gu_Ad_to_Tpa_F_FPKM , Tge_Gu_Ad_to_Tpa_F_FPKM),
  MF_box_wilcox(Tpa_Ta_Ad_to_Tpa_F_FPKM , Tge_Ta_Ad_to_Tpa_F_FPKM),
  MF_box_wilcox(Tpa_GoTe_Ad_to_Tpa_F_FPKM , Tge_GoTe_Ad_to_Tpa_F_FPKM))


Tps_MF_wilcox_df <- rbind(
  MF_box_wilcox(Tps_A_Ad_to_Tps_F_FPKM , Tdi_A_Ad_to_Tps_F_FPKM  ),
  MF_box_wilcox(Tps_B_Ad_to_Tps_F_FPKM , Tdi_B_Ad_to_Tps_F_FPKM  ),
  MF_box_wilcox(Tps_DG_Ad_to_Tps_F_FPKM , Tdi_DG_Ad_to_Tps_F_FPKM  ),
  MF_box_wilcox(Tps_FB_Ad_to_Tps_F_FPKM , Tdi_FB_Ad_to_Tps_F_FPKM  ),
  MF_box_wilcox(Tps_Fe_Ad_to_Tps_F_FPKM , Tdi_Fe_Ad_to_Tps_F_FPKM),
  MF_box_wilcox(Tps_Gu_Ad_to_Tps_F_FPKM , Tdi_Gu_Ad_to_Tps_F_FPKM),
  MF_box_wilcox(Tps_Ta_Ad_to_Tps_F_FPKM , Tdi_Ta_Ad_to_Tps_F_FPKM),
  MF_box_wilcox(Tps_GoTe_Ad_to_Tps_F_FPKM , Tdi_GoTe_Ad_to_Tps_F_FPKM))

Tce_MF_wilcox_df$wilcox_FDR <- p.adjust(Tce_MF_wilcox_df$wilcox_p, method="fdr")
Tpa_MF_wilcox_df$wilcox_FDR <- p.adjust(Tpa_MF_wilcox_df$wilcox_p, method="fdr")
Tps_MF_wilcox_df$wilcox_FDR <- p.adjust(Tps_MF_wilcox_df$wilcox_p, method="fdr")

write.csv(Tce_MF_wilcox_df, file = "Tce_MF_wilcox_df.csv", row.names = FALSE, quote = F )
write.csv(Tpa_MF_wilcox_df, file = "Tpa_MF_wilcox_df.csv", row.names = FALSE, quote = F )
write.csv(Tps_MF_wilcox_df, file = "Tps_MF_wilcox_df.csv", row.names = FALSE, quote = F )


### plot reproductive

Te_MF_df <- rbind(
  plot_MF_box_df(Tce_GoTe_Ad_to_Tce_F_FPKM , Tms_GoTe_Ad_to_Tce_F_FPKM),
  plot_MF_box_df(Tpa_GoTe_Ad_to_Tpa_F_FPKM , Tge_GoTe_Ad_to_Tpa_F_FPKM),
  plot_MF_box_df(Tps_GoTe_Ad_to_Tps_F_FPKM , Tdi_GoTe_Ad_to_Tps_F_FPKM)
)


comb_Te_MF_plot <- ggplot(Te_MF_df, aes(sp_pair, log2MF)) + 
  theme_classic() +
  geom_boxplot(aes(fill = factor(group_o)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
  coord_cartesian(ylim=c(-5,5)) +
  ylab ("log2(Male FPKM / Female FPKM)")  +
  xlab ("Species pair") + 
  scale_fill_manual(values=c("white", "darkorange","grey", "#2297E6"))   + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
  ggtitle("Te") + theme(legend.position = "none")






###################################################################################################################################################
###### plot average FPKM

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


extract_avMF = function(dat){
  
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
    as.character(eval(parse(text=paste('dat','$', 'log2MF', sep='')))),
    rep(sp, length(dat[,1])),
    rep(tiss, length(dat[,1]))
  ))
  
  colnames(df_out) <- c("gene_id", "scaf", "XA", "mid_point", "log2MF", "sp", "tiss")
  
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
  df_out$log2MF <- as.numeric(df_out$log2MF)
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




plot_avFPKM <- function(df, max_y){
  
  P1 <- ggplot(df, aes(tiss, meanFPKM)) + 
    theme_classic() +
    geom_boxplot(aes(fill = factor(sex_XA)), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T, position = position_dodge2(preserve = "single")) +
    coord_cartesian(ylim=c(0,max_y)) +
    ylab("FPKM") +
    xlab("") + 
    scale_fill_manual(values=c("F_A" = "grey", "F_X" = "firebrick2", "M_A" = "white", "M_X" = "royalblue2"), labels=c("F_A" = "Female A","F_X" =  "Female X", "M_A" = "Male A", "M_X" = "Male X")) +
    ggtitle(paste(deparse(substitute(df))))
  
  outlist = list("P1" = P1) 
  return(outlist)	
}

### with Te and Go together
plot_avFPKM2 <- function(df, max_y){
  
  P1 <- ggplot(df, aes(tiss2, meanFPKM)) + 
    theme_classic() +
    geom_boxplot(aes(fill = factor(sex_XA)), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T, position = position_dodge2(preserve = "single")) +
    coord_cartesian(ylim=c(0,max_y)) +
    ylab("FPKM") +
    xlab("") + 
    scale_fill_manual(values=c("F_A" = "grey", "F_X" = "firebrick2", "M_A" = "white", "M_X" = "royalblue2"), labels=c("F_A" = "Female A","F_X" =  "Female X", "M_A" = "Male A", "M_X" = "Male X")) +
    ggtitle(paste(deparse(substitute(df))))
  
  outlist = list("P1" = P1) 
  return(outlist)	
}




plot_avFPKM_scafs <- function(df, max_y, want_scafs, want_tiss){
  
  df_sub_1 <- subset(df, scaf %in% want_scafs)
  df_sub   <- subset(df_sub_1, tiss %in% want_tiss)
  
  print(length(df[,1]))
  print(length(df_sub_1[,1]))  
  print(length(df_sub[,1]))  
  
  P1 <- ggplot(df_sub, aes(scaf, meanFPKM)) + 
    theme_classic() +
    geom_boxplot(aes(fill = factor(sex_XA)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
    coord_cartesian(ylim=c(0,max_y)) +
    ylab("FPKM") +
    xlab("") + 
    scale_fill_manual(values=c("F_A" = "grey", "F_X" = "firebrick2", "M_A" = "white", "M_X" = "royalblue2"), labels=c("F_A" = "Female A","F_X" =  "Female X", "M_A" = "Male A", "M_X" = "Male X")) +
    ggtitle(paste(deparse(substitute(df)), want_tiss, sep = " - ")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  x_labs <- str_replace(str_replace(ggplot_build(P1)$layout$panel_params[[1]]$x$get_labels(), "Tps_LRv5b_", ""), "Tpa_LRv5a_", "")
  print(x_labs)
  
  P2 <- P1 + theme(legend.position = "none") + ggtitle(paste(strsplit(deparse(substitute(df)), "_")[[1]][1], want_tiss, sep = " - ")) +
    scale_x_discrete(labels=x_labs)
  
  outlist = list("P1" = P1, "P2" = P2) 
  return(outlist)	
}


plot_avFPKM_LGs <- function(df, max_y, want_LGs, want_tiss){
  
  df_sub_1 <- subset(df, LG %in% want_LGs)
  df_sub   <- subset(df_sub_1, tiss %in% want_tiss)
  
  print(length(df[,1]))
  print(length(df_sub_1[,1]))  
  print(length(df_sub[,1]))  
  
  P1 <- ggplot(df_sub, aes(LG, meanFPKM)) + 
    theme_classic() +
    geom_boxplot(aes(fill = factor(sex_XA)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
    coord_cartesian(ylim=c(0,max_y)) +
    ylab("FPKM") +
    xlab("") + 
    scale_fill_manual(values=c("F_A" = "grey", "F_X" = "firebrick2", "M_A" = "white", "M_X" = "royalblue2"), labels=c("F_A" = "Female A","F_X" =  "Female X", "M_A" = "Male A", "M_X" = "Male X")) +
    ggtitle(paste(deparse(substitute(df)), want_tiss, sep = " - ")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  x_labs <- str_replace(str_replace(ggplot_build(P1)$layout$panel_params[[1]]$x$get_labels(), "Tps_LRv5b_", ""), "Tpa_LRv5a_", "")
  print(x_labs)
  
  P2 <- P1 + theme(legend.position = "none") + ggtitle(paste(strsplit(deparse(substitute(df)), "_")[[1]][1], want_tiss, sep = " - ")) +
    scale_x_discrete(labels=x_labs)
  
  outlist = list("P1" = P1, "P2" = P2) 
  return(outlist)	
}


plot_MF_scafs <- function(df, want_scafs, want_tiss){
  
  df_sub_1 <- subset(df, scaf %in% want_scafs)
  df_sub   <- subset(df_sub_1, tiss %in% want_tiss)
  
  print(length(df[,1]))
  print(length(df_sub_1[,1]))  
  print(length(df_sub[,1]))  
  
  P1 <- ggplot(df_sub, aes(scaf, log2MF)) + 
    theme_classic() +
    geom_boxplot(aes(fill = factor(XA)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
    coord_cartesian(ylim=c(-4,4)) +
    ylab ("log2(Male FPKM / Female FPKM)")  +
    xlab ("Group") + 
    scale_fill_manual(values=c("white", "darkorange","grey", "#2297E6"))   + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
    ggtitle(paste(deparse(substitute(df)), want_tiss, sep = " - ")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  P1 <- P1 + theme(legend.position = "none")
  
  x_labs <- str_replace(str_replace(ggplot_build(P1)$layout$panel_params[[1]]$x$get_labels(), "Tps_LRv5b_", ""), "Tpa_LRv5a_", "")
  print(x_labs)
  P2 <- P1 + theme(legend.position = "none") + ggtitle(paste(strsplit(deparse(substitute(df)), "_")[[1]][1], want_tiss, sep = " - ")) +
    scale_x_discrete(labels=x_labs)
  
  outlist = list("P1" = P1, "P2" = P2) 
  return(outlist)	
  
}

plot_MF_LGs <- function(df, want_LGs, want_tiss){
  
  df_sub_1 <- subset(df, LG %in% want_LGs)
  df_sub   <- subset(df_sub_1, tiss %in% want_tiss)
  
  print(length(df[,1]))
  print(length(df_sub_1[,1]))  
  print(length(df_sub[,1]))  
  
  P1 <- ggplot(df_sub, aes(LG, log2MF)) + 
    theme_classic() +
    geom_boxplot(aes(fill = factor(XA)),position=position_dodge(0.6), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T) +
    coord_cartesian(ylim=c(-4,4)) +
    ylab ("log2(Male FPKM / Female FPKM)")  +
    xlab ("Group") + 
    scale_fill_manual(values=c("white", "darkorange","grey", "#2297E6"))   + geom_hline(yintercept = 0) + geom_hline(yintercept = -1, linetype = 2) +
    ggtitle(paste(deparse(substitute(df)), want_tiss, sep = " - ")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  P1 <- P1 + theme(legend.position = "none")
  x_labs <- str_replace(str_replace(ggplot_build(P1)$layout$panel_params[[1]]$x$get_labels(), "Tps_LRv5b_", ""), "Tpa_LRv5a_", "")
  print(x_labs)
  
  P2 <- P1 + theme(legend.position = "none") + ggtitle(paste(strsplit(deparse(substitute(df)), "_")[[1]][1], want_tiss, sep = " - ")) +
    scale_x_discrete(labels=x_labs)
  
  outlist = list("P1" = P1, "P2" = P2) 
  return(outlist)	
  
}



### combined

Tps_Ad_to_Tps_F_avFPKM <- as.data.frame(rbind(
  
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
  extract_avEXP_GoTe(Tps_GoTe_Ad_to_Tps_F_FPKM, "M", "Te")
))

str(Tps_Ad_to_Tps_F_avFPKM )
Tps_Ad_to_Tps_F_avFPKM$tiss <- ordered(as.factor(Tps_Ad_to_Tps_F_avFPKM $tiss), levels=c("A", "B", "DG", "FB", "Fe", "Gu", "Ta", "Go", "Te"))
Tps_Ad_to_Tps_F_avFPKM$tiss2 <- ifelse(as.character(Tps_Ad_to_Tps_F_avFPKM$tiss) == "Te", "Go", as.character(Tps_Ad_to_Tps_F_avFPKM$tiss)) 
Tps_Ad_to_Tps_F_avFPKM$tiss2 <- ordered(as.factor(Tps_Ad_to_Tps_F_avFPKM$tiss2), levels=c("A", "B", "DG", "FB", "Fe", "Gu", "Ta", "Go"))


Tdi_Ad_to_Tps_F_avFPKM <- as.data.frame(rbind(
  
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


Tdi_Ad_to_Tps_F_avFPKM$tiss  <- ordered(as.factor(Tdi_Ad_to_Tps_F_avFPKM $tiss), levels=c("A", "B", "DG", "FB", "Fe", "Gu", "Ta", "Go", "Te"))
Tdi_Ad_to_Tps_F_avFPKM$tiss2 <- ifelse(as.character(Tdi_Ad_to_Tps_F_avFPKM$tiss) == "Te", "Go", as.character(Tdi_Ad_to_Tps_F_avFPKM$tiss)) 
Tdi_Ad_to_Tps_F_avFPKM$tiss2 <- ordered(as.factor(Tdi_Ad_to_Tps_F_avFPKM$tiss2), levels=c("A", "B", "DG", "FB", "Fe", "Gu", "Ta", "Go"))


png(filename = "TpsTdi_Ad_to_Tps_F_avFPKM.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM(Tps_Ad_to_Tps_F_avFPKM , 200)$P1,
  plot_avFPKM(Tdi_Ad_to_Tps_F_avFPKM , 200)$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

pdf("TpsTdi_Ad_to_Tps_F_avFPKM.pdf", width = 8, height = 10)
plot_grid(
  plot_avFPKM(Tps_Ad_to_Tps_F_avFPKM , 200)$P1,
  plot_avFPKM(Tdi_Ad_to_Tps_F_avFPKM , 200)$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...


png(filename = "TpsTdi_Ad_to_Tps_F_avFPKM_2.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM2(Tps_Ad_to_Tps_F_avFPKM , 200)$P1,
  plot_avFPKM2(Tdi_Ad_to_Tps_F_avFPKM , 200)$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

pdf("TpsTdi_Ad_to_Tps_F_avFPKM_2.pdf", width = 8, height = 10)
plot_grid(
  plot_avFPKM2(Tps_Ad_to_Tps_F_avFPKM , 200)$P1,
  plot_avFPKM2(Tdi_Ad_to_Tps_F_avFPKM , 200)$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...




##########################################################################################################################
## Tpa
Tpa_Ad_to_Tpa_F_avFPKM <- as.data.frame(rbind(
  
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
  extract_avEXP_GoTe(Tpa_GoTe_Ad_to_Tpa_F_FPKM, "M", "Te")
))

str(Tpa_Ad_to_Tpa_F_avFPKM )
Tpa_Ad_to_Tpa_F_avFPKM$tiss <- ordered(as.factor(Tpa_Ad_to_Tpa_F_avFPKM $tiss), levels=c("A", "B", "DG", "FB", "Fe", "Gu", "Ta", "Go", "Te"))
Tpa_Ad_to_Tpa_F_avFPKM$tiss2 <- ifelse(as.character(Tpa_Ad_to_Tpa_F_avFPKM$tiss) == "Te", "Go", as.character(Tpa_Ad_to_Tpa_F_avFPKM$tiss)) 
Tpa_Ad_to_Tpa_F_avFPKM$tiss2 <- ordered(as.factor(Tpa_Ad_to_Tpa_F_avFPKM$tiss2), levels=c("A", "B", "DG", "FB", "Fe", "Gu", "Ta", "Go"))


Tge_Ad_to_Tpa_F_avFPKM <- as.data.frame(rbind(
  
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


Tge_Ad_to_Tpa_F_avFPKM$tiss  <- ordered(as.factor(Tge_Ad_to_Tpa_F_avFPKM $tiss), levels=c("A", "B", "DG", "FB", "Fe", "Gu", "Ta", "Go", "Te"))
Tge_Ad_to_Tpa_F_avFPKM$tiss2 <- ifelse(as.character(Tge_Ad_to_Tpa_F_avFPKM$tiss) == "Te", "Go", as.character(Tge_Ad_to_Tpa_F_avFPKM$tiss)) 
Tge_Ad_to_Tpa_F_avFPKM$tiss2 <- ordered(as.factor(Tge_Ad_to_Tpa_F_avFPKM$tiss2), levels=c("A", "B", "DG", "FB", "Fe", "Gu", "Ta", "Go"))


png(filename = "TpaTge_Ad_to_Tpa_F_avFPKM.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM(Tpa_Ad_to_Tpa_F_avFPKM , 200)$P1,
  plot_avFPKM(Tge_Ad_to_Tpa_F_avFPKM , 200)$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

pdf("TpaTge_Ad_to_Tpa_F_avFPKM.pdf", width = 8, height = 10)
plot_grid(
  plot_avFPKM(Tpa_Ad_to_Tpa_F_avFPKM , 200)$P1,
  plot_avFPKM(Tge_Ad_to_Tpa_F_avFPKM , 200)$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...


png(filename = "TpaTge_Ad_to_Tpa_F_avFPKM_2.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM2(Tpa_Ad_to_Tpa_F_avFPKM , 200)$P1,
  plot_avFPKM2(Tge_Ad_to_Tpa_F_avFPKM , 200)$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

pdf("TpaTge_Ad_to_Tpa_F_avFPKM_2.pdf", width = 8, height = 10)
plot_grid(
  plot_avFPKM2(Tpa_Ad_to_Tpa_F_avFPKM , 200)$P1,
  plot_avFPKM2(Tge_Ad_to_Tpa_F_avFPKM , 200)$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...





########################################################
## Tce


## combined
Tce_Ad_to_Tce_F_avFPKM <- as.data.frame(rbind(
  
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
  extract_avEXP_GoTe(Tce_GoTe_Ad_to_Tce_F_FPKM, "M", "Te")
))

str(Tce_Ad_to_Tce_F_avFPKM )
Tce_Ad_to_Tce_F_avFPKM$tiss <- ordered(as.factor(Tce_Ad_to_Tce_F_avFPKM $tiss), levels=c("A", "B", "DG", "FB", "Fe", "Gu", "Ta", "Go", "Te"))
Tce_Ad_to_Tce_F_avFPKM$tiss2 <- ifelse(as.character(Tce_Ad_to_Tce_F_avFPKM$tiss) == "Te", "Go", as.character(Tce_Ad_to_Tce_F_avFPKM$tiss)) 
Tce_Ad_to_Tce_F_avFPKM$tiss2 <- ordered(as.factor(Tce_Ad_to_Tce_F_avFPKM$tiss2), levels=c("A", "B", "DG", "FB", "Fe", "Gu", "Ta", "Go"))


Tms_Ad_to_Tce_F_avFPKM <- as.data.frame(rbind(
  
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


Tms_Ad_to_Tce_F_avFPKM$tiss  <- ordered(as.factor(Tms_Ad_to_Tce_F_avFPKM $tiss), levels=c("A", "B", "DG", "FB", "Fe", "Gu", "Ta", "Go", "Te"))
Tms_Ad_to_Tce_F_avFPKM$tiss2 <- ifelse(as.character(Tms_Ad_to_Tce_F_avFPKM$tiss) == "Te", "Go", as.character(Tms_Ad_to_Tce_F_avFPKM$tiss)) 
Tms_Ad_to_Tce_F_avFPKM$tiss2 <- ordered(as.factor(Tms_Ad_to_Tce_F_avFPKM$tiss2), levels=c("A", "B", "DG", "FB", "Fe", "Gu", "Ta", "Go"))


png(filename = "TceTms_Ad_to_Tce_F_avFPKM.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM(Tce_Ad_to_Tce_F_avFPKM , 200)$P1,
  plot_avFPKM(Tms_Ad_to_Tce_F_avFPKM , 200)$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

pdf("TceTms_Ad_to_Tce_F_avFPKM.pdf", width = 8, height = 10)
plot_grid(
  plot_avFPKM(Tce_Ad_to_Tce_F_avFPKM , 200)$P1,
  plot_avFPKM(Tms_Ad_to_Tce_F_avFPKM , 200)$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...


png(filename = "TceTms_Ad_to_Tce_F_avFPKM_2.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM2(Tce_Ad_to_Tce_F_avFPKM , 200)$P1,
  plot_avFPKM2(Tms_Ad_to_Tce_F_avFPKM , 200)$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

pdf("TceTms_Ad_to_Tce_F_avFPKM_2.pdf", width = 8, height = 10)
plot_grid(
  plot_avFPKM2(Tce_Ad_to_Tce_F_avFPKM , 200)$P1,
  plot_avFPKM2(Tms_Ad_to_Tce_F_avFPKM , 200)$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...




##################################
## look at testes expression
## This is OK as the same gene set



plot_avFPKM_testes_sex_asex <- function(df_sex, df_asex, max_y){
  df <- rbind(df_sex, df_asex)
  
  df_Te <- subset(df, df$tiss == "Te")
  df_Go <- subset(df, df$tiss == "Go")
  
  
  df_GoTe <- rbind(df_Go, df_Te)
  df_GoTe$rep_sex_XA <- paste(df_GoTe$rep_m, df_GoTe$sex_XA, sep = "_")
  df_GoTe$rep_sex_XA <- ordered(df_GoTe$rep_sex_XA, levels = c("sex_F_A", "asex_F_A","sex_F_X", "asex_F_X","sex_M_A", "asex_M_A","sex_M_X", "asex_M_X"))
  
  
  
  Wilcox_sex_F_A_asex_F_A   <- wilcox.test(subset(df_GoTe, df_GoTe$rep_sex_XA == "sex_F_A")$meanFPKM,  subset(df_GoTe, df_GoTe$rep_sex_XA == "asex_F_A")$meanFPKM,)$p.value
  Wilcox_sex_F_X_asex_F_X   <- wilcox.test(subset(df_GoTe, df_GoTe$rep_sex_XA == "sex_F_X")$meanFPKM,  subset(df_GoTe, df_GoTe$rep_sex_XA == "asex_F_X")$meanFPKM,)$p.value
  Wilcox_sex_M_A_asex_M_A   <- wilcox.test(subset(df_GoTe, df_GoTe$rep_sex_XA == "sex_M_A")$meanFPKM,  subset(df_GoTe, df_GoTe$rep_sex_XA == "asex_M_A")$meanFPKM,)$p.value
  Wilcox_sex_M_X_asex_M_X   <- wilcox.test(subset(df_GoTe, df_GoTe$rep_sex_XA == "sex_M_X")$meanFPKM,  subset(df_GoTe, df_GoTe$rep_sex_XA == "asex_M_X")$meanFPKM,)$p.value
  
  wilcox_out <- as.data.frame(cbind(
    c(paste(deparse(substitute(df_sex))), paste(deparse(substitute(df_sex))), paste(deparse(substitute(df_sex))), paste(deparse(substitute(df_sex)))),
    c("sex_F_A_asex_F_A", "sex_F_X_asex_F_X", "sex_M_A_asex_M_A", "sex_M_X_asex_M_X"),
    c(  Wilcox_sex_F_A_asex_F_A, Wilcox_sex_F_X_asex_F_X, Wilcox_sex_M_A_asex_M_A, Wilcox_sex_M_X_asex_M_X )  ))
  
  colnames(wilcox_out) <- c("Df", "comp", "wilcox_p")
  
  P1 <- ggplot(df_GoTe, aes(tiss, meanFPKM)) +
    theme_classic() +
    geom_boxplot(aes(fill = factor(rep_sex_XA)), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T, position = position_dodge2(preserve = "single")) + 
    coord_cartesian(ylim=c(0,max_y)) + 
    ylab("FPKM") +
    xlab("") + 
    scale_fill_manual(values=c("sex_F_A" = "grey", "sex_F_X" = "firebrick2", "sex_M_A" = "white", "sex_M_X" = "royalblue2", "asex_F_A" = "azure4", "asex_F_X" = "firebrick4", "asex_M_A" = "gold1", "asex_M_X" = "mediumpurple3"), 
                      labels=c("sex_F_A" = " Sexual Female A","sex_F_X" =  "Sexual Female X", "sex_M_A" = "Sexual Male A", "sex_M_X" = "Sexual Male X", "asex_F_A" = " Asexual Female A","asex_F_X" =  "Asexual Female X", "asex_M_A" = "Asexual Male A", "asex_M_X" = "Asexual Male X")) +
    ggtitle(paste(deparse(substitute(df_sex))))
  
  outlist = list("P1" = P1, "wilcox_out" = wilcox_out) 
  return(outlist)	
  
}


png(filename = "TpsTdi_Ad_to_Tps_F_avFPKM_GoTe.png", width = 8, height = 4, units = "in", bg = "white", res = 300)
plot_avFPKM_testes_sex_asex(Tps_Ad_to_Tps_F_avFPKM, Tdi_Ad_to_Tps_F_avFPKM, 200)$P1
dev.off()
getwd() ## where

png(filename = "TceTms_Ad_to_Tce_F_avFPKM_GoTe.png", width = 8, height = 4, units = "in", bg = "white", res = 300)
plot_avFPKM_testes_sex_asex(Tce_Ad_to_Tce_F_avFPKM, Tms_Ad_to_Tce_F_avFPKM, 200)$P1
dev.off()
getwd() ## where

png(filename = "TpaTge_Ad_to_Tpa_F_avFPKM_GoTe.png", width = 8, height = 4, units = "in", bg = "white", res = 300)
plot_avFPKM_testes_sex_asex(Tpa_Ad_to_Tpa_F_avFPKM, Tge_Ad_to_Tpa_F_avFPKM, 200)$P1
dev.off()
getwd() ## where

wilcox_rep_FPKM <- as.data.frame(rbind(
  plot_avFPKM_testes_sex_asex(Tce_Ad_to_Tce_F_avFPKM, Tms_Ad_to_Tce_F_avFPKM, 200)$wilcox_out,
  plot_avFPKM_testes_sex_asex(Tpa_Ad_to_Tpa_F_avFPKM, Tge_Ad_to_Tpa_F_avFPKM, 200)$wilcox_out,
  plot_avFPKM_testes_sex_asex(Tps_Ad_to_Tps_F_avFPKM, Tdi_Ad_to_Tps_F_avFPKM, 200)$wilcox_out))

wilcox_rep_FPKM$FDR <- p.adjust(wilcox_rep_FPKM$wilcox_p, method = "fdr")

write.csv(wilcox_rep_FPKM, file = "wilcox_rep_FPKM.csv")
###########################################################################
### combined

sex_avFPKM_df  <- rbind(Tce_Ad_to_Tce_F_avFPKM[c(1:10, 12)], Tpa_Ad_to_Tpa_F_avFPKM, Tps_Ad_to_Tps_F_avFPKM)
asex_avFPKM_df <- rbind(Tms_Ad_to_Tce_F_avFPKM[c(1:10, 12)], Tge_Ad_to_Tpa_F_avFPKM, Tdi_Ad_to_Tps_F_avFPKM)



plot_avFPKM_testes_sex_asex_comb <- function(df_sex, df_asex, max_y){
  df <- rbind(df_sex, df_asex)
  df_Te <- subset(df, df$tiss == "Te")
  df_Go <- subset(df, df$tiss == "Go")
  
  df_GoTe <- rbind(df_Go, df_Te)
  df_GoTe$rep_sex_XA <- paste(df_GoTe$rep_m, df_GoTe$sex_XA, sep = "_")
  df_GoTe$rep_sex_XA <- ordered(df_GoTe$rep_sex_XA, levels = c("sex_F_A", "asex_F_A","sex_F_X", "asex_F_X","sex_M_A", "asex_M_A","sex_M_X", "asex_M_X"))
  
  
  
  df_GoTe$sp_group <- ifelse(df_GoTe$sp == "Tce", "Tce_Tms",
                             ifelse(df_GoTe$sp == "Tms", "Tce_Tms",
                                    ifelse(df_GoTe$sp == "Tpa", "Tpa_Tge",
                                           ifelse(df_GoTe$sp == "Tge", "Tpa_Tge",
                                                  ifelse(df_GoTe$sp == "Tps", "Tps_Tdi",
                                                         ifelse(df_GoTe$sp == "Tdi", "Tps_Tdi","ERROR"))))))
  
  df_GoTe$sp_tiss <- paste(df_GoTe$sp_group, df_GoTe$tiss)
  
  str(df_GoTe)
  P1 <- ggplot(df_GoTe, aes(sp_tiss, meanFPKM)) +
    theme_classic() +
    geom_boxplot(aes(fill = factor(rep_sex_XA)), width = 0.5, outlier.size = 0, outlier.shape = NA, notch=T, position = position_dodge2(preserve = "single")) +
    coord_cartesian(ylim=c(0,max_y)) +
    ylab("FPKM") +
    xlab("") +
    scale_fill_manual(values=c("sex_F_A" = "grey", "sex_F_X" = "firebrick2", "sex_M_A" = "white", "sex_M_X" = "royalblue2", "asex_F_A" = "azure4", "asex_F_X" = "firebrick4", "asex_M_A" = "gold1", "asex_M_X" = "mediumpurple3"),
                      labels=c("sex_F_A" = " Sexual Female A","sex_F_X" =  "Sexual Female X", "sex_M_A" = "Sexual Male A", "sex_M_X" = "Sexual Male X", "asex_F_A" = " Asexual Female A","asex_F_X" =  "Asexual Female X", "asex_M_A" = "Asexual Male A", "asex_M_X" = "Asexual Male X")) 
  # + ggtitle(paste(deparse(substitute(df_sex))))
  
  P2 <- P1 + theme(legend.position = "none") 
  outlist = list("P1" = P1, "P2" = P2)
  return(outlist)
  
}


plot_avFPKM_testes_sex_asex_comb(sex_avFPKM_df, asex_avFPKM_df, 200)$P1
plot_avFPKM_testes_sex_asex_comb(sex_avFPKM_df, asex_avFPKM_df, 200)$P2




### plot with MF


pdf("comb_AG_MF_plot_comb_Te_MF_plot_avFPKM_testes_sex_asex_comb.pdf", width = 7, height = 6.5)
plot_grid(
  comb_Te_MF_plot,
  plot_avFPKM_testes_sex_asex_comb(sex_avFPKM_df, asex_avFPKM_df, 250)$P2, ncol = 1
)
dev.off()
getwd() ## where has my plot gone....?

pdf("comb_AG_MF_plot_comb_Te_MF_plot_avFPKM_testes_sex_asex_comb_2.pdf", width = 7, height = 6.5)
plot_grid(
  comb_Te_MF_plot,
  plot_avFPKM_testes_sex_asex_comb(sex_avFPKM_df, asex_avFPKM_df, 250)$P1, ncol = 1
)
dev.off()
getwd() ## where has my plot gone....?


pdf("comb_AG_MF_plot_comb_Te_MF_plot_avFPKM_testes_sex_asex_comb_skinny.pdf", width = 6, height = 6.5)
plot_grid(
  comb_Te_MF_plot,
  plot_avFPKM_testes_sex_asex_comb(sex_avFPKM_df, asex_avFPKM_df, 250)$P2, ncol = 1
)
dev.off()
getwd() ## where has my plot gone....?






















#################################################################################################################

want_scafs_Tps <- c(
  "Tps_LRv5b_scf1",
  "Tps_LRv5b_scf2",
  "Tps_LRv5b_scf3",
  "Tps_LRv5b_scf4",
  "Tps_LRv5b_scf5",
  "Tps_LRv5b_scf6",
  "Tps_LRv5b_scf7",
  "Tps_LRv5b_scf8",
  "Tps_LRv5b_scf9",
  "Tps_LRv5b_scf10",
  "Tps_LRv5b_scf11",
  "Tps_LRv5b_scf12")

want_scafs_Tpa <- c(
  "Tpa_LRv5a_scf1",
  "Tpa_LRv5a_scf2",
  "Tpa_LRv5a_scf3",
  "Tpa_LRv5a_scf4",
  "Tpa_LRv5a_scf5",
  "Tpa_LRv5a_scf6",
  "Tpa_LRv5a_scf7",
  "Tpa_LRv5a_scf8",
  "Tpa_LRv5a_scf9",
  "Tpa_LRv5a_scf10",
  "Tpa_LRv5a_scf11",
  "Tpa_LRv5a_scf12",
  "Tpa_LRv5a_scf13",
  "Tpa_LRv5a_scf14")


want_LGs_Tce <- c(
  "HiCLG1",
  "HiCLG2",
  "HiCLG3",
  "HiCLG4",
  "HiCLG5",
  "HiCLG6",
  "HiCLG7",
  "HiCLG8",
  "HiCLG9",
  "HiCLG10",
  "HiCLG11",
  "HiCLG12",
  "HiCLG13")


Tps_Ad_to_Tps_F_avFPKM$scaf <- ordered(Tps_Ad_to_Tps_F_avFPKM$scaf, levels = want_scafs_Tps)
Tdi_Ad_to_Tps_F_avFPKM$scaf <- ordered(Tdi_Ad_to_Tps_F_avFPKM$scaf, levels = want_scafs_Tps)

Tpa_Ad_to_Tpa_F_avFPKM$scaf <- ordered(Tpa_Ad_to_Tpa_F_avFPKM$scaf, levels = want_scafs_Tpa)
Tge_Ad_to_Tpa_F_avFPKM$scaf <- ordered(Tge_Ad_to_Tpa_F_avFPKM$scaf, levels = want_scafs_Tpa)

Tce_Ad_to_Tce_F_avFPKM$LG <- ordered(Tce_Ad_to_Tce_F_avFPKM$LG, levels = want_LGs_Tce)
Tms_Ad_to_Tce_F_avFPKM$LG <- ordered(Tms_Ad_to_Tce_F_avFPKM$LG, levels = want_LGs_Tce)


### Tps scaf

png(filename = "TpsTdi_Ad_to_Tps_F_avFPKM_scaf_A.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "A")$P1,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "A")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpsTdi_Ad_to_Tps_F_avFPKM_scaf_B.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "B")$P1,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "B")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpsTdi_Ad_to_Tps_F_avFPKM_scaf_DG.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "DG")$P1,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "DG")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpsTdi_Ad_to_Tps_F_avFPKM_scaf_FB.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "FB")$P1,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "FB")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpsTdi_Ad_to_Tps_F_avFPKM_scaf_Fe.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Fe")$P1,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Fe")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpsTdi_Ad_to_Tps_F_avFPKM_scaf_Gu.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Gu")$P1,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Gu")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpsTdi_Ad_to_Tps_F_avFPKM_scaf_Te.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Te")$P1,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Te")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...


png(filename = "TpsTdi_Ad_to_Tps_F_avFPKM_scaf_Go.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Go")$P1,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Go")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...


### Tce LG

png(filename = "TceTms_Ad_to_Tce_F_avFPKM_LG_A.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "A")$P1,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "A")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TceTms_Ad_to_Tce_F_avFPKM_LG_B.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "B")$P1,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "B")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TceTms_Ad_to_Tce_F_avFPKM_LG_DG.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "DG")$P1,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "DG")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TceTms_Ad_to_Tce_F_avFPKM_LG_FB.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "FB")$P1,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "FB")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TceTms_Ad_to_Tce_F_avFPKM_LG_Fe.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Fe")$P1,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Fe")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TceTms_Ad_to_Tce_F_avFPKM_LG_Gu.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Gu")$P1,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Gu")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TceTms_Ad_to_Tce_F_avFPKM_LG_Te.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Te")$P1,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Te")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TceTms_Ad_to_Tce_F_avFPKM_LG_Go.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Go")$P1,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Go")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

### Tpa scaf

png(filename = "TpaTge_Ad_to_Tpa_F_avFPKM_scaf_A.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "A")$P1,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "A")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpaTge_Ad_to_Tpa_F_avFPKM_scaf_B.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "B")$P1,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "B")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpaTge_Ad_to_Tpa_F_avFPKM_scaf_DG.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "DG")$P1,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "DG")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpaTge_Ad_to_Tpa_F_avFPKM_scaf_FB.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "FB")$P1,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "FB")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpaTge_Ad_to_Tpa_F_avFPKM_scaf_Fe.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Fe")$P1,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Fe")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpaTge_Ad_to_Tpa_F_avFPKM_scaf_Gu.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Gu")$P1,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Gu")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpaTge_Ad_to_Tpa_F_avFPKM_scaf_Te.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Te")$P1,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Te")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...


png(filename = "TpaTge_Ad_to_Tpa_F_avFPKM_scaf_Go.png", width = 8, height = 10, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Go")$P1,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Go")$P1, ncol = 1)
dev.off()
getwd() ## where has my plot gone...





### Tps scaf

png(filename = "TpsTdi_Ad_to_Tps_F_avFPKM_scaf_ALL.png", width = 10, height = 14, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "A")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "A")$P2,
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "B")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "B")$P2,
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "DG")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "DG")$P2,
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "FB")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "FB")$P2,
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Fe")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Fe")$P2,
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Gu")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Gu")$P2,
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Go")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Go")$P2,
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Te")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Te")$P2,
  ncol = 4)
dev.off()
getwd() ## where has my plot gone...

### Tce scaf

png(filename = "TceTms_Ad_to_Tce_F_avFPKM_LG_ALL.png", width = 10, height = 14, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "A")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "A")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "B")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "B")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "DG")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "DG")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "FB")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "FB")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Fe")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Fe")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Gu")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Gu")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Go")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Go")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Te")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Te")$P2, ncol = 4)
dev.off()
getwd() ## where has my plot gone...

### Tpa scaf

png(filename = "TpaTge_Ad_to_Tpa_F_avFPKM_scaf_ALL.png", width = 10, height = 14, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "A")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "A")$P2,
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "B")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "B")$P2,
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "DG")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "DG")$P2,
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "FB")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "FB")$P2,
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Fe")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Fe")$P2,
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Gu")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Gu")$P2,
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Go")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Go")$P2,
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Te")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Te")$P2,  
  ncol = 4)
dev.off()
getwd() ## where has my plot gone...




############ plot sex asex sp sep



### Tps scaf

png(filename = "Tps_Ad_to_Tps_F_avFPKM_scaf_ALL.png", width = 10, height = 14, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "A")$P2,
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "B")$P2,
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "DG")$P2,
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "FB")$P2,
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Fe")$P2,
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Gu")$P2,
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Go")$P2,
  plot_avFPKM_scafs(Tps_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Te")$P2,
  ncol = 3)
dev.off()
getwd() ## where has my plot gone...

png(filename = "Tdi_Ad_to_Tps_F_avFPKM_scaf_ALL.png", width = 10, height = 14, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "A")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "B")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "DG")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "FB")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Fe")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Gu")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Go")$P2,
  plot_avFPKM_scafs(Tdi_Ad_to_Tps_F_avFPKM , 200, want_scafs_Tps, "Te")$P2,
  ncol = 3)
dev.off()
getwd() ## where has my plot gone...



### Tce scaf

png(filename = "TceTms_Ad_to_Tce_F_avFPKM_LG_ALL.png", width = 10, height = 14, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "A")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "A")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "B")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "B")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "DG")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "DG")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "FB")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "FB")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Fe")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Fe")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Gu")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Gu")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Go")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Go")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Te")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Te")$P2, ncol = 3)
dev.off()
getwd() ## where has my plot gone...

## Tce scaf

png(filename = "Tce_Ad_to_Tce_F_avFPKM_LG_ALL.png", width = 10, height = 14, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "A")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "B")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "DG")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "FB")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Fe")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Gu")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Go")$P2,
  plot_avFPKM_LGs(Tce_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Te")$P2, ncol = 3)
dev.off()
getwd() ## where has my plot gone...

png(filename = "Tms_Ad_to_Tce_F_avFPKM_LG_ALL.png", width = 10, height = 14, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "A")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "B")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "DG")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "FB")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Fe")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Gu")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Go")$P2,
  plot_avFPKM_LGs(Tms_Ad_to_Tce_F_avFPKM , 200, want_LGs_Tce, "Te")$P2, ncol = 3)
dev.off()
getwd() ## where has my plot gone...



### Tpa scaf

png(filename = "Tpa_Ad_to_Tpa_F_avFPKM_scaf_ALL.png", width = 10, height = 14, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "A")$P2,
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "B")$P2,
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "DG")$P2,
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "FB")$P2,
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Fe")$P2,
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Gu")$P2,
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Go")$P2,
  plot_avFPKM_scafs(Tpa_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Te")$P2,
  ncol = 3)
dev.off()
getwd() ## where has my plot gone...


png(filename = "Tge_Ad_to_Tpa_F_avFPKM_scaf_ALL.png", width = 10, height = 14, units = "in", bg = "white", res = 300)
plot_grid(
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "A")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "B")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "DG")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "FB")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Fe")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Gu")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Go")$P2,
  plot_avFPKM_scafs(Tge_Ad_to_Tpa_F_avFPKM , 200, want_scafs_Tpa, "Te")$P2,  
  ncol = 3)
dev.off()
getwd() ## where has my plot gone...


















################################################################
## MF - scaf 




## Tps

Tps_Ad_to_Tps_F_MF <- as.data.frame(rbind(
  extract_avMF(Tps_A_Ad_to_Tps_F_FPKM),
  extract_avMF(Tps_B_Ad_to_Tps_F_FPKM),
  extract_avMF(Tps_DG_Ad_to_Tps_F_FPKM),
  extract_avMF(Tps_FB_Ad_to_Tps_F_FPKM),
  extract_avMF(Tps_Fe_Ad_to_Tps_F_FPKM),
  extract_avMF(Tps_Gu_Ad_to_Tps_F_FPKM),
  extract_avMF(Tps_Ta_Ad_to_Tps_F_FPKM),
  extract_avMF(Tps_GoTe_Ad_to_Tps_F_FPKM)  
))



## Tdi

Tdi_Ad_to_Tps_F_MF <- as.data.frame(rbind(
  extract_avMF(Tdi_A_Ad_to_Tps_F_FPKM),
  extract_avMF(Tdi_B_Ad_to_Tps_F_FPKM),
  extract_avMF(Tdi_DG_Ad_to_Tps_F_FPKM),
  extract_avMF(Tdi_FB_Ad_to_Tps_F_FPKM),
  extract_avMF(Tdi_Fe_Ad_to_Tps_F_FPKM),
  extract_avMF(Tdi_Gu_Ad_to_Tps_F_FPKM),
  extract_avMF(Tdi_Ta_Ad_to_Tps_F_FPKM),
  extract_avMF(Tdi_GoTe_Ad_to_Tps_F_FPKM)
))


## Tce

Tce_Ad_to_Tce_F_MF <- as.data.frame(rbind(
  extract_avMF(Tce_A_Ad_to_Tce_F_FPKM),
  extract_avMF(Tce_B_Ad_to_Tce_F_FPKM),
  extract_avMF(Tce_DG_Ad_to_Tce_F_FPKM),
  extract_avMF(Tce_FB_Ad_to_Tce_F_FPKM),
  extract_avMF(Tce_Fe_Ad_to_Tce_F_FPKM),
  extract_avMF(Tce_Gu_Ad_to_Tce_F_FPKM),
  extract_avMF(Tce_Ta_Ad_to_Tce_F_FPKM),
  extract_avMF(Tce_GoTe_Ad_to_Tce_F_FPKM)
))

## Tms

Tms_Ad_to_Tce_F_MF <- as.data.frame(rbind(
  
  extract_avMF(Tms_A_Ad_to_Tce_F_FPKM),
  extract_avMF(Tms_B_Ad_to_Tce_F_FPKM),
  extract_avMF(Tms_DG_Ad_to_Tce_F_FPKM),
  extract_avMF(Tms_FB_Ad_to_Tce_F_FPKM),
  extract_avMF(Tms_Fe_Ad_to_Tce_F_FPKM),
  extract_avMF(Tms_Gu_Ad_to_Tce_F_FPKM),
  extract_avMF(Tms_Ta_Ad_to_Tce_F_FPKM),
  extract_avMF(Tms_GoTe_Ad_to_Tce_F_FPKM)
))


## Tpa

Tpa_Ad_to_Tpa_F_MF <- as.data.frame(rbind(
  extract_avMF(Tpa_A_Ad_to_Tpa_F_FPKM),
  extract_avMF(Tpa_B_Ad_to_Tpa_F_FPKM),
  extract_avMF(Tpa_DG_Ad_to_Tpa_F_FPKM),
  extract_avMF(Tpa_FB_Ad_to_Tpa_F_FPKM),
  extract_avMF(Tpa_Fe_Ad_to_Tpa_F_FPKM),
  extract_avMF(Tpa_Gu_Ad_to_Tpa_F_FPKM),
  extract_avMF(Tpa_Ta_Ad_to_Tpa_F_FPKM),
  extract_avMF(Tpa_GoTe_Ad_to_Tpa_F_FPKM)
))

## Tge

Tge_Ad_to_Tpa_F_MF <- as.data.frame(rbind(
  
  extract_avMF(Tge_A_Ad_to_Tpa_F_FPKM),
  extract_avMF(Tge_B_Ad_to_Tpa_F_FPKM),
  extract_avMF(Tge_DG_Ad_to_Tpa_F_FPKM),
  extract_avMF(Tge_FB_Ad_to_Tpa_F_FPKM),
  extract_avMF(Tge_Fe_Ad_to_Tpa_F_FPKM),
  extract_avMF(Tge_Gu_Ad_to_Tpa_F_FPKM),
  extract_avMF(Tge_Ta_Ad_to_Tpa_F_FPKM),
  extract_avMF(Tge_GoTe_Ad_to_Tpa_F_FPKM)
))


Tps_Ad_to_Tps_F_MF$scaf <- ordered(as.factor(Tps_Ad_to_Tps_F_MF$scaf), levels = want_scafs_Tps)
Tdi_Ad_to_Tps_F_MF$scaf <- ordered(as.factor(Tdi_Ad_to_Tps_F_MF$scaf), levels = want_scafs_Tps)

Tce_Ad_to_Tce_F_MF$LG <- ordered(as.factor(Tce_Ad_to_Tce_F_MF$LG), levels = want_LGs_Tce)
Tms_Ad_to_Tce_F_MF$LG <- ordered(as.factor(Tms_Ad_to_Tce_F_MF$LG), levels = want_LGs_Tce)

Tpa_Ad_to_Tpa_F_MF$scaf <- ordered(as.factor(Tpa_Ad_to_Tpa_F_MF$scaf), levels = want_scafs_Tpa)
Tge_Ad_to_Tpa_F_MF$scaf <- ordered(as.factor(Tge_Ad_to_Tpa_F_MF$scaf), levels = want_scafs_Tpa)


### to_Tps
png(filename = "MF_box_to_Tps_F_MF.png", width = 10, height = 20, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "A")$P2, plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "A")$P2,
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "B")$P2, plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "B")$P2,
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "DG")$P2, plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "DG")$P2,
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "FB")$P2, plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "FB")$P2,
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "Fe")$P2, plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "Fe")$P2,
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "Gu")$P2, plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "Gu")$P2,
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "Ta")$P2, plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "Ta")$P2,
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "GoTe")$P2, plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "GoTe")$P2,
  ncol = 4)
dev.off()
getwd() ## where has my plot gone...





### to_Tce
png(filename = "MF_box_to_Tce_F_MF.png", width = 10, height = 20, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "A")$P2, plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "A")$P2,
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "B")$P2, plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "B")$P2,
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "DG")$P2, plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "DG")$P2,
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "FB")$P2, plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "FB")$P2,
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "Fe")$P2, plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "Fe")$P2,
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "Gu")$P2, plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "Gu")$P2,
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "Ta")$P2, plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "Ta")$P2,
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "GoTe")$P2, plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "GoTe")$P2,
  ncol = 4)
dev.off()
getwd() ## where has my plot gone...



### to_Tpa
png(filename = "MF_box_to_Tpa_F_MF.png", width = 10, height = 20, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "A")$P2, plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "A")$P2,
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "B")$P2, plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "B")$P2,
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "DG")$P2, plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "DG")$P2,
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "FB")$P2, plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "FB")$P2,
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "Fe")$P2, plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "Fe")$P2,
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "Gu")$P2, plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "Gu")$P2,
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "Ta")$P2, plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "Ta")$P2,
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "GoTe")$P2, plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "GoTe")$P2,
  ncol = 4)
dev.off()
getwd() ## where has my plot gone...



############################## sep

### to_Tps
png(filename = "Tps_MF_box_to_Tps_F_MF.png", width = 10, height = 15, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "A")$P2,
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "B")$P2, 
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "DG")$P2,
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "FB")$P2, 
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "Fe")$P2,
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "Gu")$P2,
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "Ta")$P2,
  plot_MF_scafs(Tps_Ad_to_Tps_F_MF, want_scafs_Tps, "GoTe")$P2,
  ncol = 3)
dev.off()
getwd() ## where has my plot gone...


png(filename = "Tdi_MF_box_to_Tps_F_MF.png", width = 10, height = 15, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "A")$P2,
  plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "B")$P2,
  plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "DG")$P2,
  plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "FB")$P2,
  plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "Fe")$P2,
  plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "Gu")$P2,
  plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "Ta")$P2,
  plot_MF_scafs(Tdi_Ad_to_Tps_F_MF, want_scafs_Tps, "GoTe")$P2,
  ncol = 3)
dev.off()
getwd() ## where has my plot gone...



### to_Tce
png(filename = "Tce_MF_box_to_Tce_F_MF.png", width = 10, height = 15, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "A")$P2,
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "B")$P2,
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "DG")$P2,
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "FB")$P2,
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "Fe")$P2,
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "Gu")$P2,
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "Ta")$P2,
  plot_MF_LGs(Tce_Ad_to_Tce_F_MF, want_LGs_Tce, "GoTe")$P2,
  ncol = 3)
dev.off()
getwd() ## where has my plot gone...

png(filename = "Tms_MF_box_to_Tce_F_MF.png", width = 10, height = 15, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "A")$P2,
  plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "B")$P2,
  plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "DG")$P2,
  plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "FB")$P2,
  plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "Fe")$P2,
  plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "Gu")$P2,
  plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "Ta")$P2,
  plot_MF_LGs(Tms_Ad_to_Tce_F_MF, want_LGs_Tce, "GoTe")$P2,
  ncol = 3)
dev.off()
getwd() ## where has my plot gone...



### to_Tpa
png(filename = "Tpa_MF_box_to_Tpa_F_MF.png", width = 10, height = 15, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "A")$P2,
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "B")$P2,
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "DG")$P2,
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "FB")$P2,
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "Fe")$P2,
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "Gu")$P2,
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "Ta")$P2,
  plot_MF_scafs(Tpa_Ad_to_Tpa_F_MF, want_scafs_Tpa, "GoTe")$P2,
  ncol = 3)
dev.off()
getwd() ## where has my plot gone...


### to_Tpa
png(filename = "Tge_MF_box_to_Tpa_F_MF.png", width = 10, height = 15, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "A")$P2,
  plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "B")$P2,
  plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "DG")$P2,
  plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "FB")$P2,
  plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "Fe")$P2,
  plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "Gu")$P2,
  plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "Ta")$P2,
  plot_MF_scafs(Tge_Ad_to_Tpa_F_MF, want_scafs_Tpa, "GoTe")$P2,
  ncol = 3)
dev.off()
getwd() ## where has my plot gone...

















































######################################################################################################################
######################################################################################################################
#### plot along chr

# heatmap

## stitch together how I want with python

write.csv(Tps_Ad_to_Tps_F_MF, file = "Tps_Ad_to_Tps_F_MF.csv", quote = F, row.names = F)
write.csv(Tdi_Ad_to_Tps_F_MF, file = "Tdi_Ad_to_Tps_F_MF.csv", quote = F, row.names = F)

Tps_Ad_to_Tps_F_avFPKM_2 <- as.data.frame(rbind(
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
  extract_avEXP_GoTe(Tps_GoTe_Ad_to_Tps_F_FPKM, "M", "Te"),
  extract_avEXP_GoTe(Tps_GoTe_Ad_to_Tps_F_FPKM, "F", "Go")
))

Tdi_Ad_to_Tps_F_avFPKM_2 <- as.data.frame(rbind(
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
  extract_avEXP_GoTe(Tdi_GoTe_Ad_to_Tps_F_FPKM, "M", "Te"),
  extract_avEXP_GoTe(Tdi_GoTe_Ad_to_Tps_F_FPKM, "F", "Go")
))

write.csv(Tps_Ad_to_Tps_F_avFPKM_2, file = "Tps_Ad_to_Tps_F_avFPKM_2.csv", quote = F, row.names = F)
write.csv(Tdi_Ad_to_Tps_F_avFPKM_2, file = "Tdi_Ad_to_Tps_F_avFPKM_2.csv", quote = F, row.names = F)


##################
write.csv(Tce_Ad_to_Tce_F_MF, file = "Tce_Ad_to_Tce_F_MF.csv", quote = F, row.names = F)
write.csv(Tms_Ad_to_Tce_F_MF, file = "Tms_Ad_to_Tce_F_MF.csv", quote = F, row.names = F)

### use the GoTe version here
Tce_Ad_to_Tce_F_avFPKM_2 <- as.data.frame(rbind(
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
  extract_avEXP_GoTe(Tce_GoTe_Ad_to_Tce_F_FPKM, "M", "Te"),
  extract_avEXP_GoTe(Tce_GoTe_Ad_to_Tce_F_FPKM, "F", "Go")
))

### use the GoTe version here
Tms_Ad_to_Tce_F_avFPKM_2 <- as.data.frame(rbind(
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
  extract_avEXP_GoTe(Tms_GoTe_Ad_to_Tce_F_FPKM, "M", "Te"),
  extract_avEXP_GoTe(Tms_GoTe_Ad_to_Tce_F_FPKM, "F", "Go")
))

write.csv(Tce_Ad_to_Tce_F_avFPKM_2, file = "Tce_Ad_to_Tce_F_avFPKM_2.csv", quote = F, row.names = F)
write.csv(Tms_Ad_to_Tce_F_avFPKM_2, file = "Tms_Ad_to_Tce_F_avFPKM_2.csv", quote = F, row.names = F)


##########################
write.csv(Tpa_Ad_to_Tpa_F_MF, file = "Tpa_Ad_to_Tpa_F_MF.csv", quote = F, row.names = F)
write.csv(Tge_Ad_to_Tpa_F_MF, file = "Tge_Ad_to_Tpa_F_MF.csv", quote = F, row.names = F)

### use the GoTe version here
Tpa_Ad_to_Tpa_F_avFPKM_2 <- as.data.frame(rbind(
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
  extract_avEXP_GoTe(Tpa_GoTe_Ad_to_Tpa_F_FPKM, "M", "Te"),
  extract_avEXP_GoTe(Tpa_GoTe_Ad_to_Tpa_F_FPKM, "F", "Go")
))

### use the GoTe version here
Tge_Ad_to_Tpa_F_avFPKM_2 <- as.data.frame(rbind(
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
  extract_avEXP_GoTe(Tge_GoTe_Ad_to_Tpa_F_FPKM, "M", "Te"),
  extract_avEXP_GoTe(Tge_GoTe_Ad_to_Tpa_F_FPKM, "F", "Go")
))

write.csv(Tpa_Ad_to_Tpa_F_avFPKM_2, file = "Tpa_Ad_to_Tpa_F_avFPKM_2.csv", quote = F, row.names = F)
write.csv(Tge_Ad_to_Tpa_F_avFPKM_2, file = "Tge_Ad_to_Tpa_F_avFPKM_2.csv", quote = F, row.names = F)








#### use 
# accessory_scripts/join_MF_for_heatmap.py   
# accessory_scripts/join_FPKM_for_heatmap.py 
# 
# python3 accessory_scripts/join_MF_for_heatmap.py   -s output/sex_ref_combSA_v2/Tps_Ad_to_Tps_F_MF.csv       -a output/sex_ref_combSA_v2/Tdi_Ad_to_Tps_F_MF.csv       -o output/sex_ref_combSA_v2/TpsTdi_Ad_to_Tps_F_MF
# python3 accessory_scripts/join_FPKM_for_heatmap.py -s output/sex_ref_combSA_v2/Tps_Ad_to_Tps_F_avFPKM_2.csv -a output/sex_ref_combSA_v2/Tdi_Ad_to_Tps_F_avFPKM_2.csv -o output/sex_ref_combSA_v2/TpsTdi_Ad_to_Tps_F_avFPKM_2
# 
# python3 accessory_scripts/join_MF_for_heatmap.py   -s output/sex_ref_combSA_v2/Tce_Ad_to_Tce_F_MF.csv       -a output/sex_ref_combSA_v2/Tms_Ad_to_Tce_F_MF.csv       -o output/sex_ref_combSA_v2/TceTms_Ad_to_Tce_F_MF
# python3 accessory_scripts/join_FPKM_for_heatmap.py -s output/sex_ref_combSA_v2/Tce_Ad_to_Tce_F_avFPKM_2.csv -a output/sex_ref_combSA_v2/Tms_Ad_to_Tce_F_avFPKM_2.csv -o output/sex_ref_combSA_v2/TceTms_Ad_to_Tce_F_avFPKM_2
# 
# python3 accessory_scripts/join_MF_for_heatmap.py   -s output/sex_ref_combSA_v2/Tpa_Ad_to_Tpa_F_MF.csv       -a output/sex_ref_combSA_v2/Tge_Ad_to_Tpa_F_MF.csv       -o output/sex_ref_combSA_v2/TpaTge_Ad_to_Tpa_F_MF
# python3 accessory_scripts/join_FPKM_for_heatmap.py -s output/sex_ref_combSA_v2/Tpa_Ad_to_Tpa_F_avFPKM_2.csv -a output/sex_ref_combSA_v2/Tge_Ad_to_Tpa_F_avFPKM_2.csv -o output/sex_ref_combSA_v2/TpaTge_Ad_to_Tpa_F_avFPKM_2
# 
# 
# 
# 


######################################################################################
###
######## read back in


##### MF
TpsTdi_Ad_to_Tps_F_MF <- read.csv("TpsTdi_Ad_to_Tps_F_MF_fheatmap.csv")

TpsTdi_Ad_to_Tps_F_MF_X <- subset(TpsTdi_Ad_to_Tps_F_MF, TpsTdi_Ad_to_Tps_F_MF$XA == "X")
TpsTdi_Ad_to_Tps_F_MF_X_for_heat = TpsTdi_Ad_to_Tps_F_MF_X[5:length((TpsTdi_Ad_to_Tps_F_MF_X[1,]))]
rownames(TpsTdi_Ad_to_Tps_F_MF_X_for_heat) <- TpsTdi_Ad_to_Tps_F_MF_X[,1]
TpsTdi_Ad_to_Tps_F_MF_X_for_heat_no_NA <- na.omit(TpsTdi_Ad_to_Tps_F_MF_X_for_heat)

TpsTdi_Ad_to_Tps_F_MF_A <- subset(TpsTdi_Ad_to_Tps_F_MF, TpsTdi_Ad_to_Tps_F_MF$XA == "A")
TpsTdi_Ad_to_Tps_F_MF_A_for_heat = TpsTdi_Ad_to_Tps_F_MF_A[5:length((TpsTdi_Ad_to_Tps_F_MF_A[1,]))]
rownames(TpsTdi_Ad_to_Tps_F_MF_A_for_heat) <- TpsTdi_Ad_to_Tps_F_MF_A[,1]
TpsTdi_Ad_to_Tps_F_MF_A_for_heat_no_NA <- na.omit(TpsTdi_Ad_to_Tps_F_MF_A_for_heat)

#### FPKM
TpsTdi_Ad_to_Tps_F_FPKM <- read.csv("TpsTdi_Ad_to_Tps_F_avFPKM_2_fheatmap.csv")

TpsTdi_Ad_to_Tps_F_FPKM_X <- subset(TpsTdi_Ad_to_Tps_F_FPKM, TpsTdi_Ad_to_Tps_F_FPKM$XA == "X")
TpsTdi_Ad_to_Tps_F_FPKM_X_for_heat = TpsTdi_Ad_to_Tps_F_FPKM_X[5:length((TpsTdi_Ad_to_Tps_F_FPKM_X[1,]))]
rownames(TpsTdi_Ad_to_Tps_F_FPKM_X_for_heat) <- TpsTdi_Ad_to_Tps_F_FPKM_X[,1]
TpsTdi_Ad_to_Tps_F_FPKM_X_for_heat_no_NA <- na.omit(TpsTdi_Ad_to_Tps_F_FPKM_X_for_heat)

TpsTdi_Ad_to_Tps_F_FPKM_A <- subset(TpsTdi_Ad_to_Tps_F_FPKM, TpsTdi_Ad_to_Tps_F_FPKM$AA == "A")
TpsTdi_Ad_to_Tps_F_FPKM_A_for_heat = TpsTdi_Ad_to_Tps_F_FPKM_A[5:length((TpsTdi_Ad_to_Tps_F_FPKM_A[1,]))]
rownames(TpsTdi_Ad_to_Tps_F_FPKM_A_for_heat) <- TpsTdi_Ad_to_Tps_F_FPKM_A[,1]
TpsTdi_Ad_to_Tps_F_FPKM_A_for_heat_no_NA <- na.omit(TpsTdi_Ad_to_Tps_F_FPKM_A_for_heat)


cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

TpsTdi_Ad_to_Tps_F_FPKM_X_for_heat_no_NA_zscore <- t(apply(TpsTdi_Ad_to_Tps_F_FPKM_X_for_heat_no_NA, 1, cal_z_score))



##### MF
TceTms_Ad_to_Tce_F_MF <- read.csv("TceTms_Ad_to_Tce_F_MF_fheatmap.csv")

TceTms_Ad_to_Tce_F_MF_X <- subset(TceTms_Ad_to_Tce_F_MF, TceTms_Ad_to_Tce_F_MF$XA == "X")
TceTms_Ad_to_Tce_F_MF_X_for_heat = TceTms_Ad_to_Tce_F_MF_X[5:length((TceTms_Ad_to_Tce_F_MF_X[1,]))]
rownames(TceTms_Ad_to_Tce_F_MF_X_for_heat) <- TceTms_Ad_to_Tce_F_MF_X[,1]
TceTms_Ad_to_Tce_F_MF_X_for_heat_no_NA <- na.omit(TceTms_Ad_to_Tce_F_MF_X_for_heat)

TceTms_Ad_to_Tce_F_MF_A <- subset(TceTms_Ad_to_Tce_F_MF, TceTms_Ad_to_Tce_F_MF$XA == "A")
TceTms_Ad_to_Tce_F_MF_A_for_heat = TceTms_Ad_to_Tce_F_MF_A[5:length((TceTms_Ad_to_Tce_F_MF_A[1,]))]
rownames(TceTms_Ad_to_Tce_F_MF_A_for_heat) <- TceTms_Ad_to_Tce_F_MF_A[,1]
TceTms_Ad_to_Tce_F_MF_A_for_heat_no_NA <- na.omit(TceTms_Ad_to_Tce_F_MF_A_for_heat)

#### FPKM
TceTms_Ad_to_Tce_F_FPKM <- read.csv("TceTms_Ad_to_Tce_F_avFPKM_2_fheatmap.csv")

TceTms_Ad_to_Tce_F_FPKM_X <- subset(TceTms_Ad_to_Tce_F_FPKM, TceTms_Ad_to_Tce_F_FPKM$XA == "X")
TceTms_Ad_to_Tce_F_FPKM_X_for_heat = TceTms_Ad_to_Tce_F_FPKM_X[5:length((TceTms_Ad_to_Tce_F_FPKM_X[1,]))]
rownames(TceTms_Ad_to_Tce_F_FPKM_X_for_heat) <- TceTms_Ad_to_Tce_F_FPKM_X[,1]
TceTms_Ad_to_Tce_F_FPKM_X_for_heat_no_NA <- na.omit(TceTms_Ad_to_Tce_F_FPKM_X_for_heat)

TceTms_Ad_to_Tce_F_FPKM_A <- subset(TceTms_Ad_to_Tce_F_FPKM, TceTms_Ad_to_Tce_F_FPKM$AA == "A")
TceTms_Ad_to_Tce_F_FPKM_A_for_heat = TceTms_Ad_to_Tce_F_FPKM_A[5:length((TceTms_Ad_to_Tce_F_FPKM_A[1,]))]
rownames(TceTms_Ad_to_Tce_F_FPKM_A_for_heat) <- TceTms_Ad_to_Tce_F_FPKM_A[,1]
TceTms_Ad_to_Tce_F_FPKM_A_for_heat_no_NA <- na.omit(TceTms_Ad_to_Tce_F_FPKM_A_for_heat)


cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

TceTms_Ad_to_Tce_F_FPKM_X_for_heat_no_NA_zscore <- t(apply(TceTms_Ad_to_Tce_F_FPKM_X_for_heat_no_NA, 1, cal_z_score))


##### MF
TpaTge_Ad_to_Tpa_F_MF <- read.csv("TpaTge_Ad_to_Tpa_F_MF_fheatmap.csv")

TpaTge_Ad_to_Tpa_F_MF_X <- subset(TpaTge_Ad_to_Tpa_F_MF, TpaTge_Ad_to_Tpa_F_MF$XA == "X")
TpaTge_Ad_to_Tpa_F_MF_X_for_heat = TpaTge_Ad_to_Tpa_F_MF_X[5:length((TpaTge_Ad_to_Tpa_F_MF_X[1,]))]
rownames(TpaTge_Ad_to_Tpa_F_MF_X_for_heat) <- TpaTge_Ad_to_Tpa_F_MF_X[,1]
TpaTge_Ad_to_Tpa_F_MF_X_for_heat_no_NA <- na.omit(TpaTge_Ad_to_Tpa_F_MF_X_for_heat)

TpaTge_Ad_to_Tpa_F_MF_A <- subset(TpaTge_Ad_to_Tpa_F_MF, TpaTge_Ad_to_Tpa_F_MF$XA == "A")
TpaTge_Ad_to_Tpa_F_MF_A_for_heat = TpaTge_Ad_to_Tpa_F_MF_A[5:length((TpaTge_Ad_to_Tpa_F_MF_A[1,]))]
rownames(TpaTge_Ad_to_Tpa_F_MF_A_for_heat) <- TpaTge_Ad_to_Tpa_F_MF_A[,1]
TpaTge_Ad_to_Tpa_F_MF_A_for_heat_no_NA <- na.omit(TpaTge_Ad_to_Tpa_F_MF_A_for_heat)

#### FPKM
TpaTge_Ad_to_Tpa_F_FPKM <- read.csv("TpaTge_Ad_to_Tpa_F_avFPKM_2_fheatmap.csv")

TpaTge_Ad_to_Tpa_F_FPKM_X <- subset(TpaTge_Ad_to_Tpa_F_FPKM, TpaTge_Ad_to_Tpa_F_FPKM$XA == "X")
TpaTge_Ad_to_Tpa_F_FPKM_X_for_heat = TpaTge_Ad_to_Tpa_F_FPKM_X[5:length((TpaTge_Ad_to_Tpa_F_FPKM_X[1,]))]
rownames(TpaTge_Ad_to_Tpa_F_FPKM_X_for_heat) <- TpaTge_Ad_to_Tpa_F_FPKM_X[,1]
TpaTge_Ad_to_Tpa_F_FPKM_X_for_heat_no_NA <- na.omit(TpaTge_Ad_to_Tpa_F_FPKM_X_for_heat)

TpaTge_Ad_to_Tpa_F_FPKM_A <- subset(TpaTge_Ad_to_Tpa_F_FPKM, TpaTge_Ad_to_Tpa_F_FPKM$AA == "A")
TpaTge_Ad_to_Tpa_F_FPKM_A_for_heat = TpaTge_Ad_to_Tpa_F_FPKM_A[5:length((TpaTge_Ad_to_Tpa_F_FPKM_A[1,]))]
rownames(TpaTge_Ad_to_Tpa_F_FPKM_A_for_heat) <- TpaTge_Ad_to_Tpa_F_FPKM_A[,1]
TpaTge_Ad_to_Tpa_F_FPKM_A_for_heat_no_NA <- na.omit(TpaTge_Ad_to_Tpa_F_FPKM_A_for_heat)


cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

TpaTge_Ad_to_Tpa_F_FPKM_X_for_heat_no_NA_zscore <- t(apply(TpaTge_Ad_to_Tpa_F_FPKM_X_for_heat_no_NA, 1, cal_z_score))





##############################################################################################################################
### heatmaps
## set breaks ## 
breaksList = c(-1114.25, -3.75, -3.25, -2.75, -2.25, -1.75, -1.25, -0.75, -0.25, 0.25,  0.75,  1.25, 1.75, 2.25, 2.75,3.25,3.75,1114.25) 


### MF
pheatmap(t(TpsTdi_Ad_to_Tps_F_MF_X_for_heat), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(TpsTdi_Ad_to_Tps_F_MF_A_for_heat), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(TpsTdi_Ad_to_Tps_F_MF_X_for_heat_no_NA), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(TpsTdi_Ad_to_Tps_F_MF_A_for_heat_no_NA), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F)

### FPKM - not great... 
pheatmap(t(log2(TpsTdi_Ad_to_Tps_F_FPKM_X_for_heat)), cluster_rows = F, cluster_cols = F, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(log2(TpsTdi_Ad_to_Tps_F_FPKM_X_for_heat_no_NA)), cluster_rows = F, cluster_cols = F, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(TpsTdi_Ad_to_Tps_F_FPKM_X_for_heat_no_NA_zscore), cluster_rows = F, cluster_cols = F, show_colnames = F,border_color = NA, legend = F)

### MF
pheatmap(t(TceTms_Ad_to_Tce_F_MF_X_for_heat), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(TceTms_Ad_to_Tce_F_MF_A_for_heat), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(TceTms_Ad_to_Tce_F_MF_X_for_heat_no_NA), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(TceTms_Ad_to_Tce_F_MF_A_for_heat_no_NA), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F)

### FPKM - not great... 
pheatmap(t(log2(TceTms_Ad_to_Tce_F_FPKM_X_for_heat)), cluster_rows = F, cluster_cols = F, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(log2(TceTms_Ad_to_Tce_F_FPKM_X_for_heat_no_NA)), cluster_rows = F, cluster_cols = F, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(TceTms_Ad_to_Tce_F_FPKM_X_for_heat_no_NA_zscore), cluster_rows = F, cluster_cols = F, show_colnames = F,border_color = NA, legend = F)


### MF
pheatmap(t(TpaTge_Ad_to_Tpa_F_MF_X_for_heat), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(TpaTge_Ad_to_Tpa_F_MF_A_for_heat), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(TpaTge_Ad_to_Tpa_F_MF_X_for_heat_no_NA), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(TpaTge_Ad_to_Tpa_F_MF_A_for_heat_no_NA), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F)

### FPKM - not great... 
pheatmap(t(log2(TpaTge_Ad_to_Tpa_F_FPKM_X_for_heat)), cluster_rows = F, cluster_cols = F, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(log2(TpaTge_Ad_to_Tpa_F_FPKM_X_for_heat_no_NA)), cluster_rows = F, cluster_cols = F, show_colnames = F,border_color = NA, legend = F)
pheatmap(t(TpaTge_Ad_to_Tpa_F_FPKM_X_for_heat_no_NA_zscore), cluster_rows = F, cluster_cols = F, show_colnames = F,border_color = NA, legend = F)



png(filename = "TpsTdi_Ad_to_Tps_F_MF_X_heatmap.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TpsTdi_Ad_to_Tps_F_MF_X_for_heat), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F, main = "TpsTdi_Ad_to_Tps_F_MF_X_for_heat")
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpsTdi_Ad_to_Tps_F_MF_A_heatmap.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TpsTdi_Ad_to_Tps_F_MF_A_for_heat), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F, main = "TpsTdi_Ad_to_Tps_F_MF_A_for_heat")
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpsTdi_Ad_to_Tps_F_MF_X_heatmap_no_NA.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TpsTdi_Ad_to_Tps_F_MF_X_for_heat_no_NA), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F, main = "TpsTdi_Ad_to_Tps_F_MF_X_for_heat_no_NA")
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpsTdi_Ad_to_Tps_F_MF_A_heatmap_no_NA.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TpsTdi_Ad_to_Tps_F_MF_A_for_heat_no_NA), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F, main = "TpsTdi_Ad_to_Tps_F_MF_A_for_heat_no_NA")
dev.off()
getwd() ## where has my plot gone...


dev.off()


png(filename = "TceTms_Ad_to_Tce_F_MF_X_heatmap.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TceTms_Ad_to_Tce_F_MF_X_for_heat), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F, main = "TceTms_Ad_to_Tce_F_MF_X_for_heat")
dev.off()
getwd() ## where has my plot gone...

png(filename = "TceTms_Ad_to_Tce_F_MF_A_heatmap.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TceTms_Ad_to_Tce_F_MF_A_for_heat), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F, main = "TceTms_Ad_to_Tce_F_MF_A_for_heat")
dev.off()
getwd() ## where has my plot gone...

png(filename = "TceTms_Ad_to_Tce_F_MF_X_heatmap_no_NA.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TceTms_Ad_to_Tce_F_MF_X_for_heat_no_NA), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F, main = "TceTms_Ad_to_Tce_F_MF_X_heatmap_no_NA")
dev.off()
getwd() ## where has my plot gone...

png(filename = "TceTms_Ad_to_Tce_F_MF_A_heatmap_no_NA.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TceTms_Ad_to_Tce_F_MF_A_for_heat_no_NA), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F, main = "TceTms_Ad_to_Tce_F_MF_A_heatmap_no_NA")
dev.off()
getwd() ## where has my plot gone...


png(filename = "TpaTge_Ad_to_Tpa_F_MF_X_heatmap.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TpaTge_Ad_to_Tpa_F_MF_X_for_heat), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F, main = "TpaTge_Ad_to_Tpa_F_MF_X_heatmap")
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpaTge_Ad_to_Tpa_F_MF_A_heatmap.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TpaTge_Ad_to_Tpa_F_MF_A_for_heat), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F, main = "TpaTge_Ad_to_Tpa_F_MF_A_heatmap")
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpaTge_Ad_to_Tpa_F_MF_X_heatmap_no_NA.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TpaTge_Ad_to_Tpa_F_MF_X_for_heat_no_NA), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F, main = "TpaTge_Ad_to_Tpa_F_MF_X_for_heat_no_NA")
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpaTge_Ad_to_Tpa_F_MF_A_heatmap_no_NA.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TpaTge_Ad_to_Tpa_F_MF_A_for_heat_no_NA), cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(c("#08103A", "#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D")))(length(breaksList)), breaks = breaksList, show_colnames = F,border_color = NA, legend = F, main = "TpaTge_Ad_to_Tpa_F_MF_A_for_heat_no_NA")
dev.off()
getwd() ## where has my plot gone...


#############################################################################################
### plot lines across.

plot_FPKM_alongchr <- function(df, chr_want, tiss_want){
  df_chr1 <- subset(df, df$tiss == tiss_want)  
  df_chr <- subset(df_chr1, df_chr1$scaf == chr_want)
  df_chr$mid_point = as.numeric(df_chr$mid_point)
  print(head(df_chr))
  print(str(df_chr))  
  ggplot(df_chr, aes(x=mid_point, log2(meanFPKM), col = sex)) +
    theme_bw() +
    geom_line(alpha=0.5) +
    ggtitle(paste("scaf", chr_want)) + 
    #scale_colour_manual(values=c("royalblue2", "red3")) +
    xlab("Position") # + ylim(0,3)
}



plot_MF_alongchr <- function(df, chr_want){
  df_chr <- subset(df, df$scaf == chr_want)
  df_chr$mid_point = as.numeric(df_chr$mid_point)
  print(head(df_chr))
  #print(str(df_chr))  
  p1 <- ggplot(df_chr, aes(x=mid_point, log2MF, col = tiss)) +
    theme_bw() +
    geom_line(alpha=0.5) +
    ggtitle(paste("scaf", chr_want)) + 
    #scale_colour_manual(values=c("royalblue2", "red3")) +
    xlab("Position") # + ylim(0,3)
  return(p1)
}


plot_MF_alongchr_tiss <- function(df, chr_want, tiss_want){
  df_chr_1 <- subset(df,       df$scaf == chr_want)
  df_chr   <- subset(df_chr_1, df_chr_1$tiss == tiss_want)
  df_chr$mid_point = as.numeric(df_chr$mid_point)
  print(head(df_chr))
  
  max_val = max(abs(df_chr_1$log2MF)) ## so same for all tissues
  print(max_val)  
  
  p1 <- ggplot(df_chr, aes(x=mid_point, log2MF)) +
    theme_bw() + geom_hline(yintercept = 0) +
    geom_line(alpha=0.8, colour = "darkred") +
    ggtitle(tiss_want) + 
    xlab("Position") + ylim(-9,9)  # ylim(-max_val,max_val) 
  return(p1)
}


### to_Tps

title_Tps_LRv5b_scf3 <- ggdraw() + draw_label("Tps - Tps_LRv5b_scf3", fontface='bold')

Tps_MFalongX_to_Tps_Tps_LRv5b_scf3 <- plot_grid(
  plot_MF_alongchr_tiss(Tps_Ad_to_Tps_F_MF, "Tps_LRv5b_scf3", "A"),
  plot_MF_alongchr_tiss(Tps_Ad_to_Tps_F_MF, "Tps_LRv5b_scf3", "B"),
  plot_MF_alongchr_tiss(Tps_Ad_to_Tps_F_MF, "Tps_LRv5b_scf3", "DG"),
  plot_MF_alongchr_tiss(Tps_Ad_to_Tps_F_MF, "Tps_LRv5b_scf3", "FB"),
  plot_MF_alongchr_tiss(Tps_Ad_to_Tps_F_MF, "Tps_LRv5b_scf3", "Fe"),
  plot_MF_alongchr_tiss(Tps_Ad_to_Tps_F_MF, "Tps_LRv5b_scf3", "Gu"),
  plot_MF_alongchr_tiss(Tps_Ad_to_Tps_F_MF, "Tps_LRv5b_scf3", "GoTe"),
  ncol = 1)

png(filename = "Tps_MFalongX_to_Tps_F.png", width = 10, height = 15, units = "in", bg = "white", res = 300)
plot_grid(title_Tps_LRv5b_scf3, Tps_MFalongX_to_Tps_Tps_LRv5b_scf3, ncol=1, rel_heights=c(0.1, 4))
dev.off()
getwd() ## where has my plot gone...


title_Tps_LRv5b_scf3 <- ggdraw() + draw_label("Tdi - Tps_LRv5b_scf3", fontface='bold')

Tdi_MFalongX_to_Tps_Tps_LRv5b_scf3 <- plot_grid(
  plot_MF_alongchr_tiss(Tdi_Ad_to_Tps_F_MF, "Tps_LRv5b_scf3", "A"),
  plot_MF_alongchr_tiss(Tdi_Ad_to_Tps_F_MF, "Tps_LRv5b_scf3", "B"),
  plot_MF_alongchr_tiss(Tdi_Ad_to_Tps_F_MF, "Tps_LRv5b_scf3", "DG"),
  plot_MF_alongchr_tiss(Tdi_Ad_to_Tps_F_MF, "Tps_LRv5b_scf3", "FB"),
  plot_MF_alongchr_tiss(Tdi_Ad_to_Tps_F_MF, "Tps_LRv5b_scf3", "Fe"),
  plot_MF_alongchr_tiss(Tdi_Ad_to_Tps_F_MF, "Tps_LRv5b_scf3", "Gu"),
  plot_MF_alongchr_tiss(Tdi_Ad_to_Tps_F_MF, "Tps_LRv5b_scf3", "GoTe"),
  ncol = 1)

png(filename = "Tdi_MFalongX_to_Tps_F.png", width = 10, height = 15, units = "in", bg = "white", res = 300)
plot_grid(title_Tps_LRv5b_scf3, Tdi_MFalongX_to_Tps_Tps_LRv5b_scf3, ncol=1, rel_heights=c(0.1, 4))
dev.off()
getwd() ## where has my plot gone...



### to_Tpa

title_Tpa_LRv5a_scf1 <- ggdraw() + draw_label("Tpa - Tpa_LRv5a_scf1", fontface='bold')

Tpa_MFalongX_to_Tpa_Tpa_LRv5a_scf1 <- plot_grid(
  plot_MF_alongchr_tiss(Tpa_Ad_to_Tpa_F_MF, "Tpa_LRv5a_scf1", "A"),
  plot_MF_alongchr_tiss(Tpa_Ad_to_Tpa_F_MF, "Tpa_LRv5a_scf1", "B"),
  plot_MF_alongchr_tiss(Tpa_Ad_to_Tpa_F_MF, "Tpa_LRv5a_scf1", "DG"),
  plot_MF_alongchr_tiss(Tpa_Ad_to_Tpa_F_MF, "Tpa_LRv5a_scf1", "FB"),
  plot_MF_alongchr_tiss(Tpa_Ad_to_Tpa_F_MF, "Tpa_LRv5a_scf1", "Fe"),
  plot_MF_alongchr_tiss(Tpa_Ad_to_Tpa_F_MF, "Tpa_LRv5a_scf1", "Gu"),
  plot_MF_alongchr_tiss(Tpa_Ad_to_Tpa_F_MF, "Tpa_LRv5a_scf1", "GoTe"),
  ncol = 1)

png(filename = "Tpa_MFalongX_to_Tpa_F.png", width = 10, height = 15, units = "in", bg = "white", res = 300)
plot_grid(title_Tpa_LRv5a_scf1, Tpa_MFalongX_to_Tpa_Tpa_LRv5a_scf1, ncol=1, rel_heights=c(0.1, 4))
dev.off()
getwd() ## where has my plot gone...

title_Tpa_LRv5a_scf1 <- ggdraw() + draw_label("Tge - Tpa_LRv5a_scf1", fontface='bold')

Tge_MFalongX_to_Tpa_Tpa_LRv5a_scf1 <- plot_grid(
  plot_MF_alongchr_tiss(Tge_Ad_to_Tpa_F_MF, "Tpa_LRv5a_scf1", "A"),
  plot_MF_alongchr_tiss(Tge_Ad_to_Tpa_F_MF, "Tpa_LRv5a_scf1", "B"),
  plot_MF_alongchr_tiss(Tge_Ad_to_Tpa_F_MF, "Tpa_LRv5a_scf1", "DG"),
  plot_MF_alongchr_tiss(Tge_Ad_to_Tpa_F_MF, "Tpa_LRv5a_scf1", "FB"),
  plot_MF_alongchr_tiss(Tge_Ad_to_Tpa_F_MF, "Tpa_LRv5a_scf1", "Fe"),
  plot_MF_alongchr_tiss(Tge_Ad_to_Tpa_F_MF, "Tpa_LRv5a_scf1", "Gu"),
  plot_MF_alongchr_tiss(Tge_Ad_to_Tpa_F_MF, "Tpa_LRv5a_scf1", "GoTe"),
  ncol = 1)

png(filename = "Tge_MFalongX_to_Tpa_F.png", width = 10, height = 15, units = "in", bg = "white", res = 300)
plot_grid(title_Tpa_LRv5a_scf1, Tge_MFalongX_to_Tpa_Tpa_LRv5a_scf1, ncol=1, rel_heights=c(0.1, 4))
dev.off()
getwd() ## where has my plot gone...


#plot_FPKM_alongchr(Tps_Ad_to_Tps_F_avFPKM, "Tps_LRv5b_scf3", "Te")



### Tce X scafs

Tce_HiCLG3 <- subset(Tce_Ad_to_Tce_F_MF, Tce_Ad_to_Tce_F_MF$LG == "HiCLG3")
Tce_HiCLG3 <- subset(Tce_Ad_to_Tce_F_MF, Tce_Ad_to_Tce_F_MF$LG == "HiCLG3")
sort(table(Tce_HiCLG3$scaf, useNA = "always"), decreasing = T)
sort(table(Tce_HiCLG3$scaf, useNA = "always"), decreasing = T)

(1661) / sum(sort(table(Tce_HiCLG3$scaf, useNA = "always"), decreasing = T))
(1661 + 851) / sum(sort(table(Tce_HiCLG3$scaf, useNA = "always"), decreasing = T))
(1661 + 851 + 598) / sum(sort(table(Tce_HiCLG3$scaf, useNA = "always"), decreasing = T))
(1661 + 851 + 598 + 583 ) / sum(sort(table(Tce_HiCLG3$scaf, useNA = "always"), decreasing = T))

# plot_MF_alongchr(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94")
# plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "GoTe")


title_Tce_LRv5a_scf94 <- ggdraw() + draw_label("Tce - Tce_LRv5a_scf94", fontface='bold')
title_Tce_LRv5a_scf101 <- ggdraw() + draw_label("Tce - Tce_LRv5a_scf101", fontface='bold')
title_Tce_LRv5a_scf96 <- ggdraw() + draw_label("Tce - Tce_LRv5a_scf96", fontface='bold')
title_Tce_LRv5a_scf93 <- ggdraw() + draw_label("Tce - Tce_LRv5a_scf93", fontface='bold')

Tce_MFalongX_to_Tce_Tce_LRv5a_scf94 <- plot_grid(
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "A"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "B"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "DG"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "FB"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "Fe"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "Gu"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "GoTe"),
  ncol = 1)


Tce_P_Tce_LRv5a_scf94  <- plot_grid(title_Tce_LRv5a_scf94, Tce_MFalongX_to_Tce_Tce_LRv5a_scf94, ncol=1, rel_heights=c(0.1, 4))


Tce_MFalongX_to_Tce_Tce_LRv5a_scf101 <- plot_grid(
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf101", "A"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf101", "B"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf101", "DG"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf101", "FB"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf101", "Fe"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf101", "Gu"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf101", "GoTe"),
  ncol = 1)

Tce_P_Tce_LRv5a_scf101 <- plot_grid(title_Tce_LRv5a_scf101, Tce_MFalongX_to_Tce_Tce_LRv5a_scf101, ncol=1, rel_heights=c(0.1, 4))


Tce_MFalongX_to_Tce_Tce_LRv5a_scf96 <- plot_grid(
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf96", "A"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf96", "B"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf96", "DG"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf96", "FB"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf96", "Fe"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf96", "Gu"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf96", "GoTe"),
  ncol = 1)

Tce_P_Tce_LRv5a_scf96 <- plot_grid(title_Tce_LRv5a_scf96, Tce_MFalongX_to_Tce_Tce_LRv5a_scf96, ncol=1, rel_heights=c(0.1, 4))


Tce_MFalongX_to_Tce_Tce_LRv5a_scf93 <- plot_grid(
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf93", "A"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf93", "B"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf93", "DG"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf93", "FB"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf93", "Fe"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf93", "Gu"),
  plot_MF_alongchr_tiss(Tce_Ad_to_Tce_F_MF, "Tce_LRv5a_scf93", "GoTe"),
  ncol = 1)

Tce_P_Tce_LRv5a_scf93 <- plot_grid(title_Tce_LRv5a_scf93, Tce_MFalongX_to_Tce_Tce_LRv5a_scf93, ncol=1, rel_heights=c(0.1, 4))


png(filename = "Tce_MFalongX_to_Tce_F.png", width = 10, height = 15, units = "in", bg = "white", res = 300)
plot_grid(Tce_P_Tce_LRv5a_scf94, Tce_P_Tce_LRv5a_scf101,  Tce_P_Tce_LRv5a_scf96, ncol=3, rel_widths =c(2, 1, 1))
dev.off()
getwd() ## where has my plot gone...



title_Tce_LRv5a_scf94 <- ggdraw() + draw_label("Tms - Tce_LRv5a_scf94", fontface='bold')
title_Tce_LRv5a_scf101 <- ggdraw() + draw_label("Tms - Tce_LRv5a_scf101", fontface='bold')
title_Tce_LRv5a_scf96 <- ggdraw() + draw_label("Tms - Tce_LRv5a_scf96", fontface='bold')
title_Tce_LRv5a_scf93 <- ggdraw() + draw_label("Tms - Tce_LRv5a_scf93", fontface='bold')

Tms_MFalongX_to_Tce_Tce_LRv5a_scf94 <- plot_grid(
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "A"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "B"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "DG"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "FB"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "Fe"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "Gu"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf94", "GoTe"),
  ncol = 1)
Tms_P_Tce_LRv5a_scf94  <- plot_grid(title_Tce_LRv5a_scf94, Tms_MFalongX_to_Tce_Tce_LRv5a_scf94, ncol=1, rel_heights=c(0.1, 4))


Tms_MFalongX_to_Tce_Tce_LRv5a_scf101 <- plot_grid(
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf101", "A"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf101", "B"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf101", "DG"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf101", "FB"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf101", "Fe"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf101", "Gu"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf101", "GoTe"),
  ncol = 1)
Tms_P_Tce_LRv5a_scf101 <- plot_grid(title_Tce_LRv5a_scf101, Tms_MFalongX_to_Tce_Tce_LRv5a_scf101, ncol=1, rel_heights=c(0.1, 4))




Tms_MFalongX_to_Tce_Tce_LRv5a_scf96 <- plot_grid(
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf96", "A"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf96", "B"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf96", "DG"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf96", "FB"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf96", "Fe"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf96", "Gu"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf96", "GoTe"),
  ncol = 1)
Tms_P_Tce_LRv5a_scf96 <- plot_grid(title_Tce_LRv5a_scf96, Tms_MFalongX_to_Tce_Tce_LRv5a_scf96, ncol=1, rel_heights=c(0.1, 4))


Tms_MFalongX_to_Tce_Tce_LRv5a_scf93 <- plot_grid(
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf93", "A"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf93", "B"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf93", "DG"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf93", "FB"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf93", "Fe"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf93", "Gu"),
  plot_MF_alongchr_tiss(Tms_Ad_to_Tce_F_MF, "Tce_LRv5a_scf93", "GoTe"),
  ncol = 1)
Tms_P_Tce_LRv5a_scf93 <- plot_grid(title_Tce_LRv5a_scf93, Tms_MFalongX_to_Tce_Tce_LRv5a_scf93, ncol=1, rel_heights=c(0.1, 4))



# png(filename = "Tms_MFalongX_to_Tce_F.png", width = 10, height = 15, units = "in", bg = "white", res = 300)
# plot_grid(Tms_P_Tce_LRv5a_scf94, Tms_P_Tce_LRv5a_scf101,  Tms_P_Tce_LRv5a_scf96, Tms_P_Tce_LRv5a_scf93, ncol=4, rel_widths =c(1, 1, 1, 1))
# dev.off()
# getwd() ## where has my plot gone...


png(filename = "Tms_MFalongX_to_Tce_F.png", width = 10, height = 15, units = "in", bg = "white", res = 300)
plot_grid(Tms_P_Tce_LRv5a_scf94, Tms_P_Tce_LRv5a_scf101,  Tms_P_Tce_LRv5a_scf96, ncol=3, rel_widths =c(2, 1, 1))
dev.off()
getwd() ## where has my plot gone...














#######################################################################################################################3
###### genes with higher expression on the X in male testes.

#######################
TpaTge_Ad_to_Tpa_F_MF_X_for_heat$Tpa_GoTe_MM  <- ifelse(TpaTge_Ad_to_Tpa_F_MF_X_for_heat$Tpa_GoTe > 0, 1, 0) 
TpaTge_Ad_to_Tpa_F_MF_X_for_heat$Tge_GoTe_MM  <- ifelse(TpaTge_Ad_to_Tpa_F_MF_X_for_heat$Tge_GoTe > 0, 1, 0) 
TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM <- subset(TpaTge_Ad_to_Tpa_F_MF_X_for_heat, TpaTge_Ad_to_Tpa_F_MF_X_for_heat$Tpa_GoTe_MM == 1 |  TpaTge_Ad_to_Tpa_F_MF_X_for_heat$Tge_GoTe_MM == 1)
TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM <- TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM[1:(length(TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM[1,])-2)]
TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM_bothsp <- subset(TpaTge_Ad_to_Tpa_F_MF_X_for_heat, TpaTge_Ad_to_Tpa_F_MF_X_for_heat$Tpa_GoTe_MM == 1 & TpaTge_Ad_to_Tpa_F_MF_X_for_heat$Tge_GoTe_MM == 1)
TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM_bothsp <- TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM_bothsp[1:(length(TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM_bothsp[1,])-2)]

write.csv(TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM, "TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM.csv")

png(filename = "TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM),        cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D","#67000D"))(length(breaksList)), breaks = breaksList, show_colnames = T,border_color = NA, legend = F, main = "TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM")
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM_bothsp.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM_bothsp), cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D","#67000D"))(length(breaksList)), breaks = breaksList, show_colnames = T,border_color = NA, legend = F, main = "TpaTge_Ad_to_Tpa_F_MF_X_for_heat_MM_bothsp")
dev.off()
getwd() ## where has my plot gone...


##########
TpsTdi_Ad_to_Tps_F_MF_X_for_heat$Tps_GoTe_MM  <- ifelse(TpsTdi_Ad_to_Tps_F_MF_X_for_heat$Tps_GoTe > 0, 1, 0) 
TpsTdi_Ad_to_Tps_F_MF_X_for_heat$Tdi_GoTe_MM  <- ifelse(TpsTdi_Ad_to_Tps_F_MF_X_for_heat$Tdi_GoTe > 0, 1, 0) 
TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM           <- subset(TpsTdi_Ad_to_Tps_F_MF_X_for_heat, TpsTdi_Ad_to_Tps_F_MF_X_for_heat$Tps_GoTe_MM == 1 |  TpsTdi_Ad_to_Tps_F_MF_X_for_heat$Tdi_GoTe_MM == 1)
TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM           <- TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM[1:(length(TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM[1,])-2)]
TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM_bothsp    <- subset(TpsTdi_Ad_to_Tps_F_MF_X_for_heat, TpsTdi_Ad_to_Tps_F_MF_X_for_heat$Tps_GoTe_MM == 1 & TpsTdi_Ad_to_Tps_F_MF_X_for_heat$Tdi_GoTe_MM == 1)
TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM_bothsp    <- TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM_bothsp[1:(length(TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM_bothsp[1,])-2)]

write.csv(TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM, "TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM.csv")

png(filename = "TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM),        cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D","#67000D"))(length(breaksList)), breaks = breaksList, show_colnames = T,border_color = NA, legend = F, main = "TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM")
dev.off()
getwd() ## where has my plot gone...

png(filename = "TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM_bothsp.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM_bothsp), cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D","#67000D"))(length(breaksList)), breaks = breaksList, show_colnames = T,border_color = NA, legend = F, main = "TpsTdi_Ad_to_Tps_F_MF_X_for_heat_MM_bothsp")
dev.off()
getwd() ## where has my plot gone...


#######################
TceTms_Ad_to_Tce_F_MF_X_for_heat$Tce_GoTe_MM  <- ifelse(TceTms_Ad_to_Tce_F_MF_X_for_heat$Tce_GoTe > 0, 1, 0) 
TceTms_Ad_to_Tce_F_MF_X_for_heat$Tms_GoTe_MM  <- ifelse(TceTms_Ad_to_Tce_F_MF_X_for_heat$Tms_GoTe > 0, 1, 0) 
TceTms_Ad_to_Tce_F_MF_X_for_heat_MM <- subset(TceTms_Ad_to_Tce_F_MF_X_for_heat, TceTms_Ad_to_Tce_F_MF_X_for_heat$Tce_GoTe_MM == 1 |  TceTms_Ad_to_Tce_F_MF_X_for_heat$Tms_GoTe_MM == 1)
TceTms_Ad_to_Tce_F_MF_X_for_heat_MM <- TceTms_Ad_to_Tce_F_MF_X_for_heat_MM[1:(length(TceTms_Ad_to_Tce_F_MF_X_for_heat_MM[1,])-2)]
TceTms_Ad_to_Tce_F_MF_X_for_heat_MM_bothsp <- subset(TceTms_Ad_to_Tce_F_MF_X_for_heat, TceTms_Ad_to_Tce_F_MF_X_for_heat$Tce_GoTe_MM == 1 & TceTms_Ad_to_Tce_F_MF_X_for_heat$Tms_GoTe_MM == 1)
TceTms_Ad_to_Tce_F_MF_X_for_heat_MM_bothsp <- TceTms_Ad_to_Tce_F_MF_X_for_heat_MM_bothsp[1:(length(TceTms_Ad_to_Tce_F_MF_X_for_heat_MM_bothsp[1,])-2)]

write.csv(TceTms_Ad_to_Tce_F_MF_X_for_heat_MM, "TceTms_Ad_to_Tce_F_MF_X_for_heat_MM.csv")

png(filename = "TceTms_Ad_to_Tce_F_MF_X_for_heat_MM.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TceTms_Ad_to_Tce_F_MF_X_for_heat_MM),        cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D","#67000D"))(length(breaksList)), breaks = breaksList, show_colnames = T,border_color = NA, legend = F, main = "TceTms_Ad_to_Tce_F_MF_X_for_heat_MM")
dev.off()
getwd() ## where has my plot gone...

png(filename = "TceTms_Ad_to_Tce_F_MF_X_for_heat_MM_bothsp.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
pheatmap(t(TceTms_Ad_to_Tce_F_MF_X_for_heat_MM_bothsp), cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D","#67000D"))(length(breaksList)), breaks = breaksList, show_colnames = T,border_color = NA, legend = F, main = "TceTms_Ad_to_Tce_F_MF_X_for_heat_MM_bothsp")
dev.off()
getwd() ## where has my plot gone...



################################################################################################################
## RANKS of X A 


### take gene intersection per tiss. calc rank of X



Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}


get_exp_ranks <- function(df_sex, df_asex, want_tiss){
  
  df_sex_tiss  <- subset(df_sex,  df_sex$tiss  == want_tiss)
  df_asex_tiss <- subset(df_asex, df_asex$tiss == want_tiss)
  
  df_sex_tiss_M   <- subset(df_sex_tiss,   df_sex_tiss$sex   == "M")
  df_asex_tiss_M  <- subset(df_asex_tiss,  df_asex_tiss$sex  == "M")
  df_sex_tiss_F   <- subset(df_sex_tiss,   df_sex_tiss$sex   == "F")
  df_asex_tiss_F  <- subset(df_asex_tiss,  df_asex_tiss$sex  == "F")
  
  gene_lists <- list(df_sex_tiss_M$gene_id, df_asex_tiss_M$gene_id, df_sex_tiss_F$gene_id, df_asex_tiss_F$gene_id)
  want_genes <- Intersect(gene_lists)
  
  print(length(df_sex_tiss_M$gene_id))
  print(length(df_asex_tiss_M$gene_id))
  print(length(df_sex_tiss_F$gene_id))
  print(length(df_asex_tiss_F$gene_id))
  print(length(want_genes))
  
  ### subset by genes intersection
  
  df_sex_tiss_M_2  <- subset(df_sex_tiss_M,  gene_id %in% want_genes)
  df_asex_tiss_M_2 <- subset(df_asex_tiss_M, gene_id %in% want_genes)
  df_sex_tiss_F_2  <- subset(df_sex_tiss_F,  gene_id %in% want_genes)
  df_asex_tiss_F_2 <- subset(df_asex_tiss_F, gene_id %in% want_genes)
  
  print(length( df_sex_tiss_M_2[,1]))
  print(length( df_asex_tiss_M_2[,1]))
  print(length( df_sex_tiss_F_2[,1]))
  print(length( df_asex_tiss_F_2[,1]))
  
  #### get ranks
  
  df_sex_tiss_M_2$exp_rank  <- rank(-df_sex_tiss_M_2$meanFPKM)
  df_asex_tiss_M_2$exp_rank <- rank(-df_asex_tiss_M_2$meanFPKM)
  df_sex_tiss_F_2$exp_rank  <- rank(-df_sex_tiss_F_2$meanFPKM)
  df_asex_tiss_F_2$exp_rank <- rank(-df_asex_tiss_F_2$meanFPKM)
  
  df_sex_tiss_M_2$exp_rank_r  <- rank(-df_sex_tiss_M_2$meanFPKM) / length(want_genes) #* 10000
  df_asex_tiss_M_2$exp_rank_r <- rank(-df_asex_tiss_M_2$meanFPKM) / length(want_genes) #* 10000
  df_sex_tiss_F_2$exp_rank_r  <- rank(-df_sex_tiss_F_2$meanFPKM) / length(want_genes) #* 10000
  df_asex_tiss_F_2$exp_rank_r <- rank(-df_asex_tiss_F_2$meanFPKM) / length(want_genes) #* 10000
  
  out_df <- as.data.frame(rbind(
    df_sex_tiss_M_2,
    df_asex_tiss_M_2,
    df_sex_tiss_F_2,
    df_asex_tiss_F_2   
  ))
  
  return(out_df)
}


get_exp_ranks_GoTe <- function(df_sex, df_asex){
  
  df_sex_Te  <- subset(df_sex,  df_sex$tiss  == "Te")
  df_sex_Go  <- subset(df_sex,  df_sex$tiss  == "Go") 
  
  df_asex_Te  <- subset(df_asex,  df_asex$tiss  == "Te")
  df_asex_Go  <- subset(df_asex,  df_asex$tiss  == "Go") 
  
  gene_lists <- list(
                     df_sex_Te$gene_id,
                     df_sex_Go$gene_id,

                     df_asex_Te$gene_id,
                     df_asex_Go$gene_id)
  want_genes <- Intersect(gene_lists)
  
  

  print(length(df_sex_Te$gene_id))
  print(length(df_sex_Go$gene_id))

  print(length(df_asex_Te$gene_id))
  print(length(df_asex_Go$gene_id))
  
  
  print(length(want_genes))
  
  ### subset by genes intersection
  

  df_sex_Te_2   <- subset(df_sex_Te,  gene_id %in% want_genes)
  df_sex_Go_2   <- subset(df_sex_Go,  gene_id %in% want_genes) 
  

  df_asex_Te_2   <- subset(df_asex_Te,  gene_id %in% want_genes)
  df_asex_Go_2   <- subset(df_asex_Go,  gene_id %in% want_genes) 

  print(length(df_sex_Te_2[,1]))
  print(length(df_sex_Go_2[,1]))
 
  print(length(df_asex_Te_2[,1]))
  print(length(df_asex_Go_2[,1]))  
  

  df_sex_Te_2$exp_rank  <- rank(-df_sex_Te_2$meanFPKM)
  df_sex_Go_2$exp_rank  <- rank(-df_sex_Go_2$meanFPKM)

  df_asex_Te_2$exp_rank  <- rank(-df_asex_Te_2$meanFPKM)
  df_asex_Go_2$exp_rank  <- rank(-df_asex_Go_2$meanFPKM)  

  df_sex_Te_2$exp_rank_r  <- rank(-df_sex_Te_2$meanFPKM) / length(want_genes) #* 10000
  df_sex_Go_2$exp_rank_r  <- rank(-df_sex_Go_2$meanFPKM) / length(want_genes) # * 10000
  df_asex_Te_2$exp_rank_r  <- rank(-df_asex_Te_2$meanFPKM) / length(want_genes) # * 10000
  df_asex_Go_2$exp_rank_r  <- rank(-df_asex_Go_2$meanFPKM) / length(want_genes) #* 10000
  
  
  out_df <- as.data.frame(rbind(

    df_sex_Te_2,
    df_sex_Go_2,

    df_asex_Te_2,
    df_asex_Go_2
  ))
  
  return(out_df)
}


TpsTdi_ranks <- as.data.frame(rbind(
  get_exp_ranks(Tps_Ad_to_Tps_F_avFPKM, Tdi_Ad_to_Tps_F_avFPKM, "A"),
  get_exp_ranks(Tps_Ad_to_Tps_F_avFPKM, Tdi_Ad_to_Tps_F_avFPKM, "B"),
  get_exp_ranks(Tps_Ad_to_Tps_F_avFPKM, Tdi_Ad_to_Tps_F_avFPKM, "DG"),
  get_exp_ranks(Tps_Ad_to_Tps_F_avFPKM, Tdi_Ad_to_Tps_F_avFPKM, "FB"),
  get_exp_ranks(Tps_Ad_to_Tps_F_avFPKM, Tdi_Ad_to_Tps_F_avFPKM, "Fe"),
  get_exp_ranks(Tps_Ad_to_Tps_F_avFPKM, Tdi_Ad_to_Tps_F_avFPKM, "Gu"),
  get_exp_ranks_GoTe(Tps_Ad_to_Tps_F_avFPKM, Tdi_Ad_to_Tps_F_avFPKM)
))
TpsTdi_ranks$rep_m <- ordered(TpsTdi_ranks$rep_m, levels = c("sex", "asex"))
str(TpsTdi_ranks)
TpsTdi_ranks_X   <- subset(TpsTdi_ranks, TpsTdi_ranks$XA == "X")
TpsTdi_ranks_X_M <- subset(TpsTdi_ranks_X, TpsTdi_ranks_X$sex == "M")
TpsTdi_ranks_X_F <- subset(TpsTdi_ranks_X, TpsTdi_ranks_X$sex == "F")



TpaTge_ranks <- as.data.frame(rbind(
  get_exp_ranks(Tpa_Ad_to_Tpa_F_avFPKM, Tge_Ad_to_Tpa_F_avFPKM, "A"),
  get_exp_ranks(Tpa_Ad_to_Tpa_F_avFPKM, Tge_Ad_to_Tpa_F_avFPKM, "B"),
  get_exp_ranks(Tpa_Ad_to_Tpa_F_avFPKM, Tge_Ad_to_Tpa_F_avFPKM, "DG"),
  get_exp_ranks(Tpa_Ad_to_Tpa_F_avFPKM, Tge_Ad_to_Tpa_F_avFPKM, "FB"),
  get_exp_ranks(Tpa_Ad_to_Tpa_F_avFPKM, Tge_Ad_to_Tpa_F_avFPKM, "Fe"),
  get_exp_ranks(Tpa_Ad_to_Tpa_F_avFPKM, Tge_Ad_to_Tpa_F_avFPKM, "Gu"),
  get_exp_ranks_GoTe(Tpa_Ad_to_Tpa_F_avFPKM, Tge_Ad_to_Tpa_F_avFPKM)
))
TpaTge_ranks$rep_m <- ordered(TpaTge_ranks$rep_m, levels = c("sex", "asex"))
str(TpaTge_ranks)
TpaTge_ranks_X   <- subset(TpaTge_ranks, TpaTge_ranks$XA == "X")
TpaTge_ranks_X_M <- subset(TpaTge_ranks_X, TpaTge_ranks_X$sex == "M")
TpaTge_ranks_X_F <- subset(TpaTge_ranks_X, TpaTge_ranks_X$sex == "F")


TceTms_ranks <- as.data.frame(rbind(
  get_exp_ranks(Tce_Ad_to_Tce_F_avFPKM, Tms_Ad_to_Tce_F_avFPKM, "A"),
  get_exp_ranks(Tce_Ad_to_Tce_F_avFPKM, Tms_Ad_to_Tce_F_avFPKM, "B"),
  get_exp_ranks(Tce_Ad_to_Tce_F_avFPKM, Tms_Ad_to_Tce_F_avFPKM, "DG"),
  get_exp_ranks(Tce_Ad_to_Tce_F_avFPKM, Tms_Ad_to_Tce_F_avFPKM, "FB"),
  get_exp_ranks(Tce_Ad_to_Tce_F_avFPKM, Tms_Ad_to_Tce_F_avFPKM, "Fe"),
  get_exp_ranks(Tce_Ad_to_Tce_F_avFPKM, Tms_Ad_to_Tce_F_avFPKM, "Gu"),
  get_exp_ranks_GoTe(Tce_Ad_to_Tce_F_avFPKM, Tms_Ad_to_Tce_F_avFPKM)
))
TceTms_ranks$rep_m <- ordered(TceTms_ranks$rep_m, levels = c("sex", "asex"))
str(TceTms_ranks)
TceTms_ranks_X   <- subset(TceTms_ranks, TceTms_ranks$XA == "X")
TceTms_ranks_X_M <- subset(TceTms_ranks_X, TceTms_ranks_X$sex == "M")
TceTms_ranks_X_F <- subset(TceTms_ranks_X, TceTms_ranks_X$sex == "F")






png(filename = "TpsTdi_ranks_X.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
plot_grid(
  ggplot(TpsTdi_ranks_X_M, aes(tiss, exp_rank)) + 
    theme_bw() +
    geom_boxplot(aes(fill = factor(rep_m))) +  scale_y_reverse() +
    scale_fill_manual(values=c( "red2", "royalblue2")) +
    ggtitle("Tps males") + theme(legend.position = "none"),
  
  ggplot(TpsTdi_ranks_X_F, aes(tiss, exp_rank)) + 
    theme_bw() +
    geom_boxplot(aes(fill = factor(rep_m))) +  scale_y_reverse() +
    scale_fill_manual(values=c( "red2", "royalblue2")) +
    ggtitle("Tps females") + theme(legend.position = "none")
)
dev.off()
getwd() ## where has my plot gone...


png(filename = "TpsTdi_rel_ranks_X.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
plot_grid(
  ggplot(TpsTdi_ranks_X_M, aes(tiss, exp_rank_r)) + 
    theme_bw() +
    geom_boxplot(aes(fill = factor(rep_m))) +  scale_y_reverse() +
    scale_fill_manual(values=c( "red2", "royalblue2")) +
    ggtitle("Tps males") + theme(legend.position = "none"),
  
  ggplot(TpsTdi_ranks_X_F, aes(tiss, exp_rank_r)) + 
    theme_bw() +
    geom_boxplot(aes(fill = factor(rep_m))) +  scale_y_reverse() +
    scale_fill_manual(values=c( "red2", "royalblue2")) +
    ggtitle("Tps females") + theme(legend.position = "none")
)
dev.off()
getwd() ## where has my plot gone...



png(filename = "TpaTge_ranks_X.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
plot_grid(
  ggplot(TpaTge_ranks_X_M, aes(tiss, exp_rank)) + 
    theme_bw() +
    geom_boxplot(aes(fill = factor(rep_m))) +  scale_y_reverse() +
    scale_fill_manual(values=c( "red2", "royalblue2")) +
    ggtitle("Tpa males") + theme(legend.position = "none"),
  
  ggplot(TpaTge_ranks_X_F, aes(tiss, exp_rank)) + 
    theme_bw() +
    geom_boxplot(aes(fill = factor(rep_m))) +  scale_y_reverse() +
    scale_fill_manual(values=c( "red2", "royalblue2")) +
    ggtitle("Tpa females") + theme(legend.position = "none")
)
dev.off()
getwd() ## where has my plot gone...


png(filename = "TpaTge_rel_ranks_X.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
plot_grid(
  ggplot(TpaTge_ranks_X_M, aes(tiss, exp_rank_r)) + 
    theme_bw() +
    geom_boxplot(aes(fill = factor(rep_m))) +  scale_y_reverse() +
    scale_fill_manual(values=c( "red2", "royalblue2")) +
    ggtitle("Tpa males") + theme(legend.position = "none"),
  
  ggplot(TpaTge_ranks_X_F, aes(tiss, exp_rank_r)) + 
    theme_bw() +
    geom_boxplot(aes(fill = factor(rep_m))) +  scale_y_reverse() +
    scale_fill_manual(values=c( "red2", "royalblue2")) +
    ggtitle("Tpa females") + theme(legend.position = "none")
)
dev.off()
getwd() ## where has my plot gone...


png(filename = "TceTms_ranks_X.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
plot_grid(
  ggplot(TceTms_ranks_X_M, aes(tiss, exp_rank)) + 
    theme_bw() +
    geom_boxplot(aes(fill = factor(rep_m))) +  scale_y_reverse() +
    scale_fill_manual(values=c( "red2", "royalblue2")) +
    ggtitle("Tce males") + theme(legend.position = "none"),
  
  ggplot(TceTms_ranks_X_F, aes(tiss, exp_rank)) + 
    theme_bw() +
    geom_boxplot(aes(fill = factor(rep_m))) +  scale_y_reverse() +
    scale_fill_manual(values=c( "red2", "royalblue2")) +
    ggtitle("Tce females") + theme(legend.position = "none")
)
dev.off()
getwd() ## where has my plot gone...

png(filename = "TceTms_rel_ranks_X.png", width = 10, height = 5, units = "in", bg = "white", res = 300)
plot_grid(
  ggplot(TceTms_ranks_X_M, aes(tiss, exp_rank_r)) + 
    theme_bw() +
    geom_boxplot(aes(fill = factor(rep_m))) +  scale_y_reverse() +
    scale_fill_manual(values=c( "red2", "royalblue2")) +
    ggtitle("Tce males") + theme(legend.position = "none"),
  
  ggplot(TceTms_ranks_X_F, aes(tiss, exp_rank_r)) + 
    theme_bw() +
    geom_boxplot(aes(fill = factor(rep_m))) +  scale_y_reverse() +
    scale_fill_manual(values=c( "red2", "royalblue2")) +
    ggtitle("Tce females") + theme(legend.position = "none")
)
dev.off()
getwd() ## where has my plot gone...












