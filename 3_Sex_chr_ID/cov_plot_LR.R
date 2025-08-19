#### plot coverage for species mapped to sexual refs 

library(ggplot2)
library(cowplot)
library(hash)
library(stringr)
library(modeest)
library(gtools)

##### CODE
### merges df by the first column
merge_df_by_col <- function(match_vec,match_rep_name){
  
  ### so merge requires that the columns NOT to match on have different names or it tries to match on those too
  ### and even if I try to specify by.y and by.x it still doesn't work
  ### so I will bang a prefix on all but the first (matching) column
  ### add first df and replace names
  
  A1 = match_vec[1]
  A1 = eval(parse(text=A1))
  match_coln <- colnames(A1)[1]
  
  new_col_n <- c()
  A1_col <- colnames(A1)
  for(el in A1_col){
    match_pre = match_rep_name[1]
    new_names <- paste(match_pre,el,sep = "")
    new_col_n = c(new_col_n,new_names )
  }
  
  colnames(A1) <- new_col_n
  colnames(A1)[1] <- match_coln
  df_merged <- A1
  
  for(i in 1:(length(match_vec))){
    curr_df_name <- match_vec[i]
    curr_df <- eval(parse(text=curr_df_name))
    match_coln <- colnames(curr_df)[1]
    curr_df_col <- colnames(curr_df)
    new_col_n <- c()
    
    for(el in curr_df_col){
      match_pre = match_rep_name[i]
      new_names <- paste(match_pre,el,sep = "")
      new_col_n = c(new_col_n,new_names )
      #print(new_names)
    }
    colnames(curr_df) <- new_col_n
    colnames(curr_df)[1] <- match_coln
    df_merged = merge(df_merged,curr_df)
  }
  return(df_merged)
}


plot_cov_along_chr <- function(df, chr_want){
  df_chr <- subset(df, df$scaf_name == chr_want)
  print(head(df_chr))
  
  p1 <- ggplot(df_chr, aes(x=mid_pos , cov)) +
    theme_bw() +
    geom_line() +
    ggtitle(paste("scaf", chr_want)) +
    xlab("Position")
  
  return(p1)	
}

peakfinder <- function(d){
  dh <- hist(d,plot=FALSE, breaks=1000)
  ins <- dh[["counts"]]
  nbins <- length(ins)
  ss <- which(rank(ins)%in%seq(from=nbins,to=nbins)) ## pick the top 3 intensities
  dh[["mids"]][ss]
}

## data Tdi 
### Notes
### good genome. X == Tps_LRv5b_scf3
### all males look right
### 2N = 24 (Tps)

setwd("/Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/angsD_LR_sliding_window")
Tdi_F_ReSeq_Di02         <- read.table("Tdi_F_ReSeq_Di02_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tdi_F_ReSeq_Di02$mid_pos <-  (Tdi_F_ReSeq_Di02$end + Tdi_F_ReSeq_Di02$start) / 2
Tdi_F_ReSeq_Di04         <- read.table("Tdi_F_ReSeq_Di04_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tdi_F_ReSeq_Di04$mid_pos <-  (Tdi_F_ReSeq_Di04$end + Tdi_F_ReSeq_Di04$start) / 2
Tdi_F_ReSeq_Di06         <- read.table("Tdi_F_ReSeq_Di06_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tdi_F_ReSeq_Di06$mid_pos <-  (Tdi_F_ReSeq_Di06$end + Tdi_F_ReSeq_Di06$start) / 2
Tdi_F_ReSeq_Di08         <- read.table("Tdi_F_ReSeq_Di08_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tdi_F_ReSeq_Di08$mid_pos <-  (Tdi_F_ReSeq_Di08$end + Tdi_F_ReSeq_Di08$start) / 2
Tdi_F_ReSeq_Di10         <- read.table("Tdi_F_ReSeq_Di10_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tdi_F_ReSeq_Di10$mid_pos <-  (Tdi_F_ReSeq_Di10$end + Tdi_F_ReSeq_Di10$start) / 2

Tdi_M_18_3997         <- read.table("Tdi_M_18-3997_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tdi_M_18_3997$mid_pos <-  (Tdi_M_18_3997$end + Tdi_M_18_3997$start) / 2
Tdi_M_18_3998         <- read.table("Tdi_M_18-3998_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tdi_M_18_3998$mid_pos <-  (Tdi_M_18_3998$end + Tdi_M_18_3998$start) / 2

### add some info
Tdi_F_ReSeq_Di02$samp_name <- rep("Di02", length(Tdi_F_ReSeq_Di02[,1]))
Tdi_F_ReSeq_Di04$samp_name <- rep("Di04", length(Tdi_F_ReSeq_Di04[,1]))
Tdi_F_ReSeq_Di06$samp_name <- rep("Di06", length(Tdi_F_ReSeq_Di06[,1]))
Tdi_F_ReSeq_Di08$samp_name <- rep("Di08", length(Tdi_F_ReSeq_Di08[,1]))
Tdi_F_ReSeq_Di10$samp_name <- rep("Di10", length(Tdi_F_ReSeq_Di10[,1]))
Tdi_M_18_3997$samp_name    <- rep("18_3997", length(Tdi_M_18_3997[,1]))
Tdi_M_18_3998$samp_name    <- rep("18_3998", length(Tdi_M_18_3998[,1]))

Tdi_F_ReSeq_Di02$sex <- rep("F", length(Tdi_F_ReSeq_Di02[,1]))
Tdi_F_ReSeq_Di04$sex <- rep("F", length(Tdi_F_ReSeq_Di04[,1]))
Tdi_F_ReSeq_Di06$sex <- rep("F", length(Tdi_F_ReSeq_Di06[,1]))
Tdi_F_ReSeq_Di08$sex <- rep("F", length(Tdi_F_ReSeq_Di08[,1]))
Tdi_F_ReSeq_Di10$sex <- rep("F", length(Tdi_F_ReSeq_Di10[,1]))
Tdi_M_18_3997$sex    <- rep("M", length(Tdi_M_18_3997[,1]))
Tdi_M_18_3998$sex    <- rep("M", length(Tdi_M_18_3998[,1]))


#############################################################################################
### chr lens

Tdi_chr_lens <- aggregate(end~scaf_name, Tdi_F_ReSeq_Di02, FUN=max)
head(Tdi_chr_lens[order(Tdi_chr_lens$end,decreasing=T),], n = 20)

### big drop after 12 (as expected)


plot_all_chr_Tdi <- function(df){
  out_plota <- plot_grid(
    plot_cov_along_chr(df, "Tps_LRv5b_scf1"),
    plot_cov_along_chr(df, "Tps_LRv5b_scf2"),
    plot_cov_along_chr(df, "Tps_LRv5b_scf3"),
    plot_cov_along_chr(df, "Tps_LRv5b_scf4"),
    plot_cov_along_chr(df, "Tps_LRv5b_scf5"),
    plot_cov_along_chr(df, "Tps_LRv5b_scf6"),
    plot_cov_along_chr(df, "Tps_LRv5b_scf7"),
    plot_cov_along_chr(df, "Tps_LRv5b_scf8"),
    plot_cov_along_chr(df, "Tps_LRv5b_scf9"),
    plot_cov_along_chr(df, "Tps_LRv5b_scf10"),
    plot_cov_along_chr(df, "Tps_LRv5b_scf11"),
    plot_cov_along_chr(df, "Tps_LRv5b_scf12"), ncol = 1)
  
  return(out_plota)
  
}

png(filename = "Tdi_F_ReSeq_Di02_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 20, units = "in", bg = "white", res = 300)
plot_all_chr_Tdi(Tdi_F_ReSeq_Di02)
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tdi_F_ReSeq_Di04_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 20, units = "in", bg = "white", res = 300)
plot_all_chr_Tdi(Tdi_F_ReSeq_Di04)
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tdi_F_ReSeq_Di06_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 20, units = "in", bg = "white", res = 300)
plot_all_chr_Tdi(Tdi_F_ReSeq_Di06)
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tdi_F_ReSeq_Di08_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 20, units = "in", bg = "white", res = 300)
plot_all_chr_Tdi(Tdi_F_ReSeq_Di08)
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tdi_F_ReSeq_Di10_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 20, units = "in", bg = "white", res = 300)
plot_all_chr_Tdi(Tdi_F_ReSeq_Di10)
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tdi_M_18_3997_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 20, units = "in", bg = "white", res = 300)
plot_all_chr_Tdi(Tdi_M_18_3997 )
dev.off()
getwd() ## where has my plot gone....

png(filename = "Tdi_M_18_3998_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 20, height = 20, units = "in", bg = "white", res = 300)
plot_all_chr_Tdi(Tdi_M_18_3998 )
dev.off()
getwd() ## where has my plot gone....



##### median cov window

want_scafs_Tdi <- c(
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




Tdi_F_ReSeq_Di02_ws     <- Tdi_F_ReSeq_Di02[Tdi_F_ReSeq_Di02$scaf_name %in% want_scafs_Tdi,]
Tdi_F_ReSeq_Di02_ws_med <- aggregate(Tdi_F_ReSeq_Di02_ws$cov, list(Tdi_F_ReSeq_Di02_ws$scaf_name), FUN=median)
colnames(Tdi_F_ReSeq_Di02_ws_med) <- c("scaf", "Tdi_F_ReSeq_Di02_med_window_cov")
Tdi_F_ReSeq_Di04_ws     <- Tdi_F_ReSeq_Di04[Tdi_F_ReSeq_Di04$scaf_name %in% want_scafs_Tdi,]
Tdi_F_ReSeq_Di04_ws_med <- aggregate(Tdi_F_ReSeq_Di04_ws$cov, list(Tdi_F_ReSeq_Di04_ws$scaf_name), FUN=median)
colnames(Tdi_F_ReSeq_Di04_ws_med) <- c("scaf", "Tdi_F_ReSeq_Di04_med_window_cov")
Tdi_F_ReSeq_Di06_ws     <- Tdi_F_ReSeq_Di06[Tdi_F_ReSeq_Di06$scaf_name %in% want_scafs_Tdi,]
Tdi_F_ReSeq_Di06_ws_med <- aggregate(Tdi_F_ReSeq_Di06_ws$cov, list(Tdi_F_ReSeq_Di06_ws$scaf_name), FUN=median)
colnames(Tdi_F_ReSeq_Di06_ws_med) <- c("scaf", "Tdi_F_ReSeq_Di06_med_window_cov")
Tdi_F_ReSeq_Di08_ws     <- Tdi_F_ReSeq_Di08[Tdi_F_ReSeq_Di08$scaf_name %in% want_scafs_Tdi,]
Tdi_F_ReSeq_Di08_ws_med <- aggregate(Tdi_F_ReSeq_Di08_ws$cov, list(Tdi_F_ReSeq_Di08_ws$scaf_name), FUN=median)
colnames(Tdi_F_ReSeq_Di08_ws_med) <- c("scaf", "Tdi_F_ReSeq_Di08_med_window_cov")
Tdi_F_ReSeq_Di10_ws     <- Tdi_F_ReSeq_Di10[Tdi_F_ReSeq_Di10$scaf_name %in% want_scafs_Tdi,]
Tdi_F_ReSeq_Di10_ws_med <- aggregate(Tdi_F_ReSeq_Di10_ws$cov, list(Tdi_F_ReSeq_Di10_ws$scaf_name), FUN=median)
colnames(Tdi_F_ReSeq_Di10_ws_med) <- c("scaf", "Tdi_F_ReSeq_Di10_med_window_cov")
Tdi_M_18_3997_ws     <- Tdi_M_18_3997[Tdi_M_18_3997$scaf_name %in% want_scafs_Tdi,]
Tdi_M_18_3997_ws_med <- aggregate(Tdi_M_18_3997_ws$cov, list(Tdi_M_18_3997_ws$scaf_name), FUN=median)
colnames(Tdi_M_18_3997_ws_med) <- c("scaf", "Tdi_M_18_3997_med_window_cov")
Tdi_M_18_3998_ws     <- Tdi_M_18_3998[Tdi_M_18_3998$scaf_name %in% want_scafs_Tdi,]
Tdi_M_18_3998_ws_med <- aggregate(Tdi_M_18_3998_ws$cov, list(Tdi_M_18_3998_ws$scaf_name), FUN=median)
colnames(Tdi_M_18_3998_ws_med) <- c("scaf", "Tdi_M_18_3998_med_window_cov")

match_dfs <- c("Tdi_F_ReSeq_Di02_ws_med", "Tdi_F_ReSeq_Di04_ws_med", "Tdi_F_ReSeq_Di06_ws_med", "Tdi_F_ReSeq_Di08_ws_med", "Tdi_F_ReSeq_Di10_ws_med", "Tdi_M_18_3997_ws_med", "Tdi_M_18_3998_ws_med")
match_name <- c("","", "", "", "", "", "")
Tdi_LRv5a_med_window_cov <- merge_df_by_col(match_dfs,match_name)
write.csv(Tdi_LRv5a_med_window_cov, "Tdi_LRv5a_med_window_cov_to_Tps.csv")


head(Tdi_F_ReSeq_Di10)

###############################################################################
### male vs female cov

#############################################################################
### get norm cov

norm_cov <- function(df){
  modal_cov <- mlv(df$cov, method = "shorth")
  df$cov_n <- df$cov / modal_cov
  return(df)
}

Tdi_F_ReSeq_Di02 <- norm_cov(Tdi_F_ReSeq_Di02)
Tdi_F_ReSeq_Di04 <- norm_cov(Tdi_F_ReSeq_Di04)
Tdi_F_ReSeq_Di06 <- norm_cov(Tdi_F_ReSeq_Di06)
Tdi_F_ReSeq_Di08 <- norm_cov(Tdi_F_ReSeq_Di08)
Tdi_F_ReSeq_Di10 <- norm_cov(Tdi_F_ReSeq_Di10)
Tdi_M_18_3997    <- norm_cov(Tdi_M_18_3997)
Tdi_M_18_3998    <- norm_cov(Tdi_M_18_3998)
Tdi_all_norm     <- rbind(Tdi_F_ReSeq_Di02, Tdi_F_ReSeq_Di04, Tdi_F_ReSeq_Di06, Tdi_F_ReSeq_Di08, Tdi_F_ReSeq_Di10, Tdi_F_ReSeq_Di02, Tdi_M_18_3997, Tdi_M_18_3998)

head(Tdi_all_norm)
tail(Tdi_all_norm)



plot_MF_cov_all <- function(df, chr_want){
  df_chr <- subset(df, df$scaf_name == chr_want)
  print(head(df_chr))
  p2 <- ggplot(df_chr, aes(x=mid_pos, cov_n, col = samp_name)) +
    theme_bw() +
    geom_line(alpha=0.5) +
    ggtitle(paste(chr_want)) + 
    #scale_colour_manual(values=c("royalblue2", "red3")) +
    xlab("Position") + ylim(0,3)
    
  p1 <- p2 + theme(legend.position="none")
  
  outlist = list("p2" = p2, "p1" = p1)
  return(outlist)	
}


plot_MF_cov_all_smooth <- function(df, chr_want){
  df_chr <- subset(df, df$scaf_name == chr_want)
  print(head(df_chr))
  p1 <- ggplot(df_chr, aes(x=mid_pos, cov_n, col = samp_name)) +
    theme_bw() +
    geom_smooth(method="auto", se=FALSE, fullrange=FALSE, level=0.95) +
    ggtitle(paste("scaf", chr_want)) + 
    #scale_colour_manual(values=c("royalblue2", "red3")) +
    xlab("Position")
  
  p2 <- p1 + theme(legend.position="none")
  
  outlist = list("p2" = p2, "p1" = p1)
  return(outlist)	
  
}


plot_MF_cov_all_Tce <- function(df, chr_want){
  df_chr <- subset(df, df$LG == chr_want)
  
  df_chr_breaks = subset(df_chr, is.na(df_chr$cov))
  print(head(df_chr))
  
  print(head(df_chr_breaks))
  p2 <- ggplot(df_chr, aes(x=LG_mid_pos, cov_n, col = samp_name)) +
    theme_bw() +
    geom_line(alpha=0.5, size = 0.5) +
    ggtitle(paste(chr_want)) + 
    #scale_colour_manual(values=c("royalblue2", "red3")) +
    xlab("Position") + ylim(0,3)  + geom_vline(xintercept = df_chr_breaks$LG_mid_pos)

  p1 <- p2 + theme(legend.position="none")
  
  outlist = list("p2" = p2, "p1" = p1)
  return(outlist)	
}



Tdi_window_legend_sep <- cowplot::get_legend(  plot_MF_cov_all(Tdi_all_norm, "Tps_LRv5b_scf1")$p2)

Tdi_window_p1 <- plot_grid(
  plot_MF_cov_all(Tdi_all_norm, "Tps_LRv5b_scf1")$p1,
  plot_MF_cov_all(Tdi_all_norm, "Tps_LRv5b_scf2")$p1,
  plot_MF_cov_all(Tdi_all_norm, "Tps_LRv5b_scf3")$p1, 
  plot_MF_cov_all(Tdi_all_norm, "Tps_LRv5b_scf4")$p1, 
  plot_MF_cov_all(Tdi_all_norm, "Tps_LRv5b_scf5")$p1, 
  plot_MF_cov_all(Tdi_all_norm, "Tps_LRv5b_scf6")$p1, 
  plot_MF_cov_all(Tdi_all_norm, "Tps_LRv5b_scf7")$p1, 
  plot_MF_cov_all(Tdi_all_norm, "Tps_LRv5b_scf8")$p1, 
  plot_MF_cov_all(Tdi_all_norm, "Tps_LRv5b_scf9")$p1, 
  plot_MF_cov_all(Tdi_all_norm, "Tps_LRv5b_scf10")$p1, 
  plot_MF_cov_all(Tdi_all_norm, "Tps_LRv5b_scf11")$p1, 
  plot_MF_cov_all(Tdi_all_norm, "Tps_LRv5b_scf12")$p1, 
  ncol = 1)


png(filename = "All_Tdi_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 10, height = 20, units = "in", bg = "white", res = 300)
plot_grid(Tdi_window_p1, Tdi_window_legend_sep, rel_widths  = c(1, 0.2))
dev.off()
getwd() ## where has my plot gone....


### scaf cov window plot

cov_plot_hist_Tdi <- function(df, max_cov) {
  df1 <- df
  samp = deparse(substitute(df))
  print(samp)
  p1 <- ggplot(df1) +
    theme_bw() +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf1"),     aes(cov, color="Tps_LRv5b_scf1"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf2"),     aes(cov, color="Tps_LRv5b_scf2"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf3"),     aes(cov, color="Tps_LRv5b_scf3"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf4"),     aes(cov, color="Tps_LRv5b_scf4"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf5"),     aes(cov, color="Tps_LRv5b_scf5"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf6"),     aes(cov, color="Tps_LRv5b_scf6"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf7"),     aes(cov, color="Tps_LRv5b_scf7"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf8"),     aes(cov, color="Tps_LRv5b_scf8"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf9"),     aes(cov, color="Tps_LRv5b_scf9"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf10"),    aes(cov, color="Tps_LRv5b_scf10"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf11"),    aes(cov, color="Tps_LRv5b_scf11"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf12"),    aes(cov, color="Tps_LRv5b_scf12"), size = 1, stat="density") +
    scale_color_manual(name = "Scaf", values = c("Tps_LRv5b_scf1" = "#1B9E77",
                                                 "Tps_LRv5b_scf2" = "#D95F02",
                                                 "Tps_LRv5b_scf3" = "red3",
                                                 "Tps_LRv5b_scf4" = "#A6761D",
                                                 "Tps_LRv5b_scf5" = "yellow2",
                                                 "Tps_LRv5b_scf6" = "#666666",
                                                 "Tps_LRv5b_scf7" = "lightblue",
                                                 "Tps_LRv5b_scf8" = "royalblue2",
                                                 "Tps_LRv5b_scf9" = "darkorchid",
                                                 "Tps_LRv5b_scf10" = "#7570B3",
                                                 "Tps_LRv5b_scf11" = "black",
                                                 "Tps_LRv5b_scf12" = "grey33")) + 
    
    xlim(c(0, max_cov)) + 
    ggtitle(samp) +
    xlab("Coverage") 
  
  p2 <- p1 + theme(legend.position="none")
  
  outlist = list("p2" = p2, "p1" = p1)
  return(outlist)	
}


Tdi_legend_sep <- cowplot::get_legend(cov_plot_hist_Tdi(Tdi_M_18_3997, 80)$p1)

Tdi_cov_plot <- plot_grid(
  cov_plot_hist_Tdi(Tdi_F_ReSeq_Di02, 30)$p2,
  cov_plot_hist_Tdi(Tdi_F_ReSeq_Di04, 30)$p2,
  cov_plot_hist_Tdi(Tdi_F_ReSeq_Di06, 30)$p2,
  cov_plot_hist_Tdi(Tdi_F_ReSeq_Di08, 30)$p2,
  cov_plot_hist_Tdi(Tdi_F_ReSeq_Di10, 30)$p2,
  cov_plot_hist_Tdi(Tdi_M_18_3997, 80)$p2,
  cov_plot_hist_Tdi(Tdi_M_18_3998, 80)$p2,
  ncol = 2)


png(filename = "Tdi_cov_plot_to_Tps.png", width  = 12, height = 15, units = "in", bg = "white", res = 300)
plot_grid(Tdi_cov_plot, Tdi_legend_sep, rel_widths  = c(1, 0.5))
dev.off()
getwd() ## where has my plot gone....?



################################################################################################################
#### class windows X A


Tdi_all_norm_c <- as.data.frame(cbind(Tdi_F_ReSeq_Di02$scaf_name, Tdi_F_ReSeq_Di02$start, Tdi_F_ReSeq_Di02$end, Tdi_F_ReSeq_Di02$mid_pos, 
                                      Tdi_F_ReSeq_Di02$cov_n, Tdi_F_ReSeq_Di04$cov_n, Tdi_F_ReSeq_Di06$cov_n, Tdi_F_ReSeq_Di08$cov_n, Tdi_F_ReSeq_Di10$cov_n, Tdi_M_18_3997$cov_n, Tdi_M_18_3998$cov_n))
colnames(Tdi_all_norm_c) <- c("scaf_name", "start", "end", "mid_pos", "Tdi_F_ReSeq_Di02_cov_n", "Tdi_F_ReSeq_Di04_cov_n", "Tdi_F_ReSeq_Di06_cov_n", "Tdi_F_ReSeq_Di08_cov_n", "Tdi_F_ReSeq_Di10_cov_n", "Tdi_M_18_3997_cov_n", "Tdi_M_18_3998_cov_n")

Tdi_all_norm_c$start <- as.numeric(Tdi_all_norm_c$start )
Tdi_all_norm_c$end   <- as.numeric(Tdi_all_norm_c$end )
Tdi_all_norm_c$mid_pos <- as.numeric(Tdi_all_norm_c$mid_pos)
Tdi_all_norm_c$Tdi_F_ReSeq_Di02_cov_n <- as.numeric(Tdi_all_norm_c$Tdi_F_ReSeq_Di02_cov_n)
Tdi_all_norm_c$Tdi_F_ReSeq_Di04_cov_n <- as.numeric(Tdi_all_norm_c$Tdi_F_ReSeq_Di04_cov_n)
Tdi_all_norm_c$Tdi_F_ReSeq_Di06_cov_n <- as.numeric(Tdi_all_norm_c$Tdi_F_ReSeq_Di06_cov_n)
Tdi_all_norm_c$Tdi_F_ReSeq_Di08_cov_n <- as.numeric(Tdi_all_norm_c$Tdi_F_ReSeq_Di08_cov_n)
Tdi_all_norm_c$Tdi_F_ReSeq_Di10_cov_n <- as.numeric(Tdi_all_norm_c$Tdi_F_ReSeq_Di10_cov_n)
Tdi_all_norm_c$Tdi_M_18_3997_cov_n    <- as.numeric(Tdi_all_norm_c$Tdi_M_18_3997_cov_n)
Tdi_all_norm_c$Tdi_M_18_3998_cov_n    <- as.numeric(Tdi_all_norm_c$Tdi_M_18_3998_cov_n)

Tdi_all_norm_c$female_cov_n = 
  (Tdi_all_norm_c$Tdi_F_ReSeq_Di02_cov_n +
     Tdi_all_norm_c$Tdi_F_ReSeq_Di04_cov_n +
     Tdi_all_norm_c$Tdi_F_ReSeq_Di06_cov_n +
     Tdi_all_norm_c$Tdi_F_ReSeq_Di08_cov_n +
     Tdi_all_norm_c$Tdi_F_ReSeq_Di10_cov_n) / 5

Tdi_all_norm_c$male_cov_n = 
  (Tdi_all_norm_c$Tdi_M_18_3997_cov_n + Tdi_all_norm_c$Tdi_M_18_3998_cov_n) / 2

Tdi_all_norm_c$MF <- log2(Tdi_all_norm_c$male_cov_n / Tdi_all_norm_c$female_cov_n)

Tdi_all_norm_c_cut <- subset(Tdi_all_norm_c, Tdi_all_norm_c$MF < -0.5) 
Tdi_Sex_chr_peak <- peakfinder(Tdi_all_norm_c_cut$MF)
Tdi_Auto_peak    <- peakfinder(Tdi_all_norm_c$MF)

Tdi_Sex_chr_peak
Tdi_Auto_peak    


Tdi_all_norm_c$XA_s <- ifelse(Tdi_all_norm_c$MF < Tdi_Sex_chr_peak + 0.1 & Tdi_all_norm_c$MF > Tdi_Sex_chr_peak - 0.1, "X", "A")
Tdi_all_norm_c$XA_l <- ifelse(Tdi_all_norm_c$MF < Tdi_Auto_peak - 0.5, "X", "A")

write.csv(Tdi_all_norm_c, "Tdi_all_norm_c_windows_100000_to_Tps.csv", row.names=FALSE)

png(filename = "Tdi_all_norm_c_windows_100000_hist_to_Tps.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
hist(Tdi_all_norm_c$MF, breaks = 400, xlim = c(-2,2))
abline(v=Tdi_Sex_chr_peak, col='red', lwd=1)
abline(v=Tdi_Auto_peak, col='blue', lwd=1)
dev.off()
getwd() ## where has my plot gone....


##################### without the scafs I WANT (i.e. the big groups)

head(Tdi_all_norm_c)
Tdi_all_norm_c_excluded <- Tdi_all_norm_c[!Tdi_all_norm_c$scaf_name %in% want_scafs_Tdi,]  

png(filename = "Tdi_all_norm_c_excluded_windows_100000_hist_to_Tps.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
hist(Tdi_all_norm_c_excluded$MF, breaks = 400, xlim = c(-2,2))
abline(v=Tdi_Sex_chr_peak, col='red', lwd=1)
abline(v=Tdi_Auto_peak, col='blue', lwd=1)
dev.off()
getwd() ## where has my plot gone....



###########################################################################################################################################################
##### plot MF by scaf dotplot


### need to get offset first

Tdi_chr_lens_s <- Tdi_chr_lens[gtools::mixedorder(Tdi_chr_lens$scaf_name), ]
head(Tdi_chr_lens_s)

offset_v <- c(0)
c_total <- 0
for(i in seq(1, length(Tdi_chr_lens_s[,1]) -1)){
  c_len <- Tdi_chr_lens_s[i,]$end
  print(c_len)
  c_total <- c_total + c_len
  offset_v <- c(offset_v, c_total)
  
}

Tdi_chr_lens_s$offset <- offset_v

## add to dict

Tdi_offset_dict <- hash()
for(i in seq(1:length(Tdi_chr_lens_s[,1]))){
  chr_n <- Tdi_chr_lens_s$scaf_name[i]
  offset_n <- Tdi_chr_lens_s$offset[i]
  Tdi_offset_dict [[chr_n]] <-  offset_n
}

## add to cov data
Tdi_chr_offset <- c()
for(i in seq(1:length(Tdi_all_norm_c[,1]))){
  scaf_n <- Tdi_all_norm_c$scaf_name[i]
  off_n <- Tdi_offset_dict[[scaf_n]]
  #### returns a NULL if key not found. Turn into NA
  if(length(off_n) == 0){
    off_n = NA
  }
  print(i)
  print(scaf_n)
  print(off_n)
  Tdi_chr_offset <- c(Tdi_chr_offset, off_n)
}

Tdi_all_norm_c <- cbind(Tdi_all_norm_c, Tdi_chr_offset)
Tdi_all_norm_c$genome_pos <- Tdi_all_norm_c$mid_pos + Tdi_all_norm_c$Tdi_chr_offset

Tdi_all_norm_c$scaf_class_1 <- ifelse(Tdi_all_norm_c$scaf_name == "Tps_LRv5b_scf1", Tdi_all_norm_c$scaf_name,
                                      ifelse(Tdi_all_norm_c$scaf_name == "Tps_LRv5b_scf2", Tdi_all_norm_c$scaf_name,
                                             ifelse(Tdi_all_norm_c$scaf_name == "Tps_LRv5b_scf3", Tdi_all_norm_c$scaf_name,
                                                    ifelse(Tdi_all_norm_c$scaf_name == "Tps_LRv5b_scf4", Tdi_all_norm_c$scaf_name,
                                                           ifelse(Tdi_all_norm_c$scaf_name == "Tps_LRv5b_scf5", Tdi_all_norm_c$scaf_name,
                                                                  ifelse(Tdi_all_norm_c$scaf_name == "Tps_LRv5b_scf6", Tdi_all_norm_c$scaf_name,
                                                                         ifelse(Tdi_all_norm_c$scaf_name == "Tps_LRv5b_scf7", Tdi_all_norm_c$scaf_name,
                                                                                ifelse(Tdi_all_norm_c$scaf_name == "Tps_LRv5b_scf8", Tdi_all_norm_c$scaf_name,
                                                                                       ifelse(Tdi_all_norm_c$scaf_name == "Tps_LRv5b_scf9", Tdi_all_norm_c$scaf_name,
                                                                                              ifelse(Tdi_all_norm_c$scaf_name == "Tps_LRv5b_scf10", Tdi_all_norm_c$scaf_name,
                                                                                                     ifelse(Tdi_all_norm_c$scaf_name == "Tps_LRv5b_scf11", Tdi_all_norm_c$scaf_name,
                                                                                                            ifelse(Tdi_all_norm_c$scaf_name == "Tps_LRv5b_scf12", Tdi_all_norm_c$scaf_name,"other"))))))))))))

Tdi_all_norm_c$scaf_class_1o <- ordered(Tdi_all_norm_c$scaf_class_1, levels= c("Tps_LRv5b_scf1", "Tps_LRv5b_scf2", "Tps_LRv5b_scf3", "Tps_LRv5b_scf4", "Tps_LRv5b_scf5", "Tps_LRv5b_scf6", "Tps_LRv5b_scf7", "Tps_LRv5b_scf8", "Tps_LRv5b_scf9", "Tps_LRv5b_scf10","Tps_LRv5b_scf11", "Tps_LRv5b_scf12", "other"))


head(Tdi_all_norm_c)
tail(Tdi_all_norm_c)


## all
ggplot(Tdi_all_norm_c, aes(x=genome_pos, MF, col = scaf_class_1o))  + geom_point(size = 0.5) + 
  theme_bw() +
  scale_color_manual(values=c("firebrick3", "black", "firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black", "firebrick3")) + ##### alter colors manually
  ylim(c(-2, 2))

## drop 'other'


Tdi_all_norm_c_LGs <- subset(Tdi_all_norm_c, Tdi_all_norm_c$scaf_class_1 != "other")



png(filename = "Tdi_cov_dotplot_LGs_to_Tps.png", width  = 12, height = 7, units = "in", bg = "white", res = 300)
plot_grid(Tdi_cov_plot, Tdi_legend_sep, rel_widths  = c(1, 0.5))
ggplot(Tdi_all_norm_c_LGs, aes(x=genome_pos, MF, col = scaf_class_1o))  + geom_point(size = 0.5) + 
  theme_bw() +
  scale_color_manual(values=c("firebrick3", "black", "firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black", "firebrick3")) + ##### alter colors manually
  ylim(c(-2, 2)) + ylab("log2(male:female coverage)") + xlab("Genome position (bp)") + theme(legend.position = "none", axis.text=element_text(size=14),  axis.title=element_text(size=14,face="bold")) + ggtitle("Tdi to Tps")
dev.off()
getwd() ## where has my plot gone....?












####################################################################################################################################################################################
###################################################################################################################################################################################
## 
####################################################################################################################################################################################
###################################################################################################################################################################################
## data Tge
### notes Tpa_LRv5a_scf1 == X chr
### 2N = 28 (Tpa)

Tge_F_CC59_A          <- read.table("Tge_F_CC59_A_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tge_F_CC59_A$mid_pos  <-  (Tge_F_CC59_A$end + Tge_F_CC59_A$start) / 2
Tge_F_CC59_C          <- read.table("Tge_F_CC59_C_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tge_F_CC59_C$mid_pos  <-  (Tge_F_CC59_C$end + Tge_F_CC59_C$start) / 2
Tge_F_CC65_B          <- read.table("Tge_F_CC65_B_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tge_F_CC65_B$mid_pos  <-  (Tge_F_CC65_B$end + Tge_F_CC65_B$start) / 2
Tge_F_CC66_A          <- read.table("Tge_F_CC66_A_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tge_F_CC66_A$mid_pos  <-  (Tge_F_CC66_A$end + Tge_F_CC66_A$start) / 2
Tge_F_CC67_A          <- read.table("Tge_F_CC67_A_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tge_F_CC67_A$mid_pos  <-  (Tge_F_CC67_A$end + Tge_F_CC67_A$start) / 2
Tge_M_10091           <- read.table("Tge_M_10091_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tge_M_10091$mid_pos   <-  (Tge_M_10091$end + Tge_M_10091$start) / 2
Tge_F_19_1002         <- read.table("Tge_F_19-1002_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tge_F_19_1002$mid_pos <-  (Tge_F_19_1002$end + Tge_F_19_1002$start) / 2
Tge_M_19_1004         <- read.table("Tge_M_19-1004_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tge_M_19_1004$mid_pos <-  (Tge_M_19_1004$end + Tge_M_19_1004$start) / 2

### add some info

Tge_F_CC59_A$samp_name  <- rep("Tge_F_CC59_A",  length(Tge_F_CC59_A[,1]))
Tge_F_CC59_C$samp_name  <- rep("Tge_F_CC59_C",  length(Tge_F_CC59_C[,1]))
Tge_F_CC65_B$samp_name  <- rep("Tge_F_CC65_B",  length(Tge_F_CC65_B[,1]))
Tge_F_CC66_A$samp_name  <- rep("Tge_F_CC66_A",  length(Tge_F_CC66_A[,1]))
Tge_F_CC67_A$samp_name  <- rep("Tge_F_CC67_A",  length(Tge_F_CC67_A[,1]))
Tge_M_10091$samp_name   <- rep("Tge_M_10091",   length(Tge_M_10091[,1]))
Tge_F_19_1002$samp_name <- rep("Tge_F_19_1002", length(Tge_F_19_1002[,1]))
Tge_M_19_1004$samp_name <- rep("Tge_M_19_1004", length(Tge_M_19_1004[,1]))


Tge_F_CC59_A$sex  <- rep("F",  length(Tge_F_CC59_A[,1]))
Tge_F_CC59_C$sex  <- rep("F",  length(Tge_F_CC59_C[,1]))
Tge_F_CC65_B$sex  <- rep("F",  length(Tge_F_CC65_B[,1]))
Tge_F_CC66_A$sex  <- rep("F",  length(Tge_F_CC66_A[,1]))
Tge_F_CC67_A$sex  <- rep("F",  length(Tge_F_CC67_A[,1]))
Tge_M_10091$sex   <- rep("M",   length(Tge_M_10091[,1]))
Tge_F_19_1002$sex <- rep("F", length(Tge_F_19_1002[,1]))
Tge_M_19_1004$sex <- rep("M", length(Tge_M_19_1004[,1]))

#############################################################################################
### chr lens

Tge_chr_lens <- aggregate(end~scaf_name, Tge_F_CC59_A, FUN=max)
head(Tge_chr_lens[order(Tge_chr_lens$end,decreasing=T),], n = 20)
### big drop after 14


want_scafs_Tge <- c(
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


Tge_F_CC59_A_ws     <- Tge_F_CC59_A[Tge_F_CC59_A$scaf_name %in% want_scafs_Tge,]
Tge_F_CC59_A_ws_med <- aggregate(Tge_F_CC59_A_ws$cov, list(Tge_F_CC59_A_ws$scaf_name), FUN=median)
colnames(Tge_F_CC59_A_ws_med) <- c("scaf", "Tge_F_CC59_A_med_window_cov")
Tge_F_CC59_C_ws     <- Tge_F_CC59_C[Tge_F_CC59_C$scaf_name %in% want_scafs_Tge,]
Tge_F_CC59_C_ws_med <- aggregate(Tge_F_CC59_C_ws$cov, list(Tge_F_CC59_C_ws$scaf_name), FUN=median)
colnames(Tge_F_CC59_C_ws_med) <- c("scaf", "Tge_F_CC59_C_med_window_cov")
Tge_F_CC65_B_ws     <- Tge_F_CC65_B[Tge_F_CC65_B$scaf_name %in% want_scafs_Tge,]
Tge_F_CC65_B_ws_med <- aggregate(Tge_F_CC65_B_ws$cov, list(Tge_F_CC65_B_ws$scaf_name), FUN=median)
colnames(Tge_F_CC65_B_ws_med) <- c("scaf", "Tge_F_CC65_B_med_window_cov")
Tge_F_CC66_A_ws     <- Tge_F_CC66_A[Tge_F_CC66_A$scaf_name %in% want_scafs_Tge,]
Tge_F_CC66_A_ws_med <- aggregate(Tge_F_CC66_A_ws$cov, list(Tge_F_CC66_A_ws$scaf_name), FUN=median)
colnames(Tge_F_CC66_A_ws_med) <- c("scaf", "Tge_F_CC66_A_med_window_cov")
Tge_F_CC67_A_ws     <- Tge_F_CC67_A[Tge_F_CC67_A$scaf_name %in% want_scafs_Tge,]
Tge_F_CC67_A_ws_med <- aggregate(Tge_F_CC67_A_ws$cov, list(Tge_F_CC67_A_ws$scaf_name), FUN=median)
colnames(Tge_F_CC67_A_ws_med) <- c("scaf", "Tge_F_CC67_A_med_window_cov")
Tge_M_10091_ws     <- Tge_M_10091[Tge_M_10091$scaf_name %in% want_scafs_Tge,]
Tge_M_10091_ws_med <- aggregate(Tge_M_10091_ws$cov, list(Tge_M_10091_ws$scaf_name), FUN=median)
colnames(Tge_M_10091_ws_med) <- c("scaf", "Tge_M_10091_med_window_cov")
Tge_F_19_1002_ws     <- Tge_F_19_1002[Tge_F_19_1002$scaf_name %in% want_scafs_Tge,]
Tge_F_19_1002_ws_med <- aggregate(Tge_F_19_1002_ws$cov, list(Tge_F_19_1002_ws$scaf_name), FUN=median)
colnames(Tge_F_19_1002_ws_med) <- c("scaf", "Tge_F_19_1002_med_window_cov")
Tge_M_19_1004_ws     <- Tge_M_19_1004[Tge_M_19_1004$scaf_name %in% want_scafs_Tge,]
Tge_M_19_1004_ws_med <- aggregate(Tge_M_19_1004_ws$cov, list(Tge_M_19_1004_ws$scaf_name), FUN=median)
colnames(Tge_M_19_1004_ws_med) <- c("scaf", "Tge_M_19_1004_med_window_cov")

match_dfs <- c("Tge_F_CC59_A_ws_med","Tge_F_CC59_C_ws_med", "Tge_F_CC65_B_ws_med", "Tge_F_CC66_A_ws_med", "Tge_F_CC67_A_ws_med", "Tge_M_10091_ws_med", "Tge_F_19_1002_ws_med" , "Tge_M_19_1004_ws_med")
match_name <- c("","", "", "", "", "", "", "")
Tge_LRv5a_med_window_cov <- merge_df_by_col(match_dfs,match_name)
write.csv(Tge_LRv5a_med_window_cov, "Tge_LRv5a_med_window_cov_to_Tpa.csv")



######
Tge_F_CC59_A <- norm_cov(Tge_F_CC59_A)
Tge_F_CC59_C <- norm_cov(Tge_F_CC59_C)
Tge_F_CC65_B <- norm_cov(Tge_F_CC65_B)
Tge_F_CC66_A <- norm_cov(Tge_F_CC66_A)
Tge_F_CC67_A <- norm_cov(Tge_F_CC67_A)
Tge_M_10091  <- norm_cov(Tge_M_10091)
Tge_F_19_1002 <- norm_cov(Tge_F_19_1002)
Tge_M_19_1004 <- norm_cov(Tge_M_19_1004)

Tge_all_norm <- rbind(
  Tge_F_CC59_A,
  Tge_F_CC59_C,
  Tge_F_CC65_B,
  Tge_F_CC66_A,
  Tge_F_CC67_A,
  Tge_M_10091,
  Tge_F_19_1002,
  Tge_M_19_1004)

Tge_window_legend_sep <- cowplot::get_legend(  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf1")$p2)

Tge_window_p1 <- plot_grid(
  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf1")$p1,
  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf2")$p1,
  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf3")$p1,
  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf4")$p1,
  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf5")$p1,
  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf6")$p1,
  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf7")$p1,
  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf8")$p1,
  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf9")$p1,
  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf10")$p1,
  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf11")$p1,
  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf12")$p1,
  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf13")$p1,
  plot_MF_cov_all(Tge_all_norm, "Tpa_LRv5a_scf14")$p1,
  ncol = 1)


png(filename = "All_Tge_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 10, height = 20, units = "in", bg = "white", res = 300)
plot_grid(Tge_window_p1, Tge_window_legend_sep, rel_widths  = c(1, 0.2))
dev.off()
getwd() ## where has my plot gone....

cov_plot_hist_Tge <- function(df, max_cov) {
  df1 <- df
  samp = deparse(substitute(df))
  print(samp)
  p1 <- ggplot(df1) +
    theme_bw() +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf1"),   aes(cov, color="Tpa_LRv5a_scf1"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf2"),   aes(cov, color="Tpa_LRv5a_scf2"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf3"),   aes(cov, color="Tpa_LRv5a_scf3"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf4"),   aes(cov, color="Tpa_LRv5a_scf4"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf5"),   aes(cov, color="Tpa_LRv5a_scf5"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf6"),   aes(cov, color="Tpa_LRv5a_scf6"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf7"),   aes(cov, color="Tpa_LRv5a_scf7"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf8"),   aes(cov, color="Tpa_LRv5a_scf8"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf9"),   aes(cov, color="Tpa_LRv5a_scf9"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf10"),   aes(cov, color="Tpa_LRv5a_scf10"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf11"),   aes(cov, color="Tpa_LRv5a_scf11"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf12"),   aes(cov, color="Tpa_LRv5a_scf12"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf13"),   aes(cov, color="Tpa_LRv5a_scf13"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf14"),   aes(cov, color="Tpa_LRv5a_scf14"), size = 1, stat="density") +
    scale_color_manual(name = "Scaf", values = c("Tpa_LRv5a_scf1"    = "red3",
                                                 "Tpa_LRv5a_scf2"    = "#1B9E77",
                                                 "Tpa_LRv5a_scf3"    = "#D95F02",
                                                 "Tpa_LRv5a_scf4"    = "#E7298A",
                                                 "Tpa_LRv5a_scf5"  = "#66A61E",
                                                 "Tpa_LRv5a_scf6"  = "#E6AB02",
                                                 "Tpa_LRv5a_scf7"    = "#A6761D",
                                                 "Tpa_LRv5a_scf8"  = "yellow2",
                                                 "Tpa_LRv5a_scf9"    = "#666666",
                                                 "Tpa_LRv5a_scf10"   = "lightblue",
                                                 "Tpa_LRv5a_scf11" = "royalblue2",
                                                 "Tpa_LRv5a_scf12"   = "darkorchid",
                                                 "Tpa_LRv5a_scf13"   = "grey",
                                                 "Tpa_LRv5a_scf14"   = "black")) + 
    
    xlim(c(0, max_cov)) + 
    ggtitle(samp) +
    xlab("Coverage") 
  
  p2 <- p1 + theme(legend.position="none")
  
  outlist = list("p2" = p2, "p1" = p1)
  return(outlist)	
}



#plot_grid(cov_plot_hist_Tge(Tge_F_CC59_A, 30)$p2, cov_plot_hist_Tge(Tge_F_CC59_A, 30)$legend_sep)

Tge_legend_sep <- cowplot::get_legend(cov_plot_hist_Tge(Tge_F_CC59_A, 30)$p1)

Tge_cov_plot <- plot_grid(
  cov_plot_hist_Tge(Tge_F_CC59_A, 30)$p2,
  cov_plot_hist_Tge(Tge_F_CC59_C, 30)$p2,
  cov_plot_hist_Tge(Tge_F_CC65_B, 30)$p2,
  cov_plot_hist_Tge(Tge_F_CC66_A, 30)$p2,
  cov_plot_hist_Tge(Tge_F_CC67_A, 30)$p2,
  cov_plot_hist_Tge(Tge_F_19_1002, 80)$p2,
  cov_plot_hist_Tge(Tge_M_10091, 20)$p2,
  cov_plot_hist_Tge(Tge_M_19_1004, 60)$p2,
  ncol = 2)


png(filename = "Tge_cov_plot_to_Tpa.png", width  = 12, height = 15, units = "in", bg = "white", res = 300)
plot_grid(Tge_cov_plot, Tge_legend_sep, rel_widths  = c(1, 0.5))
dev.off()
getwd() ## where has my plot gone....?


################################################################################################################
#### class windows X A

head(Tge_F_19_1002)



Tge_all_norm_c <- as.data.frame(cbind(Tge_F_CC59_A$scaf_name, 
                                      Tge_F_CC59_A$start, 
                                      Tge_F_CC59_A$end, 
                                      Tge_F_CC59_A$mid_pos, 
                                      Tge_F_CC59_A$cov_n, 
                                      Tge_F_CC59_C$cov_n, 
                                      Tge_F_CC65_B$cov_n, 
                                      Tge_F_CC66_A$cov_n, 
                                      Tge_F_CC67_A$cov_n, 
                                      Tge_F_19_1002$cov_n, 
                                      Tge_M_10091$cov_n, 
                                      Tge_M_19_1004$cov_n))



head(Tge_all_norm_c )



colnames(Tge_all_norm_c) <- c("scaf_name", "start", "end", "mid_pos", "Tge_F_CC59_A_cov_n", "Tge_F_CC59_C_cov_n", "Tge_F_CC65_B_cov_n", "Tge_F_CC66_A_cov_n", "Tge_F_CC67_A_cov_n", "Tge_F_19_1002_cov_n", "Tge_M_10091_cov_n", "Tge_M_19_1004_cov_n")

Tge_all_norm_c$Tge_F_CC59_A_cov_n  <- as.numeric(Tge_all_norm_c$Tge_F_CC59_A_cov_n)
Tge_all_norm_c$Tge_F_CC59_C_cov_n  <- as.numeric(Tge_all_norm_c$Tge_F_CC59_C_cov_n)
Tge_all_norm_c$Tge_F_CC65_B_cov_n  <- as.numeric(Tge_all_norm_c$Tge_F_CC65_B_cov_n)
Tge_all_norm_c$Tge_F_CC66_A_cov_n  <- as.numeric(Tge_all_norm_c$Tge_F_CC66_A_cov_n)
Tge_all_norm_c$Tge_F_CC67_A_cov_n  <- as.numeric(Tge_all_norm_c$Tge_F_CC67_A_cov_n)
Tge_all_norm_c$Tge_F_19_1002_cov_n <- as.numeric(Tge_all_norm_c$Tge_F_19_1002_cov_n)
Tge_all_norm_c$Tge_M_10091_cov_n   <- as.numeric(Tge_all_norm_c$Tge_M_10091_cov_n)
Tge_all_norm_c$Tge_M_19_1004_cov_n <- as.numeric(Tge_all_norm_c$Tge_M_19_1004_cov_n)

Tge_all_norm_c$female_cov_n = 
  (Tge_all_norm_c$Tge_F_CC59_A_cov_n + 
   Tge_all_norm_c$Tge_F_CC59_C_cov_n +
   Tge_all_norm_c$Tge_F_CC65_B_cov_n +
   Tge_all_norm_c$Tge_F_CC66_A_cov_n +
   Tge_all_norm_c$Tge_F_CC67_A_cov_n +
   Tge_all_norm_c$Tge_F_19_1002_cov_n) / 6
    
Tge_all_norm_c$male_cov_n = 
  (Tge_all_norm_c$Tge_M_10091_cov_n +
   Tge_all_norm_c$Tge_M_19_1004_cov_n) / 2

Tge_all_norm_c$MF <- log2(Tge_all_norm_c$male_cov_n / Tge_all_norm_c$female_cov_n)

Tge_all_norm_c_cut <- subset(Tge_all_norm_c, Tge_all_norm_c$MF < -0.5) 
Tge_Sex_chr_peak <- peakfinder(Tge_all_norm_c_cut$MF)
Tge_Auto_peak    <- peakfinder(Tge_all_norm_c$MF)

Tge_Sex_chr_peak
Tge_Auto_peak    


Tge_all_norm_c$XA_s <- ifelse(Tge_all_norm_c$MF < Tge_Sex_chr_peak + 0.1 & Tge_all_norm_c$MF > Tge_Sex_chr_peak - 0.1, "X", "A")
Tge_all_norm_c$XA_l <- ifelse(Tge_all_norm_c$MF < Tge_Auto_peak - 0.5, "X", "A")

write.csv(Tge_all_norm_c, "Tge_all_norm_c_windows_100000_to_Tpa.csv", row.names=FALSE)

png(filename = "Tge_all_norm_c_windows_100000_hist_to_Tpa.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
hist(Tge_all_norm_c$MF, breaks = 400, xlim = c(-2,2))
abline(v=Tge_Sex_chr_peak, col='red', lwd=1)
abline(v=Tge_Auto_peak, col='blue', lwd=1)
dev.off()
getwd() ## where has my plot gone....


##################### without the scafs I WANT (i.e. the big groups)

head(Tge_all_norm_c)
Tge_all_norm_c_excluded <- Tge_all_norm_c[!Tge_all_norm_c$scaf_name %in% want_scafs_Tge,]  

png(filename = "Tge_all_norm_c_excluded_windows_100000_hist_to_Tpa.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
hist(Tge_all_norm_c_excluded$MF, breaks = 400, xlim = c(-2,2))
abline(v=Tge_Sex_chr_peak, col='red', lwd=1)
abline(v=Tge_Auto_peak, col='blue', lwd=1)
dev.off()
getwd() ## where has my plot gone....



###########################################################################################################################################################
##### plot MF by scaf dotplot


### need to get offset first

Tge_chr_lens_s <- Tge_chr_lens[gtools::mixedorder(Tge_chr_lens$scaf_name), ]
head(Tge_chr_lens_s, n = 30)

offset_v <- c(0)
c_total <- 0
for(i in seq(1, length(Tge_chr_lens_s[,1]) -1)){
  c_len <- Tge_chr_lens_s[i,]$end
  print(c_len)
  c_total <- c_total + c_len
  offset_v <- c(offset_v, c_total)
  
}

Tge_chr_lens_s$offset <- offset_v

## add to dict

Tge_offset_dict <- hash()
for(i in seq(1:length(Tge_chr_lens_s[,1]))){
  chr_n <- Tge_chr_lens_s$scaf_name[i]
  offset_n <- Tge_chr_lens_s$offset[i]
  Tge_offset_dict [[chr_n]] <-  offset_n
}

## add to cov data
Tge_chr_offset <- c()
for(i in seq(1:length(Tge_all_norm_c[,1]))){
  scaf_n <- Tge_all_norm_c$scaf_name[i]
  off_n <- Tge_offset_dict[[scaf_n]]
  #### returns a NULL if key not found. Turn into NA
  if(length(off_n) == 0){
    off_n = NA
  }
  print(i)
  print(scaf_n)
  print(off_n)
  Tge_chr_offset <- c(Tge_chr_offset, off_n)
}

Tge_all_norm_c$mid_pos <- as.numeric(Tge_all_norm_c$mid_pos)

Tge_all_norm_c <- cbind(Tge_all_norm_c, Tge_chr_offset)
Tge_all_norm_c$genome_pos <- Tge_all_norm_c$mid_pos + Tge_all_norm_c$Tge_chr_offset

Tge_all_norm_c$scaf_class_1 <- ifelse(Tge_all_norm_c$scaf_name == "Tpa_LRv5a_scf1", Tge_all_norm_c$scaf_name,
                                      ifelse(Tge_all_norm_c$scaf_name == "Tpa_LRv5a_scf2", Tge_all_norm_c$scaf_name,
                                             ifelse(Tge_all_norm_c$scaf_name == "Tpa_LRv5a_scf3", Tge_all_norm_c$scaf_name,
                                                    ifelse(Tge_all_norm_c$scaf_name == "Tpa_LRv5a_scf4", Tge_all_norm_c$scaf_name,
                                                           ifelse(Tge_all_norm_c$scaf_name == "Tpa_LRv5a_scf5", Tge_all_norm_c$scaf_name,
                                                                  ifelse(Tge_all_norm_c$scaf_name == "Tpa_LRv5a_scf6", Tge_all_norm_c$scaf_name,
                                                                         ifelse(Tge_all_norm_c$scaf_name == "Tpa_LRv5a_scf7", Tge_all_norm_c$scaf_name,
                                                                                ifelse(Tge_all_norm_c$scaf_name == "Tpa_LRv5a_scf8", Tge_all_norm_c$scaf_name,
                                                                                       ifelse(Tge_all_norm_c$scaf_name == "Tpa_LRv5a_scf9", Tge_all_norm_c$scaf_name,
                                                                                              ifelse(Tge_all_norm_c$scaf_name == "Tpa_LRv5a_scf10", Tge_all_norm_c$scaf_name,
                                                                                                     ifelse(Tge_all_norm_c$scaf_name == "Tpa_LRv5a_scf11", Tge_all_norm_c$scaf_name,
                                                                                                            ifelse(Tge_all_norm_c$scaf_name == "Tpa_LRv5a_scf12", Tge_all_norm_c$scaf_name,
                                                                                                                   ifelse(Tge_all_norm_c$scaf_name == "Tpa_LRv5a_scf13", Tge_all_norm_c$scaf_name,
                                                                                                                          ifelse(Tge_all_norm_c$scaf_name == "Tpa_LRv5a_scf14", Tge_all_norm_c$scaf_name,"other"))))))))))))))

Tge_all_norm_c$scaf_class_1o <- ordered(Tge_all_norm_c$scaf_class_1, levels= c("Tpa_LRv5a_scf1", "Tpa_LRv5a_scf2", "Tpa_LRv5a_scf3", "Tpa_LRv5a_scf4", "Tpa_LRv5a_scf5", "Tpa_LRv5a_scf6", "Tpa_LRv5a_scf7", "Tpa_LRv5a_scf8", "Tpa_LRv5a_scf9", "Tpa_LRv5a_scf10","Tpa_LRv5a_scf11", "Tpa_LRv5a_scf12", "Tpa_LRv5a_scf13", "Tpa_LRv5a_scf14", "other"))


head(Tge_all_norm_c)
tail(Tge_all_norm_c)


## all
ggplot(Tge_all_norm_c, aes(x=genome_pos, MF, col = scaf_class_1o))  + geom_point(size = 0.5) + 
  theme_bw() +
  scale_color_manual(values=c("firebrick3", "black", "firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black", "firebrick3", "black", "firebrick3")) + ##### alter colors manually
  ylim(c(-2, 2))

## drop 'other'


Tge_all_norm_c_LGs <- subset(Tge_all_norm_c, Tge_all_norm_c$scaf_class_1 != "other")



png(filename = "Tge_cov_dotplot_LGs_to_Tpa.png", width  = 12, height = 7, units = "in", bg = "white", res = 300)
plot_grid(Tge_cov_plot, Tge_legend_sep, rel_widths  = c(1, 0.5))
ggplot(Tge_all_norm_c_LGs, aes(x=genome_pos, MF, col = scaf_class_1o))  + geom_point(size = 0.5) + 
  theme_bw() +
  scale_color_manual(values=c("firebrick3", "black", "firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black", "firebrick3","black", "firebrick3")) + ##### alter colors manually
  ylim(c(-2, 2)) + ylab("log2(male:female coverage)") + xlab("Genome position (bp)") + theme(legend.position = "none", axis.text=element_text(size=14),  axis.title=element_text(size=14,face="bold")) + ggtitle("Tge to Tpa")
dev.off()
getwd() ## where has my plot gone....?

































####################################################################################################################################################################################
###################################################################################################################################################################################
## Tms
## Notes 2N = 26 (Tce)
## Genome is in bits - the HiC was not good enough for ordering, but OK for assigning to LG as follows


Tms_F_MS_Alp03b      <- read.table("Tms_F_MS_Alp03b_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tms_F_MS_Alp03b$mid_pos  <-  (Tms_F_MS_Alp03b$end + Tms_F_MS_Alp03b$start) / 2
Tms_F_MS_Alp04b      <- read.table("Tms_F_MS_Alp04b_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tms_F_MS_Alp04b$mid_pos  <-  (Tms_F_MS_Alp04b$end + Tms_F_MS_Alp04b$start) / 2
Tms_F_ReSeq_Ms01      <- read.table("Tms_F_ReSeq_Ms01_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tms_F_ReSeq_Ms01$mid_pos  <-  (Tms_F_ReSeq_Ms01$end + Tms_F_ReSeq_Ms01$start) / 2
Tms_F_ReSeq_Ms02      <- read.table("Tms_F_ReSeq_Ms02_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tms_F_ReSeq_Ms02$mid_pos  <-  (Tms_F_ReSeq_Ms02$end + Tms_F_ReSeq_Ms02$start) / 2
Tms_F_ReSeq_Ms03      <- read.table("Tms_F_ReSeq_Ms03_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tms_F_ReSeq_Ms03$mid_pos  <-  (Tms_F_ReSeq_Ms03$end + Tms_F_ReSeq_Ms03$start) / 2
Tms_M_25_15055      <- read.table("Tms_M_25_15055_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tms_M_25_15055$mid_pos  <-  (Tms_M_25_15055$end + Tms_M_25_15055$start) / 2
Tms_M_26_15056      <- read.table("Tms_M_26_15056_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tms_M_26_15056$mid_pos  <-  (Tms_M_26_15056$end + Tms_M_26_15056$start) / 2
Tms_M_27_15057      <- read.table("Tms_M_27_15057_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tms_M_27_15057$mid_pos  <-  (Tms_M_27_15057$end + Tms_M_27_15057$start) / 2
Tms_M_34_16002      <- read.table("Tms_M_34_16002_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tms_M_34_16002$mid_pos  <-  (Tms_M_34_16002$end + Tms_M_34_16002$start) / 2


Tms_F_MS_Alp03b$LG_mid_pos  <-  (Tms_F_MS_Alp03b$LG_end + Tms_F_MS_Alp03b$LG_start) / 2
Tms_F_MS_Alp04b$LG_mid_pos  <-  (Tms_F_MS_Alp04b$LG_end + Tms_F_MS_Alp04b$LG_start) / 2
Tms_F_ReSeq_Ms01$LG_mid_pos  <-  (Tms_F_ReSeq_Ms01$LG_end + Tms_F_ReSeq_Ms01$LG_start) / 2
Tms_F_ReSeq_Ms02$LG_mid_pos  <-  (Tms_F_ReSeq_Ms02$LG_end + Tms_F_ReSeq_Ms02$LG_start) / 2
Tms_F_ReSeq_Ms03$LG_mid_pos  <-  (Tms_F_ReSeq_Ms03$LG_end + Tms_F_ReSeq_Ms03$LG_start) / 2
Tms_M_25_15055$LG_mid_pos  <-  (Tms_M_25_15055$LG_end + Tms_M_25_15055$LG_start) / 2
Tms_M_26_15056$LG_mid_pos  <-  (Tms_M_26_15056$LG_end + Tms_M_26_15056$LG_start) / 2
Tms_M_27_15057$LG_mid_pos  <-  (Tms_M_27_15057$LG_end + Tms_M_27_15057$LG_start) / 2
Tms_M_34_16002$LG_mid_pos  <-  (Tms_M_34_16002$LG_end + Tms_M_34_16002$LG_start) / 2


### add some info
Tms_F_MS_Alp03b$samp_name   <- rep("Tms_F_MS_Alp03b",  length(Tms_F_MS_Alp03b[,1]))
Tms_F_MS_Alp04b$samp_name   <- rep("Tms_F_MS_Alp04b",  length(Tms_F_MS_Alp04b[,1]))
Tms_F_ReSeq_Ms01$samp_name  <- rep("Tms_F_ReSeq_Ms01", length(Tms_F_ReSeq_Ms01[,1]))
Tms_F_ReSeq_Ms02$samp_name  <- rep("Tms_F_ReSeq_Ms02", length(Tms_F_ReSeq_Ms02[,1]))
Tms_F_ReSeq_Ms03$samp_name  <- rep("Tms_F_ReSeq_Ms03", length(Tms_F_ReSeq_Ms03[,1]))
Tms_M_25_15055$samp_name    <- rep("Tms_M_25_15055",   length(Tms_M_25_15055[,1]))
Tms_M_26_15056$samp_name    <- rep("Tms_M_26_15056",   length(Tms_M_26_15056[,1]))
Tms_M_27_15057$samp_name    <- rep("Tms_M_27_15057",   length(Tms_M_27_15057[,1]))
Tms_M_34_16002$samp_name    <- rep("Tms_M_34_16002",   length(Tms_M_34_16002[,1]))

Tms_F_MS_Alp03b$sex   <- rep("F",  length(Tms_F_MS_Alp03b[,1]))
Tms_F_MS_Alp04b$sex   <- rep("F",  length(Tms_F_MS_Alp04b[,1]))
Tms_F_ReSeq_Ms01$sex  <- rep("F", length(Tms_F_ReSeq_Ms01[,1]))
Tms_F_ReSeq_Ms02$sex  <- rep("F", length(Tms_F_ReSeq_Ms02[,1]))
Tms_F_ReSeq_Ms03$sex  <- rep("F", length(Tms_F_ReSeq_Ms03[,1]))
Tms_M_25_15055$sex    <- rep("M",   length(Tms_M_25_15055[,1]))
Tms_M_26_15056$sex    <- rep("M",   length(Tms_M_26_15056[,1]))
Tms_M_27_15057$sex    <- rep("M",   length(Tms_M_27_15057[,1]))
Tms_M_34_16002$sex    <- rep("M",   length(Tms_M_34_16002[,1]))

Tms_F_MS_Alp03b_NA <- subset(Tms_F_MS_Alp03b, is.na(Tms_F_MS_Alp03b$cov))
Tms_F_MS_Alp03b_NA$cov_n <- rep(NA, length(Tms_F_MS_Alp03b_NA[,1]))
Tms_F_MS_Alp04b_NA <- subset(Tms_F_MS_Alp04b, is.na(Tms_F_MS_Alp04b$cov))
Tms_F_MS_Alp04b_NA$cov_n <- rep(NA, length(Tms_F_MS_Alp04b_NA[,1]))
Tms_F_ReSeq_Ms01_NA <- subset(Tms_F_ReSeq_Ms01, is.na(Tms_F_ReSeq_Ms01$cov))
Tms_F_ReSeq_Ms01_NA$cov_n <- rep(NA, length(Tms_F_ReSeq_Ms01_NA[,1]))
Tms_F_ReSeq_Ms02_NA <- subset(Tms_F_ReSeq_Ms02, is.na(Tms_F_ReSeq_Ms02$cov))
Tms_F_ReSeq_Ms02_NA$cov_n <- rep(NA, length(Tms_F_ReSeq_Ms02_NA[,1]))
Tms_F_ReSeq_Ms03_NA <- subset(Tms_F_ReSeq_Ms03, is.na(Tms_F_ReSeq_Ms03$cov))
Tms_F_ReSeq_Ms03_NA$cov_n <- rep(NA, length(Tms_F_ReSeq_Ms03_NA[,1]))
Tms_M_25_15055_NA <- subset(Tms_M_25_15055, is.na(Tms_M_25_15055$cov))
Tms_M_25_15055_NA$cov_n <- rep(NA, length(Tms_M_25_15055_NA[,1]))
Tms_M_26_15056_NA <- subset(Tms_M_26_15056, is.na(Tms_M_26_15056$cov))
Tms_M_26_15056_NA$cov_n <- rep(NA, length(Tms_M_26_15056_NA[,1]))
Tms_M_27_15057_NA <- subset(Tms_M_27_15057, is.na(Tms_M_27_15057$cov))
Tms_M_27_15057_NA$cov_n <- rep(NA, length(Tms_M_27_15057_NA[,1]))
Tms_M_34_16002_NA <- subset(Tms_M_34_16002, is.na(Tms_M_34_16002$cov))
Tms_M_34_16002_NA$cov_n <- rep(NA, length(Tms_M_34_16002_NA[,1]))






#############################################################################################
### chr lens

Tms_chr_lens <- aggregate(LG_end~LG, Tms_M_34_16002, FUN=max)
head(Tms_chr_lens[order(Tms_chr_lens$LG_end,decreasing=T),], n = 20)

want_LGs_Tms <- c("HiCLG1", "HiCLG2", "HiCLG3", "HiCLG4", "HiCLG5", "HiCLG6", "HiCLG7", "HiCLG8", "HiCLG9", "HiCLG10", "HiCLG11", "HiCLG12", "HiCLG13")


Tms_F_MS_Alp03b_ws     <- na.omit(Tms_F_MS_Alp03b[Tms_F_MS_Alp03b$LG %in% want_LGs_Tms,])
Tms_F_MS_Alp03b_ws_med <- aggregate(Tms_F_MS_Alp03b_ws$cov, list(Tms_F_MS_Alp03b_ws$LG), FUN=median)
colnames(Tms_F_MS_Alp03b_ws_med) <- c("LG", "Tms_F_MS_Alp03b_med_window_cov")
Tms_F_MS_Alp04b_ws     <- na.omit(Tms_F_MS_Alp04b[Tms_F_MS_Alp04b$LG %in% want_LGs_Tms,])
Tms_F_MS_Alp04b_ws_med <- aggregate(Tms_F_MS_Alp04b_ws$cov, list(Tms_F_MS_Alp04b_ws$LG), FUN=median)
colnames(Tms_F_MS_Alp04b_ws_med) <- c("LG", "Tms_F_MS_Alp04b_med_window_cov")
Tms_F_ReSeq_Ms01_ws     <- na.omit(Tms_F_ReSeq_Ms01[Tms_F_ReSeq_Ms01$LG %in% want_LGs_Tms,])
Tms_F_ReSeq_Ms01_ws_med <- aggregate(Tms_F_ReSeq_Ms01_ws$cov, list(Tms_F_ReSeq_Ms01_ws$LG), FUN=median)
colnames(Tms_F_ReSeq_Ms01_ws_med) <- c("LG", "Tms_F_ReSeq_Ms01_med_window_cov")
Tms_F_ReSeq_Ms02_ws     <- na.omit(Tms_F_ReSeq_Ms02[Tms_F_ReSeq_Ms02$LG %in% want_LGs_Tms,])
Tms_F_ReSeq_Ms02_ws_med <- aggregate(Tms_F_ReSeq_Ms02_ws$cov, list(Tms_F_ReSeq_Ms02_ws$LG), FUN=median)
colnames(Tms_F_ReSeq_Ms02_ws_med) <- c("LG", "Tms_F_ReSeq_Ms02_med_window_cov")
Tms_F_ReSeq_Ms03_ws     <- na.omit(Tms_F_ReSeq_Ms03[Tms_F_ReSeq_Ms03$LG %in% want_LGs_Tms,])
Tms_F_ReSeq_Ms03_ws_med <- aggregate(Tms_F_ReSeq_Ms03_ws$cov, list(Tms_F_ReSeq_Ms03_ws$LG), FUN=median)
colnames(Tms_F_ReSeq_Ms03_ws_med) <- c("LG", "Tms_F_ReSeq_Ms03_med_window_cov")
Tms_M_25_15055_ws     <- na.omit(Tms_M_25_15055[Tms_M_25_15055$LG %in% want_LGs_Tms,])
Tms_M_25_15055_ws_med <- aggregate(Tms_M_25_15055_ws$cov, list(Tms_M_25_15055_ws$LG), FUN=median)
colnames(Tms_M_25_15055_ws_med) <- c("LG", "Tms_M_25_15055_med_window_cov")
Tms_M_26_15056_ws     <- na.omit(Tms_M_26_15056[Tms_M_26_15056$LG %in% want_LGs_Tms,])
Tms_M_26_15056_ws_med <- aggregate(Tms_M_26_15056_ws$cov, list(Tms_M_26_15056_ws$LG), FUN=median)
colnames(Tms_M_26_15056_ws_med) <- c("LG", "Tms_M_26_15056_med_window_cov")
Tms_M_27_15057_ws     <- na.omit(Tms_M_27_15057[Tms_M_27_15057$LG %in% want_LGs_Tms,])
Tms_M_27_15057_ws_med <- aggregate(Tms_M_27_15057_ws$cov, list(Tms_M_27_15057_ws$LG), FUN=median)
colnames(Tms_M_27_15057_ws_med) <- c("LG", "Tms_M_27_15057_med_window_cov")
Tms_M_34_16002_ws     <- na.omit(Tms_M_34_16002[Tms_M_34_16002$LG %in% want_LGs_Tms,])
Tms_M_34_16002_ws_med <- aggregate(Tms_M_34_16002_ws$cov, list(Tms_M_34_16002_ws$LG), FUN=median)
colnames(Tms_M_34_16002_ws_med) <- c("LG", "Tms_M_34_16002_med_window_cov")

match_dfs <- c("Tms_F_MS_Alp03b_ws_med", "Tms_F_MS_Alp04b_ws_med", "Tms_F_ReSeq_Ms01_ws_med", "Tms_F_ReSeq_Ms02_ws_med", "Tms_F_ReSeq_Ms03_ws_med", "Tms_M_25_15055_ws_med", "Tms_M_26_15056_ws_med", "Tms_M_27_15057_ws_med", "Tms_M_34_16002_ws_med")
match_name <- c("","", "", "", "", "", "", "")
Tms_LRv5a_med_window_cov <- merge_df_by_col(match_dfs, match_name)
write.csv(Tms_LRv5a_med_window_cov, "Tms_LRv5a_med_window_cov_to_Tce.csv")

######  DF[!is.na(DF$y),]
Tms_F_MS_Alp03b  <- norm_cov(Tms_F_MS_Alp03b[!is.na(Tms_F_MS_Alp03b$cov),])
Tms_F_MS_Alp04b  <- norm_cov(Tms_F_MS_Alp04b[!is.na(Tms_F_MS_Alp04b$cov),])
Tms_F_ReSeq_Ms01 <- norm_cov(Tms_F_ReSeq_Ms01[!is.na(Tms_F_ReSeq_Ms01$cov),])
Tms_F_ReSeq_Ms02 <- norm_cov(Tms_F_ReSeq_Ms02[!is.na(Tms_F_ReSeq_Ms02$cov),])
Tms_F_ReSeq_Ms03 <- norm_cov(Tms_F_ReSeq_Ms03[!is.na(Tms_F_ReSeq_Ms03$cov),])
Tms_M_25_15055   <- norm_cov(Tms_M_25_15055[!is.na(Tms_M_25_15055$cov),])
Tms_M_26_15056   <- norm_cov(Tms_M_26_15056[!is.na(Tms_M_26_15056$cov),])
Tms_M_27_15057   <- norm_cov(Tms_M_27_15057[!is.na(Tms_M_27_15057$cov),])
Tms_M_34_16002   <- norm_cov(Tms_M_34_16002[!is.na(Tms_M_34_16002$cov),])


Tms_F_MS_Alp03b <- rbind(Tms_F_MS_Alp03b, Tms_F_MS_Alp03b_NA)
Tms_F_MS_Alp04b <- rbind(Tms_F_MS_Alp04b, Tms_F_MS_Alp04b_NA)
Tms_F_ReSeq_Ms01 <- rbind(Tms_F_ReSeq_Ms01, Tms_F_ReSeq_Ms01_NA)
Tms_F_ReSeq_Ms02 <- rbind(Tms_F_ReSeq_Ms02, Tms_F_ReSeq_Ms02_NA)
Tms_F_ReSeq_Ms03 <- rbind(Tms_F_ReSeq_Ms03, Tms_F_ReSeq_Ms03_NA)
Tms_M_25_15055 <- rbind(Tms_M_25_15055, Tms_M_25_15055_NA)
Tms_M_26_15056 <- rbind(Tms_M_26_15056, Tms_M_26_15056_NA)
Tms_M_27_15057 <- rbind(Tms_M_27_15057, Tms_M_27_15057_NA)
Tms_M_34_16002 <- rbind(Tms_M_34_16002, Tms_M_34_16002_NA)

Tms_all_norm <- rbind(
  Tms_F_MS_Alp03b,
  Tms_F_MS_Alp04b,
  Tms_F_ReSeq_Ms01,
  Tms_F_ReSeq_Ms02,
  Tms_F_ReSeq_Ms03,
  Tms_M_25_15055,
  Tms_M_26_15056,
  Tms_M_27_15057,
  Tms_M_34_16002)


tail(Tms_all_norm)

Tms_window_legend_sep <- cowplot::get_legend(plot_MF_cov_all_Tce(Tms_all_norm, "HiCLG1")$p2)

Tms_window_p1 <- plot_grid(
  plot_MF_cov_all_Tce(Tms_all_norm, "HiCLG1")$p1,
  plot_MF_cov_all_Tce(Tms_all_norm, "HiCLG2")$p1,
  plot_MF_cov_all_Tce(Tms_all_norm, "HiCLG3")$p1,
  plot_MF_cov_all_Tce(Tms_all_norm, "HiCLG4")$p1,
  plot_MF_cov_all_Tce(Tms_all_norm, "HiCLG5")$p1,
  plot_MF_cov_all_Tce(Tms_all_norm, "HiCLG6")$p1,
  plot_MF_cov_all_Tce(Tms_all_norm, "HiCLG7")$p1,
  plot_MF_cov_all_Tce(Tms_all_norm, "HiCLG8")$p1,
  plot_MF_cov_all_Tce(Tms_all_norm, "HiCLG9")$p1,
  plot_MF_cov_all_Tce(Tms_all_norm, "HiCLG10")$p1,
  plot_MF_cov_all_Tce(Tms_all_norm, "HiCLG11")$p1,
  plot_MF_cov_all_Tce(Tms_all_norm, "HiCLG12")$p1,
  plot_MF_cov_all_Tce(Tms_all_norm, "HiCLG13")$p1,
  ncol = 1)

png(filename = "All_Tms_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 10, height = 20, units = "in", bg = "white", res = 300)
plot_grid(Tms_window_p1, Tms_window_legend_sep, rel_widths  = c(1, 0.2))
dev.off()
getwd() ## where has my plot gone....


cov_plot_hist_Tms <- function(df, max_cov) {
  df1 <- df
  samp = deparse(substitute(df))
  print(samp)
  p1 <- ggplot(df1) +
    theme_bw() +
    geom_line(data=subset(df1,LG == "HiCLG1"),   aes(cov, color="HiCLG1"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG2"), aes(cov, color="HiCLG2"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG3"),     aes(cov, color="HiCLG3"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG4"), aes(cov, color="HiCLG4"), size = 1, stat="density") +	
    geom_line(data=subset(df1,LG == "HiCLG5"), aes(cov, color="HiCLG5"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG6"),   aes(cov, color="HiCLG6"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG7"),   aes(cov, color="HiCLG7"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG8"),     aes(cov, color="HiCLG8"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG9"),     aes(cov, color="HiCLG9"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG10"),     aes(cov, color="HiCLG10"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG11"),     aes(cov, color="HiCLG11"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG12"),     aes(cov, color="HiCLG12"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG13"),     aes(cov, color="HiCLG13"), size = 1, stat="density") +
    scale_color_manual(name = "Scaf", values = c("HiCLG1" = "#1B9E77",
                                                 "HiCLG2" = "#D95F02",
                                                 "HiCLG3" = "red3",
                                                 "HiCLG4" = "#E7298A",
                                                 "HiCLG5" = "#66A61E",
                                                 "HiCLG6" = "#E6AB02",
                                                 "HiCLG7" = "#A6761D",
                                                 "HiCLG8" = "yellow2",
                                                 "HiCLG9" = "#666666",
                                                 "HiCLG10" = "lightblue",
                                                 "HiCLG11" = "royalblue2",
                                                 "HiCLG12" = "darkorchid",
                                                 "HiCLG13" = "#7570B3")) + 
    
    xlim(c(0, max_cov)) + 
    ggtitle(samp) +
    xlab("Coverage") 
  
  p2 <- p1 + theme(legend.position="none")
  
  outlist = list("p2" = p2, "p1" = p1)
  return(outlist)	
}


Tms_legend_sep <- cowplot::get_legend(cov_plot_hist_Tms(Tms_F_MS_Alp03b, 30)$p1)

Tms_cov_plot <- plot_grid(
  cov_plot_hist_Tms(Tms_F_MS_Alp03b, 30)$p2,
  cov_plot_hist_Tms(Tms_F_MS_Alp04b, 30)$p2,
  cov_plot_hist_Tms(Tms_F_ReSeq_Ms01, 30)$p2,
  cov_plot_hist_Tms(Tms_F_ReSeq_Ms02, 30)$p2,
  cov_plot_hist_Tms(Tms_F_ReSeq_Ms03, 30)$p2,
  cov_plot_hist_Tms(Tms_M_25_15055, 40)$p2,
  cov_plot_hist_Tms(Tms_M_26_15056, 40)$p2,
  cov_plot_hist_Tms(Tms_M_27_15057, 40)$p2,
  cov_plot_hist_Tms(Tms_M_34_16002, 40)$p2,
  ncol = 2)


png(filename = "Tms_to_Tce_cov_plot.png", width  = 12, height = 15, units = "in", bg = "white", res = 300)
plot_grid(Tms_cov_plot, Tms_legend_sep, rel_widths  = c(1, 0.5))
dev.off()
getwd() ## where has my plot gone....?

################################################################################################################
#### class windows X A

Tms_all_norm_c <- as.data.frame(cbind(Tms_F_MS_Alp03b$scaf_name, 
                                      Tms_F_MS_Alp03b$start, 
                                      Tms_F_MS_Alp03b$end, 
                                      Tms_F_MS_Alp03b$mid_pos, 
                                      Tms_F_MS_Alp03b$LG, 
                                      Tms_F_MS_Alp03b$cov_n, 
                                      Tms_F_MS_Alp04b$cov_n, 
                                      Tms_F_ReSeq_Ms01$cov_n, 
                                      Tms_F_ReSeq_Ms02$cov_n, 
                                      Tms_F_ReSeq_Ms03$cov_n, 
                                      Tms_M_25_15055$cov_n, 
                                      Tms_M_26_15056$cov_n, 
                                      Tms_M_27_15057$cov_n, 
                                      Tms_M_34_16002$cov_n,
                                      Tms_F_MS_Alp03b$LG_mid_pos))



head(Tms_all_norm_c, n = 30 )

colnames(Tms_all_norm_c) <- c("scaf_name", "start", "end", "mid_pos", "LG", "Tms_F_MS_Alp03b_cov_n", "Tms_F_MS_Alp04b_cov_n", "Tms_F_ReSeq_Ms01_cov_n", "Tms_F_ReSeq_Ms02_cov_n", "Tms_F_ReSeq_Ms03_cov_n", "Tms_M_25_15055_cov_n", "Tms_M_26_15056_cov_n", "Tms_M_27_15057_cov_n", "Tms_M_34_16002_cov_n", "LG_mid")


Tms_all_norm_c$Tms_F_MS_Alp03b_cov_n   <- as.numeric(Tms_all_norm_c$Tms_F_MS_Alp03b_cov_n)
Tms_all_norm_c$Tms_F_MS_Alp04b_cov_n   <- as.numeric(Tms_all_norm_c$Tms_F_MS_Alp04b_cov_n)
Tms_all_norm_c$Tms_F_ReSeq_Ms01_cov_n  <- as.numeric(Tms_all_norm_c$Tms_F_ReSeq_Ms01_cov_n)
Tms_all_norm_c$Tms_F_ReSeq_Ms02_cov_n  <- as.numeric(Tms_all_norm_c$Tms_F_ReSeq_Ms02_cov_n)
Tms_all_norm_c$Tms_F_ReSeq_Ms03_cov_n  <- as.numeric(Tms_all_norm_c$Tms_F_ReSeq_Ms03_cov_n)
Tms_all_norm_c$Tms_M_25_15055_cov_n    <- as.numeric(Tms_all_norm_c$Tms_M_25_15055_cov_n)
Tms_all_norm_c$Tms_M_26_15056_cov_n    <- as.numeric(Tms_all_norm_c$Tms_M_26_15056_cov_n)
Tms_all_norm_c$Tms_M_27_15057_cov_n    <- as.numeric(Tms_all_norm_c$Tms_M_27_15057_cov_n)
Tms_all_norm_c$Tms_M_34_16002_cov_n    <- as.numeric(Tms_all_norm_c$Tms_M_34_16002_cov_n)

Tms_all_norm_c$female_cov_n = 
  (Tms_all_norm_c$Tms_F_MS_Alp03b_cov_n +
     Tms_all_norm_c$Tms_F_MS_Alp04b_cov_n +
     Tms_all_norm_c$Tms_F_ReSeq_Ms01_cov_n + 
     Tms_all_norm_c$Tms_F_ReSeq_Ms02_cov_n +
     Tms_all_norm_c$Tms_F_ReSeq_Ms03_cov_n) / 5

Tms_all_norm_c$male_cov_n = 
  (Tms_all_norm_c$Tms_M_25_15055_cov_n +
     Tms_all_norm_c$Tms_M_26_15056_cov_n +
     Tms_all_norm_c$Tms_M_27_15057_cov_n +
     Tms_all_norm_c$Tms_M_34_16002_cov_n) / 4

Tms_all_norm_c$MF <- log2(Tms_all_norm_c$male_cov_n / Tms_all_norm_c$female_cov_n)

Tms_all_norm_c_cut <- subset(Tms_all_norm_c, Tms_all_norm_c$MF < -0.5) 
Tms_Sex_chr_peak <- peakfinder(Tms_all_norm_c_cut$MF)
Tms_Auto_peak    <- peakfinder(Tms_all_norm_c$MF)

Tms_Sex_chr_peak
Tms_Auto_peak    


Tms_all_norm_c$XA_s <- ifelse(Tms_all_norm_c$MF < Tms_Sex_chr_peak + 0.1 & Tms_all_norm_c$MF > Tms_Sex_chr_peak - 0.1, "X", "A")
Tms_all_norm_c$XA_l <- ifelse(Tms_all_norm_c$MF < Tms_Auto_peak - 0.5, "X", "A")

write.csv(Tms_all_norm_c, "Tms_all_norm_c_windows_100000_to_Tce.csv", row.names=FALSE)

png(filename = "Tms_all_norm_c_windows_100000_hist_to_Tce.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
hist(Tms_all_norm_c$MF, breaks = 400, xlim = c(-2,2))
abline(v=Tms_Sex_chr_peak, col='red', lwd=1)
abline(v=Tms_Auto_peak, col='blue', lwd=1)
dev.off()
getwd() ## where has my plot gone....


##################### without the scafs I WANT (i.e. the LGS)

head(Tms_all_norm_c)
Tms_all_norm_c_excluded <- Tms_all_norm_c[!Tms_all_norm_c$LG %in% want_LGs_Tms,]  

png(filename = "Tms_all_norm_c_excluded_windows_100000_hist_to_Tce.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
hist(Tms_all_norm_c_excluded$MF, breaks = 400, xlim = c(-2,2))
abline(v=Tms_Sex_chr_peak, col='red', lwd=1)
abline(v=Tms_Auto_peak, col='blue', lwd=1)
dev.off()
getwd() ## where has my plot gone....





###########################################################################################################################################################
##### plot MF by scaf dotplot


### need to get offset first

Tms_chr_lens_s <- Tms_chr_lens[gtools::mixedorder(Tms_chr_lens$LG), ]
head(Tms_chr_lens_s, n = 30)

offset_v <- c(0)
c_total <- 0
for(i in seq(1, length(Tms_chr_lens_s[,1]) -1)){
  c_len <- Tms_chr_lens_s[i,]$LG_end
  print(c_len)
  c_total <- c_total + c_len
  offset_v <- c(offset_v, c_total)
  
}

Tms_chr_lens_s$offset <- offset_v

## add to dict

Tms_offset_dict <- hash()
for(i in seq(1:length(Tms_chr_lens_s[,1]))){
  chr_n <- Tms_chr_lens_s$LG[i]
  offset_n <- Tms_chr_lens_s$offset[i]
  Tms_offset_dict [[chr_n]] <-  offset_n
}

## add to cov data
Tms_chr_offset <- c()
for(i in seq(1:length(Tms_all_norm_c[,1]))){
  scaf_n <- Tms_all_norm_c$LG[i]
  off_n <- Tms_offset_dict[[scaf_n]]
  #### returns a NULL if key not found. Turn into NA
  if(length(off_n) == 0){
    off_n = NA
  }
  print(i)
  print(scaf_n)
  print(off_n)
  Tms_chr_offset <- c(Tms_chr_offset, off_n)
}

Tms_all_norm_c$LG_mid <- as.numeric(Tms_all_norm_c$LG_mid)

Tms_all_norm_c <- cbind(Tms_all_norm_c, Tms_chr_offset)
Tms_all_norm_c$LG_genome_pos <- Tms_all_norm_c$LG_mid + Tms_all_norm_c$Tms_chr_offset

Tms_all_norm_c$scaf_class_1 <- ifelse(Tms_all_norm_c$LG == "HiCLG1", Tms_all_norm_c$LG,
                                      ifelse(Tms_all_norm_c$LG == "HiCLG2", Tms_all_norm_c$LG,
                                             ifelse(Tms_all_norm_c$LG == "HiCLG3", Tms_all_norm_c$LG,
                                                    ifelse(Tms_all_norm_c$LG == "HiCLG4", Tms_all_norm_c$LG,
                                                           ifelse(Tms_all_norm_c$LG == "HiCLG5", Tms_all_norm_c$LG,
                                                                  ifelse(Tms_all_norm_c$LG == "HiCLG6", Tms_all_norm_c$LG,
                                                                         ifelse(Tms_all_norm_c$LG == "HiCLG7", Tms_all_norm_c$LG,
                                                                                ifelse(Tms_all_norm_c$LG == "HiCLG8", Tms_all_norm_c$LG,
                                                                                       ifelse(Tms_all_norm_c$LG == "HiCLG9", Tms_all_norm_c$LG,
                                                                                              ifelse(Tms_all_norm_c$LG == "HiCLG10", Tms_all_norm_c$LG,
                                                                                                     ifelse(Tms_all_norm_c$LG == "HiCLG11", Tms_all_norm_c$LG,
                                                                                                            ifelse(Tms_all_norm_c$LG == "HiCLG12", Tms_all_norm_c$LG,
                                                                                                                   ifelse(Tms_all_norm_c$LG == "HiCLG13", Tms_all_norm_c$LG,"other")))))))))))))

Tms_all_norm_c$scaf_class_1o <- ordered(Tms_all_norm_c$scaf_class_1, levels= c("HiCLG1", "HiCLG2", "HiCLG3", "HiCLG4", "HiCLG5", "HiCLG6", "HiCLG7", "HiCLG8", "HiCLG9", "HiCLG10","HiCLG11", "HiCLG12", "HiCLG13", "other"))



## get breaks

Tms_LG_breaks <- subset(Tms_all_norm, is.na(Tms_all_norm$cov))
Tms_LG_breaks <- subset(Tms_LG_breaks, Tms_LG_breaks$samp_name == "Tms_M_34_16002")

## add offset

Tms_chr_offset <- c()
for(i in seq(1:length(Tms_LG_breaks[,1]))){
  scaf_n <- Tms_LG_breaks$LG[i]
  off_n <- Tms_offset_dict[[scaf_n]]
  #### returns a NULL if key not found. Turn into NA
  if(length(off_n) == 0){
    off_n = NA
  }
  print(i)
  print(scaf_n)
  print(off_n)
  Tms_chr_offset <- c(Tms_chr_offset, off_n)
}

Tms_LG_breaks$LG_mid <- as.numeric(Tms_LG_breaks$LG_mid)

Tms_LG_breaks <- cbind(Tms_LG_breaks, Tms_chr_offset)
Tms_LG_breaks$LG_genome_pos <- Tms_LG_breaks$LG_mid + Tms_LG_breaks$Tms_chr_offset


head(Tms_all_norm_c)
tail(Tms_all_norm_c)
str(Tms_all_norm_c)

png(filename = "Tms_cov_dotplot_LGs_to_Tce.png", width  = 12, height = 7, units = "in", bg = "white", res = 300) ## drops non-LG vals 
ggplot(Tms_all_norm_c, aes(x=LG_genome_pos, MF, col = scaf_class_1o))  +  geom_vline(xintercept = Tms_LG_breaks$LG_genome_pos, size = 0.1, linetype = 2) + geom_point(size = 0.5) + 
  theme_bw() +
  scale_color_manual(values=c("firebrick3", "black", "firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black", "firebrick3", "black", "firebrick3")) + ##### alter colors manually
  ylim(c(-2, 2)) + ylab("log2(male:female coverage)") + xlab("Genome position (bp)") + theme(legend.position = "none", axis.text=element_text(size=14),  axis.title=element_text(size=14,face="bold")) + ggtitle("Tms to Tce")
dev.off()
getwd() ## where has my plot gone....?































####################################################################################################################################################################################
###################################################################################################################################################################################
## Tce
## NOTES - genome fragmented - but can see clear X peak - "Tce_LRv5a_scf94"
## 2 scafs  "Tce_LRv5a_scf42" Tce_LRv5a_scf103"  - low cov in males and females - hard to seq or lots of rep content?
## Notes 2N = 26 (Tce)
## Genome is in bits - the HiC was not good enough for ordering, but OK for assigning to LG as follows

Tce_F_CC22B      <- read.table("Tce_F_CC22B_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tce_F_CC22B$mid_pos  <-  (Tce_F_CC22B$end + Tce_F_CC22B$start) / 2
Tce_F_CC22C      <- read.table("Tce_F_CC22C_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tce_F_CC22C$mid_pos  <-  (Tce_F_CC22C$end + Tce_F_CC22C$start) / 2
Tce_F_CC24B      <- read.table("Tce_F_CC24B_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tce_F_CC24B$mid_pos  <-  (Tce_F_CC24B$end + Tce_F_CC24B$start) / 2
Tce_F_CC24C      <- read.table("Tce_F_CC24C_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tce_F_CC24C$mid_pos  <-  (Tce_F_CC24C$end + Tce_F_CC24C$start) / 2
Tce_F_CC25B      <- read.table("Tce_F_CC25B_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tce_F_CC25B$mid_pos  <-  (Tce_F_CC25B$end + Tce_F_CC25B$start) / 2
Tce_M_05_HM15      <- read.table("Tce_M_05_HM15_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tce_M_05_HM15$mid_pos  <-  (Tce_M_05_HM15$end + Tce_M_05_HM15$start) / 2
Tce_M_06_HM16      <- read.table("Tce_M_06_HM16_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tce_M_06_HM16$mid_pos  <-  (Tce_M_06_HM16$end + Tce_M_06_HM16$start) / 2
Tce_M_07_HM33      <- read.table("Tce_M_07_HM33_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tce_M_07_HM33$mid_pos  <-  (Tce_M_07_HM33$end + Tce_M_07_HM33$start) / 2
Tce_M_08_HM61      <- read.table("Tce_M_08_HM61_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000_wLG.txt", sep = "\t", header = T)
Tce_M_08_HM61$mid_pos  <-  (Tce_M_08_HM61$end + Tce_M_08_HM61$start) / 2

Tce_F_CC22B$LG_mid_pos  <-  (Tce_F_CC22B$LG_end + Tce_F_CC22B$LG_start) / 2
Tce_F_CC22C$LG_mid_pos  <-  (Tce_F_CC22C$LG_end + Tce_F_CC22C$LG_start) / 2
Tce_F_CC24B$LG_mid_pos  <-  (Tce_F_CC24B$LG_end + Tce_F_CC24B$LG_start) / 2
Tce_F_CC24C$LG_mid_pos  <-  (Tce_F_CC24C$LG_end + Tce_F_CC24C$LG_start) / 2
Tce_F_CC25B$LG_mid_pos  <-  (Tce_F_CC25B$LG_end + Tce_F_CC25B$LG_start) / 2
Tce_M_05_HM15$LG_mid_pos  <-  (Tce_M_05_HM15$LG_end + Tce_M_05_HM15$LG_start) / 2
Tce_M_06_HM16$LG_mid_pos  <-  (Tce_M_06_HM16$LG_end + Tce_M_06_HM16$LG_start) / 2
Tce_M_07_HM33$LG_mid_pos  <-  (Tce_M_07_HM33$LG_end + Tce_M_07_HM33$LG_start) / 2
Tce_M_08_HM61$LG_mid_pos  <-  (Tce_M_08_HM61$LG_end + Tce_M_08_HM61$LG_start) / 2


### add some info
Tce_F_CC22B$samp_name   <- rep("Tce_F_CC22B",  length(Tce_F_CC22B[,1]))
Tce_F_CC22C$samp_name   <- rep("Tce_F_CC22C",  length(Tce_F_CC22C[,1]))
Tce_F_CC24B$samp_name  <- rep("Tce_F_CC24B", length(Tce_F_CC24B[,1]))
Tce_F_CC24C$samp_name  <- rep("Tce_F_CC24C", length(Tce_F_CC24C[,1]))
Tce_F_CC25B$samp_name  <- rep("Tce_F_CC25B", length(Tce_F_CC25B[,1]))
Tce_M_05_HM15$samp_name    <- rep("Tce_M_05_HM15",   length(Tce_M_05_HM15[,1]))
Tce_M_06_HM16$samp_name    <- rep("Tce_M_06_HM16",   length(Tce_M_06_HM16[,1]))
Tce_M_07_HM33$samp_name    <- rep("Tce_M_07_HM33",   length(Tce_M_07_HM33[,1]))
Tce_M_08_HM61$samp_name    <- rep("Tce_M_08_HM61",   length(Tce_M_08_HM61[,1]))

Tce_F_CC22B$sex   <- rep("F",  length(Tce_F_CC22B[,1]))
Tce_F_CC22C$sex   <- rep("F",  length(Tce_F_CC22C[,1]))
Tce_F_CC24B$sex  <- rep("F", length(Tce_F_CC24B[,1]))
Tce_F_CC24C$sex  <- rep("F", length(Tce_F_CC24C[,1]))
Tce_F_CC25B$sex  <- rep("F", length(Tce_F_CC25B[,1]))
Tce_M_05_HM15$sex    <- rep("M",   length(Tce_M_05_HM15[,1]))
Tce_M_06_HM16$sex    <- rep("M",   length(Tce_M_06_HM16[,1]))
Tce_M_07_HM33$sex    <- rep("M",   length(Tce_M_07_HM33[,1]))
Tce_M_08_HM61$sex    <- rep("M",   length(Tce_M_08_HM61[,1]))


Tce_F_CC22B_NA <- subset(Tce_F_CC22B, is.na(Tce_F_CC22B$cov))
Tce_F_CC22B_NA$cov_n <- rep(NA, length(Tce_F_CC22B_NA[,1]))
Tce_F_CC22C_NA <- subset(Tce_F_CC22C, is.na(Tce_F_CC22C$cov))
Tce_F_CC22C_NA$cov_n <- rep(NA, length(Tce_F_CC22C_NA[,1]))
Tce_F_CC24B_NA <- subset(Tce_F_CC24B, is.na(Tce_F_CC24B$cov))
Tce_F_CC24B_NA$cov_n <- rep(NA, length(Tce_F_CC24B_NA[,1]))
Tce_F_CC24C_NA <- subset(Tce_F_CC24C, is.na(Tce_F_CC24C$cov))
Tce_F_CC24C_NA$cov_n <- rep(NA, length(Tce_F_CC24C_NA[,1]))
Tce_F_CC25B_NA <- subset(Tce_F_CC25B, is.na(Tce_F_CC25B$cov))
Tce_F_CC25B_NA$cov_n <- rep(NA, length(Tce_F_CC25B_NA[,1]))
Tce_M_05_HM15_NA <- subset(Tce_M_05_HM15, is.na(Tce_M_05_HM15$cov))
Tce_M_05_HM15_NA$cov_n <- rep(NA, length(Tce_M_05_HM15_NA[,1]))
Tce_M_06_HM16_NA <- subset(Tce_M_06_HM16, is.na(Tce_M_06_HM16$cov))
Tce_M_06_HM16_NA$cov_n <- rep(NA, length(Tce_M_06_HM16_NA[,1]))
Tce_M_07_HM33_NA <- subset(Tce_M_07_HM33, is.na(Tce_M_07_HM33$cov))
Tce_M_07_HM33_NA$cov_n <- rep(NA, length(Tce_M_07_HM33_NA[,1]))
Tce_M_08_HM61_NA <- subset(Tce_M_08_HM61, is.na(Tce_M_08_HM61$cov))
Tce_M_08_HM61_NA$cov_n <- rep(NA, length(Tce_M_08_HM61_NA[,1]))


#############################################################################################
### chr lens


Tce_chr_lens <- aggregate(LG_end~LG, Tce_F_CC22B, FUN=max)
head(Tce_chr_lens[order(Tce_chr_lens$LG_end,decreasing=T),], n = 20)
want_LGs_Tce <- c("HiCLG1", "HiCLG2", "HiCLG3", "HiCLG4", "HiCLG5", "HiCLG6", "HiCLG7", "HiCLG8", "HiCLG9", "HiCLG10", "HiCLG11", "HiCLG12", "HiCLG13")

Tce_F_CC22B_ws     <- na.omit(Tce_F_CC22B[Tce_F_CC22B$LG %in% want_LGs_Tce,])
Tce_F_CC22B_ws_med <- aggregate(Tce_F_CC22B_ws$cov, list(Tce_F_CC22B_ws$LG), FUN=median)
colnames(Tce_F_CC22B_ws_med) <- c("scaf", "Tce_F_CC22B_med_window_cov")
Tce_F_CC22C_ws     <- na.omit(Tce_F_CC22C[Tce_F_CC22C$LG %in% want_LGs_Tce,])
Tce_F_CC22C_ws_med <- aggregate(Tce_F_CC22C_ws$cov, list(Tce_F_CC22C_ws$LG), FUN=median)
colnames(Tce_F_CC22C_ws_med) <- c("scaf", "Tce_F_CC22C_med_window_cov")
Tce_F_CC24B_ws     <- na.omit(Tce_F_CC24B[Tce_F_CC24B$LG %in% want_LGs_Tce,])
Tce_F_CC24B_ws_med <- aggregate(Tce_F_CC24B_ws$cov, list(Tce_F_CC24B_ws$LG), FUN=median)
colnames(Tce_F_CC24B_ws_med) <- c("scaf", "Tce_F_CC24B_med_window_cov")
Tce_F_CC24C_ws     <- na.omit(Tce_F_CC24C[Tce_F_CC24C$LG %in% want_LGs_Tce,])
Tce_F_CC24C_ws_med <- aggregate(Tce_F_CC24C_ws$cov, list(Tce_F_CC24C_ws$LG), FUN=median)
colnames(Tce_F_CC24C_ws_med) <- c("scaf", "Tce_F_CC24C_med_window_cov")
Tce_F_CC25B_ws     <- na.omit(Tce_F_CC25B[Tce_F_CC25B$LG %in% want_LGs_Tce,])
Tce_F_CC25B_ws_med <- aggregate(Tce_F_CC25B_ws$cov, list(Tce_F_CC25B_ws$LG), FUN=median)
colnames(Tce_F_CC25B_ws_med) <- c("scaf", "Tce_F_CC25B_med_window_cov")
Tce_M_05_HM15_ws     <- na.omit(Tce_M_05_HM15[Tce_M_05_HM15$LG %in% want_LGs_Tce,])
Tce_M_05_HM15_ws_med <- aggregate(Tce_M_05_HM15_ws$cov, list(Tce_M_05_HM15_ws$LG), FUN=median)
colnames(Tce_M_05_HM15_ws_med) <- c("scaf", "Tce_M_05_HM15_med_window_cov")
Tce_M_06_HM16_ws     <- na.omit(Tce_M_06_HM16[Tce_M_06_HM16$LG %in% want_LGs_Tce,])
Tce_M_06_HM16_ws_med <- aggregate(Tce_M_06_HM16_ws$cov, list(Tce_M_06_HM16_ws$LG), FUN=median)
colnames(Tce_M_06_HM16_ws_med) <- c("scaf", "Tce_M_06_HM16_med_window_cov")
Tce_M_07_HM33_ws     <- na.omit(Tce_M_07_HM33[Tce_M_07_HM33$LG %in% want_LGs_Tce,])
Tce_M_07_HM33_ws_med <- aggregate(Tce_M_07_HM33_ws$cov, list(Tce_M_07_HM33_ws$LG), FUN=median)
colnames(Tce_M_07_HM33_ws_med) <- c("scaf", "Tce_M_07_HM33_med_window_cov")
Tce_M_08_HM61_ws     <- na.omit(Tce_M_08_HM61[Tce_M_08_HM61$LG %in% want_LGs_Tce,])
Tce_M_08_HM61_ws_med <- aggregate(Tce_M_08_HM61_ws$cov, list(Tce_M_08_HM61_ws$LG), FUN=median)
colnames(Tce_M_08_HM61_ws_med) <- c("scaf", "Tce_M_08_HM61_med_window_cov")

match_dfs <- c("Tce_F_CC22B_ws_med", "Tce_F_CC22C_ws_med", "Tce_F_CC24B_ws_med", "Tce_F_CC24C_ws_med", "Tce_F_CC25B_ws_med", "Tce_M_05_HM15_ws_med", "Tce_M_06_HM16_ws_med", "Tce_M_07_HM33_ws_med", "Tce_M_08_HM61_ws_med")
match_name <- c("","", "", "", "", "", "", "")
Tce_LRv5a_med_window_cov <- merge_df_by_col(match_dfs,match_name)
write.csv(Tce_LRv5a_med_window_cov, "Tce_LRv5a_med_window_cov_to_Tce.csv")


######
Tce_F_CC22B  <- norm_cov(Tce_F_CC22B[!is.na(Tce_F_CC22B$cov),])
Tce_F_CC22C  <- norm_cov(Tce_F_CC22C[!is.na(Tce_F_CC22C$cov),])
Tce_F_CC24B <- norm_cov(Tce_F_CC24B[!is.na(Tce_F_CC24B$cov),])
Tce_F_CC24C <- norm_cov(Tce_F_CC24C[!is.na(Tce_F_CC24C$cov),])
Tce_F_CC25B <- norm_cov(Tce_F_CC25B[!is.na(Tce_F_CC25B$cov),])
Tce_M_05_HM15   <- norm_cov(Tce_M_05_HM15[!is.na(Tce_M_05_HM15$cov),])
Tce_M_06_HM16   <- norm_cov(Tce_M_06_HM16[!is.na(Tce_M_06_HM16$cov),])
Tce_M_07_HM33   <- norm_cov(Tce_M_07_HM33[!is.na(Tce_M_07_HM33$cov),])
Tce_M_08_HM61   <- norm_cov(Tce_M_08_HM61[!is.na(Tce_M_08_HM61$cov),])


Tce_F_CC22B <- rbind(Tce_F_CC22B, Tce_F_CC22B_NA)
Tce_F_CC22C <- rbind(Tce_F_CC22C, Tce_F_CC22C_NA)
Tce_F_CC24B <- rbind(Tce_F_CC24B, Tce_F_CC24B_NA)
Tce_F_CC24C <- rbind(Tce_F_CC24C, Tce_F_CC24C_NA)
Tce_F_CC25B <- rbind(Tce_F_CC25B, Tce_F_CC25B_NA)
Tce_M_05_HM15 <- rbind(Tce_M_05_HM15, Tce_M_05_HM15_NA)
Tce_M_06_HM16 <- rbind(Tce_M_06_HM16, Tce_M_06_HM16_NA)
Tce_M_07_HM33 <- rbind(Tce_M_07_HM33, Tce_M_07_HM33_NA)
Tce_M_08_HM61 <- rbind(Tce_M_08_HM61, Tce_M_08_HM61_NA)


Tce_all_norm <- rbind(
  Tce_F_CC22B,
  Tce_F_CC22C,
  Tce_F_CC24B,
  Tce_F_CC24C,
  Tce_F_CC25B,
  Tce_M_05_HM15,
  Tce_M_06_HM16,
  Tce_M_07_HM33,
  Tce_M_08_HM61)


tail(Tce_all_norm)


Tce_window_legend_sep <- cowplot::get_legend(  plot_MF_cov_all_Tce(Tce_all_norm, "HiCLG1")$p2)

Tce_window_p1 <- plot_grid(
  plot_MF_cov_all_Tce(Tce_all_norm, "HiCLG1")$p1,
  plot_MF_cov_all_Tce(Tce_all_norm, "HiCLG2")$p1,
  plot_MF_cov_all_Tce(Tce_all_norm, "HiCLG3")$p1,
  plot_MF_cov_all_Tce(Tce_all_norm, "HiCLG4")$p1,
  plot_MF_cov_all_Tce(Tce_all_norm, "HiCLG5")$p1,
  plot_MF_cov_all_Tce(Tce_all_norm, "HiCLG6")$p1,
  plot_MF_cov_all_Tce(Tce_all_norm, "HiCLG7")$p1,
  plot_MF_cov_all_Tce(Tce_all_norm, "HiCLG8")$p1,
  plot_MF_cov_all_Tce(Tce_all_norm, "HiCLG9")$p1,
  plot_MF_cov_all_Tce(Tce_all_norm, "HiCLG10")$p1,
  plot_MF_cov_all_Tce(Tce_all_norm, "HiCLG11")$p1,
  plot_MF_cov_all_Tce(Tce_all_norm, "HiCLG12")$p1,
  plot_MF_cov_all_Tce(Tce_all_norm, "HiCLG13")$p1,
  ncol = 1)

png(filename = "All_Tce_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 10, height = 20, units = "in", bg = "white", res = 300)
plot_grid(Tce_window_p1, Tce_window_legend_sep, rel_widths  = c(1, 0.2))
dev.off()
getwd() ## where has my plot gone....


cov_plot_hist_Tce <- function(df, max_cov) {
  df1 <- df
  samp = deparse(substitute(df))
  print(samp)
  p1 <- ggplot(df1) +
    theme_bw() +
    geom_line(data=subset(df1,LG == "HiCLG1"),   aes(cov, color="HiCLG1"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG2"), aes(cov, color="HiCLG2"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG3"),     aes(cov, color="HiCLG3"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG4"), aes(cov, color="HiCLG4"), size = 1, stat="density") +	
    geom_line(data=subset(df1,LG == "HiCLG5"), aes(cov, color="HiCLG5"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG6"),   aes(cov, color="HiCLG6"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG7"),   aes(cov, color="HiCLG7"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG8"),     aes(cov, color="HiCLG8"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG9"),     aes(cov, color="HiCLG9"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG10"),     aes(cov, color="HiCLG10"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG11"),     aes(cov, color="HiCLG11"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG12"),     aes(cov, color="HiCLG12"), size = 1, stat="density") +
    geom_line(data=subset(df1,LG == "HiCLG13"),     aes(cov, color="HiCLG13"), size = 1, stat="density") +
    scale_color_manual(name = "Scaf", values = c("HiCLG1" = "#1B9E77",
                                                 "HiCLG2" = "#D95F02",
                                                 "HiCLG3" = "red3",
                                                 "HiCLG4" = "#E7298A",
                                                 "HiCLG5" = "#66A61E",
                                                 "HiCLG6" = "#E6AB02",
                                                 "HiCLG7" = "#A6761D",
                                                 "HiCLG8" = "yellow2",
                                                 "HiCLG9" = "#666666",
                                                 "HiCLG10" = "lightblue",
                                                 "HiCLG11" = "royalblue2",
                                                 "HiCLG12" = "darkorchid",
                                                 "HiCLG13" = "#7570B3")) + 
    
    xlim(c(0, max_cov)) + 
    ggtitle(samp) +
    xlab("Coverage") 
  
  p2 <- p1 + theme(legend.position="none")
  
  outlist = list("p2" = p2, "p1" = p1)
  return(outlist)	
}


Tce_legend_sep <- cowplot::get_legend(cov_plot_hist_Tce(Tce_F_CC22B, 30)$p1)

Tce_cov_plot <- plot_grid(
  cov_plot_hist_Tce(Tce_F_CC22B, 30)$p2,
  cov_plot_hist_Tce(Tce_F_CC22C, 30)$p2,
  cov_plot_hist_Tce(Tce_F_CC24B, 30)$p2,
  cov_plot_hist_Tce(Tce_F_CC24C, 30)$p2,
  cov_plot_hist_Tce(Tce_F_CC25B, 30)$p2,
  cov_plot_hist_Tce(Tce_M_05_HM15, 40)$p2,
  cov_plot_hist_Tce(Tce_M_06_HM16, 40)$p2,
  cov_plot_hist_Tce(Tce_M_07_HM33, 40)$p2,
  cov_plot_hist_Tce(Tce_M_08_HM61, 40)$p2,
  ncol = 2)


png(filename = "Tce_to_Tce_cov_plot.png", width  = 12, height = 15, units = "in", bg = "white", res = 300)
plot_grid(Tce_cov_plot, Tce_legend_sep, rel_widths  = c(1, 0.5))
dev.off()
getwd() ## where has my plot gone....?





################################################################################################################
#### class windows X A




Tce_all_norm_c <- as.data.frame(cbind(Tce_F_CC22B$scaf_name, 
                                      Tce_F_CC22B$start, 
                                      Tce_F_CC22B$end, 
                                      Tce_F_CC22B$mid_pos, 
                                      Tce_F_CC22B$LG, 
                                      Tce_F_CC22B$cov_n, 
                                      Tce_F_CC22C$cov_n, 
                                      Tce_F_CC24B$cov_n, 
                                      Tce_F_CC24C$cov_n, 
                                      Tce_F_CC25B$cov_n, 
                                      Tce_M_05_HM15$cov_n, 
                                      Tce_M_06_HM16$cov_n, 
                                      Tce_M_07_HM33$cov_n, 
                                      Tce_M_08_HM61$cov_n,
                                      Tce_F_CC22B$LG_mid_pos
                                      ))



head(Tce_all_norm_c )

colnames(Tce_all_norm_c) <- c("scaf_name", "start", "end", "mid_pos", "LG", "Tce_F_CC22B_cov_n", "Tce_F_CC22C_cov_n", "Tce_F_CC24B_cov_n", "Tce_F_CC24C_cov_n", "Tce_F_CC25B_cov_n", "Tce_M_05_HM15_cov_n", "Tce_M_06_HM16_cov_n", "Tce_M_07_HM33_cov_n", "Tce_M_08_HM61_cov_n", "LG_mid")


Tce_all_norm_c$Tce_F_CC22B_cov_n   <- as.numeric(Tce_all_norm_c$Tce_F_CC22B_cov_n)
Tce_all_norm_c$Tce_F_CC22C_cov_n   <- as.numeric(Tce_all_norm_c$Tce_F_CC22C_cov_n)
Tce_all_norm_c$Tce_F_CC24B_cov_n  <- as.numeric(Tce_all_norm_c$Tce_F_CC24B_cov_n)
Tce_all_norm_c$Tce_F_CC24C_cov_n  <- as.numeric(Tce_all_norm_c$Tce_F_CC24C_cov_n)
Tce_all_norm_c$Tce_F_CC25B_cov_n  <- as.numeric(Tce_all_norm_c$Tce_F_CC25B_cov_n)
Tce_all_norm_c$Tce_M_05_HM15_cov_n    <- as.numeric(Tce_all_norm_c$Tce_M_05_HM15_cov_n)
Tce_all_norm_c$Tce_M_06_HM16_cov_n    <- as.numeric(Tce_all_norm_c$Tce_M_06_HM16_cov_n)
Tce_all_norm_c$Tce_M_07_HM33_cov_n    <- as.numeric(Tce_all_norm_c$Tce_M_07_HM33_cov_n)
Tce_all_norm_c$Tce_M_08_HM61_cov_n    <- as.numeric(Tce_all_norm_c$Tce_M_08_HM61_cov_n)

Tce_all_norm_c$female_cov_n = 
  (Tce_all_norm_c$Tce_F_CC22B_cov_n +
     Tce_all_norm_c$Tce_F_CC22C_cov_n +
     Tce_all_norm_c$Tce_F_CC24B_cov_n + 
     Tce_all_norm_c$Tce_F_CC24C_cov_n +
     Tce_all_norm_c$Tce_F_CC25B_cov_n) / 5

Tce_all_norm_c$male_cov_n = 
  (Tce_all_norm_c$Tce_M_05_HM15_cov_n +
     Tce_all_norm_c$Tce_M_06_HM16_cov_n +
     Tce_all_norm_c$Tce_M_07_HM33_cov_n +
     Tce_all_norm_c$Tce_M_08_HM61_cov_n) / 4

Tce_all_norm_c$MF <- log2(Tce_all_norm_c$male_cov_n / Tce_all_norm_c$female_cov_n)

Tce_all_norm_c_cut <- subset(Tce_all_norm_c, Tce_all_norm_c$MF < -0.5) 
Tce_Sex_chr_peak <- peakfinder(Tce_all_norm_c_cut$MF)
Tce_Auto_peak    <- peakfinder(Tce_all_norm_c$MF)

Tce_Sex_chr_peak
Tce_Auto_peak    


Tce_all_norm_c$XA_s <- ifelse(Tce_all_norm_c$MF < Tce_Sex_chr_peak + 0.1 & Tce_all_norm_c$MF > Tce_Sex_chr_peak - 0.1, "X", "A")
Tce_all_norm_c$XA_l <- ifelse(Tce_all_norm_c$MF < Tce_Auto_peak - 0.5, "X", "A")

write.csv(Tce_all_norm_c, "Tce_all_norm_c_windows_100000_to_Tce.csv", row.names=FALSE)

png(filename = "Tce_all_norm_c_windows_100000_hist_to_Tce.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
hist(Tce_all_norm_c$MF, breaks = 400, xlim = c(-2,2))
abline(v=Tce_Sex_chr_peak, col='red', lwd=1)
abline(v=Tce_Auto_peak, col='blue', lwd=1)
dev.off()
getwd() ## where has my plot gone....


##################### without the scafs I WANT (i.e. the LGS)

head(Tce_all_norm_c)
Tce_all_norm_c_excluded <- Tce_all_norm_c[!Tce_all_norm_c$LG %in% want_LGs_Tce,]  

png(filename = "Tce_all_norm_c_excluded_windows_100000_hist_to_Tce.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
hist(Tce_all_norm_c_excluded$MF, breaks = 400, xlim = c(-2,2))
abline(v=Tce_Sex_chr_peak, col='red', lwd=1)
abline(v=Tce_Auto_peak, col='blue', lwd=1)
dev.off()
getwd() ## where has my plot gone....








###########################################################################################################################################################
##### plot MF by scaf dotplot


### need to get offset first

Tce_chr_lens_s <- Tce_chr_lens[gtools::mixedorder(Tce_chr_lens$LG), ]
head(Tce_chr_lens_s, n = 30)

offset_v <- c(0)
c_total <- 0
for(i in seq(1, length(Tce_chr_lens_s[,1]) -1)){
  c_len <- Tce_chr_lens_s[i,]$LG_end
  print(c_len)
  c_total <- c_total + c_len
  offset_v <- c(offset_v, c_total)
  
}

Tce_chr_lens_s$offset <- offset_v

## add to dict

Tce_offset_dict <- hash()
for(i in seq(1:length(Tce_chr_lens_s[,1]))){
  chr_n <- Tce_chr_lens_s$LG[i]
  offset_n <- Tce_chr_lens_s$offset[i]
  Tce_offset_dict [[chr_n]] <-  offset_n
}

## add to cov data
Tce_chr_offset <- c()
for(i in seq(1:length(Tce_all_norm_c[,1]))){
  scaf_n <- Tce_all_norm_c$LG[i]
  off_n <- Tce_offset_dict[[scaf_n]]
  #### returns a NULL if key not found. Turn into NA
  if(length(off_n) == 0){
    off_n = NA
  }
  print(i)
  print(scaf_n)
  print(off_n)
  Tce_chr_offset <- c(Tce_chr_offset, off_n)
}

Tce_all_norm_c$LG_mid <- as.numeric(Tce_all_norm_c$LG_mid)

Tce_all_norm_c <- cbind(Tce_all_norm_c, Tce_chr_offset)
Tce_all_norm_c$LG_genome_pos <- Tce_all_norm_c$LG_mid + Tce_all_norm_c$Tce_chr_offset

Tce_all_norm_c$scaf_class_1 <- ifelse(Tce_all_norm_c$LG == "HiCLG1", Tce_all_norm_c$LG,
                                      ifelse(Tce_all_norm_c$LG == "HiCLG2", Tce_all_norm_c$LG,
                                             ifelse(Tce_all_norm_c$LG == "HiCLG3", Tce_all_norm_c$LG,
                                                    ifelse(Tce_all_norm_c$LG == "HiCLG4", Tce_all_norm_c$LG,
                                                           ifelse(Tce_all_norm_c$LG == "HiCLG5", Tce_all_norm_c$LG,
                                                                  ifelse(Tce_all_norm_c$LG == "HiCLG6", Tce_all_norm_c$LG,
                                                                         ifelse(Tce_all_norm_c$LG == "HiCLG7", Tce_all_norm_c$LG,
                                                                                ifelse(Tce_all_norm_c$LG == "HiCLG8", Tce_all_norm_c$LG,
                                                                                       ifelse(Tce_all_norm_c$LG == "HiCLG9", Tce_all_norm_c$LG,
                                                                                              ifelse(Tce_all_norm_c$LG == "HiCLG10", Tce_all_norm_c$LG,
                                                                                                     ifelse(Tce_all_norm_c$LG == "HiCLG11", Tce_all_norm_c$LG,
                                                                                                            ifelse(Tce_all_norm_c$LG == "HiCLG12", Tce_all_norm_c$LG,
                                                                                                                   ifelse(Tce_all_norm_c$LG == "HiCLG13", Tce_all_norm_c$LG,"other")))))))))))))

Tce_all_norm_c$scaf_class_1o <- ordered(Tce_all_norm_c$scaf_class_1, levels= c("HiCLG1", "HiCLG2", "HiCLG3", "HiCLG4", "HiCLG5", "HiCLG6", "HiCLG7", "HiCLG8", "HiCLG9", "HiCLG10","HiCLG11", "HiCLG12", "HiCLG13", "other"))



## get breaks

Tce_LG_breaks <- subset(Tce_all_norm, is.na(Tce_all_norm$cov))
Tce_LG_breaks <- subset(Tce_LG_breaks, Tce_LG_breaks$samp_name == "Tce_F_CC22B")

## add offset

Tce_chr_offset <- c()
for(i in seq(1:length(Tce_LG_breaks[,1]))){
  scaf_n <- Tce_LG_breaks$LG[i]
  off_n <- Tce_offset_dict[[scaf_n]]
  #### returns a NULL if key not found. Turn into NA
  if(length(off_n) == 0){
    off_n = NA
  }
  print(i)
  print(scaf_n)
  print(off_n)
  Tce_chr_offset <- c(Tce_chr_offset, off_n)
}

Tce_LG_breaks$LG_mid <- as.numeric(Tce_LG_breaks$LG_mid)

Tce_LG_breaks <- cbind(Tce_LG_breaks, Tce_chr_offset)
Tce_LG_breaks$LG_genome_pos <- Tce_LG_breaks$LG_mid + Tce_LG_breaks$Tce_chr_offset


head(Tce_all_norm_c)
tail(Tce_all_norm_c)
str(Tce_all_norm_c)

png(filename = "Tce_cov_dotplot_LGs_to_Tce.png", width  = 12, height = 7, units = "in", bg = "white", res = 300) ## drops non-LG vals 
ggplot(Tce_all_norm_c, aes(x=LG_genome_pos, MF, col = scaf_class_1o))  +  geom_vline(xintercept = Tce_LG_breaks$LG_genome_pos, size = 0.1, linetype = 2) + geom_point(size = 0.5) + 
  theme_bw() +
  scale_color_manual(values=c("firebrick3", "black", "firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black", "firebrick3", "black", "firebrick3")) + ##### alter colors manually
  ylim(c(-2, 2)) + ylab("log2(male:female coverage)") + xlab("Genome position (bp)") + theme(legend.position = "none", axis.text=element_text(size=14),  axis.title=element_text(size=14,face="bold")) + ggtitle("Tce to Tce")
dev.off()
getwd() ## where has my plot gone....?


































####################################################################################################################################################################################
###################################################################################################################################################################################
## 
####################################################################################################################################################################################
###################################################################################################################################################################################
## Tpa
## Notes - good genome - X = scf1, all males look correct
### 2N = 28 (Tpa)

Tpa_F_H54      <- read.table("Tpa_F_H54_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tpa_F_H54$mid_pos  <-  (Tpa_F_H54$end + Tpa_F_H54$start) / 2
Tpa_F_H56      <- read.table("Tpa_F_H56_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tpa_F_H56$mid_pos  <-  (Tpa_F_H56$end + Tpa_F_H56$start) / 2
Tpa_F_PA_CD      <- read.table("Tpa_F_PA_CD_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tpa_F_PA_CD$mid_pos  <-  (Tpa_F_PA_CD$end + Tpa_F_PA_CD$start) / 2
Tpa_F_PA_E      <- read.table("Tpa_F_PA_E_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tpa_F_PA_E$mid_pos  <-  (Tpa_F_PA_E$end + Tpa_F_PA_E$start) / 2
Tpa_F_Pa_AB      <- read.table("Tpa_F_Pa_AB_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tpa_F_Pa_AB$mid_pos  <-  (Tpa_F_Pa_AB$end + Tpa_F_Pa_AB$start) / 2
Tpa_M_09_Tpa      <- read.table("Tpa_M_09_Tpa_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tpa_M_09_Tpa$mid_pos  <-  (Tpa_M_09_Tpa$end + Tpa_M_09_Tpa$start) / 2
Tpa_M_10_Tpa      <- read.table("Tpa_M_10_Tpa_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tpa_M_10_Tpa$mid_pos  <-  (Tpa_M_10_Tpa$end + Tpa_M_10_Tpa$start) / 2
Tpa_M_11_Tpa      <- read.table("Tpa_M_11_Tpa_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tpa_M_11_Tpa$mid_pos  <-  (Tpa_M_11_Tpa$end + Tpa_M_11_Tpa$start) / 2
Tpa_M_12_Tpa      <- read.table("Tpa_M_12_Tpa_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tpa_M_12_Tpa$mid_pos  <-  (Tpa_M_12_Tpa$end + Tpa_M_12_Tpa$start) / 2


### add some info
Tpa_F_H54$samp_name   <- rep("Tpa_F_H54",  length(Tpa_F_H54[,1]))
Tpa_F_H56$samp_name   <- rep("Tpa_F_H56",  length(Tpa_F_H56[,1]))
Tpa_F_PA_CD$samp_name  <- rep("Tpa_F_PA_CD", length(Tpa_F_PA_CD[,1]))
Tpa_F_PA_E$samp_name  <- rep("Tpa_F_PA_E", length(Tpa_F_PA_E[,1]))
Tpa_F_Pa_AB$samp_name  <- rep("Tpa_F_Pa_AB", length(Tpa_F_Pa_AB[,1]))
Tpa_M_09_Tpa$samp_name    <- rep("Tpa_M_09_Tpa",   length(Tpa_M_09_Tpa[,1]))
Tpa_M_10_Tpa$samp_name    <- rep("Tpa_M_10_Tpa",   length(Tpa_M_10_Tpa[,1]))
Tpa_M_11_Tpa$samp_name    <- rep("Tpa_M_11_Tpa",   length(Tpa_M_11_Tpa[,1]))
Tpa_M_12_Tpa$samp_name    <- rep("Tpa_M_12_Tpa",   length(Tpa_M_12_Tpa[,1]))

Tpa_F_H54$sex   <- rep("F",  length(Tpa_F_H54[,1]))
Tpa_F_H56$sex   <- rep("F",  length(Tpa_F_H56[,1]))
Tpa_F_PA_CD$sex  <- rep("F", length(Tpa_F_PA_CD[,1]))
Tpa_F_PA_E$sex  <- rep("F", length(Tpa_F_PA_E[,1]))
Tpa_F_Pa_AB$sex  <- rep("F", length(Tpa_F_Pa_AB[,1]))
Tpa_M_09_Tpa$sex    <- rep("M",   length(Tpa_M_09_Tpa[,1]))
Tpa_M_10_Tpa$sex    <- rep("M",   length(Tpa_M_10_Tpa[,1]))
Tpa_M_11_Tpa$sex    <- rep("M",   length(Tpa_M_11_Tpa[,1]))
Tpa_M_12_Tpa$sex    <- rep("M",   length(Tpa_M_12_Tpa[,1]))


#############################################################################################
### chr lens

Tpa_chr_lens <- aggregate(end~scaf_name, Tpa_M_12_Tpa, FUN=max)
head(Tpa_chr_lens[order(Tpa_chr_lens$end,decreasing=T),], n = 20)

### big drop after Tpa_LRv5a_scf14

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

Tpa_F_H54_ws     <- Tpa_F_H54[Tpa_F_H54$scaf_name %in% want_scafs_Tpa,]
Tpa_F_H54_ws_med <- aggregate(Tpa_F_H54_ws$cov, list(Tpa_F_H54_ws$scaf_name), FUN=median)
colnames(Tpa_F_H54_ws_med) <- c("scaf", "Tpa_F_H54_med_window_cov")
Tpa_F_H56_ws     <- Tpa_F_H56[Tpa_F_H56$scaf_name %in% want_scafs_Tpa,]
Tpa_F_H56_ws_med <- aggregate(Tpa_F_H56_ws$cov, list(Tpa_F_H56_ws$scaf_name), FUN=median)
colnames(Tpa_F_H56_ws_med) <- c("scaf", "Tpa_F_H56_med_window_cov")
Tpa_F_PA_CD_ws     <- Tpa_F_PA_CD[Tpa_F_PA_CD$scaf_name %in% want_scafs_Tpa,]
Tpa_F_PA_CD_ws_med <- aggregate(Tpa_F_PA_CD_ws$cov, list(Tpa_F_PA_CD_ws$scaf_name), FUN=median)
colnames(Tpa_F_PA_CD_ws_med) <- c("scaf", "Tpa_F_PA_CD_med_window_cov")
Tpa_F_PA_E_ws     <- Tpa_F_PA_E[Tpa_F_PA_E$scaf_name %in% want_scafs_Tpa,]
Tpa_F_PA_E_ws_med <- aggregate(Tpa_F_PA_E_ws$cov, list(Tpa_F_PA_E_ws$scaf_name), FUN=median)
colnames(Tpa_F_PA_E_ws_med) <- c("scaf", "Tpa_F_PA_E_med_window_cov")
Tpa_F_Pa_AB_ws     <- Tpa_F_Pa_AB[Tpa_F_Pa_AB$scaf_name %in% want_scafs_Tpa,]
Tpa_F_Pa_AB_ws_med <- aggregate(Tpa_F_Pa_AB_ws$cov, list(Tpa_F_Pa_AB_ws$scaf_name), FUN=median)
colnames(Tpa_F_Pa_AB_ws_med) <- c("scaf", "Tpa_F_Pa_AB_med_window_cov")
Tpa_M_09_Tpa_ws     <- Tpa_M_09_Tpa[Tpa_M_09_Tpa$scaf_name %in% want_scafs_Tpa,]
Tpa_M_09_Tpa_ws_med <- aggregate(Tpa_M_09_Tpa_ws$cov, list(Tpa_M_09_Tpa_ws$scaf_name), FUN=median)
colnames(Tpa_M_09_Tpa_ws_med) <- c("scaf", "Tpa_M_09_Tpa_med_window_cov")
Tpa_M_10_Tpa_ws     <- Tpa_M_10_Tpa[Tpa_M_10_Tpa$scaf_name %in% want_scafs_Tpa,]
Tpa_M_10_Tpa_ws_med <- aggregate(Tpa_M_10_Tpa_ws$cov, list(Tpa_M_10_Tpa_ws$scaf_name), FUN=median)
colnames(Tpa_M_10_Tpa_ws_med) <- c("scaf", "Tpa_M_10_Tpa_med_window_cov")
Tpa_M_11_Tpa_ws     <- Tpa_M_11_Tpa[Tpa_M_11_Tpa$scaf_name %in% want_scafs_Tpa,]
Tpa_M_11_Tpa_ws_med <- aggregate(Tpa_M_11_Tpa_ws$cov, list(Tpa_M_11_Tpa_ws$scaf_name), FUN=median)
colnames(Tpa_M_11_Tpa_ws_med) <- c("scaf", "Tpa_M_11_Tpa_med_window_cov")
Tpa_M_12_Tpa_ws     <- Tpa_M_12_Tpa[Tpa_M_12_Tpa$scaf_name %in% want_scafs_Tpa,]
Tpa_M_12_Tpa_ws_med <- aggregate(Tpa_M_12_Tpa_ws$cov, list(Tpa_M_12_Tpa_ws$scaf_name), FUN=median)
colnames(Tpa_M_12_Tpa_ws_med) <- c("scaf", "Tpa_M_12_Tpa_med_window_cov")

match_dfs <- c("Tpa_F_H54_ws_med", "Tpa_F_H56_ws_med", "Tpa_F_PA_CD_ws_med", "Tpa_F_PA_E_ws_med", "Tpa_F_Pa_AB_ws_med", "Tpa_M_09_Tpa_ws_med", "Tpa_M_10_Tpa_ws_med", "Tpa_M_11_Tpa_ws_med", "Tpa_M_12_Tpa_ws_med")
match_name <- c("","", "", "", "", "", "", "")
Tpa_LRv5a_med_window_cov <- merge_df_by_col(match_dfs,match_name)
write.csv(Tpa_LRv5a_med_window_cov, "Tpa_LRv5a_med_window_cov_to_Tpa.csv")



######
Tpa_F_H54  <- norm_cov(Tpa_F_H54)
Tpa_F_H56  <- norm_cov(Tpa_F_H56)
Tpa_F_PA_CD <- norm_cov(Tpa_F_PA_CD)
Tpa_F_PA_E <- norm_cov(Tpa_F_PA_E)
Tpa_F_Pa_AB <- norm_cov(Tpa_F_Pa_AB)
Tpa_M_09_Tpa   <- norm_cov(Tpa_M_09_Tpa)
Tpa_M_10_Tpa   <- norm_cov(Tpa_M_10_Tpa)
Tpa_M_11_Tpa   <- norm_cov(Tpa_M_11_Tpa)
Tpa_M_12_Tpa   <- norm_cov(Tpa_M_12_Tpa)

Tpa_all_norm <- rbind(
  Tpa_F_H54,
  Tpa_F_H56,
  Tpa_F_PA_CD,
  Tpa_F_PA_E,
  Tpa_F_Pa_AB,
  Tpa_M_09_Tpa,
  Tpa_M_10_Tpa,
  Tpa_M_11_Tpa,
  Tpa_M_12_Tpa)

Tpa_window_legend_sep <- cowplot::get_legend(plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf1")$p2)


Tpa_window_p1 <- plot_grid(
  plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf1")$p1,
  plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf2")$p1,
  plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf3")$p1,
  plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf4")$p1,
  plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf5")$p1,
  plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf6")$p1,
  plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf7")$p1,
  plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf8")$p1,
  plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf9")$p1,
  plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf10")$p1,
  plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf11")$p1,
  plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf12")$p1,
  plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf13")$p1,
  plot_MF_cov_all(Tpa_all_norm, "Tpa_LRv5a_scf14")$p1,
  ncol = 1)


png(filename = "All_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 10, height = 20, units = "in", bg = "white", res = 300)
plot_grid(Tpa_window_p1, Tpa_window_legend_sep, rel_widths  = c(1, 0.2))
dev.off()
getwd() ## where has my plot gone....


cov_plot_hist_Tpa <- function(df, max_cov) {
  df1 <- df
  samp = deparse(substitute(df))
  print(samp)
  p1 <- ggplot(df1) +
    theme_bw() +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf1"),   aes(cov, color="Tpa_LRv5a_scf1"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf2"),   aes(cov, color="Tpa_LRv5a_scf2"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf3"),   aes(cov, color="Tpa_LRv5a_scf3"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf4"),   aes(cov, color="Tpa_LRv5a_scf4"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf5"),   aes(cov, color="Tpa_LRv5a_scf5"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf6"),   aes(cov, color="Tpa_LRv5a_scf6"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf7"),   aes(cov, color="Tpa_LRv5a_scf7"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf8"),   aes(cov, color="Tpa_LRv5a_scf8"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf9"),   aes(cov, color="Tpa_LRv5a_scf9"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf10"),   aes(cov, color="Tpa_LRv5a_scf10"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf11"),   aes(cov, color="Tpa_LRv5a_scf11"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf12"),   aes(cov, color="Tpa_LRv5a_scf12"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf13"),   aes(cov, color="Tpa_LRv5a_scf13"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tpa_LRv5a_scf14"),   aes(cov, color="Tpa_LRv5a_scf14"), size = 1, stat="density") +
    scale_color_manual(name = "Scaf", values = c("Tpa_LRv5a_scf1"    = "red3",
                                                 "Tpa_LRv5a_scf2"    = "#1B9E77",
                                                 "Tpa_LRv5a_scf3"    = "#D95F02",
                                                 "Tpa_LRv5a_scf4"    = "#E7298A",
                                                 "Tpa_LRv5a_scf5"  = "#66A61E",
                                                 "Tpa_LRv5a_scf6"  = "#E6AB02",
                                                 "Tpa_LRv5a_scf7"    = "#A6761D",
                                                 "Tpa_LRv5a_scf8"  = "yellow2",
                                                 "Tpa_LRv5a_scf9"    = "#666666",
                                                 "Tpa_LRv5a_scf10"   = "lightblue",
                                                 "Tpa_LRv5a_scf11" = "royalblue2",
                                                 "Tpa_LRv5a_scf12"   = "darkorchid",
                                                 "Tpa_LRv5a_scf13"   = "grey",
                                                 "Tpa_LRv5a_scf14"   = "black")) + 
    
    xlim(c(0, max_cov)) + 
    ggtitle(samp) +
    xlab("Coverage") 
  
  p2 <- p1 + theme(legend.position="none")
  
  outlist = list("p2" = p2, "p1" = p1)
  return(outlist)	
}

Tpa_legend_sep <- cowplot::get_legend(cov_plot_hist_Tpa(Tpa_F_H54, 30)$p1)

Tpa_cov_plot <- plot_grid(
  cov_plot_hist_Tpa(Tpa_F_H54, 30)$p2,
  cov_plot_hist_Tpa(Tpa_F_H56, 30)$p2,
  cov_plot_hist_Tpa(Tpa_F_PA_CD, 30)$p2,
  cov_plot_hist_Tpa(Tpa_F_PA_E, 30)$p2,
  cov_plot_hist_Tpa(Tpa_F_Pa_AB, 30)$p2,
  cov_plot_hist_Tpa(Tpa_M_09_Tpa, 40)$p2,
  cov_plot_hist_Tpa(Tpa_M_10_Tpa, 40)$p2,
  cov_plot_hist_Tpa(Tpa_M_11_Tpa, 40)$p2,
  cov_plot_hist_Tpa(Tpa_M_12_Tpa, 40)$p2,
  ncol = 2)


png(filename = "Tpa_cov_plot.png", width  = 12, height = 15, units = "in", bg = "white", res = 300)
plot_grid(Tpa_cov_plot, Tpa_legend_sep, rel_widths  = c(1, 0.5))
dev.off()
getwd() ## where has my plot gone....?



################### 
### exclude H56


Tpa_all_norm_mH56 <- rbind(
  Tpa_F_H54,
  Tpa_F_PA_CD,
  Tpa_F_PA_E,
  Tpa_F_Pa_AB,
  Tpa_M_09_Tpa,
  Tpa_M_10_Tpa,
  Tpa_M_11_Tpa,
  Tpa_M_12_Tpa)



Tpa_window_legend_sep_mH56 <- cowplot::get_legend(plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf1")$p2)

Tpa_window_p1_mH56 <- plot_grid(
  plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf1")$p1,
  plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf2")$p1,
  plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf3")$p1,
  plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf4")$p1,
  plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf5")$p1,
  plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf6")$p1,
  plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf7")$p1,
  plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf8")$p1,
  plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf9")$p1,
  plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf10")$p1,
  plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf11")$p1,
  plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf12")$p1,
  plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf13")$p1,
  plot_MF_cov_all(Tpa_all_norm_mH56, "Tpa_LRv5a_scf14")$p1,
  ncol = 1)

png(filename = "All_to_Tpa_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window_mH56.png", width = 10, height = 20, units = "in", bg = "white", res = 300)
plot_grid(Tpa_window_p1_mH56, Tpa_window_legend_sep_mH56, rel_widths  = c(1, 0.2))
dev.off()
getwd() ## where has my plot gone....

Tpa_cov_plot_mH56 <- plot_grid(
  cov_plot_hist_Tpa(Tpa_F_H54, 30)$p2,
  #cov_plot_hist_Tpa(Tpa_F_H56, 30)$p2,
  cov_plot_hist_Tpa(Tpa_F_PA_CD, 30)$p2,
  cov_plot_hist_Tpa(Tpa_F_PA_E, 30)$p2,
  cov_plot_hist_Tpa(Tpa_F_Pa_AB, 30)$p2,
  cov_plot_hist_Tpa(Tpa_M_09_Tpa, 40)$p2,
  cov_plot_hist_Tpa(Tpa_M_10_Tpa, 40)$p2,
  cov_plot_hist_Tpa(Tpa_M_11_Tpa, 40)$p2,
  cov_plot_hist_Tpa(Tpa_M_12_Tpa, 40)$p2,
  ncol = 2)


png(filename = "Tpa_cov_plot_mH56.png", width  = 12, height = 15, units = "in", bg = "white", res = 300)
plot_grid(Tpa_cov_plot_mH56, Tpa_legend_sep, rel_widths  = c(1, 0.5))
dev.off()
getwd() ## where has my plot gone....?



################################################################################################################
#### class windows X A
#### note H56 excluded

head(Tpa_F_H54)

Tpa_all_norm_c <- as.data.frame(cbind(Tpa_F_H54$scaf_name, 
                                      Tpa_F_H54$start, 
                                      Tpa_F_H54$end, 
                                      Tpa_F_H54$mid_pos, 
                                      Tpa_F_H54$cov_n, 
                                      Tpa_F_PA_CD$cov_n, 
                                      Tpa_F_PA_E$cov_n, 
                                      Tpa_F_Pa_AB$cov_n, 
                                      Tpa_M_09_Tpa$cov_n, 
                                      Tpa_M_10_Tpa$cov_n, 
                                      Tpa_M_11_Tpa$cov_n, 
                                      Tpa_M_12_Tpa$cov_n))



head(Tpa_all_norm_c )



colnames(Tpa_all_norm_c) <- c("scaf_name", "start", "end", "mid_pos", "Tpa_F_H54_cov_n", "Tpa_F_PA_CD_cov_n", "Tpa_F_PA_E_cov_n", "Tpa_F_Pa_AB_cov_n", "Tpa_M_09_Tpa_cov_n", "Tpa_M_10_Tpa_cov_n", "Tpa_M_11_Tpa_cov_n", "Tpa_M_12_Tpa_cov_n")

Tpa_all_norm_c$Tpa_F_H54_cov_n  <- as.numeric(Tpa_all_norm_c$Tpa_F_H54_cov_n)
Tpa_all_norm_c$Tpa_F_PA_CD_cov_n  <- as.numeric(Tpa_all_norm_c$Tpa_F_PA_CD_cov_n)
Tpa_all_norm_c$Tpa_F_PA_E_cov_n  <- as.numeric(Tpa_all_norm_c$Tpa_F_PA_E_cov_n)
Tpa_all_norm_c$Tpa_F_Pa_AB_cov_n  <- as.numeric(Tpa_all_norm_c$Tpa_F_Pa_AB_cov_n)
Tpa_all_norm_c$Tpa_M_09_Tpa_cov_n  <- as.numeric(Tpa_all_norm_c$Tpa_M_09_Tpa_cov_n)
Tpa_all_norm_c$Tpa_M_10_Tpa_cov_n <- as.numeric(Tpa_all_norm_c$Tpa_M_10_Tpa_cov_n)
Tpa_all_norm_c$Tpa_M_11_Tpa_cov_n   <- as.numeric(Tpa_all_norm_c$Tpa_M_11_Tpa_cov_n)
Tpa_all_norm_c$Tpa_M_12_Tpa_cov_n <- as.numeric(Tpa_all_norm_c$Tpa_M_12_Tpa_cov_n)

Tpa_all_norm_c$female_cov_n = 
  (Tpa_all_norm_c$Tpa_F_H54_cov_n + 
     Tpa_all_norm_c$Tpa_F_PA_CD_cov_n +
     Tpa_all_norm_c$Tpa_F_PA_E_cov_n +
     Tpa_all_norm_c$Tpa_F_Pa_AB_cov_n) / 4

Tpa_all_norm_c$male_cov_n = 
  (Tpa_all_norm_c$Tpa_M_11_Tpa_cov_n +
     Tpa_all_norm_c$Tpa_M_12_Tpa_cov_n+
     Tpa_all_norm_c$Tpa_M_09_Tpa_cov_n +
     Tpa_all_norm_c$Tpa_M_10_Tpa_cov_n) / 4

Tpa_all_norm_c$MF <- log2(Tpa_all_norm_c$male_cov_n / Tpa_all_norm_c$female_cov_n)

Tpa_all_norm_c_cut <- subset(Tpa_all_norm_c, Tpa_all_norm_c$MF < -0.5) 
Tpa_Sex_chr_peak <- peakfinder(Tpa_all_norm_c_cut$MF)
Tpa_Auto_peak    <- peakfinder(Tpa_all_norm_c$MF)

Tpa_Sex_chr_peak
Tpa_Auto_peak    


Tpa_all_norm_c$XA_s <- ifelse(Tpa_all_norm_c$MF < Tpa_Sex_chr_peak + 0.1 & Tpa_all_norm_c$MF > Tpa_Sex_chr_peak - 0.1, "X", "A")
Tpa_all_norm_c$XA_l <- ifelse(Tpa_all_norm_c$MF < Tpa_Auto_peak - 0.5, "X", "A")

write.csv(Tpa_all_norm_c, "Tpa_all_norm_c_windows_100000_to_Tpa.csv", row.names=FALSE)

png(filename = "Tpa_all_norm_c_windows_100000_hist_to_Tpa.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
hist(Tpa_all_norm_c$MF, breaks = 400, xlim = c(-2,2))
abline(v=Tpa_Sex_chr_peak, col='red', lwd=1)
abline(v=Tpa_Auto_peak, col='blue', lwd=1)
dev.off()
getwd() ## where has my plot gone....


##################### without the scafs I WANT (i.e. the big groups)

head(Tpa_all_norm_c)
Tpa_all_norm_c_excluded <- Tpa_all_norm_c[!Tpa_all_norm_c$scaf_name %in% want_scafs_Tpa,]  

png(filename = "Tpa_all_norm_c_excluded_windows_100000_hist_to_Tpa.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
hist(Tpa_all_norm_c_excluded$MF, breaks = 400, xlim = c(-2,2))
abline(v=Tpa_Sex_chr_peak, col='red', lwd=1)
abline(v=Tpa_Auto_peak, col='blue', lwd=1)
dev.off()
getwd() ## where has my plot gone....


###########################################################################################################################################################
##### plot MF by scaf dotplot


### need to get offset first

Tpa_chr_lens_s <- Tpa_chr_lens[gtools::mixedorder(Tpa_chr_lens$scaf_name), ]
head(Tpa_chr_lens_s, n = 30)

offset_v <- c(0)
c_total <- 0
for(i in seq(1, length(Tpa_chr_lens_s[,1]) -1)){
  c_len <- Tpa_chr_lens_s[i,]$end
  print(c_len)
  c_total <- c_total + c_len
  offset_v <- c(offset_v, c_total)
  
}

Tpa_chr_lens_s$offset <- offset_v

## add to dict

Tpa_offset_dict <- hash()
for(i in seq(1:length(Tpa_chr_lens_s[,1]))){
  chr_n <- Tpa_chr_lens_s$scaf_name[i]
  offset_n <- Tpa_chr_lens_s$offset[i]
  Tpa_offset_dict [[chr_n]] <-  offset_n
}

## add to cov data
Tpa_chr_offset <- c()
for(i in seq(1:length(Tpa_all_norm_c[,1]))){
  scaf_n <- Tpa_all_norm_c$scaf_name[i]
  off_n <- Tpa_offset_dict[[scaf_n]]
  #### returns a NULL if key not found. Turn into NA
  if(length(off_n) == 0){
    off_n = NA
  }
  print(i)
  print(scaf_n)
  print(off_n)
  Tpa_chr_offset <- c(Tpa_chr_offset, off_n)
}

Tpa_all_norm_c$mid_pos <- as.numeric(Tpa_all_norm_c$mid_pos)

Tpa_all_norm_c <- cbind(Tpa_all_norm_c, Tpa_chr_offset)
Tpa_all_norm_c$genome_pos <- Tpa_all_norm_c$mid_pos + Tpa_all_norm_c$Tpa_chr_offset

Tpa_all_norm_c$scaf_class_1 <- ifelse(Tpa_all_norm_c$scaf_name == "Tpa_LRv5a_scf1", Tpa_all_norm_c$scaf_name,
                                      ifelse(Tpa_all_norm_c$scaf_name == "Tpa_LRv5a_scf2", Tpa_all_norm_c$scaf_name,
                                             ifelse(Tpa_all_norm_c$scaf_name == "Tpa_LRv5a_scf3", Tpa_all_norm_c$scaf_name,
                                                    ifelse(Tpa_all_norm_c$scaf_name == "Tpa_LRv5a_scf4", Tpa_all_norm_c$scaf_name,
                                                           ifelse(Tpa_all_norm_c$scaf_name == "Tpa_LRv5a_scf5", Tpa_all_norm_c$scaf_name,
                                                                  ifelse(Tpa_all_norm_c$scaf_name == "Tpa_LRv5a_scf6", Tpa_all_norm_c$scaf_name,
                                                                         ifelse(Tpa_all_norm_c$scaf_name == "Tpa_LRv5a_scf7", Tpa_all_norm_c$scaf_name,
                                                                                ifelse(Tpa_all_norm_c$scaf_name == "Tpa_LRv5a_scf8", Tpa_all_norm_c$scaf_name,
                                                                                       ifelse(Tpa_all_norm_c$scaf_name == "Tpa_LRv5a_scf9", Tpa_all_norm_c$scaf_name,
                                                                                              ifelse(Tpa_all_norm_c$scaf_name == "Tpa_LRv5a_scf10", Tpa_all_norm_c$scaf_name,
                                                                                                     ifelse(Tpa_all_norm_c$scaf_name == "Tpa_LRv5a_scf11", Tpa_all_norm_c$scaf_name,
                                                                                                            ifelse(Tpa_all_norm_c$scaf_name == "Tpa_LRv5a_scf12", Tpa_all_norm_c$scaf_name,
                                                                                                                   ifelse(Tpa_all_norm_c$scaf_name == "Tpa_LRv5a_scf13", Tpa_all_norm_c$scaf_name,
                                                                                                                          ifelse(Tpa_all_norm_c$scaf_name == "Tpa_LRv5a_scf14", Tpa_all_norm_c$scaf_name,"other"))))))))))))))

Tpa_all_norm_c$scaf_class_1o <- ordered(Tpa_all_norm_c$scaf_class_1, levels= c("Tpa_LRv5a_scf1", "Tpa_LRv5a_scf2", "Tpa_LRv5a_scf3", "Tpa_LRv5a_scf4", "Tpa_LRv5a_scf5", "Tpa_LRv5a_scf6", "Tpa_LRv5a_scf7", "Tpa_LRv5a_scf8", "Tpa_LRv5a_scf9", "Tpa_LRv5a_scf10","Tpa_LRv5a_scf11", "Tpa_LRv5a_scf12", "Tpa_LRv5a_scf13", "Tpa_LRv5a_scf14", "other"))


head(Tpa_all_norm_c)
tail(Tpa_all_norm_c)


## all
ggplot(Tpa_all_norm_c, aes(x=genome_pos, MF, col = scaf_class_1o))  + geom_point(size = 0.5) + 
  theme_bw() +
  scale_color_manual(values=c("firebrick3", "black", "firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black", "firebrick3", "black", "firebrick3")) + ##### alter colors manually
  ylim(c(-2, 2))

## drop 'other'


Tpa_all_norm_c_LGs <- subset(Tpa_all_norm_c, Tpa_all_norm_c$scaf_class_1 != "other")



png(filename = "Tpa_cov_dotplot_LGs_to_Tpa.png", width  = 12, height = 7, units = "in", bg = "white", res = 300)
plot_grid(Tpa_cov_plot, Tpa_legend_sep, rel_widths  = c(1, 0.5))
ggplot(Tpa_all_norm_c_LGs, aes(x=genome_pos, MF, col = scaf_class_1o))  + geom_point(size = 0.5) + 
  theme_bw() +
  scale_color_manual(values=c("firebrick3", "black", "firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black", "firebrick3","black", "firebrick3")) + ##### alter colors manually
  ylim(c(-2, 2)) + ylab("log2(male:female coverage)") + xlab("Genome position (bp)") + theme(legend.position = "none", axis.text=element_text(size=14),  axis.title=element_text(size=14,face="bold")) + ggtitle("Tpa to Tpa")
dev.off()
getwd() ## where has my plot gone....?


dev.off()






















####################################################################################################################################################################################
###################################################################################################################################################################################
## Tps
## Notes - good genome - X = scf3, all males look correct
## ### 2N = 24 (Tps)


Tps_F_ReSeq_Ps08      <- read.table("Tps_F_ReSeq_Ps08_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tps_F_ReSeq_Ps08$mid_pos  <-  (Tps_F_ReSeq_Ps08$end + Tps_F_ReSeq_Ps08$start) / 2
Tps_F_ReSeq_Ps12      <- read.table("Tps_F_ReSeq_Ps12_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tps_F_ReSeq_Ps12$mid_pos  <-  (Tps_F_ReSeq_Ps12$end + Tps_F_ReSeq_Ps12$start) / 2
Tps_F_ReSeq_Ps14      <- read.table("Tps_F_ReSeq_Ps14_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tps_F_ReSeq_Ps14$mid_pos  <-  (Tps_F_ReSeq_Ps14$end + Tps_F_ReSeq_Ps14$start) / 2
Tps_F_ReSeq_Ps16      <- read.table("Tps_F_ReSeq_Ps16_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tps_F_ReSeq_Ps16$mid_pos  <-  (Tps_F_ReSeq_Ps16$end + Tps_F_ReSeq_Ps16$start) / 2
Tps_F_ReSeq_Ps18      <- read.table("Tps_F_ReSeq_Ps18_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tps_F_ReSeq_Ps18$mid_pos  <-  (Tps_F_ReSeq_Ps18$end + Tps_F_ReSeq_Ps18$start) / 2
Tps_M_17_HM99      <- read.table("Tps_M_17_HM99_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tps_M_17_HM99$mid_pos  <-  (Tps_M_17_HM99$end + Tps_M_17_HM99$start) / 2
Tps_M_18_HM100      <- read.table("Tps_M_18_HM100_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tps_M_18_HM100$mid_pos  <-  (Tps_M_18_HM100$end + Tps_M_18_HM100$start) / 2
Tps_M_19_HM101      <- read.table("Tps_M_19_HM101_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tps_M_19_HM101$mid_pos  <-  (Tps_M_19_HM101$end + Tps_M_19_HM101$start) / 2
Tps_M_20_15255      <- read.table("Tps_M_20_15255_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt", sep = "\t", header = T)
Tps_M_20_15255$mid_pos  <-  (Tps_M_20_15255$end + Tps_M_20_15255$start) / 2


### add some info
Tps_F_ReSeq_Ps08$samp_name   <- rep("Tps_F_ReSeq_Ps08",  length(Tps_F_ReSeq_Ps08[,1]))
Tps_F_ReSeq_Ps12$samp_name   <- rep("Tps_F_ReSeq_Ps12",  length(Tps_F_ReSeq_Ps12[,1]))
Tps_F_ReSeq_Ps14$samp_name  <- rep("Tps_F_ReSeq_Ps14", length(Tps_F_ReSeq_Ps14[,1]))
Tps_F_ReSeq_Ps16$samp_name  <- rep("Tps_F_ReSeq_Ps16", length(Tps_F_ReSeq_Ps16[,1]))
Tps_F_ReSeq_Ps18$samp_name  <- rep("Tps_F_ReSeq_Ps18", length(Tps_F_ReSeq_Ps18[,1]))
Tps_M_17_HM99$samp_name    <- rep("Tps_M_17_HM99",   length(Tps_M_17_HM99[,1]))
Tps_M_18_HM100$samp_name    <- rep("Tps_M_18_HM100",   length(Tps_M_18_HM100[,1]))
Tps_M_19_HM101$samp_name    <- rep("Tps_M_19_HM101",   length(Tps_M_19_HM101[,1]))
Tps_M_20_15255$samp_name    <- rep("Tps_M_20_15255",   length(Tps_M_20_15255[,1]))

Tps_F_ReSeq_Ps08$sex   <- rep("F",  length(Tps_F_ReSeq_Ps08[,1]))
Tps_F_ReSeq_Ps12$sex   <- rep("F",  length(Tps_F_ReSeq_Ps12[,1]))
Tps_F_ReSeq_Ps14$sex  <- rep("F", length(Tps_F_ReSeq_Ps14[,1]))
Tps_F_ReSeq_Ps16$sex  <- rep("F", length(Tps_F_ReSeq_Ps16[,1]))
Tps_F_ReSeq_Ps18$sex  <- rep("F", length(Tps_F_ReSeq_Ps18[,1]))
Tps_M_17_HM99$sex    <- rep("M",   length(Tps_M_17_HM99[,1]))
Tps_M_18_HM100$sex    <- rep("M",   length(Tps_M_18_HM100[,1]))
Tps_M_19_HM101$sex    <- rep("M",   length(Tps_M_19_HM101[,1]))
Tps_M_20_15255$sex    <- rep("M",   length(Tps_M_20_15255[,1]))


#############################################################################################
### chr lens

Tps_chr_lens <- aggregate(end~scaf_name, Tps_M_20_15255, FUN=max)
head(Tps_chr_lens[order(Tps_chr_lens$end,decreasing=T),], n = 20)
### big drop after scf12


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

Tps_F_ReSeq_Ps08_ws     <- Tps_F_ReSeq_Ps08[Tps_F_ReSeq_Ps08$scaf_name %in% want_scafs_Tps,]
Tps_F_ReSeq_Ps08_ws_med <- aggregate(Tps_F_ReSeq_Ps08_ws$cov, list(Tps_F_ReSeq_Ps08_ws$scaf_name), FUN=median)
colnames(Tps_F_ReSeq_Ps08_ws_med) <- c("scaf", "Tps_F_ReSeq_Ps08_med_window_cov")
Tps_F_ReSeq_Ps12_ws     <- Tps_F_ReSeq_Ps12[Tps_F_ReSeq_Ps12$scaf_name %in% want_scafs_Tps,]
Tps_F_ReSeq_Ps12_ws_med <- aggregate(Tps_F_ReSeq_Ps12_ws$cov, list(Tps_F_ReSeq_Ps12_ws$scaf_name), FUN=median)
colnames(Tps_F_ReSeq_Ps12_ws_med) <- c("scaf", "Tps_F_ReSeq_Ps12_med_window_cov")
Tps_F_ReSeq_Ps14_ws     <- Tps_F_ReSeq_Ps14[Tps_F_ReSeq_Ps14$scaf_name %in% want_scafs_Tps,]
Tps_F_ReSeq_Ps14_ws_med <- aggregate(Tps_F_ReSeq_Ps14_ws$cov, list(Tps_F_ReSeq_Ps14_ws$scaf_name), FUN=median)
colnames(Tps_F_ReSeq_Ps14_ws_med) <- c("scaf", "Tps_F_ReSeq_Ps14_med_window_cov")
Tps_F_ReSeq_Ps16_ws     <- Tps_F_ReSeq_Ps16[Tps_F_ReSeq_Ps16$scaf_name %in% want_scafs_Tps,]
Tps_F_ReSeq_Ps16_ws_med <- aggregate(Tps_F_ReSeq_Ps16_ws$cov, list(Tps_F_ReSeq_Ps16_ws$scaf_name), FUN=median)
colnames(Tps_F_ReSeq_Ps16_ws_med) <- c("scaf", "Tps_F_ReSeq_Ps16_med_window_cov")
Tps_F_ReSeq_Ps18_ws     <- Tps_F_ReSeq_Ps18[Tps_F_ReSeq_Ps18$scaf_name %in% want_scafs_Tps,]
Tps_F_ReSeq_Ps18_ws_med <- aggregate(Tps_F_ReSeq_Ps18_ws$cov, list(Tps_F_ReSeq_Ps18_ws$scaf_name), FUN=median)
colnames(Tps_F_ReSeq_Ps18_ws_med) <- c("scaf", "Tps_F_ReSeq_Ps18_med_window_cov")
Tps_M_17_HM99_ws     <- Tps_M_17_HM99[Tps_M_17_HM99$scaf_name %in% want_scafs_Tps,]
Tps_M_17_HM99_ws_med <- aggregate(Tps_M_17_HM99_ws$cov, list(Tps_M_17_HM99_ws$scaf_name), FUN=median)
colnames(Tps_M_17_HM99_ws_med) <- c("scaf", "Tps_M_17_HM99_med_window_cov")
Tps_M_18_HM100_ws     <- Tps_M_18_HM100[Tps_M_18_HM100$scaf_name %in% want_scafs_Tps,]
Tps_M_18_HM100_ws_med <- aggregate(Tps_M_18_HM100_ws$cov, list(Tps_M_18_HM100_ws$scaf_name), FUN=median)
colnames(Tps_M_18_HM100_ws_med) <- c("scaf", "Tps_M_18_HM100_med_window_cov")
Tps_M_19_HM101_ws     <- Tps_M_19_HM101[Tps_M_19_HM101$scaf_name %in% want_scafs_Tps,]
Tps_M_19_HM101_ws_med <- aggregate(Tps_M_19_HM101_ws$cov, list(Tps_M_19_HM101_ws$scaf_name), FUN=median)
colnames(Tps_M_19_HM101_ws_med) <- c("scaf", "Tps_M_19_HM101_med_window_cov")
Tps_M_20_15255_ws     <- Tps_M_20_15255[Tps_M_20_15255$scaf_name %in% want_scafs_Tps,]
Tps_M_20_15255_ws_med <- aggregate(Tps_M_20_15255_ws$cov, list(Tps_M_20_15255_ws$scaf_name), FUN=median)
colnames(Tps_M_20_15255_ws_med) <- c("scaf", "Tps_M_20_15255_med_window_cov")

match_dfs <- c("Tps_F_ReSeq_Ps08_ws_med", "Tps_F_ReSeq_Ps12_ws_med", "Tps_F_ReSeq_Ps14_ws_med", "Tps_F_ReSeq_Ps16_ws_med", "Tps_F_ReSeq_Ps18_ws_med", "Tps_M_17_HM99_ws_med", "Tps_M_18_HM100_ws_med", "Tps_M_19_HM101_ws_med", "Tps_M_20_15255_ws_med")
match_name <- c("","", "", "", "", "", "", "")
Tps_LRv5b_med_window_cov <- merge_df_by_col(match_dfs,match_name)
write.csv(Tps_LRv5b_med_window_cov, "Tps_LRv5b_med_window_cov_to_Tps.csv")



######
Tps_F_ReSeq_Ps08  <- norm_cov(Tps_F_ReSeq_Ps08)
Tps_F_ReSeq_Ps12  <- norm_cov(Tps_F_ReSeq_Ps12)
Tps_F_ReSeq_Ps14 <- norm_cov(Tps_F_ReSeq_Ps14)
Tps_F_ReSeq_Ps16 <- norm_cov(Tps_F_ReSeq_Ps16)
Tps_F_ReSeq_Ps18 <- norm_cov(Tps_F_ReSeq_Ps18)
Tps_M_17_HM99   <- norm_cov(Tps_M_17_HM99)
Tps_M_18_HM100   <- norm_cov(Tps_M_18_HM100)
Tps_M_19_HM101   <- norm_cov(Tps_M_19_HM101)
Tps_M_20_15255   <- norm_cov(Tps_M_20_15255)

Tps_all_norm <- rbind(
  Tps_F_ReSeq_Ps08,
  Tps_F_ReSeq_Ps12,
  Tps_F_ReSeq_Ps14,
  Tps_F_ReSeq_Ps16,
  Tps_F_ReSeq_Ps18,
  Tps_M_17_HM99,
  Tps_M_18_HM100,
  Tps_M_19_HM101,
  Tps_M_20_15255)

Tps_window_legend_sep <- cowplot::get_legend(plot_MF_cov_all(Tps_all_norm, "Tps_LRv5b_scf1")$p2)

Tps_window_p1 <- plot_grid(
  plot_MF_cov_all(Tps_all_norm, "Tps_LRv5b_scf1")$p1,
  plot_MF_cov_all(Tps_all_norm, "Tps_LRv5b_scf2")$p1,
  plot_MF_cov_all(Tps_all_norm, "Tps_LRv5b_scf3")$p1,
  plot_MF_cov_all(Tps_all_norm, "Tps_LRv5b_scf4")$p1,
  plot_MF_cov_all(Tps_all_norm, "Tps_LRv5b_scf5")$p1,
  plot_MF_cov_all(Tps_all_norm, "Tps_LRv5b_scf6")$p1,
  plot_MF_cov_all(Tps_all_norm, "Tps_LRv5b_scf7")$p1,
  plot_MF_cov_all(Tps_all_norm, "Tps_LRv5b_scf8")$p1,
  plot_MF_cov_all(Tps_all_norm, "Tps_LRv5b_scf9")$p1,
  plot_MF_cov_all(Tps_all_norm, "Tps_LRv5b_scf10")$p1,
  plot_MF_cov_all(Tps_all_norm, "Tps_LRv5b_scf11")$p1,
  plot_MF_cov_all(Tps_all_norm, "Tps_LRv5b_scf12")$p1,
  ncol = 1)

png(filename = "All_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window.png", width = 10, height = 20, units = "in", bg = "white", res = 300)
plot_grid(Tps_window_p1, Tps_window_legend_sep, rel_widths  = c(1, 0.2))
dev.off()
getwd() ## where has my plot gone....


cov_plot_hist_Tps <- function(df, max_cov) {
  df1 <- df
  samp = deparse(substitute(df))
  print(samp)
  p1 <- ggplot(df1) +
    theme_bw() +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf1"),   aes(cov, color="Tps_LRv5b_scf1"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf2"),   aes(cov, color="Tps_LRv5b_scf2"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf3"),   aes(cov, color="Tps_LRv5b_scf3"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf4"),   aes(cov, color="Tps_LRv5b_scf4"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf5"),   aes(cov, color="Tps_LRv5b_scf5"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf6"),   aes(cov, color="Tps_LRv5b_scf6"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf7"),   aes(cov, color="Tps_LRv5b_scf7"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf8"),   aes(cov, color="Tps_LRv5b_scf8"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf9"),   aes(cov, color="Tps_LRv5b_scf9"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf10"),   aes(cov, color="Tps_LRv5b_scf10"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf11"),   aes(cov, color="Tps_LRv5b_scf11"), size = 1, stat="density") +
    geom_line(data=subset(df1,scaf_name == "Tps_LRv5b_scf12"),   aes(cov, color="Tps_LRv5b_scf12"), size = 1, stat="density") +
    scale_color_manual(name = "Scaf", values = c("Tps_LRv5b_scf1"  = "#D95F02",
                                                 "Tps_LRv5b_scf2"    = "#1B9E77",
                                                 "Tps_LRv5b_scf3"    = "red3",
                                                 "Tps_LRv5b_scf4"    = "#E7298A",
                                                 "Tps_LRv5b_scf5"  = "#66A61E",
                                                 "Tps_LRv5b_scf6"  = "#E6AB02",
                                                 "Tps_LRv5b_scf7"    = "#A6761D",
                                                 "Tps_LRv5b_scf8"  = "yellow2",
                                                 "Tps_LRv5b_scf9"    = "#666666",
                                                 "Tps_LRv5b_scf10"   = "lightblue",
                                                 "Tps_LRv5b_scf11" = "royalblue2",
                                                 "Tps_LRv5b_scf12"   = "darkorchid")) + 
    
    xlim(c(0, max_cov)) + 
    ggtitle(samp) +
    xlab("Coverage") 
  
  p2 <- p1 + theme(legend.position="none")
  
  outlist = list("p2" = p2, "p1" = p1)
  return(outlist)	
}

Tps_legend_sep <- cowplot::get_legend(cov_plot_hist_Tps(Tps_F_ReSeq_Ps08, 30)$p1)

Tps_cov_plot <- plot_grid(
  cov_plot_hist_Tps(Tps_F_ReSeq_Ps08, 30)$p2,
  cov_plot_hist_Tps(Tps_F_ReSeq_Ps12, 30)$p2,
  cov_plot_hist_Tps(Tps_F_ReSeq_Ps14, 30)$p2,
  cov_plot_hist_Tps(Tps_F_ReSeq_Ps16, 30)$p2,
  cov_plot_hist_Tps(Tps_F_ReSeq_Ps18, 30)$p2,
  cov_plot_hist_Tps(Tps_M_17_HM99, 40)$p2,
  cov_plot_hist_Tps(Tps_M_18_HM100, 40)$p2,
  cov_plot_hist_Tps(Tps_M_19_HM101, 40)$p2,
  cov_plot_hist_Tps(Tps_M_20_15255, 40)$p2,
  ncol = 2)


png(filename = "Tps_cov_plot.png", width  = 12, height = 15, units = "in", bg = "white", res = 300)
plot_grid(Tps_cov_plot, Tps_legend_sep, rel_widths  = c(1, 0.5))
dev.off()
getwd() ## where has my plot gone....?



################### 
### exclude ReSeq Ps08 and Reseq_Ps12


Tps_all_norm_e0812 <- rbind(
  Tps_F_ReSeq_Ps14,
  Tps_F_ReSeq_Ps16,
  Tps_F_ReSeq_Ps18,
  Tps_M_17_HM99,
  Tps_M_18_HM100,
  Tps_M_19_HM101,
  Tps_M_20_15255)


Tps_window_legend_sep_e0812 <- cowplot::get_legend(plot_MF_cov_all(Tps_all_norm_e0812, "Tps_LRv5b_scf1")$p2)

Tps_window_p1_e0812 <- plot_grid(
  plot_MF_cov_all(Tps_all_norm_e0812, "Tps_LRv5b_scf1")$p1,
  plot_MF_cov_all(Tps_all_norm_e0812, "Tps_LRv5b_scf2")$p1,
  plot_MF_cov_all(Tps_all_norm_e0812, "Tps_LRv5b_scf3")$p1,
  plot_MF_cov_all(Tps_all_norm_e0812, "Tps_LRv5b_scf4")$p1,
  plot_MF_cov_all(Tps_all_norm_e0812, "Tps_LRv5b_scf5")$p1,
  plot_MF_cov_all(Tps_all_norm_e0812, "Tps_LRv5b_scf6")$p1,
  plot_MF_cov_all(Tps_all_norm_e0812, "Tps_LRv5b_scf7")$p1,
  plot_MF_cov_all(Tps_all_norm_e0812, "Tps_LRv5b_scf8")$p1,
  plot_MF_cov_all(Tps_all_norm_e0812, "Tps_LRv5b_scf9")$p1,
  plot_MF_cov_all(Tps_all_norm_e0812, "Tps_LRv5b_scf10")$p1,
  plot_MF_cov_all(Tps_all_norm_e0812, "Tps_LRv5b_scf11")$p1,
  plot_MF_cov_all(Tps_all_norm_e0812, "Tps_LRv5b_scf12")$p1,
  ncol = 1)

png(filename = "All_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window_e0812.png", width = 10, height = 20, units = "in", bg = "white", res = 300)
plot_grid(Tps_window_p1_e0812, Tps_window_legend_sep_e0812, rel_widths  = c(1, 0.2))
dev.off()
getwd() ## where has my plot gone....


plot_MF_cov_all_Tpsmfcolour <- function(df, chr_want){
  df_chr <- subset(df, df$scaf_name == chr_want)
  print(head(df_chr))
  ggplot(df_chr, aes(x=mid_pos, cov_n, col = samp_name)) +
    theme_bw() +
    geom_line(alpha=0.5) +
    ggtitle(paste("scaf", chr_want)) + 
    scale_colour_manual(values=c("red2", "red3", "red1", "royalblue1", "royalblue2", "royalblue3","royalblue4", "lightblue2")) +
    xlab("Position") + ylim(0,3)
}


png(filename = "All_to_Tps_LRv5b_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_100000_window_e0812_mfcolour.png", width = 20, height = 40, units = "in", bg = "white", res = 300)
plot_grid(
  plot_MF_cov_all_Tpsmfcolour(Tps_all_norm_e0812, "Tps_LRv5b_scf1"),
  plot_MF_cov_all_Tpsmfcolour(Tps_all_norm_e0812, "Tps_LRv5b_scf2"),
  plot_MF_cov_all_Tpsmfcolour(Tps_all_norm_e0812, "Tps_LRv5b_scf3"), 
  plot_MF_cov_all_Tpsmfcolour(Tps_all_norm_e0812, "Tps_LRv5b_scf4"), 
  plot_MF_cov_all_Tpsmfcolour(Tps_all_norm_e0812, "Tps_LRv5b_scf5"), 
  plot_MF_cov_all_Tpsmfcolour(Tps_all_norm_e0812, "Tps_LRv5b_scf6"), 
  plot_MF_cov_all_Tpsmfcolour(Tps_all_norm_e0812, "Tps_LRv5b_scf7"), 
  plot_MF_cov_all_Tpsmfcolour(Tps_all_norm_e0812, "Tps_LRv5b_scf8"), 
  plot_MF_cov_all_Tpsmfcolour(Tps_all_norm_e0812, "Tps_LRv5b_scf9"), 
  plot_MF_cov_all_Tpsmfcolour(Tps_all_norm_e0812, "Tps_LRv5b_scf10"), 
  plot_MF_cov_all_Tpsmfcolour(Tps_all_norm_e0812, "Tps_LRv5b_scf11"), 
  plot_MF_cov_all_Tpsmfcolour(Tps_all_norm_e0812, "Tps_LRv5b_scf12"), 
  ncol = 1)
dev.off()
getwd() ## where has my plot gone....











Tps_cov_plot_e0812 <- plot_grid(
  cov_plot_hist_Tps(Tps_F_ReSeq_Ps14, 30)$p2,
  cov_plot_hist_Tps(Tps_F_ReSeq_Ps16, 30)$p2,
  cov_plot_hist_Tps(Tps_F_ReSeq_Ps18, 30)$p2,
  cov_plot_hist_Tps(Tps_M_17_HM99, 40)$p2,
  cov_plot_hist_Tps(Tps_M_18_HM100, 40)$p2,
  cov_plot_hist_Tps(Tps_M_19_HM101, 40)$p2,
  cov_plot_hist_Tps(Tps_M_20_15255, 40)$p2,
  ncol = 2)


png(filename = "Tps_cov_plot_e0812.png", width  = 12, height = 15, units = "in", bg = "white", res = 300)
plot_grid(Tps_cov_plot_e0812, Tps_legend_sep, rel_widths  = c(1, 0.5))
dev.off()
getwd() ## where has my plot gone....?






################################################################################################################
#### class windows X A
#### note 08 and 12 excluded

head(Tps_F_ReSeq_Ps14)

Tps_all_norm_c <- as.data.frame(cbind(Tps_F_ReSeq_Ps14$scaf_name, 
                                      Tps_F_ReSeq_Ps14$start, 
                                      Tps_F_ReSeq_Ps14$end, 
                                      Tps_F_ReSeq_Ps14$mid_pos, 
                                      Tps_F_ReSeq_Ps14$cov_n, 
                                      Tps_F_ReSeq_Ps16$cov_n, 
                                      Tps_F_ReSeq_Ps18$cov_n, 
                                      Tps_M_17_HM99$cov_n, 
                                      Tps_M_18_HM100$cov_n, 
                                      Tps_M_19_HM101$cov_n, 
                                      Tps_M_20_15255$cov_n))



head(Tps_all_norm_c )



colnames(Tps_all_norm_c) <- c("scaf_name", "start", "end", "mid_pos", "Tps_F_ReSeq_Ps14_cov_n", "Tps_F_ReSeq_Ps16_cov_n", "Tps_F_ReSeq_Ps18_cov_n", "Tps_M_17_HM99_cov_n", "Tps_M_18_HM100_cov_n", "Tps_M_19_HM101_cov_n", "Tps_M_20_15255_cov_n")

Tps_all_norm_c$Tps_F_ReSeq_Ps14_cov_n  <- as.numeric(Tps_all_norm_c$Tps_F_ReSeq_Ps14_cov_n)
Tps_all_norm_c$Tps_F_ReSeq_Ps16_cov_n  <- as.numeric(Tps_all_norm_c$Tps_F_ReSeq_Ps16_cov_n)
Tps_all_norm_c$Tps_F_ReSeq_Ps18_cov_n  <- as.numeric(Tps_all_norm_c$Tps_F_ReSeq_Ps18_cov_n)
Tps_all_norm_c$Tps_M_17_HM99_cov_n  <- as.numeric(Tps_all_norm_c$Tps_M_17_HM99_cov_n)
Tps_all_norm_c$Tps_M_18_HM100_cov_n  <- as.numeric(Tps_all_norm_c$Tps_M_18_HM100_cov_n)
Tps_all_norm_c$Tps_M_19_HM101_cov_n <- as.numeric(Tps_all_norm_c$Tps_M_19_HM101_cov_n)
Tps_all_norm_c$Tps_M_20_15255_cov_n   <- as.numeric(Tps_all_norm_c$Tps_M_20_15255_cov_n)


Tps_all_norm_c$female_cov_n = 
  (Tps_all_norm_c$Tps_F_ReSeq_Ps14_cov_n + 
     Tps_all_norm_c$Tps_F_ReSeq_Ps16_cov_n +
     Tps_all_norm_c$Tps_F_ReSeq_Ps18_cov_n) / 3

Tps_all_norm_c$male_cov_n = 
  (    Tps_all_norm_c$Tps_M_17_HM99_cov_n + 
    Tps_all_norm_c$Tps_M_20_15255_cov_n +
     Tps_all_norm_c$Tps_M_18_HM100_cov_n +
     Tps_all_norm_c$Tps_M_19_HM101_cov_n) / 4

Tps_all_norm_c$MF <- log2(Tps_all_norm_c$male_cov_n / Tps_all_norm_c$female_cov_n)

Tps_all_norm_c_cut <- subset(Tps_all_norm_c, Tps_all_norm_c$MF < -0.5) 
Tps_Sex_chr_peak <- peakfinder(Tps_all_norm_c_cut$MF)
Tps_Auto_peak    <- peakfinder(Tps_all_norm_c$MF)

Tps_Sex_chr_peak
Tps_Auto_peak    


Tps_all_norm_c$XA_s <- ifelse(Tps_all_norm_c$MF < Tps_Sex_chr_peak + 0.1 & Tps_all_norm_c$MF > Tps_Sex_chr_peak - 0.1, "X", "A")
Tps_all_norm_c$XA_l <- ifelse(Tps_all_norm_c$MF < Tps_Auto_peak - 0.5, "X", "A")

write.csv(Tps_all_norm_c, "Tps_all_norm_c_windows_100000_to_Tps.csv", row.names=FALSE)

png(filename = "Tps_all_norm_c_windows_100000_hist_to_Tps.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
hist(Tps_all_norm_c$MF, breaks = 400, xlim = c(-2,2))
abline(v=Tps_Sex_chr_peak, col='red', lwd=1)
abline(v=Tps_Auto_peak, col='blue', lwd=1)
dev.off()
getwd() ## where has my plot gone....


##################### without the scafs I WANT (i.e. the big groups)

head(Tps_all_norm_c)
Tps_all_norm_c_excluded <- Tps_all_norm_c[!Tps_all_norm_c$scaf_name %in% want_scafs_Tps,]  

png(filename = "Tps_all_norm_c_excluded_windows_100000_hist_to_Tps.png", width = 7, height = 7, units = "in", bg = "white", res = 300)
hist(Tps_all_norm_c_excluded$MF, breaks = 400, xlim = c(-2,2))
abline(v=Tps_Sex_chr_peak, col='red', lwd=1)
abline(v=Tps_Auto_peak, col='blue', lwd=1)
dev.off()
getwd() ## where has my plot gone....



###########################################################################################################################################################
##### plot MF by scaf dotplot


### need to get offset first

Tps_chr_lens_s <- Tps_chr_lens[gtools::mixedorder(Tps_chr_lens$scaf_name), ]
head(Tps_chr_lens_s)

offset_v <- c(0)
c_total <- 0
for(i in seq(1, length(Tps_chr_lens_s[,1]) -1)){
  c_len <- Tps_chr_lens_s[i,]$end
  print(c_len)
  c_total <- c_total + c_len
  offset_v <- c(offset_v, c_total)
  
}

Tps_chr_lens_s$offset <- offset_v

## add to dict

Tps_offset_dict <- hash()
for(i in seq(1:length(Tps_chr_lens_s[,1]))){
  chr_n <- Tps_chr_lens_s$scaf_name[i]
  offset_n <- Tps_chr_lens_s$offset[i]
  Tps_offset_dict [[chr_n]] <-  offset_n
}

## add to cov data
Tps_chr_offset <- c()
for(i in seq(1:length(Tps_all_norm_c[,1]))){
  scaf_n <- Tps_all_norm_c$scaf_name[i]
  off_n <- Tps_offset_dict[[scaf_n]]
  #### returns a NULL if key not found. Turn into NA
  if(length(off_n) == 0){
    off_n = NA
  }
  print(i)
  print(scaf_n)
  print(off_n)
  Tps_chr_offset <- c(Tps_chr_offset, off_n)
}

Tps_all_norm_c$mid_pos <-  as.numeric(Tps_all_norm_c$mid_pos)

Tps_all_norm_c <- cbind(Tps_all_norm_c, Tps_chr_offset)
Tps_all_norm_c$genome_pos <- Tps_all_norm_c$mid_pos + Tps_all_norm_c$Tps_chr_offset

Tps_all_norm_c$scaf_class_1 <- ifelse(Tps_all_norm_c$scaf_name == "Tps_LRv5b_scf1", Tps_all_norm_c$scaf_name,
                                      ifelse(Tps_all_norm_c$scaf_name == "Tps_LRv5b_scf2", Tps_all_norm_c$scaf_name,
                                             ifelse(Tps_all_norm_c$scaf_name == "Tps_LRv5b_scf3", Tps_all_norm_c$scaf_name,
                                                    ifelse(Tps_all_norm_c$scaf_name == "Tps_LRv5b_scf4", Tps_all_norm_c$scaf_name,
                                                           ifelse(Tps_all_norm_c$scaf_name == "Tps_LRv5b_scf5", Tps_all_norm_c$scaf_name,
                                                                  ifelse(Tps_all_norm_c$scaf_name == "Tps_LRv5b_scf6", Tps_all_norm_c$scaf_name,
                                                                         ifelse(Tps_all_norm_c$scaf_name == "Tps_LRv5b_scf7", Tps_all_norm_c$scaf_name,
                                                                                ifelse(Tps_all_norm_c$scaf_name == "Tps_LRv5b_scf8", Tps_all_norm_c$scaf_name,
                                                                                       ifelse(Tps_all_norm_c$scaf_name == "Tps_LRv5b_scf9", Tps_all_norm_c$scaf_name,
                                                                                              ifelse(Tps_all_norm_c$scaf_name == "Tps_LRv5b_scf10", Tps_all_norm_c$scaf_name,
                                                                                                     ifelse(Tps_all_norm_c$scaf_name == "Tps_LRv5b_scf11", Tps_all_norm_c$scaf_name,
                                                                                                            ifelse(Tps_all_norm_c$scaf_name == "Tps_LRv5b_scf12", Tps_all_norm_c$scaf_name,"other"))))))))))))

Tps_all_norm_c$scaf_class_1o <- ordered(Tps_all_norm_c$scaf_class_1, levels= c("Tps_LRv5b_scf1", "Tps_LRv5b_scf2", "Tps_LRv5b_scf3", "Tps_LRv5b_scf4", "Tps_LRv5b_scf5", "Tps_LRv5b_scf6", "Tps_LRv5b_scf7", "Tps_LRv5b_scf8", "Tps_LRv5b_scf9", "Tps_LRv5b_scf10","Tps_LRv5b_scf11", "Tps_LRv5b_scf12", "other"))


head(Tps_all_norm_c)
tail(Tps_all_norm_c)


## all
ggplot(Tps_all_norm_c, aes(x=genome_pos, MF, col = scaf_class_1o))  + geom_point(size = 0.5) + 
  theme_bw() +
  scale_color_manual(values=c("firebrick3", "black", "firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black", "firebrick3")) + ##### alter colors manually
  ylim(c(-2, 2))

## drop 'other'


Tps_all_norm_c_LGs <- subset(Tps_all_norm_c, Tps_all_norm_c$scaf_class_1 != "other")



png(filename = "Tps_cov_dotplot_LGs_to_Tps.png", width  = 12, height = 7, units = "in", bg = "white", res = 300)
plot_grid(Tps_cov_plot, Tps_legend_sep, rel_widths  = c(1, 0.5))
ggplot(Tps_all_norm_c_LGs, aes(x=genome_pos, MF, col = scaf_class_1o))  + geom_point(size = 0.5) + 
  theme_bw() +
  scale_color_manual(values=c("firebrick3", "black", "firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black", "firebrick3")) + ##### alter colors manually
  ylim(c(-2, 2)) + ylab("log2(male:female coverage)") + xlab("Genome position (bp)") + theme(legend.position = "none", axis.text=element_text(size=14),  axis.title=element_text(size=14,face="bold")) + ggtitle("Tps to Tps")
dev.off()
getwd() ## where has my plot gone....?
dev.off()


pdf("Tps_cov_dotplot_LGs_to_Tps.pdf", width  = 12, height = 7, bg = "white")
plot_grid(Tps_cov_plot, Tps_legend_sep, rel_widths  = c(1, 0.5))
ggplot(Tps_all_norm_c_LGs, aes(x=genome_pos, MF, col = scaf_class_1o))  + geom_point(size = 0.5) + 
  theme_bw() +
  scale_color_manual(values=c("firebrick3", "black", "firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black", "firebrick3")) + ##### alter colors manually
  ylim(c(-2, 2)) + ylab("log2(male:female coverage)") + xlab("Genome position (bp)") + theme(legend.position = "none", axis.text=element_text(size=14),  axis.title=element_text(size=14,face="bold")) + ggtitle("Tps to Tps")
dev.off()
getwd() ## where has my plot gone....?
dev.off()


pdf("Tps_cov_dotplot_LGs_to_Tps.pdf", width  = 12, height = 7, bg = "white")
ggplot(Tps_all_norm_c_LGs, aes(x=genome_pos, MF, col = scaf_class_1o))  + geom_point(size = 0.5) + 
  theme_bw() +
  scale_color_manual(values=c("firebrick3", "black", "firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black","firebrick3","black", "firebrick3")) + ##### alter colors manually
  ylim(c(-2, 2)) + ylab("log2(male:female coverage)") + xlab("Genome position (bp)") + theme(legend.position = "none", axis.text=element_text(size=14),  axis.title=element_text(size=14,face="bold")) + ggtitle("Tps to Tps")
dev.off()
getwd() ## where has my plot gone....?
dev.off()

