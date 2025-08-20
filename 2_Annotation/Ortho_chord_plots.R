### chord plots for orths

library(macrosyntR)
library(ggplot2)
library(stringr)

##########################################################################################################################################
### all labels

orths <- read.delim("data/orths/Orthogroups_TpsTpaTce_onscafs_1to1.tsv", header = F)
orths_2 = cbind(orths[,3], orths[,1], orths[,2])
write.table(orths_2, "data/orths/Orthogroups_TpsTpaTce_onscafs_1to1_reorder.tsv", col.names = F, row.names =  F, quote = F, sep = "\t")

my_orthologs_with_3_sp <- load_orthologs("data/orths/Orthogroups_TpsTpaTce_onscafs_1to1_reorder.tsv",
                                         bedfiles = c("data/orths/Orthogroups_TpsTpaTce_onscafs_Tps_1to1info.bed",
                                                      "data/orths/Orthogroups_TpsTpaTce_onscafs_Tce_1to1info_LG.bed",
                                                      "data/orths/Orthogroups_TpsTpaTce_onscafs_Tpa_1to1info.bed"
                                         ))

plot_chord_diagram(my_orthologs_with_3_sp,
                   species_labels = c("Tps","Tce", "Tpa"),
                   color_by = "sp1.Chr", reorder_chromosomes = FALSE) + theme(legend.position = "none")

##########################################################################################################################################
## clean labels

### reorder orths

orths <- read.delim("data/orths/Orthogroups_TpsTpaTce_onscafs_1to1.tsv", header = F)
orths_2 = cbind(orths[,3], orths[,1], orths[,2])
write.table(orths_2, "data/orths/Orthogroups_TpsTpaTce_onscafs_1to1_reorder.tsv", col.names = F, row.names =  F, quote = F, sep = "\t")

my_orthologs_with_3_sp <- load_orthologs("data/orths/Orthogroups_TpsTpaTce_onscafs_1to1_reorder.tsv",
                                         bedfiles = c("data/orths/Orthogroups_TpsTpaTce_onscafs_Tps_c_1to1info_Scf2_Scf3_Scf6_Scf8_scf10.bed",
                                                      "data/orths/Orthogroups_TpsTpaTce_onscafs_Tce_1to1info_LG_LG11_LG12_LG13.bed",
                                                      "data/orths/Orthogroups_TpsTpaTce_onscafs_Tpa_c_1to1info_Scf10.bed"
                                         ))


my_orthologs_with_3_sp$sp1.Chr <- gsub("scf", "", my_orthologs_with_3_sp$sp1.Chr )
my_orthologs_with_3_sp$sp2.Chr <- gsub("LG", "", my_orthologs_with_3_sp$sp2.Chr )
my_orthologs_with_3_sp$sp3.Chr <- gsub("scf", "", my_orthologs_with_3_sp$sp3.Chr )

str(my_orthologs_with_3_sp)

# my_orthologs_with_3_sp$sp1.Chr <- ordered(my_orthologs_with_3_sp$sp1.Chr, levels = c("scf1", "scf2",  "scf3",  "scf4",  "scf5",  "scf6",  "scf7",  "scf8",  "scf9", "scf10", "scf11", "scf12"   ))
# my_orthologs_with_3_sp$sp2.Chr <- ordered(my_orthologs_with_3_sp$sp2.Chr, levels = c("LG6", "LG2",  "LG1", "LG3", "LG4", "LG8",  "LG5",  "LG7",  "LG9", "LG11",  "LG10",  "LG12",  "LG13"  ))
# my_orthologs_with_3_sp$sp3.Chr <- ordered(my_orthologs_with_3_sp$sp3.Chr, levels = c( "scf2", "scf10", "scf5", "scf3", "scf1", "scf4", "scf6",  "scf7", "scf8",  "scf9", "scf11" ,  "scf12",  "scf13" , "scf14" ))

my_orthologs_with_3_sp$sp1.Chr <- ifelse(my_orthologs_with_3_sp$sp1.Chr == "3", "X", my_orthologs_with_3_sp$sp1.Chr)
my_orthologs_with_3_sp$sp2.Chr <- ifelse(my_orthologs_with_3_sp$sp2.Chr == "3", "X", my_orthologs_with_3_sp$sp2.Chr)
my_orthologs_with_3_sp$sp3.Chr <- ifelse(my_orthologs_with_3_sp$sp3.Chr == "1", "X", my_orthologs_with_3_sp$sp3.Chr)

my_orthologs_with_3_sp$sp1.Chr <- ordered(my_orthologs_with_3_sp$sp1.Chr, levels = c("X" , "1", "2",  "4",  "5",  "6",  "7",  "8",  "9", "10", "11", "12"   ))
my_orthologs_with_3_sp$sp2.Chr <- ordered(my_orthologs_with_3_sp$sp2.Chr, levels = c("X" , "6", "2",  "1", "4", "8",  "5",  "7",  "9", "11",  "10",  "12",  "13"  ))
my_orthologs_with_3_sp$sp3.Chr <- ordered(my_orthologs_with_3_sp$sp3.Chr, levels = c("X" , "2", "10", "5", "3", "4", "6",  "7", "8",  "9", "11" ,  "12",  "13" , "14" ))



png(filename = "data/orths/TpsTceTpa_chord_diagram_1to1orths_alllinks_tidynames.png", width = 12, height = 3, units = "in", bg = "white", res = 300)
plot_chord_diagram(my_orthologs_with_3_sp,
                   species_labels = c("T. poppense", "T. cristinae", "T. podura"),
                   color_by = "sp1.Chr", reorder_chromosomes = FALSE, remove_non_linkage_orthologs = FALSE,  label_size = 4) + theme(legend.position = "none")

dev.off()

png(filename = "data/orths/TpsTceTpa_chord_diagram_1to1orths_siglinksonly_tidynames.png", width =14, height = 4, units = "in", bg = "white", res = 300)
plot_chord_diagram(my_orthologs_with_3_sp,
                   species_labels = c("T. poppense", "T. cristinae", "T. podura"),
                   color_by = "sp1.Chr", reorder_chromosomes = FALSE,  label_size = 4) + theme(legend.position = "none")

dev.off()



###### all Timema + bacill
##########################################################################################################################################
## clean labels

### reorder orths

orths <- read.delim("data/orths/Orthogroups_TpsTpaTceBrs_onscafs_1to1.tsv", header = F)
head(orths)
orths_2 = cbind(orths[,4], orths[,2], orths[,3], orths[,1])
write.table(orths_2, "data/orths/Orthogroups_TpsTpaTceBrs_onscafs_1to1_reorder.tsv", col.names = F, row.names =  F, quote = F, sep = "\t")

my_orthologs_with_3_sp <- load_orthologs("data/orths/Orthogroups_TpsTpaTceBrs_onscafs_1to1_reorder.tsv",
                                         bedfiles = c("data/orths/Orthogroups_TpsTpaTceBrs_onscafs_Tps_c_1to1info_Scf2_Scf3_Scf6_Scf8_scf10.bed",
                                                      "data/orths/Orthogroups_TpsTpaTceBrs_onscafs_Tce_1to1info_LG_LG11_LG12_LG13.bed",
                                                      "data/orths/Orthogroups_TpsTpaTceBrs_onscafs_Tpa_c_1to1info_Scf10.bed",
                                                      "data/orths/Orthogroups_TpsTpaTceBrs_onscafs_Brs_c_1to1info.bed"
                                         ))



### tidy labels

my_orthologs_with_3_sp$sp1.Chr <- gsub("scf", "", my_orthologs_with_3_sp$sp1.Chr )
my_orthologs_with_3_sp$sp2.Chr <- gsub("LG", "", my_orthologs_with_3_sp$sp2.Chr )
my_orthologs_with_3_sp$sp3.Chr <- gsub("scf", "", my_orthologs_with_3_sp$sp3.Chr )
my_orthologs_with_3_sp$sp4.Chr <- gsub("scf", "", my_orthologs_with_3_sp$sp4.Chr )
my_orthologs_with_3_sp$sp4.Chr <- gsub("_", ".", my_orthologs_with_3_sp$sp4.Chr )
str(my_orthologs_with_3_sp)

my_orthologs_with_3_sp$sp1.Chr <- ifelse(my_orthologs_with_3_sp$sp1.Chr == "3", "X", my_orthologs_with_3_sp$sp1.Chr)
my_orthologs_with_3_sp$sp2.Chr <- ifelse(my_orthologs_with_3_sp$sp2.Chr == "3", "X", my_orthologs_with_3_sp$sp2.Chr)
my_orthologs_with_3_sp$sp3.Chr <- ifelse(my_orthologs_with_3_sp$sp3.Chr == "1", "X", my_orthologs_with_3_sp$sp3.Chr)
my_orthologs_with_3_sp$sp4.Chr <- ifelse(my_orthologs_with_3_sp$sp4.Chr == "2", "X", my_orthologs_with_3_sp$sp4.Chr)


my_orthologs_with_3_sp$sp1.Chr <- ordered(my_orthologs_with_3_sp$sp1.Chr, levels = c("X" , "1", "2",  "4",  "5",  "6",  "7",  "8",  "9", "10", "11", "12"   ))
my_orthologs_with_3_sp$sp2.Chr <- ordered(my_orthologs_with_3_sp$sp2.Chr, levels = c("X" , "6", "2",  "1", "4", "8",  "5",  "7",  "9", "11",  "10",  "12",  "13"  ))
my_orthologs_with_3_sp$sp3.Chr <- ordered(my_orthologs_with_3_sp$sp3.Chr, levels = c("X" , "2", "10", "5", "3", "4", "6",  "7", "8",  "9", "11" ,  "12",  "13" , "14" ))
my_orthologs_with_3_sp$sp4.Chr <- ordered(my_orthologs_with_3_sp$sp4.Chr, levels = c("X" , "16", "1",  "3" , "10" ,  "11", "7" ,  "15", "14",  "17" , "12",  "8" , "13" ,"18",   "4.1", "4.2",  "9.1" ,"9.2", "5",   "6" ))


plot_chord_diagram(my_orthologs_with_3_sp,
                   species_labels = c("Tps", "Tce", "Tpa", "Brs"),
                   color_by = "sp1.Chr", reorder_chromosomes = FALSE) + theme(legend.position = "none")




png(filename = "data/orths/TpsTceTpaBrs_chord_diagram_1to1orths_alllinks_tidynames.png", width = 9, height = 5, units = "in", bg = "white", res = 300)
plot_chord_diagram(my_orthologs_with_3_sp,
                   species_labels = c("T. poppense", "T. cristinae", "T. podura", "B. rossius"),
                   color_by = "sp1.Chr", reorder_chromosomes = FALSE, remove_non_linkage_orthologs = FALSE,  label_size = 4) + theme(legend.position = "none")

dev.off()

png(filename = "data/orths/TpsTceTpaBrs_chord_diagram_1to1orths_siglinksonly_tidynames.png", width = 12, height = 6, units = "in", bg = "white", res = 300)
plot_chord_diagram(my_orthologs_with_3_sp,
                   species_labels = c("T. poppense", "T. cristinae", "T. podura", "B. rossius"),
                   color_by = "sp1.Chr", reorder_chromosomes = FALSE,  label_size = 4) + theme(legend.position = "none")

dev.off()


pdf("data/orths/TpsTceTpaBrs_chord_diagram_1to1orths_siglinksonly_tidynames.pdf", width = 12, height = 6)
plot_chord_diagram(my_orthologs_with_3_sp,
                   species_labels = c("T. poppense", "T. cristinae", "T. podura", "B. rossius"),
                   color_by = "sp1.Chr", reorder_chromosomes = FALSE,  label_size = 4) + theme(legend.position = "none")

dev.off()


