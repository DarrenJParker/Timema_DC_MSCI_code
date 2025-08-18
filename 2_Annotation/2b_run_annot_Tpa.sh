## 2a_run_annot.sh


DIR='/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot' ### main dir
SDIR='/scratch/dparker/annot' ### scratch
genome_pref="Tce_LRv5a"

#############################################################################################
### get genomes

cd $DIR
mkdir genomes

### get from DDBJ/ENA/GenBank - accession JBOIUC000000000

#############################################################################################
### Repeat identification with Repeat Modeler through Dfam TE Tools # Container v1.4 (https://github.com/Dfam-consortium/TETools)

## Info
# https://www.repeatmasker.org/RepeatModeler/

## Brůna T, Hoff KJ, Lomsadze A, Stanke M, Borodovsky M (2021). BRAKER2: automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS supported by a protein database. NAR Genom Bioinform 3: lqaa108.
### https://academic.oup.com/nargab/article/3/1/lqaa108/6066535#supplementary-data
#Repeat masking by RepeatModeler and RepeatMasker with default settings was sufficient to achieve high gene prediction accuracy in all the tested genomes except for X. tropicalis.
#BuildDatabase -engine wublast -name genome genome.fasta
#RepeatModeler -engine wublast -database genome
#RepeatMasker -engine wublast -lib genome-families.fa -xsmall genome.fasta

## maker does one round of a specif lib, then a general one from model taxa.
## https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2 ### example with boa


module load singularity/3.7.4
singularity pull docker://dfam/tetools

## build database for rep modeler
singularity exec --bind $DIR tetools_latest.sif BuildDatabase -engine ncbi -name Tpa_LRv5a  genomes/Tpa_LRv5a.fasta 


## /opt/RepeatModeler/RepeatModeler - 2.0.2
## run rep modeler to make sp rep lib # 32 cpu, 60GB ram ## took Tpa: 1-17:12:27
singularity exec --bind $DIR tetools_latest.sif RepeatModeler -engine ncbi -database Tpa_LRv5a -pa 30 -LTRStruct


#RepeatMasker version 4.1.2-p1
#No query sequence file indicated
#/opt/RepeatMasker/RepeatMasker - 4.1.2-p1
## run rep masker # 32 cpu, 60GB ram ## took Tpa 04:49:34
singularity exec --bind $DIR tetools_latest.sif RepeatMasker -engine ncbi -gff -xsmall -pa 30 -lib Tpa_LRv5a-families.fa genomes/Tpa_LRv5a.fasta 

## cp stats
cp  $DIR/genomes/*tbl $DIR/Timema_LR_genomic_code/output/repmask/


##############################################################################################################################
###########################################################################################################################
### braker2

## prots
### arthropods
wget https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz
tar -zxf odb10_arthropoda_fasta.tar.gz
cat arthropoda/Rawdata/* > arth_proteins.fasta

## includes all the insect records
#1.	level NCBI tax id
#2.	scientific name
#3.	total non-redundant count of genes in all underneath clustered species
#4.	total count of OGs built on it
#5.	total non-redundant count of species underneath
#
#6656	Arthropoda	2206003	82474	172
#6960	Hexapoda	     1959590	73149	152
#7041	Coleoptera	96027	11817	9
#7088	Lepidoptera	202146	17449	16
#7147	Diptera	     718080	38562	56
#7148	Nematocera	281317	23093	24
#7157	Culicidae	     202572	18294	17
#7164	Anopheles	     159606	16075	14
#7203	Brachycera	427316	26449	32
#7214	Drosophilidae	271870	17688	20
#7215	Drosophila	271872	17694	20


#### + Timema from prev genomes from Jaron, K. S., Parker, D. J., Anselmetti, Y., Van, P. T., Bast, J., Dumas, Z., Figuet, E., François, C. M., Hayward, K., Rossier, V., Simion, P., Robinson-Rechavi, M., Galtier, N., & Schwander, T. (2022). Convergent consequences of parthenogenesis on stick insect genomes. Science Advances, 8(8), eabg3842.

gzip -d Accessory_scripts/AllTimema_b3v08_prot.fa.gz

## join with arth prot and tidy

cat arth_proteins.fasta AllTimema_b3v08_prot.fa > arth_and_timemav8_proteins.fasta 


######################################################################################################################
########## BRAKER 2 runs
## get SIF file from https://doi.org/10.5281/zenodo.11196714
## Braker is stocastic due to augustus. Gonna run each 5 times and take the best.
# arth_proteins.fasta 
# arth_and_timemav8_proteins.fasta 

#######################################################################################################################################
### prot runs | arth+Timema |  5 times each - take best

###### 3a | Run time 10:19:01
mkdir $SDIR/Tpa_braker_prot_run_3a
cd    $SDIR/Tpa_braker_prot_run_3a
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tpa sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_prot_run_3a --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_prot_run_3a/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa


###### 3b | Run time 
mkdir $SDIR/Tpa_braker_prot_run_3b
cd    $SDIR/Tpa_braker_prot_run_3b
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tpa sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_prot_run_3b --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_prot_run_3b/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa


###### 3c | Run time 
mkdir $SDIR/Tpa_braker_prot_run_3c
cd    $SDIR/Tpa_braker_prot_run_3c
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tpa sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_prot_run_3c --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_prot_run_3c/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa


###### 3d | Run time 
mkdir $SDIR/Tpa_braker_prot_run_3d
cd    $SDIR/Tpa_braker_prot_run_3d
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tpa sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_prot_run_3d --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_prot_run_3d/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa


###### 3e | Run time 
mkdir $SDIR/Tpa_braker_prot_run_3e
cd    $SDIR/Tpa_braker_prot_run_3e
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tpa sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_prot_run_3e --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_prot_run_3e/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa


## cp back from scratch

cp -r $SDIR/Tpa_braker_prot_run_3a $DIR
cp -r $SDIR/Tpa_braker_prot_run_3b $DIR
cp -r $SDIR/Tpa_braker_prot_run_3c $DIR
cp -r $SDIR/Tpa_braker_prot_run_3d $DIR
cp -r $SDIR/Tpa_braker_prot_run_3e $DIR


awk '{print $3}' $DIR/Tpa_braker_prot_run_3a/braker.gff3 | grep "gene" | wc -l
# 47289
awk '{print $3}' $DIR/Tpa_braker_prot_run_3b/braker.gff3 | grep "gene" | wc -l
# 47259
awk '{print $3}' $DIR/Tpa_braker_prot_run_3c/braker.gff3 | grep "gene" | wc -l
# 47270
awk '{print $3}' $DIR/Tpa_braker_prot_run_3d/braker.gff3 | grep "gene" | wc -l
# 46859
awk '{print $3}' $DIR/Tpa_braker_prot_run_3e/braker.gff3 | grep "gene" | wc -l
# 45557

# tarball
cd $DIR

tar -czf Tpa_braker_prot_run_3a.tar.gz Tpa_braker_prot_run_3a
tar -czf Tpa_braker_prot_run_3b.tar.gz Tpa_braker_prot_run_3b
tar -czf Tpa_braker_prot_run_3c.tar.gz Tpa_braker_prot_run_3c
tar -czf Tpa_braker_prot_run_3d.tar.gz Tpa_braker_prot_run_3d
tar -czf Tpa_braker_prot_run_3e.tar.gz Tpa_braker_prot_run_3e


rm -r Tpa_braker_prot_run_3a/
rm -r Tpa_braker_prot_run_3b/
rm -r Tpa_braker_prot_run_3c/
rm -r Tpa_braker_prot_run_3d/
rm -r Tpa_braker_prot_run_3e/


#######################################################################################################
### map - STAR with --twopassMode Basic (to get the splice junction right)
### map invid then merge.
### map all from Tpa and Tge to Tpa

### 39 conditions (NOT using single-end GN, LG, HD)
### 36 paired, 3 single end
### should get 136 paired bams, 12 single bams
### See TpaTge_Annot_RNAseq_samples.xlsx

##########################
### Index
### Build index for star
### not using soft masked genome here - STAR ignores soft masking anyway (tested)

module load gcc
module load star/2.7.8a

mkdir $SDIR/RNAseq_mapping

STAR --runThreadN 12 \
     --runMode genomeGenerate \
     --genomeDir $SDIR/$genome_pref"_STAR" \
     --genomeFastaFiles $DIR"/genomes/"$genome_pref".fasta" \
     --genomeChrBinNbits 20 --limitGenomeGenerateRAM 79000000000


## map paired end reads


Tpa_Tge_paired_samples=(
Tpa_F_A_Ad
Tpa_F_B_Ad
Tpa_F_DG_Ad
Tpa_F_FB_Ad
Tpa_F_Fe_Ad
Tpa_F_Go_Ad
Tpa_F_Gu_Ad
Tpa_F_Ta_Ad
Tpa_M_A_Ad
Tpa_M_AG_Ad
Tpa_M_B_Ad
Tpa_M_DG_Ad
Tpa_M_FB_Ad
Tpa_M_Fe_Ad
Tpa_M_Gu_Ad
Tpa_M_Ta_Ad
Tpa_M_Te_Ad
Tpa_U_WB_Ha
Tge_F_A_Ad
Tge_F_B_Ad
Tge_F_DG_Ad
Tge_F_FB_Ad
Tge_F_Fe_Ad
Tge_F_Go_Ad
Tge_F_Gu_Ad
Tge_F_Ta_Ad
Tge_F_WB_Ha
Tge_M_A_Ad
Tge_M_AG_Ad
Tge_M_B_Ad
Tge_M_DG_Ad
Tge_M_FB_Ad
Tge_M_Fe_Ad
Tge_M_Gu_Ad
Tge_M_Ta_Ad
Tge_M_Te_Ad
)

for R1_f in /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq/trimmed_reads_paired/*_R1_*qtrimmed.fq.gz; do
R2_f=`echo $R1_f | sed 's/_R1_/_R2_/' `
out_prefix=`echo $R1_f | sed 's/.*\///' | sed 's/_R1_.*//'`
sp=`echo $R1_f | sed 's/.*\///' | sed 's/_.*//'`

for s in ${Tpa_Tge_paired_samples[@]}; do

if  [[ $out_prefix == $s* ]];
then
echo $R1_f
echo $R2_f
echo $sp
echo $genome_pref
echo $out_prefix
STAR --twopassMode Basic --genomeDir $SDIR/$genome_pref"_STAR" \
     --readFilesIn $R1_f $R2_f --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --runThreadN 12 \
     --outFileNamePrefix  $SDIR"/RNAseq_mapping/"$out_prefix"_to_"$genome_pref
fi
done
done 




### single end

mkdir $SDIR/RNAseq_mapping

## WB single end stuff
Tpa_Tge_single_samples=(
Tpa_F_WB_Ad
Tge_F_WB_Ad
Tpa_M_WB_Ad	
)

for R1_f in /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq/trimmed_reads_single/*.fq.gz; do
out_prefix=`echo $R1_f | sed 's/.*\///' | sed 's/_AandQtrimmed.*//' | sed 's/_trimmed.*//'  | sed 's/sampled.*//' | sed 's/_DJP.*//' | sed 's/_Swb.*//'  | sed 's/_Eth.*//' | sed 's/_S2e.*//'  `
sp=`echo $R1_f | sed 's/.*\///' | sed 's/_.*//'`

for s in ${Tpa_Tge_single_samples[@]}; do

if  [[ $out_prefix == $s* ]];
then
echo $R1_f
echo $sp
echo $genome_pref
echo $out_prefix
STAR --twopassMode Basic --genomeDir $SDIR/$genome_pref"_STAR" \
     --readFilesIn $R1_f --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --runThreadN 12 \
     --outFileNamePrefix  $SDIR"/RNAseq_mapping/"$out_prefix"_to_"$genome_pref
fi
done
done 




### STORE BAMs HERE /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/RNAseq_mapping

ls $SDIR"/RNAseq_mapping/"*"_to_"$genome_pref*"Aligned.sortedByCoord.out.bam" | wc -l 
# 148
mv $SDIR"/RNAseq_mapping/"*"_to_"$genome_pref*"Aligned.sortedByCoord.out.bam" /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/RNAseq_mapping

### STORE map stats HERE /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/RNAseq_mapping
mkdir /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/to_$genome_pref
cp $SDIR/RNAseq_mapping/*$genome_pref"Log.final.out" /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/to_$genome_pref
cd /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/
module load singularity/3.7.4
singularity exec -e --bind /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq,/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/ \
/work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq/multiqc_latest.sif multiqc \
to_$genome_pref/ -o to_$genome_pref"_STAR" --interactive
cp to_$genome_pref"_STAR/multiqc_report.html"  to_$genome_pref"_STAR_multiqc_report.html" 
rm -r to_$genome_pref"_STAR"
tar -czvf to_$genome_pref".tar.gz" to_$genome_pref
rm -r to_$genome_pref



###########################################################################################################################
#### Merge bams
## by sp / cond

module load gcc
module load samtools/1.12 
mkdir $SDIR"/RNAseq_mapping_merged"

Tpa_Tge_all_samples=(
Tpa_F_A_Ad
Tpa_F_B_Ad
Tpa_F_DG_Ad
Tpa_F_FB_Ad
Tpa_F_Fe_Ad
Tpa_F_Go_Ad
Tpa_F_Gu_Ad
Tpa_F_Ta_Ad
Tpa_M_A_Ad
Tpa_M_AG_Ad
Tpa_M_B_Ad
Tpa_M_DG_Ad
Tpa_M_FB_Ad
Tpa_M_Fe_Ad
Tpa_M_Gu_Ad
Tpa_M_Ta_Ad
Tpa_M_Te_Ad
Tpa_U_WB_Ha
Tge_F_A_Ad
Tge_F_B_Ad
Tge_F_DG_Ad
Tge_F_FB_Ad
Tge_F_Fe_Ad
Tge_F_Go_Ad
Tge_F_Gu_Ad
Tge_F_Ta_Ad
Tge_F_WB_Ha
Tge_M_A_Ad
Tge_M_AG_Ad
Tge_M_B_Ad
Tge_M_DG_Ad
Tge_M_FB_Ad
Tge_M_Fe_Ad
Tge_M_Gu_Ad
Tge_M_Ta_Ad
Tge_M_Te_Ad
Tpa_F_WB_Ad
Tge_F_WB_Ad
Tpa_M_WB_Ad	
)


for s in ${Tpa_Tge_all_samples[@]}; do
echo $s
samtools merge -@ 40 $SDIR"/RNAseq_mapping_merged/"$s"_to_"$genome_pref"_merged.bam" "/work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/RNAseq_mapping/"$s*"_to_"$genome_pref"Aligned.sortedByCoord.out.bam"
done

for b in $SDIR"/RNAseq_mapping_merged/"*"_to_"$genome_pref"_merged.bam"; do
echo $b
out_b=`echo $b | sed 's/.bam/_sorted.bam/'`
echo $out_b
samtools sort  -@ 40 $b -o $out_b
done

### rm merged bams
rm $SDIR"/RNAseq_mapping_merged/"*"_to_"$genome_pref"_merged.bam"
ls $SDIR"/RNAseq_mapping_merged/"*"_to_"$genome_pref"_merged_sorted.bam" -l | wc -l
# 39

#### merge for the UTR
module load gcc
module load samtools/1.12
ulimit -n 100000

samtools merge -@ 40 $SDIR"/ALL_to_"$genome_pref"_merged.bam" "/work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/RNAseq_mapping/"$s*"_to_"$genome_pref"Aligned.sortedByCoord.out.bam"
samtools sort  -@ 40 $SDIR"/ALL_to_"$genome_pref"_merged.bam" -o $SDIR"/ALL_to_"$genome_pref"_merged_sorted.bam"

rm $SDIR"/"*"_to_"$genome_pref"_merged.bam"


#######################################################################################################################################
### RNAseq runs |  5 times - take best to add UTR onto 

###### 4a 
mkdir $SDIR/Tpa_braker_RNAseq_run_4a
cd    $SDIR/Tpa_braker_RNAseq_run_4a
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_RNAseq_run_4a --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_RNAseq_run_4a/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa

###### 4b
mkdir $SDIR/Tpa_braker_RNAseq_run_4b
cd    $SDIR/Tpa_braker_RNAseq_run_4b
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_RNAseq_run_4b --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_RNAseq_run_4b/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa


###### 4c
mkdir $SDIR/Tpa_braker_RNAseq_run_4c
cd    $SDIR/Tpa_braker_RNAseq_run_4c
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_RNAseq_run_4c --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_RNAseq_run_4c/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa


###### 4d 
mkdir $SDIR/Tpa_braker_RNAseq_run_4d
cd    $SDIR/Tpa_braker_RNAseq_run_4d
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_RNAseq_run_4d --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_RNAseq_run_4d/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa



###### 4e 
mkdir $SDIR/Tpa_braker_RNAseq_run_4e
cd    $SDIR/Tpa_braker_RNAseq_run_4e
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_RNAseq_run_4e --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_RNAseq_run_4e/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa


cp -r $SDIR/Tpa_braker_RNAseq_run_4a $DIR
cp -r $SDIR/Tpa_braker_RNAseq_run_4b $DIR
cp -r $SDIR/Tpa_braker_RNAseq_run_4c $DIR
cp -r $SDIR/Tpa_braker_RNAseq_run_4d $DIR
cp -r $SDIR/Tpa_braker_RNAseq_run_4e $DIR

#8237553                braker_rnaseq_4a_Tpa              07:57:24  COMPLETED 
#8237781                braker_rnaseq_4b_Tpa              08:32:23  COMPLETED 
#8237782                braker_rnaseq_4c_Tpa              06:31:02  COMPLETED 
#8237783                braker_rnaseq_4d_Tpa              09:47:11  COMPLETED 
#8237784                braker_rnaseq_4e_Tpa              07:52:35  COMPLETED 


awk '{print $3}' $DIR/Tpa_braker_RNAseq_run_4a/braker.gff3 | grep "gene" | wc -l
# 41125
awk '{print $3}' $DIR/Tpa_braker_RNAseq_run_4b/braker.gff3 | grep "gene" | wc -l
# 39984
awk '{print $3}' $DIR/Tpa_braker_RNAseq_run_4c/braker.gff3 | grep "gene" | wc -l
# 43474
awk '{print $3}' $DIR/Tpa_braker_RNAseq_run_4d/braker.gff3 | grep "gene" | wc -l
# 40195
awk '{print $3}' $DIR/Tpa_braker_RNAseq_run_4e/braker.gff3 | grep "gene" | wc -l
# 38070



#########################################################################################################################################################################################################
### UTR

########################################## bigmem 
## add the UTR ## this will delete the bam.
### use modified gushr bigmem ## 230GB mem 48 cores

cd $DIR
git clone https://github.com/Gaius-Augustus/GUSHR.git
mv GUSHR GUSHR_DJP_bigmem
cp $DIR/accessory_scripts/GUSHR_DJP_bigmem/gushr.py GUSHR_DJP_bigmem/

## make a config file to point braker to it
gushr_config_DJP_bigmem.txt 
GUSHR_PATH=/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/GUSHR_DJP_bigmem

### cp braker dicts
cp -r $DIR/Tpa_braker_RNAseq_run_4a $SDIR/Tpa_braker_RNAseq_run_4a_utrBMAll
cp -r $DIR/Tpa_braker_RNAseq_run_4b $SDIR/Tpa_braker_RNAseq_run_4b_utrBMAll
cp -r $DIR/Tpa_braker_RNAseq_run_4c $SDIR/Tpa_braker_RNAseq_run_4c_utrBMAll
cp -r $DIR/Tpa_braker_RNAseq_run_4d $SDIR/Tpa_braker_RNAseq_run_4d_utrBMAll
cp -r $DIR/Tpa_braker_RNAseq_run_4e $SDIR/Tpa_braker_RNAseq_run_4e_utrBMAll


### cp merged bam
nohup cp $SDIR/ALL_to_Tpa_LRv5a_merged_sorted.bam $SDIR/Tpa_braker_RNAseq_run_4a_utrBMAll &
nohup cp $SDIR/ALL_to_Tpa_LRv5a_merged_sorted.bam $SDIR/Tpa_braker_RNAseq_run_4b_utrBMAll &
nohup cp $SDIR/ALL_to_Tpa_LRv5a_merged_sorted.bam $SDIR/Tpa_braker_RNAseq_run_4c_utrBMAll &
nohup cp $SDIR/ALL_to_Tpa_LRv5a_merged_sorted.bam $SDIR/Tpa_braker_RNAseq_run_4d_utrBMAll &
nohup cp $SDIR/ALL_to_Tpa_LRv5a_merged_sorted.bam $SDIR/Tpa_braker_RNAseq_run_4e_utrBMAll &


module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_RNAseq_run_4a_utrBMAll \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tpa_braker_RNAseq_run_4a_utrBMAll/ALL_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_RNAseq_run_4a_utrBMAll/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tpa_braker_RNAseq_run_4a_utrBMAll/augustus.hints.gtf --skipAllTraining
module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_RNAseq_run_4b_utrBMAll \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tpa_braker_RNAseq_run_4b_utrBMAll/ALL_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_RNAseq_run_4b_utrBMAll/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tpa_braker_RNAseq_run_4b_utrBMAll/augustus.hints.gtf --skipAllTraining

module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_RNAseq_run_4c_utrBMAll \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tpa_braker_RNAseq_run_4c_utrBMAll/ALL_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_RNAseq_run_4c_utrBMAll/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tpa_braker_RNAseq_run_4c_utrBMAll/augustus.hints.gtf --skipAllTraining

module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_RNAseq_run_4d_utrBMAll \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tpa_braker_RNAseq_run_4d_utrBMAll/ALL_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_RNAseq_run_4d_utrBMAll/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tpa_braker_RNAseq_run_4d_utrBMAll/augustus.hints.gtf --skipAllTraining

module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tpa_braker_RNAseq_run_4e_utrBMAll \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tpa_braker_RNAseq_run_4e_utrBMAll/ALL_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tpa_braker_RNAseq_run_4e_utrBMAll/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tpa --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tpa_braker_RNAseq_run_4e_utrBMAll/augustus.hints.gtf --skipAllTraining



cp -r $SDIR/Tpa_braker_RNAseq_run_4a_utrBMAll $DIR
cp -r $SDIR/Tpa_braker_RNAseq_run_4b_utrBMAll $DIR
cp -r $SDIR/Tpa_braker_RNAseq_run_4c_utrBMAll $DIR
cp -r $SDIR/Tpa_braker_RNAseq_run_4d_utrBMAll $DIR
cp -r $SDIR/Tpa_braker_RNAseq_run_4e_utrBMAll $DIR



############################################################
### merge prot + RNASeq + UTR
# Tpa_braker_prot_run_3a, Tpa_braker_RNAseq_run_4c


mkdir /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/annotations_v2/
mkdir /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/annotations_v2/orig_merge

tar -zxf Tpa_braker_prot_run_3a.tar.gz
./TSEBRA/bin/fix_gtf_ids.py --gtf Tpa_braker_prot_run_3a/braker.gtf   --out Tpa_braker_prot_run_3a/braker_fixed.gtf
./TSEBRA/bin/fix_gtf_ids.py --gtf Tpa_braker_RNAseq_run_4c/braker.gtf --out Tpa_braker_RNAseq_run_4c/braker_fixed.gtf
./TSEBRA/bin/fix_gtf_ids.py --gtf Tpa_braker_RNAseq_run_4c_utrBMAll/braker_utr.gtf --out Tpa_braker_RNAseq_run_4c_utrBMAll/braker_utr_fixed.gtf

#### no-UTR version pref_braker1.cfg (RNAseq evi more weight)
./TSEBRA/bin/tsebra.py -g Tpa_braker_prot_run_3a/braker_fixed.gtf,Tpa_braker_RNAseq_run_4c/braker_fixed.gtf \
-c TSEBRA/config/pref_braker1.cfg \
-e Tpa_braker_prot_run_3a/hintsfile.gff,Tpa_braker_RNAseq_run_4c/hintsfile.gff \
-o Tpa_braker_prot_run_3a_Tpa_braker_RNAseq_run_4c_combined_2.gtf

#### UTR version pref_braker1.cfg (RNAseq evi more weight)
./TSEBRA/bin/tsebra.py -g Tpa_braker_prot_run_3a/braker_fixed.gtf,Tpa_braker_RNAseq_run_4c/braker_fixed.gtf,Tpa_braker_RNAseq_run_4c_utrBMAll/braker_utr_fixed.gtf \
-c TSEBRA/config/pref_braker1.cfg \
-e Tpa_braker_prot_run_3a/hintsfile.gff,Tpa_braker_RNAseq_run_4c/hintsfile.gff,Tpa_braker_RNAseq_run_4c_utrBMAll/hintsfile.gff \
-o Tpa_braker_prot_run_3a_Tpa_braker_RNAseq_run_4cUTR_combined_2.gtf

### add to v2 
mv Tpa_braker_prot_run_3a_Tpa_braker_RNAseq_run_4c_combined_2.gtf    $DIR/annotations_v2/orig_merge
mv Tpa_braker_prot_run_3a_Tpa_braker_RNAseq_run_4cUTR_combined_2.gtf $DIR/annotations_v2/orig_merge

python3 ~/Gen_BioInf/fix_braker_gtf.py -i $DIR/annotations_v2/orig_merge/Tpa_braker_prot_run_3a_Tpa_braker_RNAseq_run_4c_combined_2.gtf -d 6 -p Tpa -G 
python3 ~/Gen_BioInf/fix_braker_gtf.py -i $DIR/annotations_v2/orig_merge/Tpa_braker_prot_run_3a_Tpa_braker_RNAseq_run_4cUTR_combined_2.gtf -d 6 -p Tpa -G 
python3 ~/Gen_BioInf/UTR_GUSHR_fix.py -g $DIR/annotations_v2/orig_merge/Tpa_braker_prot_run_3a_Tpa_braker_RNAseq_run_4c_combined_2_FBGgi.gff -u $DIR/annotations_v2/orig_merge/Tpa_braker_prot_run_3a_Tpa_braker_RNAseq_run_4cUTR_combined_2_FBGgi.gff -m 1000 -G
cp $DIR/annotations_v2/orig_merge/Tpa_braker_prot_run_3a_Tpa_braker_RNAseq_run_4c_combined_2_FBGgiUGF1000M0.gff $DIR/annotations_v2
mv $DIR/annotations_v2/Tpa_braker_prot_run_3a_Tpa_braker_RNAseq_run_4c_combined_2_FBGgiUGF1000M0.gff $DIR/annotations_v2/Tpa_braker_prot_run_3a_RNAseq_run_4c_combined_2_FBGgiUGF1000M0.gff 

### run BUSCO on the annotation

module load gcc
module load bedtools2/2.29.2

for i in $DIR/annotations_v2/Tpa*_combined_2_FBGgiUGF1000M0.gff; do
echo $i
out_pref=`echo $i | sed 's/.gff.*//'`
echo $out_pref
awk '$3 == "gene"' $i > $out_pref"_genes.gff"
cut -f 1,4,5,9 $out_pref"_genes.gff"        | sort -k1,3 -u > $out_pref"_genes_coord.temp" ## as some rev orient genes
cut -f 4       $out_pref"_genes_coord.temp" | sed 's/.*ID=//' | sed 's/;.*//' > $out_pref"_genes_names.temp"
cut -f 1,2,3   $out_pref"_genes_coord.temp" > $out_pref"_genes_coord_2.temp" 
paste $out_pref"_genes_coord_2.temp" $out_pref"_genes_names.temp" > $out_pref"_genes.bed"
bedtools getfasta -fi genomes/$genome_pref.fasta -bed  $out_pref"_genes.bed" -name >  $out_pref"_genes.fa"
python3 ~/DJP_py_scripts/fasta_len_0.2.py $out_pref"_genes.fa"
rm $out_pref"_genes.gff"
rm $out_pref"_genes_coord.temp"
rm $out_pref"_genes_names.temp"
rm $out_pref"_genes_coord_2.temp"
rm $out_pref"_genes.bed"
done


### run BUSCO
module load singularity/3.7.4
for f in $DIR/annotations_v2/Tpa*_genes.fa; do
out_f=`echo $f | sed 's/.fa/_BUSCO/' `
echo $f
echo $out_f
singularity exec --bind $DIR /work/FAC/FBM/DEE/tschwand/default/dparker/busco_v5.3.2_cv1.sif busco -c 30 -m genome -i $f -o $out_f -l insecta_odb10 --offline --download_path $DIR/busco_downloads/ 
done

mkdir annotations_v2/BUSCO
mv    work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/annotations_v2/* annotations_v2/BUSCO/

##### add missing busco

python3 ~/Gen_BioInf/Add_missing_busco.py -m annotations_v2/BUSCO/Tpa_braker_*_combined_2_FBGgiUGF1000M0_genes_BUSCO/run_insecta_odb10/missing_busco_list.tsv -b BUSCO/$genome_pref/ -g annotations_v2/Tpa_*_combined_2_FBGgiUGF1000M0.gff 
#N missing BUSCO genes:
#29
#Number of complete busco found in genome busco run:
#19


