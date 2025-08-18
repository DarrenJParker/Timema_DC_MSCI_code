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

## calc in the /scratch

cd $SDIR

module load singularity/3.7.4
#singularity pull docker://dfam/tetools

## build database for rep modeler
singularity exec --bind $DIR,$SDIR $DIR/tetools_latest.sif BuildDatabase -engine ncbi -name $genome_pref $DIR/genomes/$genome_pref.fasta


## /opt/RepeatModeler/RepeatModeler - 2.0.2
## run rep modeler to make sp rep lib # 32 cpu, 80GB ram ## took: 2-02:22:38
singularity exec --bind $DIR,$SDIR $DIR/tetools_latest.sif RepeatModeler -engine ncbi -database $genome_pref -pa 30 -LTRStruct


#RepeatMasker version 4.1.2-p1
#No query sequence file indicated
#/opt/RepeatMasker/RepeatMasker - 4.1.2-p1
## run rep masker # 32 cpu, 80GB ram ## took 07:01:29

genome_file=`echo $DIR/genomes/$genome_pref".fasta"`
singularity exec --bind $DIR,$SDIR $DIR/tetools_latest.sif RepeatMasker -engine ncbi -gff -xsmall -pa 30 -lib $SDIR/$genome_pref"-families.fa" $genome_file

## cp stats
cp  $DIR/genomes/*tbl $DIR/Timema_LR_genomic_code/output/repmask/



######################################################################################################################
########## BRAKER 2 run
## get SIF file from https://doi.org/10.5281/zenodo.11196714
## Braker is stocastic due to augustus. Gonna rune each 5 times and take the best.
## prots from 2b_run_annot_Tpa.sh
# arth_proteins.fasta 
# arth_and_timemav8_proteins.fasta 

#######################################################################################################################################
### prot runs | arth+Timema |  5 times each - take best

###### 3a | Run time 10:09:17
mkdir $SDIR/Tce_braker_prot_run_3a
cd    $SDIR/Tce_braker_prot_run_3a
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tce sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_prot_run_3a --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_prot_run_3a/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce


###### 3b | Run time 
mkdir $SDIR/Tce_braker_prot_run_3b
cd    $SDIR/Tce_braker_prot_run_3b
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tce sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_prot_run_3b --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_prot_run_3b/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce


###### 3c | Run time 
mkdir $SDIR/Tce_braker_prot_run_3c
cd    $SDIR/Tce_braker_prot_run_3c
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tce sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_prot_run_3c --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_prot_run_3c/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce


###### 3d | Run time 
mkdir $SDIR/Tce_braker_prot_run_3d
cd    $SDIR/Tce_braker_prot_run_3d
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tce sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_prot_run_3d --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_prot_run_3d/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce


###### 3e | Run time 
mkdir $SDIR/Tce_braker_prot_run_3e
cd    $SDIR/Tce_braker_prot_run_3e
git clone https://github.com/Gaius-Augustus/Augustus.git ### as using Tce sep from RNASeq
cd $DIR

module load singularity/3.7.4
singularity exec --bind $DIR,$SDIR braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_prot_run_3e --genome $DIR/genomes/$genome_pref.fasta.masked --gff3 \
--prot_seq=$DIR/arth_and_timemav8_proteins.fasta --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_prot_run_3e/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce


## cp back from scratch

cp -r $SDIR/Tce_braker_prot_run_3a $DIR
cp -r $SDIR/Tce_braker_prot_run_3b $DIR
cp -r $SDIR/Tce_braker_prot_run_3c $DIR
cp -r $SDIR/Tce_braker_prot_run_3d $DIR
cp -r $SDIR/Tce_braker_prot_run_3e $DIR


awk '{print $3}' $DIR/Tce_braker_prot_run_3a/braker.gff3 | grep "gene" | wc -l
# 39911
awk '{print $3}' $DIR/Tce_braker_prot_run_3b/braker.gff3 | grep "gene" | wc -l
# 38760
awk '{print $3}' $DIR/Tce_braker_prot_run_3c/braker.gff3 | grep "gene" | wc -l
# 39718
awk '{print $3}' $DIR/Tce_braker_prot_run_3d/braker.gff3 | grep "gene" | wc -l
# 39693
awk '{print $3}' $DIR/Tce_braker_prot_run_3e/braker.gff3 | grep "gene" | wc -l
# 40345


# tarball 
cd $DIR


tar -czf Tce_braker_prot_run_3a.tar.gz Tce_braker_prot_run_3a
tar -czf Tce_braker_prot_run_3b.tar.gz Tce_braker_prot_run_3b
tar -czf Tce_braker_prot_run_3c.tar.gz Tce_braker_prot_run_3c
tar -czf Tce_braker_prot_run_3d.tar.gz Tce_braker_prot_run_3d
tar -czf Tce_braker_prot_run_3e.tar.gz Tce_braker_prot_run_3e

rm -r Tce_braker_prot_run_3a/
rm -r Tce_braker_prot_run_3b/
rm -r Tce_braker_prot_run_3c/
rm -r Tce_braker_prot_run_3d/
rm -r Tce_braker_prot_run_3e/


#######################################################################################################
### map - STAR with --twopassMode Basic (to get the splice junction right)
### map invid then merge.
### map all from Tce and Tms to Tce

### 41 conditions (NOT using single-end GN, LG, HD)
### 38 paired, 3 single end
### should get 194 paired bams, 12 single bams
### See TceTms_Annot_RNAseq_samples.xlsx


Tce_Tms_paired_samples=(
Tce_F_A_Ad
Tce_F_A_Ju
Tce_M_A_Ad
Tce_M_A_Ju
Tce_M_AG_Ad
Tce_F_B_Ad
Tce_M_B_Ad
Tce_F_DG_Ad
Tce_M_DG_Ad
Tce_F_FB_Ad
Tce_M_FB_Ad
Tce_F_Fe_Ad
Tce_M_Fe_Ad
Tce_F_Go_Ad
Tce_F_Gu_Ad
Tce_M_Gu_Ad
Tce_F_Ta_Ad
Tce_M_Ta_Ad
Tce_M_Te_Ad
Tce_U_WB_Ha
Tms_F_A_Ad
Tms_M_A_Ad
Tms_M_AG_Ad
Tms_F_B_Ad
Tms_M_B_Ad
Tms_F_DG_Ad
Tms_M_DG_Ad
Tms_F_FB_Ad
Tms_M_FB_Ad
Tms_F_Fe_Ad
Tms_M_Fe_Ad
Tms_F_Go_Ad
Tms_F_Gu_Ad
Tms_M_Gu_Ad
Tms_F_Ta_Ad
Tms_M_Ta_Ad
Tms_M_Te_Ad
Tms_F_WB_Ha
)

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

for R1_f in /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq/trimmed_reads_paired/*_R1_*qtrimmed.fq.gz; do
R2_f=`echo $R1_f | sed 's/_R1_/_R2_/' `
out_prefix=`echo $R1_f | sed 's/.*\///' | sed 's/_R1_.*//'`
sp=`echo $R1_f | sed 's/.*\///' | sed 's/_.*//'`

for s in ${Tce_Tms_paired_samples[@]}; do

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

Tce_Tms_single_samples=(
Tms_F_WB_Ad
Tce_F_WB_Ad
Tce_M_WB_Ad
)

for R1_f in /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq/trimmed_reads_single/*.fq.gz; do
out_prefix=`echo $R1_f | sed 's/.*\///' | sed 's/_AandQtrimmed.*//' | sed 's/_trimmed.*//'  | sed 's/sampled.*//' | sed 's/_DJP.*//' | sed 's/_Swb.*//'  | sed 's/_Eth.*//' | sed 's/_S2e.*//'  `
sp=`echo $R1_f | sed 's/.*\///' | sed 's/_.*//'`

for s in ${Tce_Tms_single_samples[@]}; do

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
# 206
mv $SDIR"/RNAseq_mapping/"*"_to_"$genome_pref*"Aligned.sortedByCoord.out.bam" /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/RNAseq_mapping


### STORE map stats HERE /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/
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




############################################################################################################################
#### Merge bams
## by sp / cond

module load gcc
module load samtools/1.12 
mkdir $SDIR"/RNAseq_mapping_merged"

Tce_Tms_all_samples=(
Tce_F_A_Ad
Tce_F_A_Ju
Tce_M_A_Ad
Tce_M_A_Ju
Tce_M_AG_Ad
Tce_F_B_Ad
Tce_M_B_Ad
Tce_F_DG_Ad
Tce_M_DG_Ad
Tce_F_FB_Ad
Tce_M_FB_Ad
Tce_F_Fe_Ad
Tce_M_Fe_Ad
Tce_F_Go_Ad
Tce_F_Gu_Ad
Tce_M_Gu_Ad
Tce_F_Ta_Ad
Tce_M_Ta_Ad
Tce_M_Te_Ad
Tce_U_WB_Ha
Tms_F_A_Ad
Tms_M_A_Ad
Tms_M_AG_Ad
Tms_F_B_Ad
Tms_M_B_Ad
Tms_F_DG_Ad
Tms_M_DG_Ad
Tms_F_FB_Ad
Tms_M_FB_Ad
Tms_F_Fe_Ad
Tms_M_Fe_Ad
Tms_F_Go_Ad
Tms_F_Gu_Ad
Tms_M_Gu_Ad
Tms_F_Ta_Ad
Tms_M_Ta_Ad
Tms_M_Te_Ad
Tms_F_WB_Ha
Tms_F_WB_Ad
Tce_F_WB_Ad
Tce_M_WB_Ad
)


for s in ${Tce_Tms_all_samples[@]}; do
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
# 41

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
mkdir $SDIR/Tce_braker_RNAseq_run_4a
cd    $SDIR/Tce_braker_RNAseq_run_4a
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_RNAseq_run_4a --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_RNAseq_run_4a/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce

###### 4b
mkdir $SDIR/Tce_braker_RNAseq_run_4b
cd    $SDIR/Tce_braker_RNAseq_run_4b
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_RNAseq_run_4b --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_RNAseq_run_4b/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce

###### 4c
mkdir $SDIR/Tce_braker_RNAseq_run_4c
cd    $SDIR/Tce_braker_RNAseq_run_4c
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_RNAseq_run_4c --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_RNAseq_run_4c/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce

###### 4d
mkdir $SDIR/Tce_braker_RNAseq_run_4d
cd    $SDIR/Tce_braker_RNAseq_run_4d
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_RNAseq_run_4d --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_RNAseq_run_4d/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce

###### 4e
mkdir $SDIR/Tce_braker_RNAseq_run_4e
cd    $SDIR/Tce_braker_RNAseq_run_4e
git clone https://github.com/Gaius-Augustus/Augustus.git ### as doing sep RNAseq and prot runs

module load singularity/3.7.4
bam_list=`ls -1p $SDIR/RNAseq_mapping_merged/*_to_$genome_pref"_merged_sorted.bam" | xargs echo | sed 's/ /,/g'`
singularity exec --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_RNAseq_run_4e --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 \
--bam $bam_list  --softmasking --cores 30 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_RNAseq_run_4e/Augustus/config/ \
--GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce


cp -r $SDIR/Tce_braker_RNAseq_run_4a $DIR
cp -r $SDIR/Tce_braker_RNAseq_run_4b $DIR
cp -r $SDIR/Tce_braker_RNAseq_run_4c $DIR
cp -r $SDIR/Tce_braker_RNAseq_run_4d $DIR
cp -r $SDIR/Tce_braker_RNAseq_run_4e $DIR

#8081276                braker_rnaseq_4a_Tce              10:06:36  COMPLETED 
#8083302                braker_rnaseq_4b_Tce              08:45:28  COMPLETED 
#8083303                braker_rnaseq_4c_Tce              09:41:13  COMPLETED 
#8083308                braker_rnaseq_4d_Tce              09:43:52  COMPLETED 
#8083316                braker_rnaseq_4e_Tce              07:32:05  COMPLETED 

awk '{print $3}' $DIR/Tce_braker_RNAseq_run_4a/braker.gff3 | grep "gene" | wc -l
# 49055
awk '{print $3}' $DIR/Tce_braker_RNAseq_run_4b/braker.gff3 | grep "gene" | wc -l
# 47804
awk '{print $3}' $DIR/Tce_braker_RNAseq_run_4c/braker.gff3 | grep "gene" | wc -l
# 48077
awk '{print $3}' $DIR/Tce_braker_RNAseq_run_4d/braker.gff3 | grep "gene" | wc -l
# 49053
awk '{print $3}' $DIR/Tce_braker_RNAseq_run_4e/braker.gff3 | grep "gene" | wc -l
# 49226


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
cp -r $DIR/Tce_braker_RNAseq_run_4a $SDIR/Tce_braker_RNAseq_run_4a_utrBMAll
cp -r $DIR/Tce_braker_RNAseq_run_4b $SDIR/Tce_braker_RNAseq_run_4b_utrBMAll
cp -r $DIR/Tce_braker_RNAseq_run_4c $SDIR/Tce_braker_RNAseq_run_4c_utrBMAll
cp -r $DIR/Tce_braker_RNAseq_run_4d $SDIR/Tce_braker_RNAseq_run_4d_utrBMAll
cp -r $DIR/Tce_braker_RNAseq_run_4e $SDIR/Tce_braker_RNAseq_run_4e_utrBMAll


### cp merged bam
nohup cp $SDIR/ALL_to_Tce_LRv5a_merged_sorted.bam $SDIR/Tce_braker_RNAseq_run_4a_utrBMAll &
nohup cp $SDIR/ALL_to_Tce_LRv5a_merged_sorted.bam $SDIR/Tce_braker_RNAseq_run_4b_utrBMAll &
nohup cp $SDIR/ALL_to_Tce_LRv5a_merged_sorted.bam $SDIR/Tce_braker_RNAseq_run_4c_utrBMAll &
nohup cp $SDIR/ALL_to_Tce_LRv5a_merged_sorted.bam $SDIR/Tce_braker_RNAseq_run_4d_utrBMAll &
nohup cp $SDIR/ALL_to_Tce_LRv5a_merged_sorted.bam $SDIR/Tce_braker_RNAseq_run_4e_utrBMAll &


module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_RNAseq_run_4a_utrBMAll \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tce_braker_RNAseq_run_4a_utrBMAll/ALL_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_RNAseq_run_4a_utrBMAll/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tce_braker_RNAseq_run_4a_utrBMAll/augustus.hints.gtf --skipAllTraining
module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_RNAseq_run_4b_utrBMAll \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tce_braker_RNAseq_run_4b_utrBMAll/ALL_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_RNAseq_run_4b_utrBMAll/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tce_braker_RNAseq_run_4b_utrBMAll/augustus.hints.gtf --skipAllTraining

module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_RNAseq_run_4c_utrBMAll \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tce_braker_RNAseq_run_4c_utrBMAll/ALL_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_RNAseq_run_4c_utrBMAll/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tce_braker_RNAseq_run_4c_utrBMAll/augustus.hints.gtf --skipAllTraining

module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_RNAseq_run_4d_utrBMAll \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tce_braker_RNAseq_run_4d_utrBMAll/ALL_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_RNAseq_run_4d_utrBMAll/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tce_braker_RNAseq_run_4d_utrBMAll/augustus.hints.gtf --skipAllTraining

module load singularity/3.7.4
singularity exec --cleanenv --env-file $DIR/gushr_config_DJP_bigmem.txt --bind $DIR,$SDIR $DIR/braker-2.1.6_augustus-3.4.0_ubuntu-20.04.sif braker.pl --workingdir=$SDIR/Tce_braker_RNAseq_run_4e_utrBMAll \
			     --genome $DIR/genomes/$genome_pref".fasta.masked" --gff3 --bam $SDIR/Tce_braker_RNAseq_run_4e_utrBMAll/ALL_to_$genome_pref"_merged_sorted.bam" --softmasking --cores 46 \
				 --AUGUSTUS_CONFIG_PATH=$SDIR/Tce_braker_RNAseq_run_4e_utrBMAll/Augustus/config/ --GENEMARK_PATH=$DIR/gmes_linux_64 --species=Tce --addUTR=on \
				 --AUGUSTUS_hints_preds=$SDIR/Tce_braker_RNAseq_run_4e_utrBMAll/augustus.hints.gtf --skipAllTraining



cp -r $SDIR/Tce_braker_RNAseq_run_4a_utrBMAll $DIR
cp -r $SDIR/Tce_braker_RNAseq_run_4b_utrBMAll $DIR
cp -r $SDIR/Tce_braker_RNAseq_run_4c_utrBMAll $DIR
cp -r $SDIR/Tce_braker_RNAseq_run_4d_utrBMAll $DIR
cp -r $SDIR/Tce_braker_RNAseq_run_4e_utrBMAll $DIR

#8236105       braker_rnaseq_4a_utrBMAll_Tce              20:37:04  COMPLETED 
#8236128       braker_rnaseq_4b_utrBMAll_Tce              19:26:43  COMPLETED 
#8236129       braker_rnaseq_4c_utrBMAll_Tce              19:13:14  COMPLETED 
#8236130       braker_rnaseq_4d_utrBMAll_Tce              19:16:42  COMPLETED 
#8236131       braker_rnaseq_4e_utrBMAll_Tce              19:07:10  COMPLETED 



############################################################
### merge prot + RNASeq + UTR
# Tce_braker_prot_run_3e, Tce_braker_RNAseq_run_4e

mkdir /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/annotations_v2/
mkdir /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/annotations_v2/orig_merge

tar -zxf Tce_braker_prot_run_3e.tar.gz
./TSEBRA/bin/fix_gtf_ids.py --gtf Tce_braker_prot_run_3e/braker.gtf   --out Tce_braker_prot_run_3e/braker_fixed.gtf
./TSEBRA/bin/fix_gtf_ids.py --gtf Tce_braker_RNAseq_run_4e/braker.gtf --out Tce_braker_RNAseq_run_4e/braker_fixed.gtf
./TSEBRA/bin/fix_gtf_ids.py --gtf Tce_braker_RNAseq_run_4e_utrBMAll/braker_utr.gtf --out Tce_braker_RNAseq_run_4e_utrBMAll/braker_utr_fixed.gtf

#### no-UTR version pref_braker1.cfg (RNAseq evi more weight)
./TSEBRA/bin/tsebra.py -g Tce_braker_prot_run_3e/braker_fixed.gtf,Tce_braker_RNAseq_run_4e/braker_fixed.gtf \
-c TSEBRA/config/pref_braker1.cfg \
-e Tce_braker_prot_run_3e/hintsfile.gff,Tce_braker_RNAseq_run_4e/hintsfile.gff \
-o Tce_braker_prot_run_3e_Tce_braker_RNAseq_run_4e_combined_2.gtf

#### UTR version pref_braker1.cfg (RNAseq evi more weight)
./TSEBRA/bin/tsebra.py -g Tce_braker_prot_run_3e/braker_fixed.gtf,Tce_braker_RNAseq_run_4e/braker_fixed.gtf,Tce_braker_RNAseq_run_4e_utrBMAll/braker_utr_fixed.gtf \
-c TSEBRA/config/pref_braker1.cfg \
-e Tce_braker_prot_run_3e/hintsfile.gff,Tce_braker_RNAseq_run_4e/hintsfile.gff,Tce_braker_RNAseq_run_4e_utrBMAll/hintsfile.gff \
-o Tce_braker_prot_run_3e_Tce_braker_RNAseq_run_4eUTR_combined_2.gtf

### add to v2 
mv Tce_braker_prot_run_3e_Tce_braker_RNAseq_run_4e_combined_2.gtf    $DIR/annotations_v2/orig_merge
mv Tce_braker_prot_run_3e_Tce_braker_RNAseq_run_4eUTR_combined_2.gtf $DIR/annotations_v2/orig_merge

python3 Accessory_scripts/fix_braker_gtf.py -i $DIR/annotations_v2/orig_merge/Tce_braker_prot_run_3e_Tce_braker_RNAseq_run_4e_combined_2.gtf -d 6 -p Tce -G 
python3 Accessory_scripts/fix_braker_gtf.py -i $DIR/annotations_v2/orig_merge/Tce_braker_prot_run_3e_Tce_braker_RNAseq_run_4eUTR_combined_2.gtf -d 6 -p Tce -G 
python3 Accessory_scripts/UTR_GUSHR_fix.py -g $DIR/annotations_v2/orig_merge/Tce_braker_prot_run_3e_Tce_braker_RNAseq_run_4e_combined_2_FBGgi.gff -u $DIR/annotations_v2/orig_merge/Tce_braker_prot_run_3e_Tce_braker_RNAseq_run_4eUTR_combined_2_FBGgi.gff -m 1000 -G
cp $DIR/annotations_v2/orig_merge/Tce_braker_prot_run_3e_Tce_braker_RNAseq_run_4e_combined_2_FBGgiUGF1000M0.gff $DIR/annotations_v2
mv $DIR/annotations_v2/Tce_braker_prot_run_3e_Tce_braker_RNAseq_run_4e_combined_2_FBGgiUGF1000M0.gff $DIR/annotations_v2/Tce_braker_prot_run_3e_RNAseq_run_4e_combined_2_FBGgiUGF1000M0.gff 


### run BUSCO on the annotation

module load gcc
module load bedtools2/2.29.2

for i in $DIR/annotations_v2/Tce*_combined_2_FBGgiUGF1000M0.gff; do
echo $i
out_pref=`echo $i | sed 's/.gff.*//'`
echo $out_pref
awk '$3 == "gene"' $i > $out_pref"_genes.gff"
cut -f 1,4,5,9 $out_pref"_genes.gff"        | sort -k1,3 -u > $out_pref"_genes_coord.temp" ## as some rev orient genes
cut -f 4       $out_pref"_genes_coord.temp" | sed 's/.*ID=//' | sed 's/;.*//' > $out_pref"_genes_names.temp"
cut -f 1,2,3   $out_pref"_genes_coord.temp" > $out_pref"_genes_coord_2.temp" 
paste $out_pref"_genes_coord_2.temp" $out_pref"_genes_names.temp" > $out_pref"_genes.bed"
bedtools getfasta -fi genomes/$genome_pref.fasta -bed  $out_pref"_genes.bed" -name >  $out_pref"_genes.fa"
python3 Accessory_scripts/fasta_len_0.2.py $out_pref"_genes.fa"
rm $out_pref"_genes.gff"
rm $out_pref"_genes_coord.temp"
rm $out_pref"_genes_names.temp"
rm $out_pref"_genes_coord_2.temp"
rm $out_pref"_genes.bed"
done

### run BUSCO
module load singularity/3.7.4
for f in $DIR/annotations_v2/Tce*_genes.fa; do
out_f=`echo $f | sed 's/.fa/_BUSCO/' `
echo $f
echo $out_f
singularity exec --bind $DIR /work/FAC/FBM/DEE/tschwand/default/dparker/busco_v5.3.2_cv1.sif busco -c 30 -m genome -i $f -o $out_f -l insecta_odb10 --offline --download_path $DIR/busco_downloads/ 
done

mkdir annotations_v2/BUSCO
mv    work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/annotations_v2/* annotations_v2/BUSCO/

##### add missing busco

python3 Accessory_scripts/Add_missing_busco.py -m annotations_v2/BUSCO/Tce_braker_*_combined_2_FBGgiUGF1000M0_genes_BUSCO/run_insecta_odb10/missing_busco_list.tsv -b BUSCO/$genome_pref/ -g annotations_v2/Tce_*_combined_2_FBGgiUGF1000M0.gff 
#N missing BUSCO genes:
#87
#Number of complete busco found in genome busco run:
#69
























