#### mapping RNAseq to LR genomes

### stratagy. basically the same mapping proceedure as for braker annot, but with the mtDNA. Doing for all samples but only using matching samples for expression analyses (see R scripts)
### will store here /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/RNAseq_mapping
### and mapping stats: /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/


DIR='/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot'
SDIR='/scratch/dparker/annot'
genome_pref="Tpa_LRv5a_mtDNAv350"

#######################################################################################################
### map - STAR with --twopassMode Basic (to get the splice junction right)
### map all from Tpa and Tge to Tpa

### 39 conditions (NOT using single-end GN, LG, HD)
### 36 paired, 3 single end
### should get 136 paired bams, 12 single bams
### See TpaTge_Annot_RNAseq_samples.xlsx

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

##########################
### Index
### Build index for star

dcsrsoft use old # Switch to the 2021 stack
module load gcc
module load star/2.7.8a

mkdir $SDIR/RNAseq_mapping

STAR --runThreadN 12 \
     --runMode genomeGenerate \
     --genomeDir $SDIR/$genome_pref"_STAR" \
     --genomeFastaFiles $DIR"/genomes/"$genome_pref".fasta" \
     --genomeChrBinNbits 20 --limitGenomeGenerateRAM 79000000000





for R1_f in /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq/trimmed_*/*_R1_*qtrimmed.fq.gz; do
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

## WB single end stuff
Tpa_Tge_single_samples=(
Tpa_F_WB_Ad
Tge_F_WB_Ad
Tpa_M_WB_Ad	
)

for R1_f in /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq/trimmedSE*/*.fq.gz; do
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
     --outFileNamePrefix "/work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/RNAseq_mapping/"$out_prefix"_to_"$genome_pref
fi
done
done 


ls "/work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/RNAseq_mapping/"*"_to_"$genome_pref*"Aligned.sortedByCoord.out.bam" | wc -l
# 148

### STORE map stats HERE /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/
mkdir /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/to_$genome_pref
cp /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/RNAseq_mapping/*$genome_pref"Log.final.out" /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/to_$genome_pref
cd /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/
module load  singularity/3.8.5
singularity exec -e --bind /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq,/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/STAR_mapping_stats/ \
/work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/READS/RNAseq/multiqc_latest.sif multiqc \
to_$genome_pref/ -o to_$genome_pref"_STAR" --interactive
cp to_$genome_pref"_STAR/multiqc_report.html"  to_$genome_pref"_STAR_multiqc_report.html" 
rm -r to_$genome_pref"_STAR"
tar -czvf to_$genome_pref".tar.gz" to_$genome_pref
rm -r to_$genome_pref


################################################################################################################
#### get annot ready for htseq

mkdir /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/gff_exons
cd    /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/gff_exons
mkdir REFS

cp /nas/FAC/FBM/DEE/tschwand/timema_lr_genomes/D2c/timema_genomes/v2/*gff REFS
cp /nas/FAC/FBM/DEE/tschwand/timema_lr_genomes/D2c/timema_genomes/v2/*fasta REFS

module load gcc cufflinks/2.2.1

for gff in ./REFS/*_mtDNAv350_v2.2.gff; do
echo $gff
gff_1=`echo $gff   | sed 's/.gff/_wExon.gff/'` 
gff_2=`echo $gff_1 | sed 's/.gff/_wExon.gff/'`
gff_3=`echo $gff_2 | sed 's/.gff/_wExon.gff/'`
new_gff=`echo $gff | sed 's/.gff/_allwExon.gff/'`
fasta=`echo $gff | sed 's/_v2.*/.fasta/'`
out_trans=`echo $new_gff | sed 's/.gff/_trans.fa/'`
python3 ~/Gen_BioInf/gff_add_exon_to_feature.py -i $gff    -f rRNA
python3 ~/Gen_BioInf/gff_add_exon_to_feature.py -i $gff_1  -f ncRNA
python3 ~/Gen_BioInf/gff_add_exon_to_feature.py -i $gff_2  -f tRNA
cp $gff_3 $new_gff
#rm $gff $gff_1 $gff_2 $gff_3
gffread $new_gff -g $fasta -w $out_trans
grep ">" $out_trans | sed 's/.*gene=//' | sed 's/ .*//' | sort | uniq | wc -l
done


### count reads - STAR - HTSeq

mkdir /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/STAR_HTseq_out
cd    /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/STAR_HTseq_out

# get gffs (exon)
cp /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/gff_exons/REFS/*v2.2_allwExon.gff . 

module load gcc htseq/0.11.2

for f in "/work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/RNAseq_mapping/"replace*"_to_"$genome_pref*"Aligned.sortedByCoord.out.bam" ; do
gff=`echo $genome_pref"_v2.2_allwExon.gff"`
out_name=`echo $f | sed 's/Aligned.sortedByCoord.out.bam.*//' | sed 's/.*\///'`
echo $f
echo $gff
echo $out_name
htseq-count --order=pos --type=exon --idattr=gene_id --stranded=reverse --format=bam $f $gff > $out_name"_HTseq.txt"
done




