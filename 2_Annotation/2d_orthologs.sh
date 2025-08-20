## othologs

## stored here: data/gene_fasta
## unzip to use

###################################################################################################################################################
## Timema only 
##### filter for genes ON nice chromosomes -

## test cp /Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/fasta/Tps_LRv5b_aa_alllongest.fa.gz ./files

python3 Accessory_scripts/filter_isoforms_by_scaf.py \
-c data/read_counts/to_Tps_H2E.counts_wGLSL.csv \
-s scaf -i files/Tps_LRv5b_aa_alllongest.fa 
#scaf
#389
#{'Tps_LRv5b_scf10', 'Tps_LRv5b_scf8', 'Tps_LRv5b_scf4', 'Tps_LRv5b_scf9', 'Tps_LRv5b_scf5', 'Tps_LRv5b_scf6', 'Tps_LRv5b_scf2', 'Tps_LRv5b_scf1', 'Tps_LRv5b_scf7', 'Tps_LRv5b_scf11', 'Tps_LRv5b_scf3', 'Tps_LRv5b_scf12'}
#12
#34810
#37905
#
#Unwrapping fasta file
#Read 37390 sequences from files/Tps_LRv5b_aa_alllongest.fa
#number of seqs kept
#34383
#Note - excluded:
#112 ncRNAs
#255 tRNAs
#60 rRNAs
#that ARE on the wanted scafs


python3 Accessory_scripts/filter_isoforms_by_scaf.py \
-c data/read_counts/to_Tce_H2E.counts_wGLSLLG.csv \
-s LG -i files/Tce_LRv5a_aa_alllongest.fa 
#LG
#211
#{'HiCLG11', 'HiCLG8', 'HiCLG5', 'HiCLG12', 'HiCLG3', 'HiCLG13', 'HiCLG6', 'HiCLG10', 'HiCLG2', 'HiCLG4', 'HiCLG1', 'HiCLG7', 'HiCLG9'}
#13
#30155
#32202
#
#Unwrapping fasta file
#Read 31704 sequences from files/Tce_LRv5a_aa_alllongest.fa
#number of seqs kept
#29786
#Note - excluded:
#84 ncRNAs
#252 tRNAs
#33 rRNAs
#that ARE on the wanted scafs


python3 Accessory_scripts/filter_isoforms_by_scaf.py \
-c data/read_counts/to_Tpa_H2E.counts_wGLSL.csv \
-s scaf -i files/Tpa_LRv5a_aa_alllongest.fa 
#scaf
#151
#{'Tpa_LRv5a_scf13', 'Tpa_LRv5a_scf9', 'Tpa_LRv5a_scf5', 'Tpa_LRv5a_scf12', 'Tpa_LRv5a_scf3', 'Tpa_LRv5a_scf8', 'Tpa_LRv5a_scf6', 'Tpa_LRv5a_scf11', 'Tpa_LRv5a_scf4', 'Tpa_LRv5a_scf10', 'Tpa_LRv5a_scf14', 'Tpa_LRv5a_scf2', 'Tpa_LRv5a_scf7', 'Tpa_LRv5a_scf1'}
#14
#32271
#33719
#
#Unwrapping fasta file
#Read 33210 sequences from files/Tpa_LRv5a_aa_alllongest.fa
#number of seqs kept
#31838
#Note - excluded:
#112 ncRNAs
#278 tRNAs
#43 rRNAs
#that ARE on the wanted scafs


for f in *_onscafs.fa ; do
python3 Accessory_scripts/fasta_prot_trim_stop.py -f $f -R
done

mkdir  prot_seqs_onscafs
cp *_onscafs_aans.fa  prot_seqs_onscafs


#dcsrsoft use vitalit
#module add Bioinformatics/Software/vital-it
#module load Phylogeny/OrthoFinder/2.3.8 
#orthofinder -f  /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths/prot_seqs_onscafs

DIR='/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths'
module load gcc singularity
#singularity  pull docker://davidemms/orthofinder ## OrthoFinder version 2.5.5 Copyright (C) 2014 David Emms
singularity exec --bind $DIR $DIR/orthofinder_latest.sif orthofinder -f  /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths/prot_seqs_onscafs -t 10


## store
#data/orths/Orthogroups_TpsTpaTce_onscafs.txt 

#################################################################################################################################
#### add Bacillus to call orths

python3 Accessory_scripts/filter_isoforms_by_scaf_bacillus.py \
-c /work/FAC/FBM/DEE/tschwand/rnaseq_erc/dparker/Asex_males/output/orthologs/Brsri_droso_annot_gene.gtf -i Brsri_v3_aa_longest.fa 

#{'Brsri_v3_scf16', 'Brsri_v3_scf9_1', 'Brsri_v3_scf12', 'Brsri_v3_scf2', 'Brsri_v3_scf4_1', 'Brsri_v3_scf7', 'Brsri_v3_scf11', 'Brsri_v3_scf17', 'Brsri_v3_scf18', 'Brsri_v3_scf4_2', 'Brsri_v3_scf13', 'Brsri_v3_scf10', 'Brsri_v3_scf8', 'Brsri_v3_scf15', 'Brsri_v3_scf9_2', 'Brsri_v3_scf1', 'Brsri_v3_scf6', 'Brsri_v3_scf14', 'Brsri_v3_scf5', 'Brsri_v3_scf3'}
#20
#23540
#24058
#
#Unwrapping fasta file
#Read 24056 sequences from Brsri_v3_aa_longest.fa
#these genes should be in the output! ERROR!
#G23403
#these genes should be in the output! ERROR!
#G43363
#number of seqs kept
# both G23403 and G43363 have no CDS.

python3 Accessory_scripts/fasta_prot_trim_stop.py -f Brsri_v3_aa_longest_onscafs.fa -R


#dcsrsoft use vitalit
#module add Bioinformatics/Software/vital-it
#module load Phylogeny/OrthoFinder/2.3.8 
#orthofinder -f  /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths/prot_seqs_onscafs_wBacillus


DIR='/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths'
module load gcc singularity
#singularity  pull docker://davidemms/orthofinder ## OrthoFinder version 2.5.5 Copyright (C) 2014 David Emms
singularity exec --bind $DIR $DIR/orthofinder_latest.sif orthofinder -f  /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths/prot_seqs_onscafs_wBacillus -t 10


### store

#data/orths/Orthogroups_TpsTpaTceBrs_onscafs.txt 

#### plot

##  Ortho_chord_plots.R

## then invert coords as needed

python3 Accessory_scripts/invert_coords.py -i output/orthologs/Orthogroups_TpsTpaTceBrs_onscafs_Tps_c_1to1info.bed -c scf2
python3 Accessory_scripts/invert_coords.py -i output/orthologs/Orthogroups_TpsTpaTceBrs_onscafs_Tps_c_1to1info_scf2.bed -c scf3
python3 Accessory_scripts/invert_coords.py -i output/orthologs/Orthogroups_TpsTpaTceBrs_onscafs_Tps_c_1to1info_scf2_scf3.bed -c scf6
python3 Accessory_scripts/invert_coords.py -i output/orthologs/Orthogroups_TpsTpaTceBrs_onscafs_Tps_c_1to1info_scf2_scf3_scf6.bed -c scf8
python3 Accessory_scripts/invert_coords.py -i output/orthologs/Orthogroups_TpsTpaTceBrs_onscafs_Tps_c_1to1info_scf2_scf3_scf6_scf8.bed -c scf10
python3 Accessory_scripts/invert_coords.py -i output/orthologs/Orthogroups_TpsTpaTceBrs_onscafs_Tce_1to1info_LG.bed  -c LG11
python3 Accessory_scripts/invert_coords.py -i output/orthologs/Orthogroups_TpsTpaTceBrs_onscafs_Tce_1to1info_LG_LG11.bed  -c LG12
python3 Accessory_scripts/invert_coords.py -i output/orthologs/Orthogroups_TpsTpaTceBrs_onscafs_Tce_1to1info_LG_LG11_LG12.bed  -c LG13
python3 Accessory_scripts/invert_coords.py -i output/orthologs/Orthogroups_TpsTpaTceBrs_onscafs_Tpa_c_1to1info.bed  -c Scf10

















































































############# blast to Tch - 
### Tch genome - Should have been published with Nosil P, Villoutreix R, de Carvalho CF, Feder JL, Parchman TL, Gompert Z (2020). Ecology shapes epistasis in a genotype–phenotype–fitness map for stick insect colour. Nature Ecology & Evolution 4: 1673–1684., but they didn't do it
### Zach sent me a copy - stored here: /nas/FAC/FBM/DEE/tschwand/timema_lr_genomes/D2c/timema_chumash_29Feb2020_N4ago.fasta.gz
### THIS GENOME has unresoved haplotypes... NOT USING - SEE 3b_Tch_Xchr.sh


cd /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths/
cp /nas/FAC/FBM/DEE/tschwand/timema_lr_genomes/D2c/timema_chumash_29Feb2020_N4ago.fasta.gz .
gzip -d timema_chumash_29Feb2020_N4ago.fasta.gz 

module load gcc        
module load blast-plus

makeblastdb -in timema_chumash_29Feb2020_N4ago.fasta -out Tch_db  -dbtype nucl
tblastn -query prot_seqs/Tce_LRv5a_aa_alllongest_aans.fa -db Tch_db  -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out Tce_LRv5a_aa_alllongest_to_Tch_blast_out.csv
tblastn -query prot_seqs/Tpa_LRv5a_aa_alllongest_aans.fa -db Tch_db  -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out Tpa_LRv5a_aa_alllongest_to_Tch_blast_out.csv
tblastn -query prot_seqs/Tps_LRv5b_aa_alllongest_aans.fa -db Tch_db  -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out Tps_LRv5b_aa_alllongest_to_Tch_blast_out.csv


## split for blasting
python3 ~/DJP_py_scripts/get_missed_blasts_at_end_of_file.py  prot_seqs/Tce_LRv5a_aa_alllongest_aans.fa  Tce_LRv5a_aa_alllongest_to_Tch_blast_out.csv temp
python3 ~/DJP_py_scripts/get_missed_blasts_at_end_of_file.py  prot_seqs/Tpa_LRv5a_aa_alllongest_aans.fa  Tpa_LRv5a_aa_alllongest_to_Tch_blast_out.csv temp
python3 ~/DJP_py_scripts/get_missed_blasts_at_end_of_file.py  prot_seqs/Tps_LRv5b_aa_alllongest_aans.fa  Tps_LRv5b_aa_alllongest_to_Tch_blast_out.csv temp

python3 ~/DJP_py_scripts/get_seq_with_no_blast_hit.py -b Tce_LRv5a_aa_alllongest_to_Tch_blast_out_wNA.csv -f prot_seqs/Tce_LRv5a_aa_alllongest_aans.fa -o R2_Tce_LRv5a_aa_alllongest_aans -n 2000
python3 ~/DJP_py_scripts/get_seq_with_no_blast_hit.py -b Tpa_LRv5a_aa_alllongest_to_Tch_blast_out_wNA.csv -f prot_seqs/Tpa_LRv5a_aa_alllongest_aans.fa -o R2_Tpa_LRv5a_aa_alllongest_aans -n 2000
python3 ~/DJP_py_scripts/get_seq_with_no_blast_hit.py -b Tps_LRv5b_aa_alllongest_to_Tch_blast_out_wNA.csv -f prot_seqs/Tps_LRv5b_aa_alllongest_aans.fa -o R2_Tps_LRv5b_aa_alllongest_aans -n 2000




tblastn -query R2_Tce_LRv5a_aa_alllongest_aans0_s.fasta -db Tch_db  -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out R2_Tce_LRv5a_aa_alllongest_aans0_to_Tch_blast_out.csv
tblastn -query R2_Tpa_LRv5a_aa_alllongest_aans0_s.fasta -db Tch_db  -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out R2_Tpa_LRv5a_aa_alllongest_aans0_to_Tch_blast_out.csv	
tblastn -query R2_Tps_LRv5b_aa_alllongest_aans0_s.fasta -db Tch_db  -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out R2_Tps_LRv5b_aa_alllongest_aans0_to_Tch_blast_out.csv	



sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans1/g' TceTch_R2_0_blast.sh > TceTch_R2_1_blast.sh 
sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans2/g' TceTch_R2_0_blast.sh > TceTch_R2_2_blast.sh
sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans3/g' TceTch_R2_0_blast.sh > TceTch_R2_3_blast.sh
sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans4/g' TceTch_R2_0_blast.sh > TceTch_R2_4_blast.sh
sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans5/g' TceTch_R2_0_blast.sh > TceTch_R2_5_blast.sh
sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans6/g' TceTch_R2_0_blast.sh > TceTch_R2_6_blast.sh
sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans7/g' TceTch_R2_0_blast.sh > TceTch_R2_7_blast.sh
sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans8/g' TceTch_R2_0_blast.sh > TceTch_R2_8_blast.sh
sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans9/g' TceTch_R2_0_blast.sh > TceTch_R2_9_blast.sh
sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans10/g' TceTch_R2_0_blast.sh > TceTch_R2_10_blast.sh
sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans11/g' TceTch_R2_0_blast.sh > TceTch_R2_11_blast.sh
sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans12/g' TceTch_R2_0_blast.sh > TceTch_R2_12_blast.sh
sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans13/g' TceTch_R2_0_blast.sh > TceTch_R2_13_blast.sh
sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans14/g' TceTch_R2_0_blast.sh > TceTch_R2_14_blast.sh
sed 's/R2_Tce_LRv5a_aa_alllongest_aans0/R2_Tce_LRv5a_aa_alllongest_aans15/g' TceTch_R2_0_blast.sh > TceTch_R2_15_blast.sh

sbatch TceTch_R2_0_blast.sh
sbatch TceTch_R2_1_blast.sh
sbatch TceTch_R2_2_blast.sh
sbatch TceTch_R2_3_blast.sh
sbatch TceTch_R2_4_blast.sh
sbatch TceTch_R2_5_blast.sh
sbatch TceTch_R2_6_blast.sh
sbatch TceTch_R2_7_blast.sh
sbatch TceTch_R2_8_blast.sh
sbatch TceTch_R2_9_blast.sh
sbatch TceTch_R2_10_blast.sh
sbatch TceTch_R2_11_blast.sh
sbatch TceTch_R2_12_blast.sh
sbatch TceTch_R2_13_blast.sh
sbatch TceTch_R2_14_blast.sh
sbatch TceTch_R2_15_blast.sh


sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans1/g' TpsTch_R2_0_blast.sh > TpsTch_R2_1_blast.sh 
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans2/g' TpsTch_R2_0_blast.sh > TpsTch_R2_2_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans3/g' TpsTch_R2_0_blast.sh > TpsTch_R2_3_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans4/g' TpsTch_R2_0_blast.sh > TpsTch_R2_4_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans5/g' TpsTch_R2_0_blast.sh > TpsTch_R2_5_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans6/g' TpsTch_R2_0_blast.sh > TpsTch_R2_6_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans7/g' TpsTch_R2_0_blast.sh > TpsTch_R2_7_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans8/g' TpsTch_R2_0_blast.sh > TpsTch_R2_8_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans9/g' TpsTch_R2_0_blast.sh > TpsTch_R2_9_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans10/g' TpsTch_R2_0_blast.sh > TpsTch_R2_10_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans11/g' TpsTch_R2_0_blast.sh > TpsTch_R2_11_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans12/g' TpsTch_R2_0_blast.sh > TpsTch_R2_12_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans13/g' TpsTch_R2_0_blast.sh > TpsTch_R2_13_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans14/g' TpsTch_R2_0_blast.sh > TpsTch_R2_14_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans15/g' TpsTch_R2_0_blast.sh > TpsTch_R2_15_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans16/g' TpsTch_R2_0_blast.sh > TpsTch_R2_16_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans17/g' TpsTch_R2_0_blast.sh > TpsTch_R2_17_blast.sh
sed 's/R2_Tps_LRv5b_aa_alllongest_aans0/R2_Tps_LRv5b_aa_alllongest_aans18/g' TpsTch_R2_0_blast.sh > TpsTch_R2_18_blast.sh


sbatch TpsTch_R2_0_blast.sh
sbatch TpsTch_R2_1_blast.sh
sbatch TpsTch_R2_2_blast.sh
sbatch TpsTch_R2_3_blast.sh
sbatch TpsTch_R2_4_blast.sh
sbatch TpsTch_R2_5_blast.sh
sbatch TpsTch_R2_6_blast.sh
sbatch TpsTch_R2_7_blast.sh
sbatch TpsTch_R2_8_blast.sh
sbatch TpsTch_R2_9_blast.sh
sbatch TpsTch_R2_10_blast.sh
sbatch TpsTch_R2_11_blast.sh
sbatch TpsTch_R2_12_blast.sh
sbatch TpsTch_R2_13_blast.sh
sbatch TpsTch_R2_14_blast.sh
sbatch TpsTch_R2_15_blast.sh
sbatch TpsTch_R2_16_blast.sh
sbatch TpsTch_R2_17_blast.sh
sbatch TpsTch_R2_18_blast.sh



sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans1/g' TpaTch_R2_0_blast.sh > TpaTch_R2_1_blast.sh 
sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans2/g' TpaTch_R2_0_blast.sh > TpaTch_R2_2_blast.sh
sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans3/g' TpaTch_R2_0_blast.sh > TpaTch_R2_3_blast.sh
sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans4/g' TpaTch_R2_0_blast.sh > TpaTch_R2_4_blast.sh
sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans5/g' TpaTch_R2_0_blast.sh > TpaTch_R2_5_blast.sh
sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans6/g' TpaTch_R2_0_blast.sh > TpaTch_R2_6_blast.sh
sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans7/g' TpaTch_R2_0_blast.sh > TpaTch_R2_7_blast.sh
sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans8/g' TpaTch_R2_0_blast.sh > TpaTch_R2_8_blast.sh
sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans9/g' TpaTch_R2_0_blast.sh > TpaTch_R2_9_blast.sh
sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans10/g' TpaTch_R2_0_blast.sh > TpaTch_R2_10_blast.sh
sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans11/g' TpaTch_R2_0_blast.sh > TpaTch_R2_11_blast.sh
sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans12/g' TpaTch_R2_0_blast.sh > TpaTch_R2_12_blast.sh
sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans13/g' TpaTch_R2_0_blast.sh > TpaTch_R2_13_blast.sh
sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans14/g' TpaTch_R2_0_blast.sh > TpaTch_R2_14_blast.sh
sed 's/R2_Tpa_LRv5a_aa_alllongest_aans0/R2_Tpa_LRv5a_aa_alllongest_aans15/g' TpaTch_R2_0_blast.sh > TpaTch_R2_15_blast.sh

sbatch TpaTch_R2_0_blast.sh
sbatch TpaTch_R2_1_blast.sh
sbatch TpaTch_R2_2_blast.sh
sbatch TpaTch_R2_3_blast.sh
sbatch TpaTch_R2_4_blast.sh
sbatch TpaTch_R2_5_blast.sh
sbatch TpaTch_R2_6_blast.sh
sbatch TpaTch_R2_7_blast.sh
sbatch TpaTch_R2_8_blast.sh
sbatch TpaTch_R2_9_blast.sh
sbatch TpaTch_R2_10_blast.sh
sbatch TpaTch_R2_11_blast.sh
sbatch TpaTch_R2_12_blast.sh
sbatch TpaTch_R2_13_blast.sh
sbatch TpaTch_R2_14_blast.sh
sbatch TpaTch_R2_15_blast.sh


#### collect up




for f in ./R2*_Tch_blast_out.csv ; do
echo $f
fasta_f=`echo $f | sed 's/_to_Tch_blast_out.csv/_s.fasta/'`
echo $fasta_f
python3 ~/DJP_py_scripts/get_missed_blasts_at_end_of_file.py  $fasta_f $f temp
done


#(base) [dparker@curnagl orths]$ cat  *Tce_*_wNA.csv | cut -f 1 -d ',' | sort | uniq | wc -l
#14150
#(base) [dparker@curnagl orths]$ cat  *Tps_*_wNA.csv | cut -f 1 -d ',' | sort | uniq | wc -l
#7531
#(base) [dparker@curnagl orths]$ cat  *Tpa_*_wNA.csv | cut -f 1 -d ',' | sort | uniq | wc -l
#22705


cat  *Tce_*_wNA.csv > R1R2_Tce_LRv5a_aa_alllongest_to_Tch_blast_out_wNA.csv
cat  *Tpa_*_wNA.csv > R1R2_Tpa_LRv5a_aa_alllongest_to_Tch_blast_out_wNA.csv
cat  *Tps_*_wNA.csv > R1R2_Tps_LRv5b_aa_alllongest_to_Tch_blast_out_wNA.csv



python3 ~/DJP_py_scripts/get_seq_with_no_blast_hit.py -b R1R2_Tce_LRv5a_aa_alllongest_to_Tch_blast_out_wNA.csv -f prot_seqs/Tce_LRv5a_aa_alllongest_aans.fa -o R3_Tce_LRv5a_aa_alllongest_aans -n 500
python3 ~/DJP_py_scripts/get_seq_with_no_blast_hit.py -b R1R2_Tpa_LRv5a_aa_alllongest_to_Tch_blast_out_wNA.csv -f prot_seqs/Tpa_LRv5a_aa_alllongest_aans.fa -o R3_Tpa_LRv5a_aa_alllongest_aans -n 500
python3 ~/DJP_py_scripts/get_seq_with_no_blast_hit.py -b R1R2_Tps_LRv5b_aa_alllongest_to_Tch_blast_out_wNA.csv -f prot_seqs/Tps_LRv5b_aa_alllongest_aans.fa -o R3_Tps_LRv5b_aa_alllongest_aans -n 500


tblastn -query R3_Tce_LRv5a_aa_alllongest_aans0_s.fasta -db Tch_db  -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out R3_Tce_LRv5a_aa_alllongest_aans0_to_Tch_blast_out.csv
tblastn -query R3_Tpa_LRv5a_aa_alllongest_aans0_s.fasta -db Tch_db  -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out R3_Tpa_LRv5a_aa_alllongest_aans0_to_Tch_blast_out.csv
tblastn -query R3_Tps_LRv5b_aa_alllongest_aans0_s.fasta -db Tch_db  -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out R3_Tps_LRv5b_aa_alllongest_aans0_to_Tch_blast_out.csv




sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans1/g' TceTch_R3_0_blast.sh > TceTch_R3_1_blast.sh 
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans2/g' TceTch_R3_0_blast.sh > TceTch_R3_2_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans3/g' TceTch_R3_0_blast.sh > TceTch_R3_3_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans4/g' TceTch_R3_0_blast.sh > TceTch_R3_4_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans5/g' TceTch_R3_0_blast.sh > TceTch_R3_5_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans6/g' TceTch_R3_0_blast.sh > TceTch_R3_6_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans7/g' TceTch_R3_0_blast.sh > TceTch_R3_7_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans8/g' TceTch_R3_0_blast.sh > TceTch_R3_8_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans9/g' TceTch_R3_0_blast.sh > TceTch_R3_9_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans10/g' TceTch_R3_0_blast.sh > TceTch_R3_10_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans11/g' TceTch_R3_0_blast.sh > TceTch_R3_11_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans12/g' TceTch_R3_0_blast.sh > TceTch_R3_12_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans13/g' TceTch_R3_0_blast.sh > TceTch_R3_13_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans14/g' TceTch_R3_0_blast.sh > TceTch_R3_14_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans15/g' TceTch_R3_0_blast.sh > TceTch_R3_15_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans16/g' TceTch_R3_0_blast.sh > TceTch_R3_16_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans17/g' TceTch_R3_0_blast.sh > TceTch_R3_17_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans18/g' TceTch_R3_0_blast.sh > TceTch_R3_18_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans19/g' TceTch_R3_0_blast.sh > TceTch_R3_19_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans20/g' TceTch_R3_0_blast.sh > TceTch_R3_20_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans21/g' TceTch_R3_0_blast.sh > TceTch_R3_21_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans22/g' TceTch_R3_0_blast.sh > TceTch_R3_22_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans23/g' TceTch_R3_0_blast.sh > TceTch_R3_23_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans24/g' TceTch_R3_0_blast.sh > TceTch_R3_24_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans25/g' TceTch_R3_0_blast.sh > TceTch_R3_25_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans26/g' TceTch_R3_0_blast.sh > TceTch_R3_26_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans27/g' TceTch_R3_0_blast.sh > TceTch_R3_27_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans28/g' TceTch_R3_0_blast.sh > TceTch_R3_28_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans29/g' TceTch_R3_0_blast.sh > TceTch_R3_29_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans30/g' TceTch_R3_0_blast.sh > TceTch_R3_30_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans31/g' TceTch_R3_0_blast.sh > TceTch_R3_31_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans32/g' TceTch_R3_0_blast.sh > TceTch_R3_32_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans33/g' TceTch_R3_0_blast.sh > TceTch_R3_33_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans34/g' TceTch_R3_0_blast.sh > TceTch_R3_34_blast.sh
sed 's/R3_Tce_LRv5a_aa_alllongest_aans0/R3_Tce_LRv5a_aa_alllongest_aans35/g' TceTch_R3_0_blast.sh > TceTch_R3_35_blast.sh


sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans1/g' TpaTch_R3_0_blast.sh > TpaTch_R3_1_blast.sh 
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans2/g' TpaTch_R3_0_blast.sh > TpaTch_R3_2_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans3/g' TpaTch_R3_0_blast.sh > TpaTch_R3_3_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans4/g' TpaTch_R3_0_blast.sh > TpaTch_R3_4_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans5/g' TpaTch_R3_0_blast.sh > TpaTch_R3_5_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans6/g' TpaTch_R3_0_blast.sh > TpaTch_R3_6_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans7/g' TpaTch_R3_0_blast.sh > TpaTch_R3_7_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans8/g' TpaTch_R3_0_blast.sh > TpaTch_R3_8_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans9/g' TpaTch_R3_0_blast.sh > TpaTch_R3_9_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans10/g' TpaTch_R3_0_blast.sh > TpaTch_R3_10_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans11/g' TpaTch_R3_0_blast.sh > TpaTch_R3_11_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans12/g' TpaTch_R3_0_blast.sh > TpaTch_R3_12_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans13/g' TpaTch_R3_0_blast.sh > TpaTch_R3_13_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans14/g' TpaTch_R3_0_blast.sh > TpaTch_R3_14_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans15/g' TpaTch_R3_0_blast.sh > TpaTch_R3_15_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans16/g' TpaTch_R3_0_blast.sh > TpaTch_R3_16_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans17/g' TpaTch_R3_0_blast.sh > TpaTch_R3_17_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans18/g' TpaTch_R3_0_blast.sh > TpaTch_R3_18_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans19/g' TpaTch_R3_0_blast.sh > TpaTch_R3_19_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans20/g' TpaTch_R3_0_blast.sh > TpaTch_R3_20_blast.sh
sed 's/R3_Tpa_LRv5a_aa_alllongest_aans0/R3_Tpa_LRv5a_aa_alllongest_aans21/g' TpaTch_R3_0_blast.sh > TpaTch_R3_21_blast.sh


sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans1/g' TpsTch_R3_0_blast.sh > TpsTch_R3_1_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans2/g' TpsTch_R3_0_blast.sh > TpsTch_R3_2_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans3/g' TpsTch_R3_0_blast.sh > TpsTch_R3_3_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans4/g' TpsTch_R3_0_blast.sh > TpsTch_R3_4_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans5/g' TpsTch_R3_0_blast.sh > TpsTch_R3_5_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans6/g' TpsTch_R3_0_blast.sh > TpsTch_R3_6_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans7/g' TpsTch_R3_0_blast.sh > TpsTch_R3_7_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans8/g' TpsTch_R3_0_blast.sh > TpsTch_R3_8_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans9/g' TpsTch_R3_0_blast.sh > TpsTch_R3_9_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans10/g' TpsTch_R3_0_blast.sh > TpsTch_R3_10_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans11/g' TpsTch_R3_0_blast.sh > TpsTch_R3_11_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans12/g' TpsTch_R3_0_blast.sh > TpsTch_R3_12_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans13/g' TpsTch_R3_0_blast.sh > TpsTch_R3_13_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans14/g' TpsTch_R3_0_blast.sh > TpsTch_R3_14_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans15/g' TpsTch_R3_0_blast.sh > TpsTch_R3_15_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans16/g' TpsTch_R3_0_blast.sh > TpsTch_R3_16_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans17/g' TpsTch_R3_0_blast.sh > TpsTch_R3_17_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans18/g' TpsTch_R3_0_blast.sh > TpsTch_R3_18_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans19/g' TpsTch_R3_0_blast.sh > TpsTch_R3_19_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans20/g' TpsTch_R3_0_blast.sh > TpsTch_R3_20_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans21/g' TpsTch_R3_0_blast.sh > TpsTch_R3_21_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans22/g' TpsTch_R3_0_blast.sh > TpsTch_R3_22_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans23/g' TpsTch_R3_0_blast.sh > TpsTch_R3_23_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans24/g' TpsTch_R3_0_blast.sh > TpsTch_R3_24_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans25/g' TpsTch_R3_0_blast.sh > TpsTch_R3_25_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans26/g' TpsTch_R3_0_blast.sh > TpsTch_R3_26_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans27/g' TpsTch_R3_0_blast.sh > TpsTch_R3_27_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans28/g' TpsTch_R3_0_blast.sh > TpsTch_R3_28_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans29/g' TpsTch_R3_0_blast.sh > TpsTch_R3_29_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans30/g' TpsTch_R3_0_blast.sh > TpsTch_R3_30_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans31/g' TpsTch_R3_0_blast.sh > TpsTch_R3_31_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans32/g' TpsTch_R3_0_blast.sh > TpsTch_R3_32_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans33/g' TpsTch_R3_0_blast.sh > TpsTch_R3_33_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans34/g' TpsTch_R3_0_blast.sh > TpsTch_R3_34_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans35/g' TpsTch_R3_0_blast.sh > TpsTch_R3_35_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans36/g' TpsTch_R3_0_blast.sh > TpsTch_R3_36_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans37/g' TpsTch_R3_0_blast.sh > TpsTch_R3_37_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans38/g' TpsTch_R3_0_blast.sh > TpsTch_R3_38_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans39/g' TpsTch_R3_0_blast.sh > TpsTch_R3_39_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans40/g' TpsTch_R3_0_blast.sh > TpsTch_R3_40_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans41/g' TpsTch_R3_0_blast.sh > TpsTch_R3_41_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans42/g' TpsTch_R3_0_blast.sh > TpsTch_R3_42_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans43/g' TpsTch_R3_0_blast.sh > TpsTch_R3_43_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans44/g' TpsTch_R3_0_blast.sh > TpsTch_R3_44_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans45/g' TpsTch_R3_0_blast.sh > TpsTch_R3_45_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans46/g' TpsTch_R3_0_blast.sh > TpsTch_R3_46_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans47/g' TpsTch_R3_0_blast.sh > TpsTch_R3_47_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans48/g' TpsTch_R3_0_blast.sh > TpsTch_R3_48_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans49/g' TpsTch_R3_0_blast.sh > TpsTch_R3_49_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans50/g' TpsTch_R3_0_blast.sh > TpsTch_R3_50_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans51/g' TpsTch_R3_0_blast.sh > TpsTch_R3_51_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans52/g' TpsTch_R3_0_blast.sh > TpsTch_R3_52_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans53/g' TpsTch_R3_0_blast.sh > TpsTch_R3_53_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans54/g' TpsTch_R3_0_blast.sh > TpsTch_R3_54_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans55/g' TpsTch_R3_0_blast.sh > TpsTch_R3_55_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans56/g' TpsTch_R3_0_blast.sh > TpsTch_R3_56_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans57/g' TpsTch_R3_0_blast.sh > TpsTch_R3_57_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans58/g' TpsTch_R3_0_blast.sh > TpsTch_R3_58_blast.sh
sed 's/R3_Tps_LRv5b_aa_alllongest_aans0/R3_Tps_LRv5b_aa_alllongest_aans59/g' TpsTch_R3_0_blast.sh > TpsTch_R3_59_blast.sh



sbatch TceTch_R3_1_blast.sh
sbatch TceTch_R3_2_blast.sh
sbatch TceTch_R3_3_blast.sh
sbatch TceTch_R3_4_blast.sh
sbatch TceTch_R3_5_blast.sh
sbatch TceTch_R3_6_blast.sh
sbatch TceTch_R3_7_blast.sh
sbatch TceTch_R3_8_blast.sh
sbatch TceTch_R3_9_blast.sh
sbatch TceTch_R3_10_blast.sh
sbatch TceTch_R3_11_blast.sh
sbatch TceTch_R3_12_blast.sh
sbatch TceTch_R3_13_blast.sh
sbatch TceTch_R3_14_blast.sh
sbatch TceTch_R3_15_blast.sh
sbatch TceTch_R3_16_blast.sh
sbatch TceTch_R3_17_blast.sh
sbatch TceTch_R3_18_blast.sh
sbatch TceTch_R3_19_blast.sh
sbatch TceTch_R3_20_blast.sh
sbatch TceTch_R3_21_blast.sh
sbatch TceTch_R3_22_blast.sh
sbatch TceTch_R3_23_blast.sh
sbatch TceTch_R3_24_blast.sh
sbatch TceTch_R3_25_blast.sh
sbatch TceTch_R3_26_blast.sh
sbatch TceTch_R3_27_blast.sh
sbatch TceTch_R3_28_blast.sh
sbatch TceTch_R3_29_blast.sh
sbatch TceTch_R3_30_blast.sh
sbatch TceTch_R3_31_blast.sh
sbatch TceTch_R3_32_blast.sh
sbatch TceTch_R3_33_blast.sh
sbatch TceTch_R3_34_blast.sh
sbatch TceTch_R3_35_blast.sh



sbatch TpaTch_R3_1_blast.sh
sbatch TpaTch_R3_2_blast.sh
sbatch TpaTch_R3_3_blast.sh
sbatch TpaTch_R3_4_blast.sh
sbatch TpaTch_R3_5_blast.sh
sbatch TpaTch_R3_6_blast.sh
sbatch TpaTch_R3_7_blast.sh
sbatch TpaTch_R3_8_blast.sh
sbatch TpaTch_R3_9_blast.sh
sbatch TpaTch_R3_10_blast.sh
sbatch TpaTch_R3_11_blast.sh
sbatch TpaTch_R3_12_blast.sh
sbatch TpaTch_R3_13_blast.sh
sbatch TpaTch_R3_14_blast.sh
sbatch TpaTch_R3_15_blast.sh
sbatch TpaTch_R3_16_blast.sh
sbatch TpaTch_R3_17_blast.sh
sbatch TpaTch_R3_18_blast.sh
sbatch TpaTch_R3_19_blast.sh
sbatch TpaTch_R3_20_blast.sh
sbatch TpaTch_R3_21_blast.sh

sbatch TpsTch_R3_1_blast.sh
sbatch TpsTch_R3_2_blast.sh
sbatch TpsTch_R3_3_blast.sh
sbatch TpsTch_R3_4_blast.sh
sbatch TpsTch_R3_5_blast.sh
sbatch TpsTch_R3_6_blast.sh
sbatch TpsTch_R3_7_blast.sh
sbatch TpsTch_R3_8_blast.sh
sbatch TpsTch_R3_9_blast.sh
sbatch TpsTch_R3_10_blast.sh
sbatch TpsTch_R3_11_blast.sh
sbatch TpsTch_R3_12_blast.sh
sbatch TpsTch_R3_13_blast.sh
sbatch TpsTch_R3_14_blast.sh
sbatch TpsTch_R3_15_blast.sh
sbatch TpsTch_R3_16_blast.sh
sbatch TpsTch_R3_17_blast.sh
sbatch TpsTch_R3_18_blast.sh
sbatch TpsTch_R3_19_blast.sh
sbatch TpsTch_R3_20_blast.sh
sbatch TpsTch_R3_21_blast.sh
sbatch TpsTch_R3_22_blast.sh
sbatch TpsTch_R3_23_blast.sh
sbatch TpsTch_R3_24_blast.sh
sbatch TpsTch_R3_25_blast.sh
sbatch TpsTch_R3_26_blast.sh
sbatch TpsTch_R3_27_blast.sh
sbatch TpsTch_R3_28_blast.sh
sbatch TpsTch_R3_29_blast.sh
sbatch TpsTch_R3_30_blast.sh
sbatch TpsTch_R3_31_blast.sh
sbatch TpsTch_R3_32_blast.sh
sbatch TpsTch_R3_33_blast.sh
sbatch TpsTch_R3_34_blast.sh
sbatch TpsTch_R3_35_blast.sh
sbatch TpsTch_R3_36_blast.sh
sbatch TpsTch_R3_37_blast.sh
sbatch TpsTch_R3_38_blast.sh
sbatch TpsTch_R3_39_blast.sh
sbatch TpsTch_R3_40_blast.sh
sbatch TpsTch_R3_41_blast.sh
sbatch TpsTch_R3_42_blast.sh
sbatch TpsTch_R3_43_blast.sh
sbatch TpsTch_R3_44_blast.sh
sbatch TpsTch_R3_45_blast.sh
sbatch TpsTch_R3_46_blast.sh
sbatch TpsTch_R3_47_blast.sh
sbatch TpsTch_R3_48_blast.sh
sbatch TpsTch_R3_49_blast.sh
sbatch TpsTch_R3_50_blast.sh
sbatch TpsTch_R3_51_blast.sh
sbatch TpsTch_R3_52_blast.sh
sbatch TpsTch_R3_53_blast.sh
sbatch TpsTch_R3_54_blast.sh
sbatch TpsTch_R3_55_blast.sh
sbatch TpsTch_R3_56_blast.sh
sbatch TpsTch_R3_57_blast.sh
sbatch TpsTch_R3_58_blast.sh
sbatch TpsTch_R3_59_blast.sh

#
#
#cut -f 1 -d "," X_soft_HOG_longest_aa_to_Brsri_v3_blast_out.csv | sort | uniq | wc -l
##209
#
#python3 ~/Documents/Gen_BioInf/blast_get_most_likely_hit.py -i X_soft_HOG_longest_aa_to_Brsri_v3_blast_out.csv -o X_soft_HOG_longest_aa_to_Brsri_v3_blast_out_e1e-20_q80 -f X_soft_HOG_longest_aa.fa -e 1e-20 -v 10 -j 5 -q 80
#python3 ~/Documents/Gen_BioInf/blast_get_most_likely_hit.py -i X_soft_HOG_longest_aa_to_Brsri_v3_blast_out.csv -o X_soft_HOG_longest_aa_to_Brsri_v3_blast_out_e1e-20_q50 -f X_soft_HOG_longest_aa.fa -e 1e-20 -v 10 -j 5 -q 50
#
#
# drp22jhz$ python3 Accessory_scripts/X_orths.py -x data/read_counts/ -r to_Tce_H2E.counts_wGLSLLG.csv,to_Tps_H2E.counts_wGLSL.csv,to_Tpa_H2E.counts_wGLSL.csv -g output/orthologs/Orthogroups.txt 
#




































#### OLD
#### get from Dryococelus australis; Stuart OP, Cleave R, Magrath MJL, Mikheyev AS (2023). Genome of the Lord Howe Island Stick Insect Reveals a Highly Conserved Phasmid X Chromosome. Genome Biol Evol 15.
#
#mkdir -p /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths
#cd       /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths
#
#curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_029891345.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCA_029891345.1.zip" -H "Accept: application/zip"
#unzip GCA_029891345.1.zip
#
#
#### get prots - just doing Tce, Tpa, Tps and Dryococelus australis
#
#mkdir files
#
#cp  /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/fasta/Tps_LRv5b_aa_alllongest.fa.gz ./files
#cp  /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/fasta/Tpa_LRv5a_aa_alllongest.fa.gz ./files
#cp  /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/Timema_LR_genomic_code/output/fasta/Tce_LRv5a_aa_alllongest.fa.gz ./files
#cp ncbi_dataset/data/GCA_029891345.1/protein.faa  ./files
#
#cd /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths/ncbi_dataset/data/GCA_029891345.1
#
#module load gcc cufflinks/2.2.1 
#gffread genomic.gff  -g  GCA_029891345.1_Daus_2.0_genomic.fna  -y Daus_prot.fa
#python3 ~/DJP_py_scripts/Take_longest_isoform_from_gffread.py -i  Daus_prot.fa -o _longest.fa
#
#grep ">" Daus_prot.fa_longest.fa | wc -l 
## 33111
#grep ">" Daus_prot.fa | wc -l 
## 33111
#
#### ^^^ not needed
#
#cd       /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths
#cp  /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths/ncbi_dataset/data/GCA_029891345.1/Daus_prot.fa_longest.fa files/
#gzip -d files/*gz
#
#
#
#
#######################################
#####
#
#for f in /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths/files/* ; do
#	
#python3 ~/DJP_py_scripts/fasta_prot_trim_stop.py -f $f -R
#done
#
#
#mkdir  /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths/prot_seqs
#cp /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths/files/*_aans.fa  /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths/prot_seqs
#
#
#dcsrsoft use vitalit
#module add Bioinformatics/Software/vital-it
#module load Phylogeny/OrthoFinder/2.3.8 
#orthofinder -f  /work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/orths/prot_seqs
#
#
#






#
#
#
#
#
#
#
#
#
#
###### run orthologer
#module load singularity/3.8.5 
#singularity  pull docker://ezlabgva/orthologer:v3.0.0
#for x in $(ls files/*.fa); do echo "+$(basename $x .fa) $x"; done > mydata.txt
#
#
#docker run -u $(id -u) -v $(pwd):/odbwork ezlabgva/orthologer:v3.0.0  ./orthologer.sh  manage -f mydata.txt
#
#
#singularity exec --bind $pwd orthologer_v3.0.0.sif ./orthologer.sh --workingdir=$pwd manage -f mydata.txt
#
#singularity exec  --bind $pwd -u $(id -u) -v $(pwd):/odbwork ezlabgva/orthologer:v3.0.0 setup_odb.sh 
#
#export PATH=${LEM_ORTHO_ROOT}/workflow/scripts:$PATH












