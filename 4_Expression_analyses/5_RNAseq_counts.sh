#### Store results from 4a_*

mv  STAR_HTseq_out_exon/*to_Tce_LRv5a_mtDNAv350*txt  STAR_HTseq_out_exon/to_Tce_LRv5a_mtDNAv350
mv  STAR_HTseq_out_exon/*to_Tpa_LRv5a_mtDNAv350*txt   STAR_HTseq_out_exon/to_Tpa_LRv5a_mtDNAv350
mv  STAR_HTseq_out_exon/*to_Tps_LRv5b_mtDNAv350*txt   STAR_HTseq_out_exon/to_Tps_LRv5b_mtDNAv350

### Join

mkdir data/read_counts

python3 Accessory_scripts/HTSeq_to_edgeR.py -i  /Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/STAR_HTseq_out_exon/to_Tbi_LRv4a_mtDNAv350/     -o data/read_counts/to_Tbi -e .txt
python3 Accessory_scripts/HTSeq_to_edgeR.py -i  /Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/STAR_HTseq_out_exon/to_Tce_LRv5a_mtDNAv350/     -o data/read_counts/to_Tce -e .txt
python3 Accessory_scripts/HTSeq_to_edgeR.py -i  /Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/STAR_HTseq_out_exon/to_Tcm_LRv5a_mtDNAv350/     -o data/read_counts/to_Tcm -e .txt
python3 Accessory_scripts/HTSeq_to_edgeR.py -i  /Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/STAR_HTseq_out_exon/to_Tpa_LRv5a_mtDNAv350/     -o data/read_counts/to_Tpa -e .txt
python3 Accessory_scripts/HTSeq_to_edgeR.py -i  /Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/STAR_HTseq_out_exon/to_Tps_LRv5b_mtDNAv350/     -o data/read_counts/to_Tps -e .txt
python3 Accessory_scripts/HTSeq_to_edgeR.py -i  /Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/STAR_HTseq_out_exon/to_Tms_LRv5a_mtDNAv350/     -o data/read_counts/to_Tms -e .txt
python3 Accessory_scripts/HTSeq_to_edgeR.py -i  /Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/STAR_HTseq_out_exon/to_Tsi_LRv5b_mtDNAv350/     -o data/read_counts/to_Tsi -e .txt
python3 Accessory_scripts/HTSeq_to_edgeR.py -i  /Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/STAR_HTseq_out_exon/to_Tge_A_LRv5a_mtDNAv350/   -o data/read_counts/to_Tge_A -e .txt
python3 Accessory_scripts/HTSeq_to_edgeR.py -i  /Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/STAR_HTseq_out_exon/to_Tge_H20_LRv5a_mtDNAv350/ -o data/read_counts/to_Tge_H20 -e .txt
python3 Accessory_scripts/HTSeq_to_edgeR.py -i  /Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/STAR_HTseq_out_exon/to_Tdi_LRv5a_mtDNAv350/     -o data/read_counts/to_Tdi -e .txt


### add total exon len
### GET GFFS WITH EXONS FOR All features 

cp /Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/gffs/*gz .
gzip -d *gz

for gff in ./*_mtDNAv350_v2.1.gff; do
echo $gff
gff_1=`echo $gff   | sed 's/.gff/_wExon.gff/'` 
gff_2=`echo $gff_1 | sed 's/.gff/_wExon.gff/'`
gff_3=`echo $gff_2 | sed 's/.gff/_wExon.gff/'`
new_gff=`echo $gff | sed 's/.gff/_allwExon.gff/'`
fasta=`echo $gff | sed 's/_v2.*/.fasta/'`
out_trans=`echo $new_gff | sed 's/.gff/_trans.fa/'`
python3 Accessory_scripts/gff_add_exon_to_feature.py -i $gff    -f rRNA
python3 Accessory_scripts/gff_add_exon_to_feature.py -i $gff_1  -f ncRNA
python3 Accessory_scripts/gff_add_exon_to_feature.py -i $gff_2  -f tRNA
cp $gff_3 $new_gff
done

### add info to gene exp files

for f in ./*_mtDNAv350_v2.1_allwExon.gff ; do
	sp=`echo $f | sed 's/.*\///' | sed 's/_LRv.*//'`
	echo $sp
	python3 Accessory_scripts/gff_feature_lengths_and_scaf.py -i $f -o $sp -g data/read_counts/to_$sp"_H2E.counts.csv" 
done


### add gene position

for f in ./*_mtDNAv350_v2.1_allwExon.gff ; do
sp=`echo $f | sed 's/.*\///' | sed 's/_LRv.*//'`
echo $sp
python3 Accessory_scripts/add_gene_position.py -i $f -r "data/read_counts/to_"$sp"_H2E.counts_wGLS.csv"
done

rm *gff
rm *_lengths.csv

python3 Accessory_scripts/add_HiC_LG_to_Tce_counts.py -i data/read_counts/to_Tce_H2E.counts_wGLSL.csv 

### saved
## data/read_counts/

### also for 6_EdgeR_dosage.R: join_MF_for_heatmap.py


python3 Accessory_scripts/join_MF_for_heatmap.py   -s output/sex_ref_combSA_v2/Tps_Ad_to_Tps_F_MF.csv       -a output/sex_ref_combSA_v2/Tdi_Ad_to_Tps_F_MF.csv       -o output/sex_ref_combSA_v2/TpsTdi_Ad_to_Tps_F_MF
python3 Accessory_scripts/join_FPKM_for_heatmap.py -s output/sex_ref_combSA_v2/Tps_Ad_to_Tps_F_avFPKM_2.csv -a output/sex_ref_combSA_v2/Tdi_Ad_to_Tps_F_avFPKM_2.csv -o output/sex_ref_combSA_v2/TpsTdi_Ad_to_Tps_F_avFPKM_2

python3 Accessory_scripts/join_MF_for_heatmap.py   -s output/sex_ref_combSA_v2/Tce_Ad_to_Tce_F_MF.csv       -a output/sex_ref_combSA_v2/Tms_Ad_to_Tce_F_MF.csv       -o output/sex_ref_combSA_v2/TceTms_Ad_to_Tce_F_MF
python3 Accessory_scripts/join_FPKM_for_heatmap.py -s output/sex_ref_combSA_v2/Tce_Ad_to_Tce_F_avFPKM_2.csv -a output/sex_ref_combSA_v2/Tms_Ad_to_Tce_F_avFPKM_2.csv -o output/sex_ref_combSA_v2/TceTms_Ad_to_Tce_F_avFPKM_2

python3 Accessory_scripts/join_MF_for_heatmap.py   -s output/sex_ref_combSA_v2/Tpa_Ad_to_Tpa_F_MF.csv       -a output/sex_ref_combSA_v2/Tge_Ad_to_Tpa_F_MF.csv       -o output/sex_ref_combSA_v2/TpaTge_Ad_to_Tpa_F_MF
python3 Accessory_scripts/join_FPKM_for_heatmap.py -s output/sex_ref_combSA_v2/Tpa_Ad_to_Tpa_F_avFPKM_2.csv -a output/sex_ref_combSA_v2/Tge_Ad_to_Tpa_F_avFPKM_2.csv -o output/sex_ref_combSA_v2/TpaTge_Ad_to_Tpa_F_avFPKM_2






