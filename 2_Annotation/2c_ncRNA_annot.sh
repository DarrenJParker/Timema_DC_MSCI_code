#

DIR='/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot'
SDIR='/scratch/dparker/annot'


##############################################################################
### Infernal
##  https://docs.rfam.org/en/latest/genome-annotation.html

## set up in default
wget eddylab.org/infernal/infernal-1.1.2.tar.gz
tar xf infernal-1.1.2.tar.gz
cd infernal-1.1.2
./configure 
make
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
./src/cmpress Rfam.cm



cd $DIR
mkdir infernal_runs
cd $DIR/infernal_runs

### workout size
for f in  $DIR/genomes/*.fasta; do
echo $f
/work/FAC/FBM/DEE/tschwand/default/dparker/infernal-1.1.2/easel/miniapps/esl-seqstat $f
done


# Because we want millions of nucleotides on both strands, we multiply 1339820256 by 2, and divide by 1,000,000 = 2679.64051
# needed for accurate eval


##/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/genomes/Tce_LRv5a.fasta	1211418897	2422.837794
##/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/annot/genomes/Tpa_LRv5a.fasta 	1145684399	2291.368798

/work/FAC/FBM/DEE/tschwand/default/dparker/infernal-1.1.2/src/cmscan -Z 2422.837794 --cut_ga --rfam --nohmmonly --tblout Tce_LRv5a-genome.tblout --fmt 2 --cpu 40 \
--clanin /work/FAC/FBM/DEE/tschwand/default/dparker/infernal-1.1.2/Rfam.clanin /work/FAC/FBM/DEE/tschwand/default/dparker/infernal-1.1.2/Rfam.cm $DIR/genomes/Tce_LRv5a.fasta > Tce_LRv5a-genome.cmscan

/work/FAC/FBM/DEE/tschwand/default/dparker/infernal-1.1.2/src/cmscan -Z 2291.368798 --cut_ga --rfam --nohmmonly --tblout Tpa_LRv5a-genome.tblout --fmt 2 --cpu 40 \
--clanin /work/FAC/FBM/DEE/tschwand/default/dparker/infernal-1.1.2/Rfam.clanin /work/FAC/FBM/DEE/tschwand/default/dparker/infernal-1.1.2/Rfam.cm $DIR/genomes/Tpa_LRv5a.fasta > Tpa_LRv5a-genome.cmscan



### filter overlapping output

for f in infernal_runs/T*-genome.tblout; do
echo $f
out_f=`echo $f | sed 's/-genome.tblout/-genome.deoverlapped.tblout/'`
echo $out_f
grep -v " = " $f > $out_f
done


#ncRNA derived pseudogenes pose the biggest problem for eukaryotic genome annotation using Rfam/Infernal. Many genomes contain repeat elements that are derived from a non-coding RNA gene, sometimes in huge copy number. For example, Alu repeats in human are evolutionarily related to SRP RNA, and the active B2 SINE in mouse is recently derived from a tRNA.
#In addition, specific RNA genes appear to have undergone massive pseudogene expansions in certain genomes. For example, searching the human genome using the Rfam U6 family yields over 1000 hits, all with very high score. These are not “false positives” in the sequence analysis sense, because they are closely related by sequence to the real U6 genes, but they completely overwhelm the small number (only 10s) of expected real U6 genes.
#At present we don’t have computational methods to distinguish the real genes from the pseudogenes (of course the standard protein coding gene tricks - in frame stop codons and the like - are useless).
# The sensible and precedented method for ncRNA annotation in large vertebrate genomes is to annotate the easy-to-identify RNAs, such as tRNAs and rRNAs,
# and then trust only hits with very high sequence identity (>95% over >95% of the sequence length) to an experimentally verified real gene. tRNAscan-SE has a very nice method for detecting tRNA pseudogenes.
#

### so don't trust the ncRNA (other than rRNA and tRNA)?

### need to filter output and make gff.
### filter eval.
### store full output.

for f in infernal_runs/*-genome.deoverlapped.tblout; do
prefix=`echo $f | sed 's/-genome.deoverlapped.tblout.*//' | sed 's/.*\///' | sed 's/_LR.*//'`
echo $f
echo $prefix
python3 ~/Gen_BioInf/infernal_to_gff.py -i $f -g annotations_v2/$prefix"_braker_prot_run"*_combined_2_FBGgiUGF1000M0_wB.gff ### default evalue (1e-10)
done


#Tce: Number of infernal records kept: 473
#Tpa: Number of infernal records kept: 484



### store 
mkdir Timema_LR_genomic_code/output/infernal
cp infernal_runs/* Timema_LR_genomic_code/output/infernal
cd Timema_LR_genomic_code/output/
tar -czvf infernal.tar.gz infernal/
rm -r infernal





