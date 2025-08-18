## 2a_prep_RNAseq_for_annot.sh

##########################################################################################################################
### get reads

## TceTms: 2019 WB single end stuff, WB dev stuff, fasteris adults, Old antennae 
## TpaTge: 2019 WB single end stuff, WB dev stuff, fasteris adults

mkdir raw_reads_paired
mkdir raw_reads_single

###################################################################################################
#### Trim paired reads


mkdir trimmed_reads_paired
cp Accessory_scripts/Trimmomatic_adapters/AllIllumina-PEadapters.fa .

module load gcc
module load trimmomatic/0.39
cd   trimmed_reads_paired
for i in raw_reads_paired/*_R1*.fq.gz ; do
foo1=`echo $i`
basename=`echo $foo1 | sed 's/_R1.*//' | sed 's/.*\///'`
infileR1=`echo $foo1`
infileR2=`echo $foo1 | sed 's/_R1_CC.fq.gz/_R2_CC.fq.gz/'`
outfileR1=`echo "./"$basename"_R1_qtrimmed.fq.gz"`
outfileR2=`echo "./"$basename"_R2_qtrimmed.fq.gz"`
outfileR1_UP=`echo "./"$basename"_R1_qtrimmed_UNPAIRED.fq.gz"`
outfileR2_UP=`echo "./"$basename"_R2_qtrimmed_UNPAIRED.fq.gz"`

echo $infileR1
echo $infileR2
echo $outfileR1
echo $outfileR1_UP
echo $outfileR2
echo $outfileR2_UP

trimmomatic PE -threads 30 $infileR1 $infileR2 $outfileR1 $outfileR1_UP $outfileR2 $outfileR2_UP ILLUMINACLIP:AllIllumina-PEadapters.fa:3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:80 CROP:100

done

rm *_UNPAIRED.fq.gz


###################################################################################################
#### Trim single reads

mkdir trimmed_reads_single
cp Accessory_scripts/Trimmomatic_adapters/AllIllumina-PEadapters.fa .

module load gcc
module load trimmomatic/0.39
cd   trimmed_reads_single
for i in raw_reads_single/*_R1*.fq.gz ; do
foo1=`echo $i`
basename=`echo $foo1 | sed 's/_R1_.*//' | sed 's/.*\///'`
infileR1=`echo $foo1`
outfileR1=`echo "./"$basename"_R1_qtrimmed.fq.gz"`

echo $infileR1
echo $outfileR1

trimmomatic SE -threads 30 $infileR1 $infileR2 $outfileR1 $outfileR1_UP $outfileR2 $outfileR2_UP ILLUMINACLIP:AllIllumina-PEadapters.fa:3:25:6 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:80 CROP:100

done



