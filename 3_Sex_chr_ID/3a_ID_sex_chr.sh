#### ID sex chr

DIR='/work/FAC/FBM/DEE/tschwand/timema_lr_genomes/dparker/sex_chr'
SDIR='/scratch/dparker/sex_chr'
mkdir -p $DIR
mkdir -p $SDIR


#########################################################################################################################
### get references

## use the ones with mtDNA
## genomes: DDBJ/ENA/GenBank under the accessions JBOIUC000000000 (T. cristinae), JBOIUD000000000 (Tpa), JBJCJD000000000 (Tps)

mkdir -p $DIR/REFS
cd $DIR/REFS


### map sexual species
#########################################################################################################################
#### get reads

### reads stored /nas/FAC/FBM/DEE/tschwand/sex_chromosomes/D2c/raw_reads/

### sex  species males   (Tce, Tpa, Tps) - SRA under the following project PRJNA725673
### sex  species females (Tce, Tpa, Tps) - SRA under the following project PRJNA670663
### asex species males (Tge, Tms, Tdi)   - SRA under the following project PRJNA1293928, PRJNA808673
### asex species females (Tge, Tms, Tdi) - SRA under the following project PRJNA670663


cd /nas/FAC/FBM/DEE/tschwand/sex_chromosomes/D2c/raw_reads/Tge/

########## clean and trim

mkdir -p $SDIR/temp_reads_for_cc
cd $SDIR

# rename 
for file in $DIR/READS/raw_reads/*.gz; do

	echo $file
	i=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/_L._R._.*//' `
	lane_N=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/_R._.*//' | sed 's/.*_//'`
	R_N=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/.*_L._R*//' | sed 's/_.*//'`	
	file_number=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/.*_L._R._//' | sed 's/_.*//' `
	Group_N=`echo $file | sed 's/.fastq.gz.*//' | sed 's/.*\///' | sed 's/.*_L._R._..._//' | sed 's/_//g'  `
	echo $i
	echo $lane_N	
	echo $R_N
	echo $file_number
	echo $Group_N

	sp=""

	if [ $i = "13_Tbi" ]; then
		   sp="Tbi_M"
	elif [ $i = "14_Tbi" ]; then
		   sp="Tbi_M"
	elif [ $i = "15_Tbi" ]; then
		   sp="Tbi_M"
	elif [ $i = "16_Tbi" ]; then
		   sp="Tbi_M"
	
	elif [ $i = "CC86B" ]; then
			sp="Tbi_F"
	elif [ $i = "CC86C" ]; then
		   sp="Tbi_F"
	elif [ $i = "CC87B" ]; then
		   sp="Tbi_F"
	elif [ $i = "CC87C" ]; then
		   sp="Tbi_F"
	elif [ $i = "CC88B" ]; then
		   sp="Tbi_F"


	elif [ $i = "05_HM15" ]; then
		   sp="Tce_M"
	elif [ $i = "06_HM16" ]; then
		   sp="Tce_M"
	elif [ $i = "07_HM33" ]; then
		   sp="Tce_M"
	elif [ $i = "08_HM61" ]; then
		   sp="Tce_M"

	elif [ $i = "CC22B" ]; then
		   sp="Tce_F"
	elif [ $i = "CC22C" ]; then
		   sp="Tce_F"
	elif [ $i = "CC24B" ]; then
		   sp="Tce_F"
	elif [ $i = "CC24C" ]; then
		   sp="Tce_F"
	elif [ $i = "CC25B" ]; then
		   sp="Tce_F"

	elif [ $i = "01_HM148" ]; then
			sp="Tcm_M"
	elif [ $i = "02_HM149" ]; then
		   sp="Tcm_M"
	elif [ $i = "03_HM150" ]; then
		   sp="Tcm_M"
	elif [ $i = "04_HM151" ]; then
		   sp="Tcm_M"

	elif [ $i = "HM217" ]; then
		   sp="Tcm_F"
	elif [ $i = "HM218" ]; then
		   sp="Tcm_F"
	elif [ $i = "HM219" ]; then
		   sp="Tcm_F"
	elif [ $i = "HM220" ]; then
		   sp="Tcm_F"
	elif [ $i = "HM221" ]; then
		   sp="Tcm_F"
										
	elif [ $i = "25_15055" ]; then
		   sp="Tms_M"
	elif [ $i = "26_15056" ]; then
		   sp="Tms_M"
	elif [ $i = "27_15057" ]; then
		   sp="Tms_M"
	elif [ $i = "34_16002" ]; then
		   sp="Tms_M"

	elif [ $i = "ReSeq_Ms01" ]; then
		   sp="Tms_F"
	elif [ $i = "ReSeq_Ms02" ]; then
		   sp="Tms_F"
	elif [ $i = "ReSeq_Ms03" ]; then
		   sp="Tms_F"
	elif [ $i = "MS_Alp03b" ]; then
		   sp="Tms_F"
	elif [ $i = "MS_Alp04b" ]; then
		   sp="Tms_F"

	elif [ $i = "09_Tpa" ]; then
		   sp="Tpa_M"
	elif [ $i = "10_Tpa" ]; then
		   sp="Tpa_M"		   
	elif [ $i = "11_Tpa" ]; then
		   sp="Tpa_M"		   
	elif [ $i = "12_Tpa" ]; then
		   sp="Tpa_M"		   

	elif [ $i = "Pa_AB" ]; then
		   sp="Tpa_F"
	elif [ $i = "PA_CD" ]; then
		   sp="Tpa_F"
	elif [ $i = "PA_E" ]; then
		   sp="Tpa_F"
	elif [ $i = "H54" ]; then
		   sp="Tpa_F"
	elif [ $i = "H56" ]; then
		   sp="Tpa_F"


	elif [ $i = "17_HM99" ]; then
		   sp="Tps_M"	
	elif [ $i = "18_HM100" ]; then
		   sp="Tps_M"
	elif [ $i = "19_HM101" ]; then
		   sp="Tps_M"
	elif [ $i = "20_15255" ]; then
		   sp="Tps_M"	

	elif [ $i = "ReSeq_Ps14" ]; then
		   sp="Tps_F"
	elif [ $i = "ReSeq_Ps16" ]; then
		   sp="Tps_F"
	elif [ $i = "ReSeq_Ps18" ]; then
		   sp="Tps_F"
	elif [ $i = "ReSeq_Ps08" ]; then
		   sp="Tps_F"
	elif [ $i = "ReSeq_Ps12" ]; then
		   sp="Tps_F"


	elif [ $i = "10091" ]; then
		   sp="Tge_M"
	elif [ $i = "19-1002" ]; then
		   sp="Tge_F"
	elif [ $i = "19-1004" ]; then
		   sp="Tge_M"

	elif [ $i = "CC59_A" ]; then
		   sp="Tge_F"
	elif [ $i = "CC59_C" ]; then
		   sp="Tge_F"
	elif [ $i = "CC65_B" ]; then
		   sp="Tge_F"
	elif [ $i = "CC66_A" ]; then
		   sp="Tge_F"
	elif [ $i = "CC67_A" ]; then
		   sp="Tge_F"

	elif [ $i = "18-3997" ]; then
		   sp="Tdi_M"
	elif [ $i = "18-3998" ]; then
		   sp="Tdi_M"

	elif [ $i = "ReSeq_Di02" ]; then
		   sp="Tdi_F"
	elif [ $i = "ReSeq_Di04" ]; then
		   sp="Tdi_F"
	elif [ $i = "ReSeq_Di06" ]; then
		   sp="Tdi_F"
	elif [ $i = "ReSeq_Di08" ]; then
		   sp="Tdi_F"
	elif [ $i = "ReSeq_Di10" ]; then
		   sp="Tdi_F"
		   

	elif [ $i = "ReSeq_S14" ]; then
		   sp="Tsi_F"
	elif [ $i = "ReSeq_Si01" ]; then
		   sp="Tsi_F"
	elif [ $i = "ReSeq_Si03" ]; then
		   sp="Tsi_F"
	elif [ $i = "ReSeq_Si16" ]; then
		   sp="Tsi_F"
	elif [ $i = "ReSeq_Si18" ]; then
		   sp="Tsi_F"

	elif [ $i = "ReSeq_Te07" ]; then
		   sp="Tte_F"
	elif [ $i = "ReSeq_Te08" ]; then
		   sp="Tte_F"
	elif [ $i = "ReSeq_Te09" ]; then
		   sp="Tte_F"
	elif [ $i = "ReSeq_Te10" ]; then
		   sp="Tte_F"
	elif [ $i = "ReSeq_Te11" ]; then
		   sp="Tte_F"
		   
	else
			foo5="nope"
	fi

	new_name=`echo $sp"_"$i"_R"$R_N"_G"$Group_N$lane_N"_"$file_number".fq.gz"` 
	echo $new_name
	echo ""
	cp $file $SDIR"/temp_reads_for_cc/"$new_name
done


#### cat and clean by read group

cd $SDIR
mkdir $SDIR/cat_clean_by_read_group

for i in temp_reads_for_cc/*_R1_* ; do
	foo_R1=`echo $i`
	foo_R2=`echo $i | sed 's/_R1_/_R2_/g'`
	base_out_name=`echo $i | sed 's/.*\///' | sed 's/R1.*//' `
	read_group=`echo $i | sed 's/.*\///' | sed 's/.*R1_//' | sed 's/_.*//'`	
	
	out_R1=`echo "cat_clean_by_read_group/"$base_out_name$read_group"_R1_CC.fq"  `
	out_R2=`echo "cat_clean_by_read_group/"$base_out_name$read_group"_R2_CC.fq"  `
	echo $foo_R1
	echo $foo_R2
	echo $base_out_name
	echo $read_group
	echo $out_R1
	echo $out_R2
	echo ""
	zcat $foo_R1  | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' >> $out_R1
	zcat $foo_R2  | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' >> $out_R2
done


### zip

gzip $SDIR/cat_clean_by_read_group/*fq

rm -r $DIR/READS/raw_reads
rm -r $SDIR/temp_reads_for_cc 


#####################################################################################################
#####################################################################################################
### read trimming


mkdir  $DIR/READS/trimmed_reads_by_RG
cd     $DIR/READS/trimmed_reads_by_RG

cp Accessory_scripts/AllIllumina-PEadapters.fa .

module load gcc
module load trimmomatic/0.39

for i in $SDIR/cat_clean_by_read_group/*_R1_CC.fq.gz ; do
        foo1=`echo $i`
		basename=`echo $foo1 | sed 's/_R1_CC.fq.gz.*//' | sed 's/.*\///'`
        infileR1=`echo $foo1`
        infileR2=`echo $foo1 | sed 's/_R1_CC.fq.gz/_R2_CC.fq.gz/'`
        outfileR1=`echo "./"$basename"_R1_qtrimmed.fq.gz"`
        outfileR2=`echo "./"$basename"_R2_qtrimmed.fq.gz"`
        outfileR1_UP=`echo "./"$basename"_R1_qtrimmed_UNPAIRED.fq.gz"`
        outfileR2_UP=`echo "./"$basename"_R2_qtrimmed_UNPAIRED.fq.gz"`
        
        echo $infileR1
        echo $infileR2
		echo $basename
        echo $outfileR1
        echo $outfileR1_UP
        echo $outfileR2
        echo $outfileR2_UP
        trimmomatic PE -threads 30 $infileR1 $infileR2 $outfileR1 $outfileR1_UP $outfileR2 $outfileR2_UP ILLUMINACLIP:AllIllumina-PEadapters.fa:3:25:6 LEADING:9 TRAILING:9 SLIDINGWINDOW:4:15 MINLEN:90
done


## tidy up

rm $DIR/READS/trimmed_reads_by_RG/*UNPAIRED.fq.gz
rm $SDIR/cat_clean_by_read_group/*fq.gz


#####################################################################################################
#####################################################################################################
### map reads

### prep refs

module load gcc
module load bwa/0.7.17 
module load samtools/1.15.1 
 
cd   $DIR/REFS

for f in ./*fasta; do
echo $f
bwa index $f
done


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### map reads as paired reads with BWA # gave 40GB RAM
### removes multi + supp reads

### Tce

mkdir $SDIR/BWA_out
mkdir $SDIR/BWA_out/mapped_as_paired
mkdir $SDIR/BWA_out/flagstat_out_paired

module load gcc
module load bwa/0.7.17
module load samtools/1.15

read_dir="$DIR/READS/trimmed_reads_by_RG"
ref_dir="$DIR/REFS"
map_out_dir="$SDIR/BWA_out/mapped_as_paired"
flag_out_dir="$SDIR/BWA_out/flagstat_out_paired"
mapqfilt="30"
ref_v="Tce_LRv5a_mtDNAv350"

for i in $read_dir/Tce*_R1_qtrimmed.fq.gz; do
	read_file_name_R1=`echo $i`
	read_file_name_R2=`echo $i | sed 's/_R1_/_R2_/' `		
    base_read_name=`echo $i | sed 's/.fq.gz.*//' |  sed 's/_R1_qtrimmed//' | sed 's/.*\///'`
    base_read_name3=`echo $i | sed 's/_qtrimmed.*//' | sed 's/.*\///' | sed 's/_G.*//'`
	badnames=`echo $base_read_name"_badnames.txt"`
    infile=`echo $i`
    ref_fa=`echo $ref_dir"/"$ref_v".fasta"`
    outsam=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA.sam"`
	outbam=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a.bam"`
	outbam_sorted=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a_sorted.bam"`
	flagstat_out_sam=`echo $flag_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_flagstat_out.txt"`
	flagstat_out_bam=`echo $flag_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a_flagstat_out.txt"`
	IFS='_' read -r -a sp_want_list <<< "$base_read_name"
	readgroup=`echo ${sp_want_list[-1]}`
	readgroup_txt=`echo "@RG\tID:"$readgroup"\tSM:"$base_read_name3"\tPL:ILLUMINA\tLB:"$base_read_name3`
		
		
	echo $read_file_name_R1
	echo $read_file_name_R2
    echo $base_read_name
	echo $base_read_name3
    echo $ref_fa
    echo $outsam
	echo $outbam
	echo $outbam_sorted
	echo $flagstat_out_sam
	echo $flagstat_out_bam
	echo $readgroup
	echo $readgroup_txt
	echo $badnames
	echo ""

	## map
    bwa mem -t 16 -R $readgroup_txt $ref_fa $read_file_name_R1 $read_file_name_R2 > $outsam
        
    #flagstat
	samtools flagstat $outsam > $flagstat_out_sam

    # filter ## filter both reads out to avoid broken flags
	samtools view -S  $outsam | fgrep -e 'XA:Z:' -e 'SA:Z:' | cut -f 1 > $badnames
	samtools view -Sh $outsam | fgrep -vf $badnames | samtools view -Shbq $mapqfilt - > $outbam
		
	# sort bam
	samtools sort $outbam > $outbam_sorted
		
	#flagstat
	samtools flagstat $outbam_sorted > $flagstat_out_bam
		
	#tidyup
	rm $outsam
	rm $outbam
        
done

rm *_badnames.txt


### Tpa

mkdir $SDIR/BWA_out
mkdir $SDIR/BWA_out/mapped_as_paired
mkdir $SDIR/BWA_out/flagstat_out_paired

module load gcc
module load bwa/0.7.17
module load samtools/1.15

read_dir="$DIR/READS/trimmed_reads_by_RG"
ref_dir="$DIR/REFS"
map_out_dir="$SDIR/BWA_out/mapped_as_paired"
flag_out_dir="$SDIR/BWA_out/flagstat_out_paired"
mapqfilt="30"
ref_v="Tpa_LRv5a_mtDNAv350"

for i in $read_dir/Tpa*_R1_qtrimmed.fq.gz; do
	read_file_name_R1=`echo $i`
	read_file_name_R2=`echo $i | sed 's/_R1_/_R2_/' `		
    base_read_name=`echo $i | sed 's/.fq.gz.*//' |  sed 's/_R1_qtrimmed//' | sed 's/.*\///'`
    base_read_name3=`echo $i | sed 's/_qtrimmed.*//' | sed 's/.*\///' | sed 's/_G.*//'`
	badnames=`echo $base_read_name"_badnames.txt"`
    infile=`echo $i`
    ref_fa=`echo $ref_dir"/"$ref_v".fasta"`
    outsam=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA.sam"`
	outbam=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a.bam"`
	outbam_sorted=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a_sorted.bam"`
	flagstat_out_sam=`echo $flag_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_flagstat_out.txt"`
	flagstat_out_bam=`echo $flag_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a_flagstat_out.txt"`
	IFS='_' read -r -a sp_want_list <<< "$base_read_name"
	readgroup=`echo ${sp_want_list[-1]}`
	readgroup_txt=`echo "@RG\tID:"$readgroup"\tSM:"$base_read_name3"\tPL:ILLUMINA\tLB:"$base_read_name3`
		
		
	echo $read_file_name_R1
	echo $read_file_name_R2
    echo $base_read_name
	echo $base_read_name3
    echo $ref_fa
    echo $outsam
	echo $outbam
	echo $outbam_sorted
	echo $flagstat_out_sam
	echo $flagstat_out_bam
	echo $readgroup
	echo $readgroup_txt
	echo $badnames
	echo ""

	## map
    bwa mem -t 16 -R $readgroup_txt $ref_fa $read_file_name_R1 $read_file_name_R2 > $outsam
        
    #flagstat
	samtools flagstat $outsam > $flagstat_out_sam

    # filter ## filter both reads out to avoid broken flags
	samtools view -S  $outsam | fgrep -e 'XA:Z:' -e 'SA:Z:' | cut -f 1 > $badnames
	samtools view -Sh $outsam | fgrep -vf $badnames | samtools view -Shbq $mapqfilt - > $outbam
		
	# sort bam
	samtools sort $outbam > $outbam_sorted
		
	#flagstat
	samtools flagstat $outbam_sorted > $flagstat_out_bam
		
	#tidyup
	rm $outsam
	rm $outbam
        
done

rm *_badnames.txt



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### map reads as paired reads with BWA # gave 40GB RAM
### removes multi + supp reads

### Tps

mkdir $SDIR/BWA_out
mkdir $SDIR/BWA_out/mapped_as_paired
mkdir $SDIR/BWA_out/flagstat_out_paired

module load gcc
module load bwa/0.7.17
module load samtools/1.15

read_dir="$DIR/READS/trimmed_reads_by_RG"
ref_dir="$DIR/REFS"
map_out_dir="$SDIR/BWA_out/mapped_as_paired"
flag_out_dir="$SDIR/BWA_out/flagstat_out_paired"
mapqfilt="30"
ref_v="Tps_LRv5b_mtDNAv350"

for i in $read_dir/Tps*_R1_qtrimmed.fq.gz; do
	read_file_name_R1=`echo $i`
	read_file_name_R2=`echo $i | sed 's/_R1_/_R2_/' `		
    base_read_name=`echo $i | sed 's/.fq.gz.*//' |  sed 's/_R1_qtrimmed//' | sed 's/.*\///'`
    base_read_name3=`echo $i | sed 's/_qtrimmed.*//' | sed 's/.*\///' | sed 's/_G.*//'`
	badnames=`echo $base_read_name"_badnames.txt"`
    infile=`echo $i`
    ref_fa=`echo $ref_dir"/"$ref_v".fasta"`
    outsam=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA.sam"`
	outbam=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a.bam"`
	outbam_sorted=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a_sorted.bam"`
	flagstat_out_sam=`echo $flag_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_flagstat_out.txt"`
	flagstat_out_bam=`echo $flag_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a_flagstat_out.txt"`
	IFS='_' read -r -a sp_want_list <<< "$base_read_name"
	readgroup=`echo ${sp_want_list[-1]}`
	readgroup_txt=`echo "@RG\tID:"$readgroup"\tSM:"$base_read_name3"\tPL:ILLUMINA\tLB:"$base_read_name3`
		
		
	echo $read_file_name_R1
	echo $read_file_name_R2
    echo $base_read_name
	echo $base_read_name3
    echo $ref_fa
    echo $outsam
	echo $outbam
	echo $outbam_sorted
	echo $flagstat_out_sam
	echo $flagstat_out_bam
	echo $readgroup
	echo $readgroup_txt
	echo $badnames
	echo ""

	## map
    bwa mem -t 16 -R $readgroup_txt $ref_fa $read_file_name_R1 $read_file_name_R2 > $outsam
        
    #flagstat
	samtools flagstat $outsam > $flagstat_out_sam

    # filter ## filter both reads out to avoid broken flags
	samtools view -S  $outsam | fgrep -e 'XA:Z:' -e 'SA:Z:' | cut -f 1 > $badnames
	samtools view -Sh $outsam | fgrep -vf $badnames | samtools view -Shbq $mapqfilt - > $outbam
		
	# sort bam
	samtools sort $outbam > $outbam_sorted
		
	#flagstat
	samtools flagstat $outbam_sorted > $flagstat_out_bam
		
	#tidyup
	rm $outsam
	rm $outbam
        
done

rm *_badnames.txt





#########################################################################################################################
##### merge bams by samp


mkdir $SDIR/BWA_out/mapped_as_paired_merged

Tce_samples=(
Tce_F_CC22B
Tce_F_CC22C
Tce_F_CC24B
Tce_F_CC24C
Tce_F_CC25B
Tce_M_05_HM15
Tce_M_06_HM16
Tce_M_07_HM33
Tce_M_08_HM61
)

Tpa_samples=(
Tpa_F_H54
Tpa_F_H56
Tpa_F_Pa_AB
Tpa_F_PA_CD
Tpa_F_PA_E
Tpa_M_09_Tpa
Tpa_M_10_Tpa
Tpa_M_11_Tpa
Tpa_M_12_Tpa
)

Tps_samples=(
Tps_F_ReSeq_Ps08
Tps_F_ReSeq_Ps12
Tps_F_ReSeq_Ps14
Tps_F_ReSeq_Ps16
Tps_F_ReSeq_Ps18
Tps_M_17_HM99
Tps_M_18_HM100
Tps_M_19_HM101
Tps_M_20_15255
)


genome_pref="Tce_LRv5a_mtDNAv350"
for s in ${Tce_samples[@]}; do
echo $s
samtools merge -@ 40 $SDIR"/BWA_out/mapped_as_paired_merged/"$s"_to_"$genome_pref"_pe_BWA_mapqfilt_30a_sorted.bam" $SDIR"/BWA_out/mapped_as_paired/"$s*"_to_"$genome_pref*".bam"
done

genome_pref="Tpa_LRv5a_mtDNAv350"
for s in ${Tpa_samples[@]}; do
echo $s
samtools merge -@ 40 $SDIR"/BWA_out/mapped_as_paired_merged/"$s"_to_"$genome_pref"_pe_BWA_mapqfilt_30a_sorted.bam" $SDIR"/BWA_out/mapped_as_paired/"$s*"_to_"$genome_pref*".bam"
done

genome_pref="Tps_LRv5b_mtDNAv350"
for s in ${Tps_samples[@]}; do
echo $s
samtools merge -@ 40 $SDIR"/BWA_out/mapped_as_paired_merged/"$s"_to_"$genome_pref"_pe_BWA_mapqfilt_30a_sorted.bam" $SDIR"/BWA_out/mapped_as_paired/"$s*"_to_"$genome_pref*".bam"
done

#########################################################################################################################
##### remove PCR duplicates  (Not just mark, as I don't think angsD pays attention to this) (actually it should, but I also need to calc cov.)
###  60GB RAM


module load gcc
module load bwa/0.7.17
module load samtools/1.15.1
module load picard/2.26.2

for i in  $SDIR/BWA_out/mapped_as_paired_merged/*_mapqfilt_30a_sorted.bam; do
    outbam=`echo $i | sed 's/_mapqfilt_30a_sorted.bam/_mapqfilt_30aDR_sorted.bam/'`
        flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`
        metric_file=`echo $outbam | sed 's/.bam/_metric.txt/'`

        echo $i
        echo $outbam
        echo $metric_file
        echo $flagstat_out_bam

        picard MarkDuplicates REMOVE_DUPLICATES=true \
        INPUT=$i \
    OUTPUT=$outbam \
    METRICS_FILE=$metric_file

        samtools flagstat $outbam > $flagstat_out_bam

done




###############################################################################################################################
### First indel realignment ##

module load gcc
module load bwa/0.7.17
module load samtools/1.15.1
module load picard/2.26.2
module load gatk/3.8.1 


## index and dict of fasta

for f in $DIR/REFS/*_LR*_mtDNAv350.fasta  ; do 
    dict_name=`echo $f | sed 's/.fasta/.dict/'`
    samtools faidx $f
	picard CreateSequenceDictionary R= $f O= $dict_name
done


# indel realignment Tce

ref_v="Tce_LRv5a_mtDNAv350"
for i in  $SDIR/BWA_out/mapped_as_paired_merged/*_to_Tce_*_mapqfilt_30aDR_sorted.bam; do
    outbam=`echo $i | sed 's/_mapqfilt_30aDR_sorted.bam/_mapqfilt_30aDRra_sorted.bam/'`
	flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`
	interval_file=`echo $outbam | sed 's/.bam/.intervals/'`
	ref_fa=`echo $DIR"/REFS/"$ref_v".fasta"`
	
	echo $i
	echo $outbam
	echo $flagstat_out_bam
	echo $interval_file
	echo $ref_fa
	echo ""

	# index bam
	samtools index $i
	
	## make target intervals list
	gatk -T RealignerTargetCreator -R $ref_fa -I $i -o $interval_file
	
	## realign
	gatk -T IndelRealigner -R $ref_fa --targetIntervals $interval_file -I $i -o $outbam
	samtools flagstat $outbam > $flagstat_out_bam	

done


# indel realignment Tpa

ref_v="Tpa_LRv5a_mtDNAv350"
for i in  $SDIR/BWA_out/mapped_as_paired_merged/*_to_Tpa_*_mapqfilt_30aDR_sorted.bam; do
    outbam=`echo $i | sed 's/_mapqfilt_30aDR_sorted.bam/_mapqfilt_30aDRra_sorted.bam/'`
	flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`
	interval_file=`echo $outbam | sed 's/.bam/.intervals/'`
	ref_fa=`echo $DIR"/REFS/"$ref_v".fasta"`
	
	echo $i
	echo $outbam
	echo $flagstat_out_bam
	echo $interval_file
	echo $ref_fa
	echo ""

	# index bam
	samtools index $i
	
	## make target intervals list
	gatk -T RealignerTargetCreator -R $ref_fa -I $i -o $interval_file
	
	## realign
	gatk -T IndelRealigner -R $ref_fa --targetIntervals $interval_file -I $i -o $outbam
	samtools flagstat $outbam > $flagstat_out_bam	

done

# indel realignment Tps

ref_v="Tps_LRv5b_mtDNAv350"
for i in  $SDIR/BWA_out/mapped_as_paired_merged/*_to_Tps_*_mapqfilt_30aDR_sorted.bam; do
    outbam=`echo $i | sed 's/_mapqfilt_30aDR_sorted.bam/_mapqfilt_30aDRra_sorted.bam/'`
	flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`
	interval_file=`echo $outbam | sed 's/.bam/.intervals/'`
	ref_fa=`echo $DIR"/REFS/"$ref_v".fasta"`
	
	echo $i
	echo $outbam
	echo $flagstat_out_bam
	echo $interval_file
	echo $ref_fa
	echo ""

	# index bam
	samtools index $i
	
	## make target intervals list
	gatk -T RealignerTargetCreator -R $ref_fa -I $i -o $interval_file
	
	## realign
	gatk -T IndelRealigner -R $ref_fa --targetIntervals $interval_file -I $i -o $outbam
	samtools flagstat $outbam > $flagstat_out_bam	

done



### STORE out of scratch

mkdir -p $DIR/BWA_out/mapped_as_paired_merged/

cp $SDIR/BWA_out/mapped_as_paired_merged/*_pe_BWA_mapqfilt_30aDRra_sorted.bam $DIR/BWA_out/mapped_as_paired_merged/

### also take the flagstats

cp -r $SDIR/BWA_out/flagstat_out_paired $DIR/BWA_out/
cp $SDIR/BWA_out/mapped_as_paired_merged/*_flagstat_out.txt $DIR/BWA_out/flagstat_out_paired



### store bams /nas/FAC/FBM/DEE/tschwand/timema_lr_genomes/D2c/timema_genomes/BAMs

nohup cp *_30aDRra_sorted.bam /nas/FAC/FBM/DEE/tschwand/timema_lr_genomes/D2c/timema_genomes/BAMs &


##################################################################################################################################################################################3
##### coverage - genome hist 


module load gcc
module load bedtools2 # bedtools v2.30.0

for s in $DIR/BWA_out/mapped_as_paired_merged_recip/*_30aDRra_sorted.bam; do
outfile=`echo $s | sed 's/_sorted.bam/_cov.txt/'`
echo $s
echo $outfile
genomeCoverageBed -ibam $s > $outfile
done

#### plot coverage

cd $DIR

for i in $DIR/BWA_out/mapped_as_paired_merged/*_cov.txt; do
	echo $i
	python3 Accessory_scripts/genomeCoverageBed_tidier_wholegenomecov.py -i $i
done

module load gcc r/4.2.1 

for i in $DIR/BWA_out/mapped_as_paired_merged/*_genomecov.txt; do
echo $i
Rscript Accessory_scripts/plot_genome_cov.R $i
done

## store plots and cov ests

mkdir ../annot/Timema_LR_genomic_code/output/coverage
cp BWA_out/mapped_as_paired_merged/*.txtcovest* ../annot/Timema_LR_genomic_code/output/coverage


##################################################################################################################################################################################3
##### coverage - per base


module load gcc
module load bedtools2 # bedtools v2.30.0

for s in $DIR/BWA_out/mapped_as_paired_merged/*_30aDRra_sorted.bam; do
outfile=`echo $s | sed 's/_sorted.bam/_covBEDGRAPHbga.txt/'`
echo $s
echo $outfile
genomeCoverageBed -ibam $s -bga > $outfile
done



### per window
## ~30 mins to run per samp

cd /work/FAC/FBM/DEE/tschwand/asex_sinergia/dparker/mapping_to_LR_genomes/BWA_out/mapped_as_paired_merged

#### 10000 window
for s in $DIR/BWA_out/mapped_as_paired_merged/*_30aDRra_covBEDGRAPHbga.txt; do
outfile=`echo $s | sed 's/_covBEDGRAPHbga.txt//'`
echo $s
echo $outfile
python3 Accessory_scripts/bedgraph_cov_to_windows.py -i $s -w 10000 -o $outfile
done

#### 100000 window
for s in $DIR/BWA_out/mapped_as_paired_merged/*_30aDRra_covBEDGRAPHbga.txt; do
outfile=`echo $s | sed 's/_covBEDGRAPHbga.txt//'`
echo $s
echo $outfile
python3 Accessory_scripts/bedgraph_cov_to_windows.py -i $s -w 100000 -o $outfile
done


###########################################################################################################################################
####  map asex to sex ref

### Tms

mkdir $SDIR/BWA_out
mkdir $SDIR/BWA_out/mapped_as_paired
mkdir $SDIR/BWA_out/flagstat_out_paired

module load gcc
module load bwa/0.7.17
module load samtools/1.15

read_dir="$DIR/READS/trimmed_reads_by_RG"
ref_dir="$DIR/REFS"
map_out_dir="$SDIR/BWA_out/mapped_as_paired"
flag_out_dir="$SDIR/BWA_out/flagstat_out_paired"
mapqfilt="30"
ref_v="Tce_LRv5a_mtDNAv350"

for i in $read_dir/Tms*_R1_qtrimmed.fq.gz; do
	read_file_name_R1=`echo $i`
	read_file_name_R2=`echo $i | sed 's/_R1_/_R2_/' `		
    base_read_name=`echo $i | sed 's/.fq.gz.*//' |  sed 's/_R1_qtrimmed//' | sed 's/.*\///'`
    base_read_name3=`echo $i | sed 's/_qtrimmed.*//' | sed 's/.*\///' | sed 's/_G.*//'`
	badnames=`echo $base_read_name"_badnames.txt"`
    infile=`echo $i`
    ref_fa=`echo $ref_dir"/"$ref_v".fasta"`
    outsam=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA.sam"`
	outbam=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a.bam"`
	outbam_sorted=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a_sorted.bam"`
	flagstat_out_sam=`echo $flag_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_flagstat_out.txt"`
	flagstat_out_bam=`echo $flag_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a_flagstat_out.txt"`
	IFS='_' read -r -a sp_want_list <<< "$base_read_name"
	readgroup=`echo ${sp_want_list[-1]}`
	readgroup_txt=`echo "@RG\tID:"$readgroup"\tSM:"$base_read_name3"\tPL:ILLUMINA\tLB:"$base_read_name3`
		
		
	echo $read_file_name_R1
	echo $read_file_name_R2
    echo $base_read_name
	echo $base_read_name3
    echo $ref_fa
    echo $outsam
	echo $outbam
	echo $outbam_sorted
	echo $flagstat_out_sam
	echo $flagstat_out_bam
	echo $readgroup
	echo $readgroup_txt
	echo $badnames
	echo ""

	## map
    bwa mem -t 16 -R $readgroup_txt $ref_fa $read_file_name_R1 $read_file_name_R2 > $outsam
        
    #flagstat
	samtools flagstat $outsam > $flagstat_out_sam

    # filter ## filter both reads out to avoid broken flags
	samtools view -S  $outsam | fgrep -e 'XA:Z:' -e 'SA:Z:' | cut -f 1 > $badnames
	samtools view -Sh $outsam | fgrep -vf $badnames | samtools view -Shbq $mapqfilt - > $outbam
		
	# sort bam
	samtools sort $outbam > $outbam_sorted
		
	#flagstat
	samtools flagstat $outbam_sorted > $flagstat_out_bam
		
	#tidyup
	rm $outsam
	rm $outbam
        
done

rm *_badnames.txt


### Tge

mkdir $SDIR/BWA_out
mkdir $SDIR/BWA_out/mapped_as_paired
mkdir $SDIR/BWA_out/flagstat_out_paired

module load gcc
module load bwa/0.7.17
module load samtools/1.15

read_dir="$DIR/READS/trimmed_reads_by_RG"
ref_dir="$DIR/REFS"
map_out_dir="$SDIR/BWA_out/mapped_as_paired"
flag_out_dir="$SDIR/BWA_out/flagstat_out_paired"
mapqfilt="30"
ref_v="Tpa_LRv5a_mtDNAv350"

for i in $read_dir/Tge*_R1_qtrimmed.fq.gz; do
	read_file_name_R1=`echo $i`
	read_file_name_R2=`echo $i | sed 's/_R1_/_R2_/' `		
    base_read_name=`echo $i | sed 's/.fq.gz.*//' |  sed 's/_R1_qtrimmed//' | sed 's/.*\///'`
    base_read_name3=`echo $i | sed 's/_qtrimmed.*//' | sed 's/.*\///' | sed 's/_G.*//'`
	badnames=`echo $base_read_name"_badnames.txt"`
    infile=`echo $i`
    ref_fa=`echo $ref_dir"/"$ref_v".fasta"`
    outsam=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA.sam"`
	outbam=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a.bam"`
	outbam_sorted=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a_sorted.bam"`
	flagstat_out_sam=`echo $flag_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_flagstat_out.txt"`
	flagstat_out_bam=`echo $flag_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a_flagstat_out.txt"`
	IFS='_' read -r -a sp_want_list <<< "$base_read_name"
	readgroup=`echo ${sp_want_list[-1]}`
	readgroup_txt=`echo "@RG\tID:"$readgroup"\tSM:"$base_read_name3"\tPL:ILLUMINA\tLB:"$base_read_name3`
		
		
	echo $read_file_name_R1
	echo $read_file_name_R2
    echo $base_read_name
	echo $base_read_name3
    echo $ref_fa
    echo $outsam
	echo $outbam
	echo $outbam_sorted
	echo $flagstat_out_sam
	echo $flagstat_out_bam
	echo $readgroup
	echo $readgroup_txt
	echo $badnames
	echo ""

	## map
    bwa mem -t 16 -R $readgroup_txt $ref_fa $read_file_name_R1 $read_file_name_R2 > $outsam
        
    #flagstat
	samtools flagstat $outsam > $flagstat_out_sam

    # filter ## filter both reads out to avoid broken flags
	samtools view -S  $outsam | fgrep -e 'XA:Z:' -e 'SA:Z:' | cut -f 1 > $badnames
	samtools view -Sh $outsam | fgrep -vf $badnames | samtools view -Shbq $mapqfilt - > $outbam
		
	# sort bam
	samtools sort $outbam > $outbam_sorted
		
	#flagstat
	samtools flagstat $outbam_sorted > $flagstat_out_bam
		
	#tidyup
	rm $outsam
	rm $outbam
        
done

rm *_badnames.txt


### Tdi

mkdir $SDIR/BWA_out
mkdir $SDIR/BWA_out/mapped_as_paired
mkdir $SDIR/BWA_out/flagstat_out_paired

module load gcc
module load bwa/0.7.17
module load samtools/1.15

read_dir="$DIR/READS/trimmed_reads_by_RG"
ref_dir="$DIR/REFS"
map_out_dir="$SDIR/BWA_out/mapped_as_paired"
flag_out_dir="$SDIR/BWA_out/flagstat_out_paired"
mapqfilt="30"
ref_v="Tps_LRv5b_mtDNAv350"

for i in $read_dir/Tdi*_R1_qtrimmed.fq.gz; do
	read_file_name_R1=`echo $i`
	read_file_name_R2=`echo $i | sed 's/_R1_/_R2_/' `		
    base_read_name=`echo $i | sed 's/.fq.gz.*//' |  sed 's/_R1_qtrimmed//' | sed 's/.*\///'`
    base_read_name3=`echo $i | sed 's/_qtrimmed.*//' | sed 's/.*\///' | sed 's/_G.*//'`
	badnames=`echo $base_read_name"_badnames.txt"`
    infile=`echo $i`
    ref_fa=`echo $ref_dir"/"$ref_v".fasta"`
    outsam=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA.sam"`
	outbam=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a.bam"`
	outbam_sorted=`echo $map_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a_sorted.bam"`
	flagstat_out_sam=`echo $flag_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_flagstat_out.txt"`
	flagstat_out_bam=`echo $flag_out_dir"/"$base_read_name"_to_"$ref_v"_pe_BWA_mapqfilt_"$mapqfilt"a_flagstat_out.txt"`
	IFS='_' read -r -a sp_want_list <<< "$base_read_name"
	readgroup=`echo ${sp_want_list[-1]}`
	readgroup_txt=`echo "@RG\tID:"$readgroup"\tSM:"$base_read_name3"\tPL:ILLUMINA\tLB:"$base_read_name3`
		
		
	echo $read_file_name_R1
	echo $read_file_name_R2
    echo $base_read_name
	echo $base_read_name3
    echo $ref_fa
    echo $outsam
	echo $outbam
	echo $outbam_sorted
	echo $flagstat_out_sam
	echo $flagstat_out_bam
	echo $readgroup
	echo $readgroup_txt
	echo $badnames
	echo ""

	## map
    bwa mem -t 16 -R $readgroup_txt $ref_fa $read_file_name_R1 $read_file_name_R2 > $outsam
        
    #flagstat
	samtools flagstat $outsam > $flagstat_out_sam

    # filter ## filter both reads out to avoid broken flags
	samtools view -S  $outsam | fgrep -e 'XA:Z:' -e 'SA:Z:' | cut -f 1 > $badnames
	samtools view -Sh $outsam | fgrep -vf $badnames | samtools view -Shbq $mapqfilt - > $outbam
		
	# sort bam
	samtools sort $outbam > $outbam_sorted
		
	#flagstat
	samtools flagstat $outbam_sorted > $flagstat_out_bam
		
	#tidyup
	rm $outsam
	rm $outbam
        
done

rm *_badnames.txt





#### merge bams

mkdir $SDIR/BWA_out/mapped_as_paired_merged

Tdi_samples=(
Tdi_F_ReSeq_Di02
Tdi_F_ReSeq_Di04
Tdi_F_ReSeq_Di06
Tdi_F_ReSeq_Di08
Tdi_F_ReSeq_Di10
Tdi_M_18-3997
Tdi_M_18-3998
)

Tge_samples=(
Tge_F_CC59_A
Tge_F_CC59_C
Tge_F_CC65_B
Tge_F_CC66_A
Tge_F_CC67_A
Tge_F_19-1002
Tge_M_10091
Tge_M_19-1004
)

Tms_samples=(
Tms_F_MS_Alp03b
Tms_F_MS_Alp04b
Tms_F_ReSeq_Ms01
Tms_F_ReSeq_Ms02
Tms_F_ReSeq_Ms03
Tms_M_25_15055
Tms_M_26_15056
Tms_M_27_15057
Tms_M_34_16002
)

genome_pref="Tce_LRv5a_mtDNAv350"
for s in ${Tms_samples[@]}; do
echo $s
samtools merge -@ 40 $SDIR"/BWA_out/mapped_as_paired_merged/"$s"_to_"$genome_pref"_pe_BWA_mapqfilt_30a_sorted.bam" $SDIR"/BWA_out/mapped_as_paired/"$s*"_to_"$genome_pref*".bam"
done

genome_pref="Tpa_LRv5a_mtDNAv350"
for s in ${Tge_samples[@]}; do
echo $s
samtools merge -@ 40 $SDIR"/BWA_out/mapped_as_paired_merged/"$s"_to_"$genome_pref"_pe_BWA_mapqfilt_30a_sorted.bam" $SDIR"/BWA_out/mapped_as_paired/"$s*"_to_"$genome_pref*".bam"
done

genome_pref="Tps_LRv5b_mtDNAv350"
for s in ${Tdi_samples[@]}; do
echo $s
samtools merge -@ 40 $SDIR"/BWA_out/mapped_as_paired_merged/"$s"_to_"$genome_pref"_pe_BWA_mapqfilt_30a_sorted.bam" $SDIR"/BWA_out/mapped_as_paired/"$s*"_to_"$genome_pref*".bam"
done







#########################################################################################################################
##### remove PCR duplicates  (Not just mark, as I don't think angsD pays attention to this) (actually it should, but I also need to calc cov.)
###  60GB RAM


module load gcc
module load bwa/0.7.17
module load samtools/1.15.1
module load picard/2.26.2

for i in  $SDIR/BWA_out/mapped_as_paired_merged/*_mapqfilt_30a_sorted.bam; do
    outbam=`echo $i | sed 's/_mapqfilt_30a_sorted.bam/_mapqfilt_30aDR_sorted.bam/'`
        flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`
        metric_file=`echo $outbam | sed 's/.bam/_metric.txt/'`

        echo $i
        echo $outbam
        echo $metric_file
        echo $flagstat_out_bam

        picard MarkDuplicates REMOVE_DUPLICATES=true \
        INPUT=$i \
    OUTPUT=$outbam \
    METRICS_FILE=$metric_file

        samtools flagstat $outbam > $flagstat_out_bam

done




###############################################################################################################################
### First indel realignment ##

module load gcc
module load bwa/0.7.17
module load samtools/1.15.1
module load picard/2.26.2
module load gatk/3.8.1 


## index and dict of fasta

for f in $DIR/REFS/*_LR*_mtDNAv350.fasta  ; do 
    dict_name=`echo $f | sed 's/.fasta/.dict/'`
    samtools faidx $f
	picard CreateSequenceDictionary R= $f O= $dict_name
done


# indel realignment 

ref_v="Tbi_LRv4a_mtDNAv350"
for i in  $SDIR/BWA_out/mapped_as_paired_merged/*_to_Tbi_*_mapqfilt_30aDR_sorted.bam; do
    outbam=`echo $i | sed 's/_mapqfilt_30aDR_sorted.bam/_mapqfilt_30aDRra_sorted.bam/'`
	flagstat_out_bam=`echo $outbam | sed 's/.bam/_flagstat_out.txt/'`
	interval_file=`echo $outbam | sed 's/.bam/.intervals/'`
	ref_fa=`echo $DIR"/REFS/"$ref_v".fasta"`
	
	echo $i
	echo $outbam
	echo $flagstat_out_bam
	echo $interval_file
	echo $ref_fa
	echo ""

	# index bam
	samtools index $i
	
	## make target intervals list
	gatk -T RealignerTargetCreator -R $ref_fa -I $i -o $interval_file
	
	## realign
	gatk -T IndelRealigner -R $ref_fa --targetIntervals $interval_file -I $i -o $outbam
	samtools flagstat $outbam > $flagstat_out_bam	

done


#### changed the reads and ref with sed

sed 's/_to_Tbi_/_to_Tce_/g'     indelR_Tbi_.sh | sed 's/Tbi_LRv4a_mtDNAv350/Tce_LRv5a_mtDNAv350/g'     > indelR_Tce_.sh
sed 's/_to_Tbi_/_to_Tpa_/g'     indelR_Tbi_.sh | sed 's/Tbi_LRv4a_mtDNAv350/Tpa_LRv5a_mtDNAv350/g'     > indelR_Tpa_.sh
sed 's/_to_Tbi_/_to_Tps_/g'     indelR_Tbi_.sh | sed 's/Tbi_LRv4a_mtDNAv350/Tps_LRv5b_mtDNAv350/g'     > indelR_Tps_.sh


### STORE out of scratch

mkdir -p $DIR/BWA_out/mapped_as_paired_merged_recip

cp $SDIR/BWA_out/mapped_as_paired_merged/*_pe_BWA_mapqfilt_30aDRra_sorted.bam $DIR/BWA_out/mapped_as_paired_merged_recip/

### also take the flagstats

cp $SDIR/BWA_out/flagstat_out_paired/*_flagstat_out.txt $DIR/BWA_out/flagstat_out_paired
cp $SDIR/BWA_out/mapped_as_paired_merged/*_flagstat_out.txt $DIR/BWA_out/flagstat_out_paired



### store bams /nas/FAC/FBM/DEE/tschwand/timema_lr_genomes/D2c/timema_genomes/BAMs

cd $DIR/BWA_out/mapped_as_paired_merged_recip
nohup cp *_30aDRra_sorted.bam /nas/FAC/FBM/DEE/tschwand/timema_lr_genomes/D2c/timema_genomes/BAMs &



##################################################################################################################################################################################3
##### coverage - genome hist 


module load gcc
module load bedtools2 # bedtools v2.30.0

for s in $DIR/BWA_out/mapped_as_paired_merged_recip/*_30aDRra_sorted.bam; do
outfile=`echo $s | sed 's/_sorted.bam/_cov.txt/'`
echo $s
echo $outfile
genomeCoverageBed -ibam $s > $outfile
done




#### plot coverage

cd $DIR

for i in $DIR/BWA_out/mapped_as_paired_merged_recip/*_cov.txt; do
echo $i
python3 Accessory_scripts/genomeCoverageBed_tidier_wholegenomecov.py -i $i
done

module load gcc r/4.2.1 

for i in $DIR/BWA_out/mapped_as_paired_merged_recip/*_genomecov.txt; do
echo $i
Rscript Accessory_scripts/plot_genome_cov.R $i
done


## store plots and cov ests

mkdir ../annot/Timema_LR_genomic_code/output/coverage
cp BWA_out/mapped_as_paired_merged_recip/*.txtcovest* ../annot/Timema_LR_genomic_code/output/coverage


##################################################################################################################################################################################3
##### coverage - per base


module load gcc
module load bedtools2 # bedtools v2.30.0

for s in $DIR/BWA_out/mapped_as_paired_merged_recip/*_30aDRra_sorted.bam; do
outfile=`echo $s | sed 's/_sorted.bam/_covBEDGRAPHbga.txt/'`
echo $s
echo $outfile
genomeCoverageBed -ibam $s -bga > $outfile
done



### per window
## ~30 mins to run per samp

#### 10000 window
for s in $DIR/BWA_out/mapped_as_paired_merged_recip/*_30aDRra_covBEDGRAPHbga.txt; do
outfile=`echo $s | sed 's/_covBEDGRAPHbga.txt//'`
echo $s
echo $outfile
python3  Accessory_scripts/bedgraph_cov_to_windows.py -i $s -w 10000 -o $outfile
done

#### 100000 window
for s in $DIR/BWA_out/mapped_as_paired_merged_recip/*_30aDRra_covBEDGRAPHbga.txt; do
outfile=`echo $s | sed 's/_covBEDGRAPHbga.txt//'`
echo $s
echo $outfile
python3  Accessory_scripts/bedgraph_cov_to_windows.py -i $s -w 100000 -o $outfile
done


### tidy up to_Tce (add the HiC --> LG info) for plotting
### info here:
### renamed by size
#OLD	Size	new
#LG3	159.323	HiCLG1
#LG2	152.505	HiCLG2
#LG8	136.351	HiCLG3
#LG13	85.394	HiCLG4
#LG7	77.854	HiCLG5
#LG1	76.422	HiCLG6
#LG5	76.182	HiCLG7
#LG9	76.098	HiCLG8
#LG10	73.404	HiCLG9
#LG6	67.818	HiCLG10
#LG4	65.285	HiCLG11
#LG11	64.343	HiCLG12
#LG12	41.071	HiCLG13


for f in /Users/drp22jhz/Documents/University/Lausanne/Timema_LR_genomes/Timema_LR_genomic_code/output/angsD_LR_sliding_window/*_to_Tce_LRv5a_mtDNAv350_pe_BWA_mapqfilt_30aDRra_cov_window_100000.txt; do
echo $f
python3 Accessory_scripts/add_HiC_LG_to_Tce.py -i $f
done


## plot
## cov_plot_LR.R

