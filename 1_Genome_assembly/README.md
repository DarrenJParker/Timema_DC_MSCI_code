## Pipeline description for genome assembly.

Genome assemblies have been deposited at DDBJ/ENA/GenBank under the accessions JBOIUC000000000 (T. cristinae (Tce)) and JBOIUD000000000 (T. podura (Tpa)).

## Table of contents

1. [Oxford Nanopore reads filtration](#1)
2. [Contig assembly](#2)
3. [Contamination removal](#3)
4. [Haplotype purging](#4)



#### <a name="1"></a>1) Oxford Nanopore reads filtration. 

Reads for Tpa available in the European Nucleotide Archive under the following project PRJEB94285.

Using [Filtlong](https://github.com/rrwick/Filtlong).

```
module load filtlong/0.2.0 

filtlong --min_length 1000 --keep_percent 90 --target_bases 69050000000 Tce_nanopore_raw_reads.fastq.gz | gzip > Tce_nanopore_filtlong50x.fastq.gz
filtlong --min_length 1000 --keep_percent 90 --target_bases 69050000000 Tpa_nanopore_raw_reads.fastq.gz | gzip > Tpa_nanopore_filtlong50x.fastq.gz
```

#### <a name="2"></a>2) Contig assembly

2.1) Nanopore reads assembly

Using [Flye](https://github.com/fenderglass/Flye).

```
source flye/2.8.1.sh

flye --nano-raw  Tce_nanopore_filtlong50x.fastq --out-dir Tce_flye --genome-size 1.381g
flye --nano-raw  Tpa_nanopore_filtlong50x.fastq --out-dir Tpa_flye --genome-size 1.381g
```

2.2) First step of polishing using filtered Oxford Nanopore reads and Racon. 

Using [Minimap2](https://github.com/lh3/minimap2) and [Racon](https://github.com/lbcb-sci/racon).

```
# Map the reads against the contig assembly generated previously

module load minimap2/2.19

minimap2 -c -x map-ont Tce_flye_contig.fasta Tce_nanopore_filtlong50x.fastq > Tce_flye_contig_filtlong50x.paf
minimap2 -c -x map-ont Tpa_flye_contig.fasta Tpa_nanopore_filtlong50x.fastq > Tpa_flye_contig_filtlong50x.paf

# Use Racon for polishing

source racon/1.4.3.sh

racon Tce_nanopore_filtlong50x.fastq Tce_flye_contig_filtlong50x.paf Tce_flye_contig.fasta > Tce_flye_racon.fasta
racon Tpa_nanopore_filtlong50x.fastq Tpa_flye_contig_filtlong50x.paf Tpa_flye_contig.fasta > Tpa_flye_racon.fasta
```

2.3) Second step of polishing using filtered Illumina short reads and three rounds of Pilon. 

Using [BWA](https://github.com/lh3/bwa), [Samtools](https://www.htslib.org/) and [Pilon](https://github.com/broadinstitute/pilon).

```
module load bwa/0.7.17
module load samtools/1.12

# Round 1

bwa mem Tce_flye_racon Tce_illumina_R1.fastq.gz Tce_illumina_R2.fastq.gz | samtools view -bS - | samtools sort - > Tce_flye_racon_pilon1.sorted.bam
java -jar pilon-1.23.jar --genome Tce_flye_racon --frags Tce_flye_racon_pilon1.sorted.bam --output Tce_flye_racon_pilon1 --diploid

bwa mem Tpa_flye_racon Tpa_illumina_R1.fastq.gz Tpa_illumina_R2.fastq.gz | samtools view -bS - | samtools sort - > Tpa_flye_racon_pilon1.sorted.bam
java -jar pilon-1.23.jar --genome Tpa_flye_racon --frags Tpa_flye_racon_pilon1.sorted.bam --output Tpa_flye_racon_pilon1 --diploid


# Round 2

bwa mem Tce_flye_racon_pilon1 Tce_illumina_R1.fastq.gz Tce_illumina_R2.fastq.gz | samtools view -bS - | samtools sort - > Tce_flye_racon_pilon2.sorted.bam
java -jar pilon-1.23.jar --genome Tce_flye_racon_pilon1 --frags Tce_flye_racon_pilon2.sorted.bam --output Tce_flye_racon_pilon2 --diploid

bwa mem Tpa_flye_racon_pilon1 Tpa_illumina_R1.fastq.gz Tpa_illumina_R2.fastq.gz | samtools view -bS - | samtools sort - > Tpa_flye_racon_pilon2.sorted.bam
java -jar pilon-1.23.jar --genome Tpa_flye_racon_pilon1 --frags Tpa_flye_racon_pilon2.sorted.bam --output Tpa_flye_racon_pilon2 --diploid


# Round 3

bwa mem Tce_flye_racon_pilon2 Tce_illumina_R1.fastq.gz Tce_illumina_R2.fastq.gz | samtools view -bS - | samtools sort - > Tce_flye_racon_pilon3.sorted.bam
java -jar pilon-1.23.jar --genome Tce_flye_racon_pilon2 --frags Tce_flye_racon_pilon3.sorted.bam --output Tce_flye_racon_pilon3 --diploid

bwa mem Tpa_flye_racon_pilon2 Tpa_illumina_R1.fastq.gz Tpa_illumina_R2.fastq.gz | samtools view -bS - | samtools sort - > Tpa_flye_racon_pilon3.sorted.bam
java -jar pilon-1.23.jar --genome Tpa_flye_racon_pilon2 --frags Tpa_flye_racon_pilon3.sorted.bam --output Tpa_flye_racon_pilon3 --diploid


```

2.4) Only for *Timema cristinae* (Tce)

```
source canu/2.0.sh

canu -p Tce -d Tce genomeSize=1381m correctedErrorRate=0.105 -pacbio-raw ../Tce_pacbio_raw_reads.fastq
```


#### <a name="3"></a>3) Contamination removal

Using [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [BlobTools](https://github.com/DRL/blobtools) and [BBMap](https://sourceforge.net/projects/bbmap/).

3.1) Get the coverage of each contigs from the Flye output in `contigs_stats.txt`.

3.2) Blastn against nt:

```
module add Blast/ncbi-blast/2.7.1+
 
blastn -query Tce_flye_racon_pilon3.fasta -db /nt \
-outfmt '6 qseqid staxids bitscore evalue std sscinames sskingdoms stitle' \
-max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -out Tce_flye_racon_pilon3.vs.nt.max10.1e25.blastn.out

blastn -query Tpa_flye_racon_pilon3.fasta -db /nt \
-outfmt '6 qseqid staxids bitscore evalue std sscinames sskingdoms stitle' \
-max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -out Tpa_flye_racon_pilon3.vs.nt.max10.1e25.blastn.out
```

3.3) Run BlobTools:

```
source blobtools/1.1.1.sh

blobtools create -i Tce_flye_racon_pilon3.fasta -t Tce_flye_racon_pilon3.vs.nt.max10.1e25.blastn.out \
--nodes data/nodes.dmp --names data/names.dmp \
-c coverage.txt -x bestsumorder -o Tce_flye_racon_pilon3
blobtools blobplot -i Tce_flye_racon_pilon3.blobDB.json --sort count --hist count -x bestsumorder
blobtools view -i Tce_flye_racon_pilon3.blobDB.json --hits --rank all -x bestsumorder

blobtools create -i Tpa_flye_racon_pilon3.fasta -t Tpa_flye_racon_pilon3.vs.nt.max10.1e25.blastn.out \
--nodes data/nodes.dmp --names data/names.dmp \
-c coverage.txt -x bestsumorder -o Tpa_flye_racon_pilon3
blobtools blobplot -i Tpa_flye_racon_pilon3.blobDB.json --sort count --hist count -x bestsumorder
blobtools view -i Tpa_flye_racon_pilon3.blobDB.json --hits --rank all -x bestsumorder


```

3.4) Filter out contigs without hits to metazoans

```
python contamination_filtration.py -s contamination_identification -i1 Tce_flye_racon_pilon3.blobDB.table.txt
python contamination_filtration.py -s contamination_identification -i1 Tpa_flye_racon_pilon3.blobDB.table.txt

module add UHTS/Analysis/BBMap/37.82

filterbyname.sh in=Tce_flye_racon_pilon3.fasta names=contaminant.txt out=Tce_flye_racon_pilon3_blobtools.fasta include=f 
filterbyname.sh in=Tpa_flye_racon_pilon3.fasta names=contaminant.txt out=Tpa_flye_racon_pilon3_blobtools.fasta include=f 
```

#### <a name="4"></a>4) Haplotype purging

Using [Minimap2](https://github.com/lh3/minimap2), [Samtools](https://www.htslib.org/) and [Purge Haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/).

4.1) Map the filtered Oxford Nanopore reads against the decontaminated genome

```
module load minimap2/2.19

minimap2 -ax map-ont Tce_flye_racon_pilon3_blobtools.fasta Tce_nanopore_filtlong50x.fastq --secondary=no | samtools sort -o Tce_map.sorted.bam
minimap2 -ax map-ont Tpa_flye_racon_pilon3_blobtools.fasta Tpa_nanopore_filtlong50x.fastq --secondary=no | samtools sort -o Tpa_map.sorted.bam
```

4.2) Run Purge Haplotigs

```
source purge_haplotigs/1.1.1.sh

purge_haplotigs cov -i Tce_map.sorted.bam.gencov -l 3 -m 25 -h 190 -j 101
purge_haplotigs purge -g Tce_flye_racon_pilon3_blobtools.fasta -c coverage_stats.csv -o Tce_flye_racon_pilon3_blobtools_purgehap

purge_haplotigs cov -i Tpa_map.sorted.bam.gencov -l 3 -m 45 -h 195 -j 101 
purge_haplotigs purge -g Tpa_flye_racon_pilon3_blobtools.fasta -c coverage_stats.csv -o Tpa_flye_racon_pilon3_blobtools_purgehap

```

#### <a name="5"></a>5) Hi-C scaffolding

Using [BWA](https://github.com/lh3/minimap2), [Juicer](https://www.htslib.org/) and [3D-DNA](https://github.com/aidenlab/3d-dna).

5.1) Map Hi-C reads

```
# Create a BWA reference

module add UHTS/Aligner/bwa/0.7.17;

bwa index Tce_flye_racon_pilon3_blobtools_purgehap.fasta 

# Generate a restriction sites file 

python juicer/generate_site_positions.py Sau3AI draft Tce_flye_racon_pilon3_blobtools_purgehap.fasta


# Create chrom.sizes file

awk 'BEGIN{OFS="\t"}{print $1, $NF}' draft_Sau3AI.txt > Tce_flye_racon_pilon3_blobtools_purgehap.chrom.sizes


# Run Juicer

source juicer/1.6.sh

./scripts/juicer.sh -g Tce -z Tce_flye_racon_pilon3_blobtools_purgehap.fasta -y draft_Sau3AI.txt -p Tce_flye_racon_pilon3_blobtools_purgehap.chrom.sizes -D Tce_juicer

# Create a BWA reference

module add UHTS/Aligner/bwa/0.7.17;

bwa index Tpa_flye_racon_pilon3_blobtools_purgehap.fasta

# Generate a restriction sites file 

python juicer/generate_site_positions.py Sau3AI draft Tpa_flye_racon_pilon3_blobtools_purgehap.fasta


# Create chrom.sizes file

awk 'BEGIN{OFS="\t"}{print $1, $NF}' draft_Sau3AI.txt > Tpa_flye_racon_pilon3_blobtools_purgehap.chrom.sizes

# Run Juicer

source juicer/1.6.sh

./scripts/juicer.sh -g Tpa -z Tpa_flye_racon_pilon3_blobtools_purgehap.fasta -y draft_Sau3AI.txt -p Tpa_flye_racon_pilon3_blobtools_purgehap.chrom.sizes -D Tpa_juicer

```


5.2) Perform scaffolding

```
source 3d-dna/v180922.sh

3d-dna/run-asm-pipeline.sh --editor-coarse-resolution 250000 Tce_flye_racon_pilon3_blobtools_purgehap.fasta merged_nodups.txt
3d-dna/run-asm-pipeline.sh --editor-coarse-resolution 50000  Tpa_flye_racon_pilon3_blobtools_purgehap.fasta merged_nodups.txt
```

