#!/bin/bash
#
#SBATCH --job-name=geno_proc
#SBATCH --cpus-per-task=8
#SBATCH --mem=90GB
#SBATCH --time=65:00:00
#SBATCH --output=./pipeline.out
#SBATCH --error=./pipeline.err
#SBATCH --partition=skylake_0096
#SBATCH --qos=skylake_0096
#SBATCH --export=NONE

###############################
###############################
###### GENOME PROCESSING ######
## now a streamlined version ##
###############################
###############################

####################
## 0) pre-processing

## first steps to do: 
# create a directory following SR style
# mkdir Result/Sequences/$species_name/$individual_name
# go to that directory
# create a list with the accession IDs and the name SraAccList.txt, e.g.
# echo "ERRXXXX
# ERRXXXXX" > SraAccList.txt
# in this directory, you should also run this pipeline script


####################
#### RUNNING IT ####
####################

## define the individual here
ind=$(pwd | tr "\/" "\n" | tail -1 )

module load --auto miniconda3 samtools/1.14-gcc-12.2.0-mdhu6gw 
module unload openjdk
conda activate Software/miniconda3/envs/java8
SPARK_JAVA_OPTS+=" -Dspark.local.dir=tmp,tmp2 -Dhadoop.tmp.dir=tmp2"
SPARK_LOCAL_DIRS="tmp2,tmp"
export SPARK_JAVA_OPTS
export SPARK_LOCAL_DIRS


####################
## 1) SRA processing

if [[ `ls -1 *.fastq 2>/dev/null | wc -l ` -gt 0 ]]; then
gzip ./*.fastq
fi

if [[ `ls -1 *.fastq.gz 2>/dev/null | wc -l ` -eq 0 ]]; then
echo "SRA download"
#download SRA files
Software/sratoolkit.3.0.6-centos_linux64/bin/prefetch $(<SraAccList.txt) --max-size u -O ./
#SRA to fastq conversion
Software/sratoolkit.3.0.6-centos_linux64/bin/fasterq-dump  $(<SraAccList.txt) 
#validate SRA files
Software/sratoolkit.3.0.6-centos_linux64/bin/vdb-validate $(<SraAccList.txt) 
#count the reads
gzip ./*.fastq
ls -l ./*.fastq.gz | awk '{ print $9 }' | while read l; do echo -ne $l"\t"; awk 'END{print NR/4}' $l >> ./readcount.txt; done
fi

##########################
## 2) FASTQC quality check

# fastQC
if [[ ! -e ./fastqc-output ]]; then
mkdir fastqc-output
echo "FASTQC"
for infile in ./*.fastq.gz; do
/usr/bin/time Software/FastQC/fastqc -j Software/jdk-17.0.2/bin/java -t 8 -o ./fastqc-output $infile
done
fi

######################
## 3) adapter trimming

#Trimmomatic
if [[ ! -e ./filtered ]]; then mkdir filtered;echo "TRIMMING"; fi

while read NAME in; do
if [[ `ls -1 ./filtered/$NAME* 2>/dev/null | wc -l ` -eq 0 ]]; then
echo "trimming $NAME"

# paired end
if [[ -f ${NAME}_1.fastq.gz ]]; then
echo "paired end"
/usr/bin/time java -jar Software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 -phred33 ${NAME}_1.fastq.gz ${NAME}_2.fastq.gz \
./filtered/${NAME}_1p.fastq.gz ./filtered/${NAME}_1un.fastq.gz \
./filtered/${NAME}_2p.fastq.gz ./filtered/${NAME}_2un.fastq.gz \
ILLUMINACLIP:Software/Trimmomatic-0.39/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fi

# single end
if [[ -f ${NAME}.fastq.gz ]]; then
echo "single end"
/usr/bin/time java -jar Software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 8 -phred33 ${NAME}.fastq.gz \
./filtered/${NAME}_p.fastq.gz \
ILLUMINACLIP:Software/Trimmomatic-0.39/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fi

fi
done < SraAccList.txt

######################
## 4) mapping
echo "MAPPING"

if [[ ! -e ./mapping ]]; then
mkdir mapping
mkdir mapping/removeunmapped
mkdir mapping/bamreadgroup
mkdir mapping/markduplicationspark
fi

## mapping genome with BWA 
## directly pipe to bam & sorting & remove unmapped

while read NAME in; do
# check if output file exists
if [[ ! -f ./mapping/removeunmapped/${NAME}.bam ]] ; then
echo "mapping $NAME"

# check if paired-end fastq exists, and if so, do the mapping
if [[ -f ./filtered/${NAME}_1p.fastq.gz ]] ; then
if [[ ! -f ./mapping/removeunmapped/${NAME}_a.bam ]] ; then
echo "paired end"
/usr/bin/time Software/bwa-0.7.16a/bwa mem -M -t 8 Result/ref_genome/humanHG38  ./filtered/${NAME}_1p.fastq.gz ./filtered/${NAME}_2p.fastq.gz | samtools sort -@ 8 | samtools view -@ 8 -b -F 4 -o ./mapping/removeunmapped/${NAME}_a.bam
fi
fi

# check if single-end fastq exists, and if so, do the mapping
if [[ -f ./filtered/${NAME}_p.fastq.gz ]]; then
if [[ ! -f ./mapping/removeunmapped/${NAME}_b.bam ]]; then
echo "single end"
/usr/bin/time Software/bwa-0.7.16a/bwa mem -M -t 8 Result/ref_genome/humanHG38  ./filtered/${NAME}_p.fastq.gz | samtools sort -@ 8 | samtools view -@ 8 -b -F 4 -o ./mapping/removeunmapped/${NAME}_b.bam 
fi
fi

# check if both paired and single end bams exist, and if so, do the merging; if only paired end exists, define this as main bam
if [[ -f ./mapping/removeunmapped/${NAME}_a.bam ]]; then
echo "merge"
if [[ -f ./mapping/removeunmapped/${NAME}_b.bam ]]; then
/usr/bin/time samtools merge -@ 8 -o ./mapping/removeunmapped/${NAME}.bam ./mapping/removeunmapped/${NAME}_b.bam ./mapping/removeunmapped/${NAME}_a.bam; else mv ./mapping/removeunmapped/${NAME}_a.bam ./mapping/removeunmapped/${NAME}.bam
fi
fi

# check if only single end bam exists, and if so, define it as main bam
if [[ -f ./mapping/removeunmapped/${NAME}_b.bam ]]; then
mv ./mapping/removeunmapped/${NAME}_b.bam ./mapping/removeunmapped/${NAME}.bam
fi

# check everything went ok
if [[ ! -f ./mapping/removeunmapped/${NAME}.bam ]]; then echo "ERROR AT MAPPING $NAME";exit;fi
fi
done < SraAccList.txt


######################################
## 5) add Read Group & mark duplicates
echo "ADD READ GROUP"

# add read group
while read NAME in; do
if [[ ! -f ./mapping/bamreadgroup/${NAME}.bam ]]; then
echo "RG $NAME"
/usr/bin/time java -jar Software/picard-2.21.4/picard.jar AddOrReplaceReadGroups \
I=./mapping/removeunmapped/${NAME}.bam O=./mapping/bamreadgroup/${NAME}.bam \
RGLB=lib3 \
RGPL=ILLUMINA \
RGPU=unit3 \
RGSM=${ind}
fi

if [[ ! -f ./mapping/bamreadgroup/${NAME}.bam ]]; then echo "ERROR AT ADDRG $NAME";exit;fi

done < SraAccList.txt

## mark duplication

while read NAME in; do
if [[ ! -f ./mapping/markduplicationspark/${NAME}.bam ]]; then
echo "MARKDUPLICATES $NAME"
rm -r ./mapping/markduplicationspark/${NAME}.*
/usr/bin/time Software/gatk-4.1.4.0/gatk MarkDuplicatesSpark -I ./mapping/bamreadgroup/${NAME}.bam  \
-O ./mapping/markduplicationspark/${NAME}.bam \
-M ./mapping/markduplicationspark/marked_dup_metrics.$NAME.txt 

if [[ ! -f ./mapping/markduplicationspark/${NAME}.bam ]]; then echo "ERROR AT MARKDUP $NAME";exit;fi
fi

done < SraAccList.txt

######################################
## 6) merging

mkdir mapping/merge
echo "MERGING"

if [[ ! -f ./merge/$ind.merge.cram ]]; then
/usr/bin/time samtools merge -f -c -O CRAM --reference Result/ref_genome/HG38-ucsc/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --write-index -@ 8 ./mapping/markduplicationspark/*.bam -o ./mapping/merge/$ind.merge.cram
fi

echo "DONE"

exit


## RUN
#sbatch Great_Ape_genomes/processing/basic_processing.sh


