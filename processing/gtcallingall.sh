#!/bin/bash
#
#SBATCH --job-name=snpcalling
#SBATCH --cpus-per-task=4
#SBATCH -N 1
#SBATCH --ntasks=23
#SBATCH --mem=736GB
#SBATCH --time=72:00:00
#SBATCH --partition=skylake_0768 
#SBATCH --qos=skylake_0768 

ind=$(pwd | tr "\/" "\n" | tail -1 )
echo $ind
module load --auto samtools/1.14-gcc-12.2.0-mdhu6gw anaconda3/2022.05-gcc-12.2.0-oqiw76n
SPARK_JAVA_OPTS+=" -Dspark.local.dir=tmp,tmp2 -Dhadoop.tmp.dir=tmp2"
export SPARK_JAVA_OPTS
echo "start"
if [[ ! -e ./snpcalling ]]; then
mkdir snpcalling
fi

for CHR in {1..24};do
if [[ ${CHR} == 23 ]]; then CHR="X"; fi
if [[ ${CHR} == 24 ]]; then CHR="Y"; fi
Software/gatk-4.1.4.0/gatk --java-options "-Xmx32g" HaplotypeCaller \
-R Result/ref_genome/HG38-ucsc/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
-I ./mapping/merge/${ind}.merge.cram \
-O  ./snpcalling/${ind}_${CHR}.g.vcf.gz --native-pair-hmm-threads 4 -L chr${CHR} \
-ERC GVCF &
done

wait

echo "done"

