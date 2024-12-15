#!/bin/bash
#
#SBATCH --job-name=g_stats
#SBATCH --cpus-per-task=1
#SBATCH -N 1
#SBATCH --time=6:00:00
#SBATCH --output=./snpcalling/g_stat_%j.out
#SBATCH --error=./snpcalling/g_stat_%j.err
#SBATCH --partition=skylake_0096
#SBATCH --qos=skylake_0096
#SBATCH --export=NONE

module load --auto samtools/1.14-gcc-12.2.0-mdhu6gw anaconda3/2022.05-gcc-12.2.0-oqiw76n

ind=$(pwd | tr "\/" "\n" | tail -1 )
echo $ind

if [[ ! -f ./merge/$ind.merge.cram.crai ]]; then
echo "creating index"
samtools index ./mapping/merge/${ind}.merge.cram
fi

## mosdepth
module load --auto miniconda3
echo "mosdepth"
mosdepth -x -n -f Result/ref_genome/HG38-ucsc/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta mapping/merge/mosdepth ./mapping/merge/${ind}.merge.cram

## samtools
echo "samtools stats"
samtools stats  ./mapping/merge/${ind}.merge.cram > ./mapping/merge/samtoolsstats_${ind}.txt 

## bcftools
echo "bcftools stats"
for CHR in  {1..24}; do
if [[ ${CHR} == 23 ]]; then CHR="X"; fi
if [[ ${CHR} == 24 ]]; then CHR="Y"; fi
echo $CHR
Software/bcftools/bcftools-1.16/bcftools stats  ./snpcalling/${ind}_${CHR}.g.vcf.gz > ./snpcalling/bcftoolsstats_${ind}_${CHR}.txt 
done

echo "done"

exit


