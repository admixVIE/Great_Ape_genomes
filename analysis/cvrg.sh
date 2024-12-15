#!/bin/bash
#
#SBATCH --job-name=sustats
#SBATCH --cpus-per-task=4
#SBATCH --mem 25G
#SBATCH --ntasks=1
#SBATCH --time=15:00:00
#SBATCH --output=logs/su_stat_%j.out
#SBATCH --error=logs/su_stat_%j.err
#SBATCH --partition=skylake_0096
#SBATCH --qos=skylake_0096
#SBATCH --array=197,173,30,42

# method should be "lastpos" if checking the last position for integrity of the file
# method should be "centcov" if calculating the cumulative coverage distribution

method="lastpos"

basedir="summary_stats/"
ID=$SLURM_ARRAY_TASK_ID

if [[ $method == "lastpos" ]] ; then

file=$(sed -n ${ID}p ${basedir}/pathlist.txt)
ind=$(echo $file | tr "\/" "\n" | tail -1 )

if [[ -f $file/snpcalling/lastpos.txt ]]; then rm $file/snpcalling/lastpos.txt; fi

for chr in {1..24};do
echo $chr
if [[ ${chr} == 23 ]]; then chr="X"; fi
if [[ ${chr} == 24 ]]; then chr="Y"; fi
zcat $file/snpcalling/"$ind"_"$chr".g.vcf.gz | tail -n 1 >> $file/snpcalling/lastpos.txt &
done  
wait

exit
fi


if [[ $method == "centcov" ]] ; then

module load --auto r/4.2.2-gcc-12.2.0-ekejlri bcftools

Rscript ./Great_Ape/genome_pipeline/analysis/central_coverage.R ${ID}

exit
fi

