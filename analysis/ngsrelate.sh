#!/bin/bash
#
#SBATCH --job-name=ngsr
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=500GB
#SBATCH --time=72:00:00  
#SBATCH --output=./logs/ng_%j.out
#SBATCH --error=./logs/ng_%j.err
#SBATCH --array=2-3

ID=$SLURM_ARRAY_TASK_ID
basedir=greatapes/

module load --auto ngsrelate

cd $basedir

if [[ ${ID} == 1 ]]; then spec="gorilla"; fi
if [[ ${ID} == 2 ]]; then spec="pan"; fi
if [[ ${ID} == 3 ]]; then spec="pongo"; fi
ngsRelate -h "${spec}"_all/allchr.vcf.gz -T GT -O plink/ngsr/"${spec}"_all.res -c 1

echo "done"

exit

