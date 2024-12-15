#!/bin/bash
#
#SBATCH --job-name=snpcalling
#SBATCH --cpus-per-task=1
#SBATCH -n 1
#SBATCH --output=./anc_%j.out
#SBATCH --error=./anc_%j.err
#SBATCH --mem=10GB
#SBATCH --time=1:30:00
#SBATCH --partition=skylake_0096
#SBATCH --qos=skylake_0096
#SBATCH --array=1-23

cd Result/ImportDB/{spec} 

ID=$SLURM_ARRAY_TASK_ID

ind=$(pwd | tr "\/" "\n" | tail -1 )
echo $ind
cd ..
dir=$(pwd)
cd ~/logs

module load --auto r/4.2.2-gcc-12.2.0-ekejlri
Rscript ~/Great_Ape/genome_pipeline/processing/create_ancestral.R $ID $ind $dir

echo "done"

exit
