#!/bin/bash
#
#SBATCH --job-name=ImportDB1
#SBATCH --cpus-per-task=8
#SBATCH --partition=skylake_0384
#SBATCH --qos=skylake_0384
#SBATCH -N 1
#SBATCH --mem=220GB
#SBATCH --time=1-03:45:00
#SBATCH --output=log/%x_%a.out
#SBATCH --error=log/%x_%a.err
#SBATCH --array=1-23
if [[ ${SLURM_ARRAY_TASK_ID} == 23 ]]; then SLURM_ARRAY_TASK_ID="X"; fi

module load python

/usr/bin/time Software/gatk-4.1.4.0/gatk --java-options "-Xmx80g -Xms80g" \
       GenomicsDBImport --genomicsdb-workspace-path Result/ImportDB/{species}/my_database_chr${SLURM_ARRAY_TASK_ID} \
       -L chr${SLURM_ARRAY_TASK_ID} \
       --sample-name-map Result/ImportDB/{species}/${SLURM_ARRAY_TASK_ID}.list \
       --tmp-dir Result/ImportDB/{species}/log --reader-threads 8


