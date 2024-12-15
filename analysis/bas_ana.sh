#!/bin/bash
#
#SBATCH --job-name=bas_ana
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mem=25GB
#SBATCH --time=20:00:00  
#SBATCH --output=./logs/ba_%j.out
#SBATCH --error=./logs/ba_%j.err
###SBATCH --array=1-70

ID=$SLURM_ARRAY_TASK_ID

samp=$(sed -n ${ID}p ./Great_Ape/genome_pipeline/files/captive.txt)

method="roh"
# famcor = correct fam files for f-stats
# fstatc = precompute f-statistics for a given individual
# fstatbest = best fit f-statistics all individuals
# hucon = human contamination test
# caga = rareCAGA preparation (liftover) & geolocalization
# roh = RoH calling

DIR=${TMPDIR}

if [[ $method == "hucon" ]] ; then DIR="Result/ref_genome/HG38-ucsc/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"; module load --auto r/4.2.2-gcc-12.2.0-ekejlri; fi

if [[ $method =~ ^("famcor"|"fstatc"|"fstatbest")$ ]] ; then module load R/4.2.3; fi

if [[ $method == "caga" ]] ; then
  module load R/4.2.3 bcftools conda; conda activate geo ; conda activate gdal;conda activate proj; LD_LIBRARY_PATH=$LD_LIBRARY_PATH
  export BCFTOOLS_PLUGINS=programs/bcftools-plugin && bcftools +liftover
  
  # liftover of genotypes from hg38 to hg19
  echo "A
  B" | cat - ./Great_Ape/genome_pipeline/files/captive_chimps.txt > greatapes/newlist.txt
  bcftools +$BCFTOOLS_PLUGINS/liftover.so -Ou greatapes/pan_all/chr21.vcf.gz -- -s ./GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -f mk_data/nea_contest/hg19.fa -c ./GRCh38/hg38ToHg19.over.chain.gz | bcftools sort | bcftools view -S greatapes/newlist.txt -a -A | bcftools view -m 2 -M 2 -Ov -o greatapes/analysis/chr21.lifted.vcf.gz -W
  
  cd rareCAGA
  Rscript fulltest.R "greatapes/analysis/chr21.lifted.vcf.gz captive"
  echo "done"
  exit
  fi

if [[ $method == "roh" ]] ; then
  module load bcftools/1.21
  for spec in gorilla pan pongo; do
    echo "rohs in" $spec
    for chr in {1..22};do
      bcftools roh -O rz --AF-tag AF -G30 -I --threads 8 -o greatapes/analysis/"${spec}"."${chr}".rohs.txt greatapes/"${spec}"_all/chr"${chr}".vcf.gz
      done
    done
  echo "done"
  exit
  fi

# the main script to run for some of these steps
Rscript ./Great_Ape/genome_pipeline/analysis/basic_analysis.R ${samp} ${DIR} ${method}

if [[ $method == "fstatc" ]]; then rm -r $TMPDIR/${samp}; fi

echo "done"

exit

