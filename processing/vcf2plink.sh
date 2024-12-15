#!/bin/bash
#
#SBATCH --job-name=plink
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --mem=20GB
#SBATCH --time=5:00:00
#SBATCH --output=~/logs/pp_%j.out
#SBATCH --error=~/logs/pp_%j.err

## initial preparation: convert vcf files to eigenstrat:
module load conda bcftools samtools htslib
conda activate plink2-2.00a5

#################################### create slim version of all-chromosome vcf file ####################################
for spec in pan gorilla pongo; do
  echo $spec
  while read chr in;do 
    echo $spec/chr$chr.vcf.gz >> plink/tmp/$spec.vcflist.txt
  done < plink/chroms.txt
  bcftools concat -n -f plink/tmp/$spec.vcflist.txt | bcftools annotate -x 'INFO' | bcftools annotate -Oz -x 'FORMAT' --write-index -o $spec/allchr.vcf.gz
done
  
## get variable positions
for spec in pan gorilla pongo; do
  echo $spec
  bcftools view -H $spec/allchr.vcf.gz | cut -f 1-2 | bgzip > $spec/variable_pos.txt.gz
done



#################################### create all-individual plink files ####################################
for spec in pan gorilla pongo; do
  echo $spec
  
  if [[ -f plink/mergelist.txt ]]; then rm plink/mergelist.txt; fi 
  
  while read chr in;do 
    echo $chr
    plink2 --vcf $spec/chr$chr.vcf.gz --max-alleles 2 --snps-only --make-pgen --maf 0.00 --out plink/tmp/$chr
    echo plink/tmp/$chr >> plink/mergelist.txt
  done < plink/chroms.txt
  echo "merge"
  plink2 --pmerge-list plink/mergelist.txt --make-pgen --out plink/$spec
done
echo "done"



