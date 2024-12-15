#!/bin/bash
#
#SBATCH --job-name=filt
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=3
#SBATCH --mem=150GB
#SBATCH --time=23:00:00
#SBATCH --output=./logs/filt_%j.out
#SBATCH --error=./logs/filt_%j.err
#SBATCH --array=1

## initial preparation: convert vcf files to eigenstrat:
module load conda bcftools/1.20  samtools htslib
module load java/17.0.6
basedir=""
cd $basedir
chrom=$SLURM_ARRAY_TASK_ID
if [[ ${chrom} == 23 ]]; then chrom="X"; fi

# method should be "filt" if doing the basic filtering procedure
# method should be "check" if testing the number of SNPs
# method should be "final" if generating the final reliable SNP set
# method should be "test" if generating bcftools stats test output

method="filt"

if [[ $method == "filt" ]] ; then
# within central 98% for each individual
# allele imbalance per individual
# mapability 35

for spec in pan_all gorilla_all pongo_all; do
  echo $spec
    lw=3;hw=4; if [[ ${chrom} == "X" ]]; then lw=5;hw=6; fi
    allctr=""
    allalb=""
    while read ind in; do
      ctf=$(awk -v var="$ind" '$2==var' $basedir/analysis/coverage_cutoffs.txt)
      ctl=$(echo ${ctf} | cut -d " " -f $lw ); ctu=$(echo ${ctf} | cut -d " " -f $hw )
      if [[ ${ind} == "C" ]]; then ind="C_a"; fi
      if [[ ${ind} == "Suma" ]]; then ind="Sum"; fi
      samp=$( bcftools query -l $basedir/$spec/chr$chrom.vcf.gz | nl -v 0 | awk -v var="$ind" '$2==var' | cut -f 1 | sed 's/ //g' )
      ctr="'(FMT/DP[$samp]<$ctl | FMT/DP[$samp]>$ctu)'"
      alb="'(FMT/AD[$samp:0]/FMT/DP[$samp]<0.15 | FMT/AD[$samp:1]/FMT/DP[$samp]<0.15) & FMT/GT[$samp]==\"het\"'"
      allctr="$allctr | bcftools filter -e $ctr -S ."
      allalb="$allalb | bcftools filter -e $alb -S ."
    done < $spec.individuals
    allctr=$(echo $allctr | sed "s/| //")
    allalb=$(echo $allalb | sed "s/| //")
  
    bcftools view $basedir/$spec/chr$chrom.vcf.gz | eval "$allctr" | eval $allalb | java -jar ~/jvarkit.jar vcfbigwig -T mapability -B k36.Umap.MultiTrackMappability.bw | bcftools view -Oz -W -o $basedir/$spec/chr$chrom.filteranno.vcf.gz &
  done
wait

exit
fi


if [[ $method == "check" ]] ; then

## check that the ratio of original SNPs and filtered SNPs is approximately equal (i.e. no data loss)
for spec in pan_all gorilla_all pongo_all; do
  echo $spec
  for chrom in {1..22} X; do
  v1=$(zcat $basedir/$spec/chr$chrom.filteranno.vcf.gz | wc -l )
  v2=$(zcat $basedir/$spec/chr$chrom.vcf.gz | wc -l )
  res=$(echo "scale=5 ; $v1 / $v2" | bc )
  echo "$chrom $res" >> ~/logs/$spec.check.txt &
  done
  wait
done

exit
fi



if [[ $method == "final" ]] ; then

## generate the final filtered SNP set (analogous to vcf2plink.sh)
for spec in pan_all gorilla_all pongo_all; do
  echo $spec
  rm plink/tmp/$spec.filt.vcflist.txt
  rm $basedir/$spec/allchr.filtered.vcf.gz
  for chrom in {1..22} X; do
    echo $basedir/$spec/chr$chrom.filteranno.vcf.gz >> plink/tmp/$spec.filt.vcflist.txt
  done
  bcftools concat --threads 4 -n -f plink/tmp/$spec.filt.vcflist.txt | bcftools view --threads 4 -V indels | bcftools view --threads 4 -U -v snps | bcftools view --threads 4 -a | bcftools view --threads 4 -e 'GT=="./."' |  bcftools view --threads 4 -m2 -M2 | egrep "#|mapability=1" | bcftools annotate --threads 4 -x 'INFO' | bcftools annotate --threads 4 -Oz -x 'FORMAT' -W -o $basedir/$spec/allchr.filtered.vcf.gz &
done
wait

exit
fi


if [[ $method == "test" ]] ; then

## generate bcftools stats output
for chrom in {1..22} X; do
  echo $chrom
  for spec in pan_all gorilla_all pongo_all; do
    bcftools stats -s - $basedir/$spec/allchr.filtered.vcf.gz > $basedir/analysis/$spec.bstata.txt &
    bcftools stats -s - $basedir/$spec/chr$chrom.filteranno.vcf.gz > $basedir/analysis/$spec.$chrom.bstatf.txt &
    bcftools stats -s - $basedir/$spec/chr$chrom.vcf.gz > $basedir/analysis/$spec.$chrom.bstatr.txt &
  done
  wait
done

exit
fi

