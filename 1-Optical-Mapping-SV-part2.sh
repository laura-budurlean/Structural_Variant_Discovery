#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=#set value
#SBATCH -t #set value
#SBATCH --mem=#set value
#SBATCH --partition=#set value
#SBATCH --job-name=1-Optical-Mapping-part2
#SBATCH -o /scripts/Bionano/script.pack/out_err_files/.%J.%N.out
#SBATCH -e /scripts/Bionano/script.pack/out_err_files/.%J.%N.err
#SBATCH --array=1-# #set value for number of samples


#modules needed

#module load bedtools
#module load python
#module load R

array=(sample1 sample2) #set sample names
sample=${array[$SLURM_ARRAY_TASK_ID-1]}


#NOTE: Make a folder for each sample you have and put the Bionano SMAP or VCF file in it and rename the file with the corresponding sample name. Example: .../sample3/sample3.vcf or sample3.smap

#Set the directory for the path to pack of scripts, input, output, and name of samples

path=/scripts/Bionano/script.pack  #set the path to the pack of scripts
InDir=/path/to/Bionano_RareVariantData	#set the directory that contains the folder with your sample names
OutDir=/gpfs/Labs/IPM/Project_ALL/Bionano_RareVariantData/scripts_output  #set the directory for outputing processed Bionno SV calls
SamInDir=/gpfs/Labs/IPM/Project_ALL/Bionano_RareVariantData/${sample}  #set sample directories

#Continue from part1


#6 extract duplication
grep duplication $SamInDir/${sample}.smap | awk '{print "chr"$3, int($7), int($8), int($8-$7), "dup"}' OFS="\t" | sort | uniq | sed 's/chr23/chrX/g'| sed 's/chr24/chrY/g' > $OutDir/${sample}.duplication.bed
 
#Further SV formating to make Bionano SVs compatible with SVs from WGS (orientations) 
cd $OutDir

#Seperate large intra-chr alteration by orientation:
awk '$1==$3' ${sample}.bed-non-recurrent-translocation-from-raw-unfiltered-list.txt | awk '$5=="5to3"'| cut -f1-2,4,6 | awk '{print $1,$2,$3,$3-$2,$4}'| sed 's/ /\t/g' > ./${sample}.intra.5-3.bed
 
awk '$1==$3' ${sample}.bed-non-recurrent-translocation-from-raw-unfiltered-list.txt | awk '$5=="5to5" || $5=="3to3"'| cut -f1-2,4 | awk '{print $1,$2,$3,$3-$2,$4}'| sed 's/ /\t/g' > ./${sample}.intra.5-5and3-3.bed
 
awk '$1==$3' ${sample}.bed-non-recurrent-translocation-from-raw-unfiltered-list.txt | awk '$5=="3to5"'| cut -f1-2,4,6 | awk '{print $1,$2,$3,$3-$2,$4}'| sed 's/ /\t/g' > ./${sample}.intra.3-5.bed
 
 
#inter-chr TL
awk '$1!=$3' ${sample}.bed-non-recurrent-translocation-from-raw-unfiltered-list.txt | cut -f1-5| sort | uniq  > ./${sample}.inter.txt
 
#add deletion
cut -f1-5 ${sample}.Bionano.del.uniq.fi50.rmGap.rm1MBrec.bed >> ./${sample}.intra.5-3.bed

#add insertion
 cut -f1-5 ${sample}.Bionano.ins.uniq.fi50.confidence >> ./${sample}.intra.3-5.bed
 
#inversion
grep -i partial ${sample}.inversion.txt | cut -f1-2 | grep -v NA| sed 's/:/\t/g'| sed 's/chr23/chrX/g'|sed 's/chr24/chrY/g' >> ./${sample}.intra.5-5and3-3.bed.1
grep -i partial ${sample}.inversion.txt | cut -f1,3 | grep -v NA| sed 's/:/\t/g'| sed 's/chr23/chrX/g'|sed 's/chr24/chrY/g' >> ./${sample}.intra.5-5and3-3.bed.1
grep -i paired ${sample}.inversion.txt | cut -f1,3 | grep -v NA| sed 's/:/\t/g'| sed 's/chr23/chrX/g' |sed 's/chr24/chrY/g' >> ./${sample}.intra.5-5and3-3.bed.1
 
cat ${sample}.intra.5-5and3-3.bed.1 | sort | uniq | awk '{print $1,$2,$3,$3-$2,$4}'| sed 's/ /\t/g'  > ${sample}.inversion.bed
rm ${sample}.intra.5-5and3-3.bed.1
cat ${sample}.intra.5-5and3-3.bed ${sample}.inversion.bed > ${sample}.5533.bed #inversion file
rm ${sample}.inversion.bed ${sample}.intra.5-5and3-3.bed ${sample}.Bionano.del.uniq.fi50.rmGap.rm1MBrec.bed ${sample}.Bionano.ins.uniq.fi50.confidence  ${sample}.inversion.txt # ${sample}.bed-non-recurrent-translocation-from-raw-unfiltered-list.txt  
 
#remove trailing tabs
  for file in $OutDir/*
    do
    sed -i 's/[\t]*$//' $file
  done
