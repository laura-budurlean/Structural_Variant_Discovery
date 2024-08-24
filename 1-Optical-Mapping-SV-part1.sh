#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=#set value
#SBATCH -t #set value
#SBATCH --mem=#set value
#SBATCH --partition=#set value
#SBATCH --job-name=1-Optical-Mapping part1
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

cd $InDir


#1. Extract SVs: extract TLs from VCF or SMAP file. Bionano VCF file: You can use either the VCF or SMAP for this pipeline. Here we will be using the unfiltered SMAP file in the step below.

#grep SVTYPE=BND $SamInDir/${sample}.vcf | awk '{match($8, /.*CHR2=(.*);END=(.*);SVLEN=(.*);CT=(.*);EXPERIMENT=(.*)/, arr);print $1, $2, arr[1], arr[2], arr[4],$10}' OFS='\t' > $OutDir/translocation.${sample}.bed

grep trans $SamInDir/${sample}.smap | awk -F'\t' -v OFS='\t' '$9>0.05 {print "chr"$3, $7, "chr"$4, $8, $24, $26, $41, $44}' > $OutDir/translocation.${sample}.bed


#Getting the proper columns from the VCF file and obtain a BED file with proper orientation definitions. I.e. 3to5 5to3 3to3 etc. For this run this python script. Ignore the mv awk mv parts below...

python $path/convertTLs.py $OutDir/translocation.${sample}.bed $OutDir/fixedTLs.${sample}.bed
 

#2 extract deletion
cd $SamInDir
$path/bin/extractBionano-deletion $sample
mv ${sample}.Bionano.del.uniq.fi50.rmGap.rm1MBrec.bed $OutDir


#3 extract insertion
$path/bin/extractBionano-insertion $sample
mv ${sample}.Bionano.ins.uniq.fi50.confidence $OutDir


#4 extract inversion
$path/bin/extractINV-transform-bnd.awk $SamInDir/${sample} > $OutDir/${sample}.inversion.txt
cd $path


#5 remove false positive translocations 
#Run the Bionano-remove-recurrent.sh script. Then run 1-Optical-Mapping-SV-part2.sh, the continuation of this process.


