#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=#set value
#SBATCH -t #set value
#SBATCH --mem=#set value
#SBATCH --partition= #set value
#SBATCH --job-name=remove-recurrent
#SBATCH -o /scripts/Bionano/script.pack/out_err_files/.%J.%N.out
#SBATCH -e /scripts/Bionano/script.pack/out_err_files/.%J.%N.err



#Remove recurrent translocations


cd /directory/to/Bionano_RareVariantData

for sample in #list samples

do
	#Run R-script 
Rscript --vanilla /scripts/Bionano/script.pack/remove-recurrent.R /Bionano_RareVariantData/scripts_output/fixedTLs.${sample}.bed /Bionano_RareVariantData/scripts_output 

mv /Bionano_RareVariantData/scripts_output/.bed-non-recurrent-translocation-from-raw-unfiltered-list.txt /gpfs/Labs/IPM/Project_ALL/Bionano_RareVariantData/scripts_output/${sample}.bed-non-recurrent-translocation-from-raw-unfiltered-list.txt

done
