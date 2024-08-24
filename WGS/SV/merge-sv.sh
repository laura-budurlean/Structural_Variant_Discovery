#!/bin/bash
#SBATCH --job-name=merge-sv
#SBATCH -e /Proc_Results_hg38/err_files/.%J.%N.err
#SBATCH -o /Proc_Results_hg38/err_files/.%J.%N.out
#SBATCH -p #set value
#SBATCH -N 1
#SBATCH --ntasks-per-node=#set value
#SBATCH -t #set value
#SBATCH --mem=#set value
#SBATCH --partition=#set value


#this script merges lumpy and delly SV calls from the WGS data
path=/Proc_Results_hg38
cd $path

CUR_DIR=$PWD

#modules needed

#bedtools= /ri/shared/modules/bedtools/2.27.1/bin/
#module load bedtools
#module load python


for sample in #list your samples

do 
 
	OUT_DIR=merged_delly_lumpy_sv_calls
	if [ ! -d $OUT_DIR ]
	then
		mkdir $OUT_DIR
	fi
    #remove trailing tabs in WGS lumpy delly filtered files
  for vcf in $path/${sample}/SV_Delly_filtered/*.vcf
    do
    sed -i 's/[\t]*$//' $vcf
  done


#Intersecting DEL/DUP/INV from Lumpy only if reciprocally overlapped with Delly
	for svType in DEL DUP INV
	do
		bedtools intersect -a $path/${sample}/SV_Delly_filtered/${sample}.lumpy.${svType}.filtered.vcf -b $path/${sample}/SV_Delly_filtered/${sample}.delly.${svType}.filtered.vcf -f 0.5 -F 0.5 -wo | cut -f1-7 | sort | uniq > $path/merged_delly_lumpy_sv_calls/${sample}.${svType}.merged.bed
       	svNum=`cut -f1-3 $path/merged_delly_lumpy_sv_calls/${sample}.${svType}.merged.bed | sort | uniq | wc -l`
       	echo $sample $svType $svNum >> $path/sv-summary-0.5-del-dup-inv.txt
	done
	

	#For insertions, those will only come from Delly.
	cp $path/${sample}/SV_Delly_filtered/${sample}.delly.INS.filtered.vcf $path/merged_delly_lumpy_sv_calls/${sample}.delly.INS.filtered.bed
	svNum=`cat $path/${sample}/SV_Delly_filtered/${sample}.delly.INS.filtered.vcf | cut -f1-3 | sort | uniq | wc -l`
	echo $sample "INS" $svNum >> $path/sv-summary-0.5-insertions.txt
	

	#Inter-TRA from Lumpy if closer than 50bp to Delly calls
	buffer=50

	#Compare DELLY and Lumpy BND with 50 bp vicinity. Merging translocations that are withing 50bp of each other in Delly/Lumpy
	/scripts/WGS/SV/merge-BND.py $path/${sample}/SV_Delly_filtered/${sample}.lumpy.inter.filtered.vcf $path/${sample}/SV_Delly_filtered/${sample}.delly.inter.filtered.vcf $buffer > $path/merged_delly_lumpy_sv_calls/${sample}.inter.merged.bed
	svNum=`cut -f1-4 $path/merged_delly_lumpy_sv_calls/${sample}.inter.merged.bed | sort | uniq | wc -l`
	echo $sample "Inter-TRA" $svNum >> $path/sv-summary-0.5-Inter-Translocations.txt

	#remove DEL overlapped with Gap +- 50bp and Centromere +- 1kb
	cd $path/merged_delly_lumpy_sv_calls
	bedtools intersect -a <(cut -f1-3 $sample.DEL.merged.bed | sort | uniq) -b /scripts/WGS/SV/hg38.gap.50bp.flanking.goldenpath -f 0.1 -v > 1
	bedtools intersect -a 1 -b /scripts/WGS/SV/hg38.centromere.1kb.bed -v | sort -k1,1V -k2,2g -k3,3g > $sample.DEL.nogap.nocentromere.bed
	rm 1
 
done

#remove trailing tabs from filtered and merged lumpy delly calls
  for bed in /Proc_Results_hg38/merged_delly_lumpy_sv_calls/*.bed
    do
      sed -i 's/[\t]*$//' $bed
  done
