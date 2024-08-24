#!/bin/bash
#SBATCH --job-name=2-compare
#SBATCH -e /compare_WGS_OM/err_files/%A_%a.err
#SBATCH -o /compare_WGS_OM/err_files/.%J.%N.out
#SBATCH -p #set value for partition
#SBATCH -N 1
#SBATCH --ntasks-per-node= #set value
#SBATCH -t #set value
#SBATCH --mem= #set value
#SBATCH --partition=#set value 
#SBATCH --array=#-# #set value for number of samples

#NOTE BEFORE STARTING: set up the wgs input directory, which contains folders with names: cnv, ins, inv, del, dup

#NOTE BEFORE STARTING: in the following lines, set the input directory of OM data, which is the output from 1-Optical-Mapping-SV.sh scripts; the input directory for WGS data, and the path to output

path=/scripts/Bionano/script.pack #path to the scripts
OMinp=/gpfs/Labs/IPM/Project_ALL/Bionano_RareVariantData/scripts_output #set optical mapping input
WGSinp=/WGS_data/Proc_Results_hg38/merged_delly_lumpy_sv_calls #set WGS input path for SVs
FreeCNVs=/path/to/FreeCNV/output #wgs CNV path
output=/path/to/compare_WGS_OM  #set output directory for this script
bedtools=/gpfs/ri/shared/modules7/bedtools/2.30.0/bin #bedtools path

#remove decimals from the OM input files in scripts_output. Bionano sometimes adds decimals to coordinate locations which will give errors if not fixed.
for file in $OMinp/*.bed 
  do
    sed -i 's/\([0-9]*\)\.[0-9]*/\1/g' $file 
done


#START
#set sample array
array=(sample1 sample2)
sample=${array[$SLURM_ARRAY_TASK_ID-1]}

echo $sample

#I. deletion
	# wgs union and CNV loss

cd $output

grep loss $FreeCNV/${sample}.bam_CNVs_with_chr.bed | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $3-$2}' > $output/${sample}.cnv.loss.bed

#intersect copy loss w deletions
bedtools intersect -a $output/${sample}.cnv.loss.bed -b $WGSinp/${sample}.DEL.nogap.nocentromere.bed -f 0.8 -F 0.8 -v > $output/${sample}.cnv.loss.bed.specific
 
cat $WGSinp/${sample}.DEL.nogap.nocentromere.bed $output/${sample}.cnv.loss.bed.specific | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $3-$2}'| sort | uniq > $output/${sample}.wgs-del-cnv.bed #wgs deletions and wgs copy losses, not confidence filtered yet

rm ${sample}.cnv.loss.bed ${sample}.cnv.loss.bed.specific
	
# intersect deletion
bedtools intersect -a $OMinp/${sample}.intra.5-3.bed -b $output/${sample}.wgs-del-cnv.bed -F 0.5 -wa -wb | awk '{a=(($4-$9)>=0?($4-$9):($9-$4));if(a <= $9*0.3){print $0}}' | sort | uniq > $output/${sample}.del.confidence.bed #wgs deletions and wgs copy losses with confidence filters applied

bedtools intersect -a $OMinp/${sample}.intra.5-3.bed -b $output/${sample}.del.confidence.bed -f 1 -F 1 -v | awk '{print $1,$2,$3,$4,"optical"}' OFS="\t" > $output/${sample}.saphyr.specific.deletion.bed

cat $output/${sample}.del.confidence.bed | awk 'BEGIN{OFS="\t"}{print $6,$7,$8,$9,"both"}' | sort | uniq > $output/${sample}.wgs.del.overlap

bedtools intersect -a $output/${sample}.wgs-del-cnv.bed  -b $output/${sample}.wgs.del.overlap -f 1 -F 1 -v | sort | uniq | awk '{print $0, "wgs"}' OFS="\t" > $output/${sample}.wgs.specific.deletion.bed


#deletion union
cat $output/${sample}.wgs.del.overlap $output/${sample}.saphyr.specific.deletion.bed $output/${sample}.wgs.specific.deletion.bed | sort | uniq > $output/${sample}.del.union.loss.bed
rm $output/${sample}.wgs.del.overlap 


#inversion, intersect inversion
bedtools intersect -a $OMinp/${sample}.5533.bed -b $WGSinp/${sample}.INV.merged.bed -F 0.1 -wa -wb | cut -f1-7 | awk 'BEGIN{OFS="\t"}{print $0,$7-$6}' |  awk '{a=(($4-$8)>=0?($4-$8):($8-$4));if(a <= $8*0.3){print $0}}' | sort | uniq > $output/${sample}.inv.confidence.bed

#inversion union
bedtools intersect -a $OMinp/${sample}.5533.bed -b $output/${sample}.inv.confidence.bed -v -f 1 -F 1 | awk '{print $0, "optical"}' OFS="\t" > $output/${sample}.saphyr.specific.inv.bed
cat $output/${sample}.inv.confidence.bed | cut -f5-8 | awk '{print $0, "both"}' OFS="\t" > $output/${sample}.inv.overlap

bedtools intersect -a $WGSinp/${sample}.INV.merged.bed -b $output/${sample}.inv.overlap -v -f 1 -F 1  | awk '{print $1, $2, $3, $3-$2, "wgs"}' OFS="\t" > $output/${sample}.wgs.specific.inv.bed

cat $output/${sample}.saphyr.specific.inv.bed $output/${sample}.inv.overlap $output/${sample}.wgs.specific.inv.bed | sort | uniq | sed 's/ //g' | sed 's/\t\t/\t/g' > $output/${sample}.inv.union.bed
rm $output/${sample}.inv.overlap

#duplication and gain, wgs
cut -f1-3 $WGSinp/${sample}.DUP.merged.bed | awk 'BEGIN{OFS="\t"}{print $0, $3-$2, "dup"}' > $output/${sample}.dup.bed

cd $output

grep gain $FreeCNV/${sample}.bam_CNVs_with_chr.bed | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $3-$2}' > $output/${sample}.cnv.gain.bed

bedtools intersect -a $output/${sample}.cnv.gain.bed -b $output/${sample}.dup.bed -f 0.7 -F 0.7 -v | awk '{print $0, "gain"}' OFS="\t" > $output/${sample}.cnv.gain.specific.bed

cat $output/${sample}.cnv.gain.specific.bed ${sample}.dup.bed | awk '{print $0, "wgs"}' OFS="\t" |sort | uniq > $output/${sample}.wgs.gain.dup.union.bed

rm $output/${sample}.cnv.gain.specific.bed $output/${sample}.cnv.gain.bed ${sample}.dup.bed
	
#compare wgs and bionano
bedtools intersect -b $output/${sample}.wgs.gain.dup.union.bed -a $OMinp/${sample}.duplication.bed -F 0.5 -wa -wb | awk '{a=(($4-$9)>=0?($4-$9):($9-$4));if(a <= $9*0.3){print $0}}' > $output/${sample}.dup.overlap.bionano

cut -f6-10 $output/${sample}.dup.overlap.bionano | awk '{print $0, "both"}' OFS="\t" > $output/${sample}.dup.overlap.wgs

bedtools intersect -a $OMinp/${sample}.duplication.bed -b $output/${sample}.dup.overlap.bionano -f 1 -F 1 -v | awk '{print $0, "optical"}' OFS="\t" > $output/${sample}.dup.bionano.specific
cat $output/${sample}.dup.overlap.wgs $output/${sample}.dup.bionano.specific $output/${sample}.wgs.gain.dup.union.bed | awk '!a[$1$2$3]++' > $output/${sample}.gain.dup.union.bed.temp
	
#compare to bionano insertion    
bedtools intersect -a $output/${sample}.gain.dup.union.bed.temp -b $OMinp/${sample}.intra.3-5.bed -wa -wb | awk '{a=(($4-$10)>=0?($4-$10):($10-$4));if(a <= $10*0.3){print $0}}' OFS="\t" | grep dup | sed 's/wgs/both/g' > $output/${sample}.wgs-dup.om-ins.bed

bedtools intersect -a $output/${sample}.gain.dup.union.bed.temp -b $output/${sample}.wgs-dup.om-ins.bed -f 1 -F 1 -v > $output/${sample}.gain.dup.union.bed.temp.temp
cat $output/${sample}.wgs-dup.om-ins.bed $output/${sample}.gain.dup.union.bed.temp.temp | cut -f1-6| sort | uniq > $output/${sample}.gain.dup.union.bed 

#rm DUP from OM insertion 
cut -f7-11 $output/${sample}.wgs-dup.om-ins.bed > $output/${sample}.om-ins-dup.bed
                                 
bedtools intersect -a $OMinp/${sample}.intra.3-5.bed -b $output/${sample}.om-ins-dup.bed -f 1 -F 1 -v > $OMinp/${sample}.intra.3-5.bed.mod


for mod2 in /Bionano_RareVariantData/scripts_output/*.mod
    do
      sed -i 's/[\t]*$//' $mod2
  done

# insertion
cut -f1-3,6 $WGSinp/${sample}.delly.INS.filtered.bed | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, length($4)}' > $output/${sample}.wgs.ins.bed

# intersect insertion   
bedtools intersect -a $OMinp/${sample}.intra.3-5.bed.mod -b $output/${sample}.wgs.ins.bed -wa -wb |  awk '{a=(($4-$9)>=0?($4-$9):($9-$4));if(a <= $9*0.3){print $0}}' | sort | uniq > $output/${sample}.confidence.insertion.bed 

bedtools intersect -a $OMinp/${sample}.intra.3-5.bed.mod -b $output/${sample}.confidence.insertion.bed -f 1 -F 1 -v | cut -f1-4 | awk '{print $0, "optical"}' OFS="\t" > $output/${sample}.saphyr.specific.ins.bed
cut -f6-9 $output/${sample}.confidence.insertion.bed | awk '{print $0, "both"}' > $output/${sample}.wgs.ins.overlap

bedtools intersect -a $output/${sample}.wgs.ins.bed -b $output/${sample}.wgs.ins.overlap -v -f 1 -F 1 | awk '{print $0, "wgs"}' OFS="\t" | sort | uniq > $output/${sample}.wgs.specific.ins.bed

# insertion union
cd $output
cat $output/${sample}.saphyr.specific.ins.bed $output/${sample}.wgs.ins.overlap $output/${sample}.wgs.specific.ins.bed > $output/${sample}.ins.union.bed #an output of this script
rm $output/${sample}.wgs.ins.bed $output/${sample}.wgs.ins.overlap

# compare TL
python /Bionano/script.pack/compare-Bionano-wgs-tl.py $OMinp/${sample}.inter.txt $WGSinp/${sample}.inter.merged.bed $output

#remove trailing tabs 
for bed in $output/*.bed
    do
      sed -i 's/[\t]*$//' $bed 
  done

for mod in $output/*.mod
    do
      sed -i 's/[\t]*$//' $mod
  done
