#!/bin/bash
#SBATCH --job-name=3-stratification
#SBATCH -e /compare_WGS_OM/err_files/%A_%a.err
#SBATCH -o /compare_WGS_OM/err_files/.%J.%N.out
#SBATCH -p #set value for partition
#SBATCH -N 1
#SBATCH --ntasks-per-node=#set value
#SBATCH -t #set value
#SBATCH --mem=#set value
#SBATCH --partition=#set value
#SBATCH --array=1-# #set values for number of samples


path=/scripts/Bionano/script.pack #set the path for the scripts
input=/gpfs/Labs/IPM/Project_ALL/compare_WGS_OM  #set the input of combined SV calls from WGS and OM, which is by default the output from script 2-compare.sh

#make a directory for somatic calls
mkdir $input/somatic


array=(sample1 sample2) #input sample names
sample=${array[$SLURM_ARRAY_TASK_ID-1]}

echo $sample
F=0.5
cd $input

#Now we remove DGV calls

#deletion and copy loss
#germline and somatic
bedtools intersect -a ${sample}.del.union.loss.bed -b $path/filter/DGV.deletion.loss.hg38.bed -wa -wb -F $F |awk '{a=(($4-$9)>=0?($4-$9):($9-$4));if(a <= $4*0.3){print $0}}'| cut -f1-5 | sort | uniq  > somatic/${sample}.del.loss.germline.bed 

bedtools intersect -a ${sample}.del.union.loss.bed -b somatic/${sample}.del.loss.germline.bed -f 1 -F 1 -v | sort | uniq > somatic/${sample}.del.loss.somatic.bed.1 

cd $input/somatic

bedtools intersect -a ${sample}.del.loss.somatic.bed.1 -b $path/filter/cross-sample-same.del.loss.somatic.bed -f 1 -F 1 | sort | uniq >> ${sample}.del.loss.germline.bed

bedtools intersect -a ${sample}.del.loss.somatic.bed.1 -b $path/filter/cross-sample-same.del.loss.somatic.bed -f 1 -F 1 -v >  ${sample}.del.loss.somatic.bed.2 ##

bedtools intersect -a ${sample}.del.loss.somatic.bed.2 -b $path/filter/5normal-deletion.bed -wa -wb  | awk '{a=(($4-$10)>=0?($4-$10):($10-$4));if(a <= $10*0.3){print $0}}'| cut -f1-5 > ${sample}.loss.norm.bed

bedtools intersect -a ${sample}.del.loss.somatic.bed.2 -b ${sample}.loss.norm.bed -f 1 -F 1 -v > ${sample}.del.loss.somatic.bed.3 ##
rm ${sample}.loss.norm.bed

bedtools intersect -a ${sample}.del.loss.somatic.bed.3 -b $path/filter/bionano-del-db.bed -f 0.5 -F 0.5 -wa -wb |  awk '{a=(($4-$9)>=0?($4-$9):($9-$4));if(a <= $9*0.3){print $0}}' OFS='-' | cut -f1-5 | sort | uniq >> ${sample}.del.loss.germline.bed

bedtools intersect -a ${sample}.del.loss.somatic.bed.3 -b ${sample}.del.loss.germline.bed -f 1 -F 1 -v > ${sample}.del.loss.somatic.bed.4 #the final output of deletions and copy loss??

cd $input

#Duplication and CNV gain
#germline and somtaic
bedtools intersect -a ${sample}.gain.dup.union.bed -b $path/filter/DGV.gain.dup.hg38.bed.HF  -f 0.7 -F 0.7  -wa | sort | uniq > somatic/${sample}.cnv.gain.germline.bed #HF
bedtools intersect -a ${sample}.gain.dup.union.bed -b somatic/${sample}.cnv.gain.germline.bed -f 1 -F 1 -v | sort | uniq  > somatic/${sample}.cnv.gain.somatic.bed.1 ##
bedtools intersect -a somatic/${sample}.cnv.gain.somatic.bed.1 -b $path/filter/bionano-dup-db.bed -wa -f 0.5 -F 0.5 | sort | uniq >> somatic/${sample}.cnv.gain.germline.bed
bedtools intersect -a somatic/${sample}.cnv.gain.somatic.bed.1 -b $path/filter/bionano-dup-db.bed  -f 0.5 -F 0.5 -v | sort | uniq > somatic/${sample}.cnv.gain.somatic.bed.2

cd $input/somatic
bedtools intersect -a ${sample}.cnv.gain.somatic.bed.2 -b $path/filter/cross-sample-same.cnv.gain.somatic.bed -f 1 -F 1  | sort | uniq >> ${sample}.cnv.gain.germline.bed
bedtools intersect -a ${sample}.cnv.gain.somatic.bed.2 -b $path/filter/cross-sample-same.cnv.gain.somatic.bed -f 1 -F 1 -v > ${sample}.cnv.gain.somatic.bed.4
cd $input

#insertion
#germline and somtaic
bedtools intersect -a ${sample}.ins.union.bed -b $path/filter/DGV.dup-ins.bed  -f 1 -F 1 -wa | cut -f1-5|  sort | uniq > somatic/${sample}.ins.germline.bed
bedtools intersect -a ${sample}.ins.union.bed -b $path/filter/DGV.dup-ins.bed -wa -wb  | awk '{a=(($4-$12)>=0?($4-$12):($12-$4));if(a <= $12*0.3){print $0}}' | cut -f1-5 | sort | uniq >> somatic/${sample}.ins.germline.bed
bedtools intersect -a ${sample}.ins.union.bed -b somatic/${sample}.ins.germline.bed -f 1 -F 1 -v | sort |uniq > somatic/${sample}.ins.somatic.bed.1

cd $input/somatic
bedtools intersect -a ${sample}.ins.somatic.bed.1 -b $path/filter/cross-sample-same.ins.somatic.bed -f 1 -F 1 | sort | uniq >> ${sample}.ins.germline.bed
bedtools intersect -a ${sample}.ins.somatic.bed.1 -b $path/filter/cross-sample-same.ins.somatic.bed -f 1 -F 1 -v > ${sample}.ins.somatic.bed.2 ##
bedtools intersect -a ${sample}.ins.somatic.bed.2 -b $path/filter/5normal-insertion.bed -wa -wb  | awk '{a=(($4-$10)>=0?($4-$10):($10-$4));if(a <= $10*0.3){print $0}}'| cut -f1-5 > ${sample}.ins.norm.bed
bedtools intersect -a ${sample}.ins.somatic.bed.2 -b ${sample}.ins.norm.bed -f 1 -F 1 -v > ${sample}.ins.somatic.bed.3 ##
rm ${sample}.ins.norm.bed
bedtools intersect -a ${sample}.ins.somatic.bed.3 -b $path/filter/bionano-ins-db.bed -f 0.5 -F 0.5 -wa -wb |  awk '{a=(($4-$9)>=0?($4-$9):($9-$4));if(a <= $9*0.3){print $0}}' OFS='-' | cut -f1-5 | sort | uniq >> ${sample}.ins.germline.bed
bedtools intersect -a ${sample}.ins.somatic.bed.3 -b ${sample}.ins.germline.bed -f 1 -F 1 -v > ${sample}.ins.somatic.bed.4 
cd $input

#inversion
#germline and somtaic
bedtools intersect -a ${sample}.inv.union.bed -b $path/filter/DGV.inversion.bed -wa -wb | awk '{a=(($4-$12)>=0?($4-$12):($12-$4));if(a <= $12*0.3){print $0}}' | cut -f1-5| sort | uniq > somatic/${sample}.inv.germline.bed
bedtools intersect -a ${sample}.inv.union.bed -b somatic/${sample}.inv.germline.bed -f 1 -F 1 -v | sort | uniq > somatic/${sample}.inv.somatic.bed.1 ##

cd $input/somatic
bedtools intersect -a ${sample}.inv.somatic.bed.1 -b $path/filter/cross-sample-same.inv.somatic.bed -f 1 -F 1 | sort | uniq >> ${sample}.inv.germline.bed
bedtools intersect -a ${sample}.inv.somatic.bed.1 -b $path/filter/cross-sample-same.inv.somatic.bed -f 1 -F 1 -v > ${sample}.inv.somatic.bed.2 ##
bedtools intersect -a ${sample}.inv.somatic.bed.2 -b $path/filter/bionano-inv-db.bed -wa -wb | awk '$2==$8 || $2==$7 || $3==$8 || $3==$7' | cut -f1-5| sort | uniq >> ${sample}.inv.germline.bed
bedtools intersect -a ${sample}.inv.somatic.bed.2 -b ${sample}.inv.germline.bed -f 1 -F 1 -v > ${sample}.inv.somatic.bed.4 
rm *.1 *.2  *.3
