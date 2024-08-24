#!/bin/bash
#SBATCH --job-name=filter-sv
#SBATCH -e /Proc_Results_hg38/err_files/.%J.%N.err
#SBATCH -o /Proc_Results_hg38/err_files/.%J.%N.out
#SBATCH -N 1
#SBATCH --ntasks-per-node=#set value
#SBATCH -t #set value
#SBATCH --mem=#set value
#SBATCH --partition=#set value


#Combining and filtering Delly and Lumpy SV files for WGS data

WGSinp=/path/to/wgs/files
cd $WGSinp 

#sample vcf files

for cell in #list your samples
do

mkdir $WGSinp/${cell}/SV_Delly_filtered
outdir=$WGSinp/${cell}/SV_Delly_filtered

#DELLY part: filter for delly DEL/DUP/INV. SVLEN >= 50bp;PASS;auto+X
  for svType in DEL DUP INV
  do
		grep -v "^#" $WGSinp/${cell}/${cell}.delly.${svType}.vcf | grep -v chrY | gawk '$7 == "PASS"{start=$2;match($8, /.*;END=([0-9]*);.*/, end);match($8, /.*;PE=([0-9]*);.*/, PE);if(PE[1]==""){PE[1]=0};match($8, /.*;SR=([0-9]*);.*/, SR);if(SR[1]==""){SR[1]=0};len=(end[1]-start);if(len >= 50 || len <= -50) {if(end[1] < start){t=start;start=end[1];end[1]=t};split($10, arr, ":");print $1, start, end[1], PE[1], SR[1], arr[1], arr[3], arr[4], arr[8]}}' OFS='\t' | sort -k1,1V -k2,2g -k3,3g > $outdir/${cell}.delly.${svType}.filtered.vcf
		grep -v "^#" $WGSinp/${cell}/${cell}.delly.${svType}.vcf | grep -v chrY | gawk '$7 == "PASS"{start=$2;match($8, /.*;END=([0-9]*);.*/, end);len=(end[1]-start);if(len < 50 && len > -50) {if(end[1] < start){t=start;start=end[1];end[1]=t};split($10, arr, ":");print $1, start, end[1], arr[1], arr[3], arr[4], arr[8]}}' OFS='\t' | sort -k1,1V -k2,2g -k3,3g > $outdir/${cell}.delly.${svType}.lessthan50bp.filtered.vcf
	
#Delly INS. SVLEN >= 50bp;PASS;auto+X, old Delly
	  grep -v "^#" $WGSinp/${cell}/${cell}.delly.INS.vcf | grep -v chrY | gawk '$7 == "PASS"{match($8, /.*;INSLEN=([0-9]*);.*/, l);len=l[1];if(len >= 50) {match($8, /.*;PE=([0-9]*);.*/, PE);if(PE[1]==""){PE[1]=0};match($8, /.*;SR=([0-9]*);.*/, SR);if(SR[1]==""){SR[1]=0};split($10, arr, ":");print $1, $2, $2+1, PE[1], SR[1], $5, arr[1], arr[3], arr[4]}}' OFS='\t' | sort -k1,1V -k2,2g -k3,3g > $outdir/${cell}.delly.INS.filtered.vcf
  	grep -v "^#"$WGSinp/${cell}/${cell}.delly.INS.vcf | grep -v chrY | gawk '$7 == "PASS"{match($8, /.*;INSLEN=([0-9]*);.*/, l);len=l[1];if(len < 50) {match($8, /.*;PE=([0-9]*);.*/, PE);if(PE[1]==""){PE[1]=0};match($8, /.*;SR=([0-9]*);.*/, SR);if(SR[1]==""){SR[1]=0};split($10, arr, ":");print $1, $2, $2+1, PE[1], SR[1], $5, arr[1], arr[3], arr[4]}}' OFS='\t' | sort -k1,1V -k2,2g -k3,3g > $outdir/${cell}.delly.INS.lessthan50bp.filtered.vcf

#Delly BND only inter is called by delly, for delly, chrA is always > chrB, so need to switch

	  grep -v "^#" $WGSinp/${cell}/${cell}.delly.BND.vcf | grep -v chrY | gawk '$7 == "PASS"{ \
		match($8, /.*;CHR2=(.*);END=([0-9]*).*PE=([0-9]*);.*/, end); \
		match($8, /.*;SR=([0-9]*);.*/, SR); \
		if(SR[1]==""){SR[1]=0};split($10, arr, ":"); \
		match($8, /.*;CT=([0-9]*)to([0-9]*);.*/, symbol); \
		s=symbol[1]""symbol[2]; \
		if(s=="53"){s="N[["} \
		else if(s=="35"){s="]]N"} \
		else if(s=="55"){s="N]]"} \
		else {s=="[[N"} \
		print end[1], end[2], $1, $2, s, end[3], SR[1], arr[1], arr[3], arr[4]}' OFS='\t' | sort -k1,1V -k2,2g -k3,3V -k4,4g > $outdir/${cell}.delly.inter.vcf
  done

	
#LUMPY PART: filter for Lumpy DEL/DUP/INV. SVLEN >= 50bp;auto+X
	for svType in DEL DUP INV
	do
		grep -v "^#" $WGSinp/${cell}/${cell}.lumpy.sv.vcf | grep -v chrY | grep SVTYPE=$svType | gawk '{match($8, /.*;SVLEN=(-?[0-9]*);END=([0-9]*);.*/, end);if(end[1] >= 50 || end[1] <= -50) {start=$2;if(end[2] < start){t=start;start=end[2];end[2]=t};split($10, arr, ":");print $1, start, end[2], arr[3], arr[4], arr[1], arr[5], $6, arr[8]}}' OFS='\t' | sort -k1,1V -k2,2g -k3,3g > $outdir/${cell}.lumpy.${svType}.filtered.vcf
		grep -v "^#" $WGSinp/${cell}/${cell}.lumpy.sv.vcf | grep -v chrY | grep SVTYPE=$svType | gawk '{match($8, /.*;SVLEN=(-?[0-9]*);END=([0-9]*);.*/, end);if(end[1] < 50 && end[1] > -50) {start=$2;if(end[2] < start){t=start;start=end[2];end[2]=t};split($10, arr, ":");print $1, start, end[2], arr[3], arr[4],  arr[1], arr[5], $6, arr[8]}}' OFS='\t' | sort -k1,1V -k2,2g -k3,3g > $outdir/${cell}.lumpy.${svType}.lessthan50bp.filtered.vcf
	done

	PE_SR=10
	SR=2

	
	for svType in DEL DUP INV
	do                                                   
		awk -v pe_sr=$PE_SR -v sr=$SR '$4+$5>=pe_sr && $5>=sr' $outdir/$cell.delly.${svType}.filtered.vcf > $outdir/1
		rm $outdir/$cell.delly.${svType}.filtered.vcf 
		mv $outdir/1 $outdir/$cell.delly.${svType}.filtered.vcf
   
		awk -v pe_sr=$PE_SR -v sr=$SR '$4+$5>=pe_sr && $5>=sr' $outdir/$cell.lumpy.${svType}.filtered.vcf > $outdir/2
		rm $outdir/$cell.lumpy.${svType}.filtered.vcf 
		mv $outdir/2 $outdir/$cell.lumpy.${svType}.filtered.vcf
	done

#LUMPY BND. inter, make sure chrA < chrB in the bed file, orientation changed with chromosome
	grep -v "^#" $WGSinp/${cell}/${cell}.lumpy.sv.vcf | grep -v chrY | grep SVTYPE=BND | gawk '{ \
	split($1, pos1, "r");split($10, arr, ":");match($5, /.*chr(.*):([0-9]*).*/, pos2); \
	if(pos2[1]!=pos1[2]){ \
		match($5, /(.*)chr(.*):([0-9]*)(.*)/, symbol); \
		s=symbol[1]""symbol[4]; \
		gsub(/[AGCT]/, "N", s); \
		if(pos2[1] < pos1[2]){ \
			t=pos1[2];pos1[2]=pos2[1];pos2[1]=t; \
			t=$2;$2=pos2[2];pos2[2]=t; \
			if(s=="]]N"){s="N[["} \
			else if(s=="N[["){s="]]N"} \
			else {} \
		} \
		print "chr"pos1[2], $2, "chr"pos2[1], pos2[2], s, arr[3], arr[4], arr[1], arr[5], $6, arr[8]}}' OFS='\t' | sort -k1,1V -k2,2g -k3,3V -k4,4g | uniq > $outdir/${cell}.lumpy.inter.vcf

#LUMPY BND. intra, make sure posA < posB in the bed file, orientation changed with chromosome
		grep -v "^#" $WGSinp/${cell}/${cell}.lumpy.sv.vcf | grep -v chrY | grep SVTYPE=BND | gawk '{ \
			split($1, pos1, "r");split($10, arr, ":"); \
			match($5, /.*chr(.*):([0-9]*).*/, pos2); \
			if(pos2[1]==pos1[2]){ \
				match($5, /(.*)chr(.*):([0-9]*)(.*)/, symbol); \
				s=symbol[1]""symbol[4];gsub(/[AGCT]/, "N", s); \
				if(pos2[2] < $2){ \
					t=pos2[2];pos2[2]=$2;$2=t; \
					if(s=="]]N"){s="N[["} \
					else if(s=="N[["){s="]]N"} \
					else {} \
				} \
				print $1, $2, "chr"pos2[1], pos2[2], s, arr[3], arr[4], arr[1], arr[5], $6, arr[8]}}' OFS='\t' | sort -k1,1V -k2,2g -k3,3V -k4,4g | uniq > $outdir/${cell}.lumpy.intra.vcf

#Filter inter BND, called by Lumpy and Delly. PE or SR > a threshold
		PE_SR=10
		SR=2
		
	awk -v pe_sr=$PE_SR -v sr=$SR '$6+$7>=pe_sr && $7>=sr' $outdir/$cell.delly.inter.vcf > $outdir/$cell.delly.inter.filtered.vcf
	awk -v pe_sr=$PE_SR -v sr=$SR '$6+$7>=pe_sr && $7>=sr' $outdir/$cell.lumpy.inter.vcf > $outdir/$cell.lumpy.inter.filtered.vcf

      ### filter intra BND, called only by Lumpy, length >= 50bp,	
	awk -v pe_sr=$PE_SR -v sr=$SR '($4-$2 >= 50) && ($10 > 100) && ($6+$7>=pe_sr) && ($7>=sr)' $outdir/${cell}.lumpy.intra.vcf > $outdir/${cell}.lumpy.intra.filtered.bed
	awk '$4-$2 > 1000000' $outdir/${cell}.lumpy.intra.filtered.bed > $outdir/${cell}.lumpy.intra.filtered.1Mb.bed
	awk '$4-$2 <= 1000000' $outdir/${cell}.lumpy.intra.filtered.bed > $outdir/${cell}.lumpy.intra.filtered.lt.1Mb.bed
	
	cd ../
done
