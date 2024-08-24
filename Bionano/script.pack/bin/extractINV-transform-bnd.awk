#!/bin/awk -f
$1 !~ /^#/ && $10 ~ /^inversion/ {
start1=$7;end1=$8;
id1=$1;
getline line2;split(line2, a);
start2=a[7];end2=a[8];
id2=a[1];
if (end2!=-1) {
	if (end1 > end2) {
		bnd_1_left=start2;
		bnd_1_right=start1;
		ori_1="5to5";
		bnd_2_left=end2;
		bnd_2_right=end1;
		ori_2="3to3";
		
	} else {
		bnd_1_left=start1;
                bnd_1_right=start2;
		ori_1="5to5"; 
                bnd_2_left=end1;
                bnd_2_right=end2;
		ori_2="3to3";
	}
printf ("chr%s\t%d:%d\t%d:%d\t%s\n", $3, bnd_1_left, bnd_1_right, bnd_2_left, bnd_2_right, "Paired_inversion");
} else {
	if (start2 > end1) {
		bnd_1_left=start1;
		bnd_1_right=start2;
		ori1="5to5"
printf ("chr%s\t%d:%d\tNA\t%s\n", $3, bnd_1_left, bnd_1_right,  "Partial_inversion");

	} else {
		bnd_1_left=start2;
		bnd_1_right=start1;
		ori="3to3"
printf("chr%s\tNA\t%d:%d\t%s\n",  $3, bnd_1_left, bnd_1_right, "Partial_inversion");

	}
}
#a[$3":"int(bnd-1-left)"-"int(bnd_1_right)]="chr"$3"\t"int(bnd_1_left"\t"int(bnd_1_right);
#a[$3":"int(bnd-1-left)"-"int(bnd_2_right)]="chr"$3"\t"int(bnd_2_left"\t"int(bnd_1_right);
}
#END {
#	for (i in a) print a[i];
#}
