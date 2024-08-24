#!/usr/bin/env python

from __future__ import print_function
"""
@author: dearxxj
"""

import sys
sys.path.append("/.../python/3.7.1/lib/python3.7/site-packages") #set this path as needed
import pandas as pd

if(len(sys.argv) != 4): 
	print("Usage: python compare-Bionano-wgs-tl.py Bionano.bed wgs.bed out_path")
	exit(1)

Bionano_file=sys.argv[1]
wgs_file=sys.argv[2]
out_path=sys.argv[3]
cells = [Bionano_file.split('/')[-1].split('.')[0]]

CHR_ORDER={"chr1":1, "chr2":2, "chr3":3, "chr4":4, "chr5":6, "chr6":6, 
           "chr7":7, "chr8":8, "chr9":9, "chr10":10, "chr11":11, "chr12":12,
           "chr13":13, "chr14":14, "chr15":15, "chr16":16, "chr17":17, "chr18":18, 
           "chr19":19, "chr20":20, "chr21":21, "chr22":22, "chrX":23, "chrY":24}

def swap_hic_tr(x):
    if (CHR_ORDER[x[1]] > CHR_ORDER[x[5]]) or (CHR_ORDER[x[1]] == CHR_ORDER[x[5]] and x[2] > x[6]):
        (x[1], x[5], x[2], x[6], x[3], x[7], x[4], x[8]) = (x[5], x[1], x[6], x[2], x[7], x[3], x[8], x[4])
    return x

def swap_kyr_tr(x):
    if (CHR_ORDER[x[1]] > CHR_ORDER[x[4]]) or (CHR_ORDER[x[1]] == CHR_ORDER[x[4]] and x[2] > x[5]):
        (x[1], x[4], x[2], x[5], x[3], x[6]) = (x[4], x[1], x[5], x[2], x[6], x[3])
    return x

def swap_fus_tr(x):
    if (CHR_ORDER[x[0]] > CHR_ORDER[x[2]]) or (CHR_ORDER[x[0]] == CHR_ORDER[x[2]] and x[1] > x[3]):
        (x[1], x[3], x[0], x[2]) = (x[3], x[1], x[2], x[0])
        x[4] = x[4][::-1]
    return x

def as_hic_ori(x):
    return {'5to3':'+-', '5to5':'++', '3to5':'-+', '3to3':'--',
            'N[[':'+-', 'N]]':'++', ']]N':'-+', '[[N':'--'}.get(x, None)  ## if not found, use default value, safer than dict

for cell in cells:
    bnd_irs = pd.read_table(Bionano_file, header=None)
    bnd_wgs = pd.read_table(wgs_file, header=None)
    res = pd.DataFrame(columns=['chrA', 'chrB', 'Orientation', \
                                'Bionano_pos_A', 'Bionano_pos_B', 'WGS_pos_A', \
                                'WGS_pos_B'], index=range(5000))
    
    nrow_Bionano = bnd_irs.shape[0]
    nrow_wgs = bnd_wgs.shape[0]
    
    ## fill in Bionano breakpints
    for i in range(nrow_Bionano):
        res.iloc[i, range(5)] = list(bnd_irs.iloc[i, [0, 2]]) + [as_hic_ori(bnd_irs.iloc[i, 4])]+ [bnd_irs.iloc[i, 1], bnd_irs.iloc[i, 3]]
    
    irs_wgs_cutoff = 500000

    total_row = nrow_Bionano
    ## compare and append wgs breakpoints
    total_row_Bionano = nrow_Bionano ## record total_row to where Bionano records stop
    for i in range(nrow_wgs):
        found = False
        for j in range(nrow_Bionano):
            ## consider orientation or not, 
            if res.loc[j, 'Orientation'] == as_hic_ori(bnd_wgs.iloc[i, 4]) and \
               res.loc[j, 'chrA'] == bnd_wgs.iloc[i, 0] and res.loc[j, 'chrB'] == bnd_wgs.iloc[i, 2] and \
               (res.loc[j, 'Bionano_pos_A']-irs_wgs_cutoff) <= bnd_wgs.iloc[i, 1] and (res.loc[j, 'Bionano_pos_A']+irs_wgs_cutoff) >= bnd_wgs.iloc[i, 1] and \
               (res.loc[j, 'Bionano_pos_B']-irs_wgs_cutoff) <= bnd_wgs.iloc[i, 3] and (res.loc[j, 'Bionano_pos_B']+irs_wgs_cutoff) >= bnd_wgs.iloc[i, 3]:
                   res.loc[j, ['WGS_pos_A', 'WGS_pos_B']] = list(bnd_wgs.iloc[i, [1, 3]])
                   found = True
        if not found:
            res.loc[total_row, ['chrA', 'chrB', 'Orientation', 'WGS_pos_A', 'WGS_pos_B']] = list(bnd_wgs.iloc[i, [0, 2]]) + [as_hic_ori(bnd_wgs.iloc[i, 4])] + list(bnd_wgs.iloc[i, [1, 3]])
            total_row = total_row + 1


    res = res.iloc[range(total_row), :]
    res.to_csv(out_path +'/' + cell + '.tl.2way.txt', sep = '\t', index = False, na_rep = '.')



