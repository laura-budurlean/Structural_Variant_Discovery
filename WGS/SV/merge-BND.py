#!/bin/env python

# compare two BND find to find out BNDs that are within a given distance.
import sys

file_1=sys.argv[1]
file_2=sys.argv[2]
buffer=int(sys.argv[3])

f_1 = open(file_1, 'r')
f_2 = open(file_2, 'r')

res = []
lines_2 = f_2.readlines() ## store file_2 into the memory

for line_1 in f_1:
  pos_1 = line_1.strip()
  chr_1 = pos_1.split('\t')[0:5]
  chr_1[1] = int(chr_1[1]);chr_1[3] = int(chr_1[3])
  for line_2 in lines_2:
    pos_2 = line_2.strip()
    chr_2 = pos_2.split('\t')[0:5]
    chr_2[1] = int(chr_2[1]);chr_2[3] = int(chr_2[3])
## translocation orientation must be the same to be included as final
    if (chr_1[4] != chr_2[4]):
      continue
    if (chr_1[0] != chr_2[0]) or ((chr_1[0] == chr_2[0]) and ((chr_1[1] - buffer > chr_2[1]) or chr_1[1] + buffer < chr_2[1])):
      continue
    if (chr_1[0] == chr_2[0]) and (chr_1[1] + buffer >= chr_2[1]) and (chr_1[1] - buffer <= chr_2[1]):
      if (chr_1[2] == chr_2[2]) and (chr_1[3] + buffer >= chr_2[3]) and (chr_1[3] - buffer <= chr_2[3]):
        res.append(pos_1 + '\t' + pos_2)
      
for records in res:
  print(records)

f_1.close
f_2.close
