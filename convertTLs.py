#!/usr/bin/env/python

# usage:
# python scriptname.py <inputfile> <outputfile>
# 		mv <outputfile> <inputfile>


import re
import sys


outfile = open(sys.argv[2], "w")


with open(sys.argv[1], "r") as infile:
	for line in infile:	# iterate through each line in the file
		col = line.split() # split each line into the columns
		replacement = ""	
		if len(col) > 4:	# to make sure we don't run into issues with NA 
			orientation = col[4] # the orientation is in the 5th column (vim starts from 0)

			# change the orientation
			if orientation == "+/+":
				replacement = "3to5"
			elif orientation == "+/-":
				replacement = "3to3"
			elif orientation == "-/+":
				replacement = "5to5"
			elif orientation == "-/-":
				replacement = "5to3"
	
			# this writes out the first 4 columns tab delimited, 
			# then inserts the replacement orientation, & finally
			# the rest of the columns afterwards
			outfile.write("\t".join(col[:4]) + "\t" + replacement + "\t" + "\t".join(col[5:]) + "\n")
		

		else:	# if it's an NA 
			outfile.write(line)


outfile.close()
