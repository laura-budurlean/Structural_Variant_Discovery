#!/usr/bin/env Rscript


# accept command line arguments
sys.argv = commandArgs(trailingOnly=TRUE)

if (length(sys.argv) != 2) {
	stop("Usage: Rscript --vanilla remove-recurrent.R input_file output_path\n");
}

# the irys recurrent translocation list
recurrent.trans <- read.csv("/scripts/Bionano/script.pack/Bionano-recurrent-translocation-1Mb.txt", header = FALSE, stringsAsFactors = FALSE, skip = 1, sep = "\t")

# the recently generated translocation bed file that we will be using
our.trans <- read.table(sys.argv[1], comment.char = '#', stringsAsFactors = FALSE, header = FALSE, sep = '\t')

# the comparison threshold
threshold <- 5e5

# number of rows in the irys recurrent translocation list
recurrent.row <- nrow(recurrent.trans)

# number of rows in our bed file of translocations
our.row <- nrow(our.trans)

# generating a summary vector that lists whether or not a value is recurrent. Default: TRUE
recurrent.totals <- rep(T, our.row)

# iterate through each row and test it against the recurrent translocation list
for (i in 1:our.row) {
	for (j in 1:recurrent.row) {
		
		# check if orientation, chrA, chrB, and pos1 & pos2 match 
		if ((our.trans[i, 'V5'] == recurrent.trans[j, 'V5']) &&
			(our.trans[i, 'V1'] == recurrent.trans[j, 'V1']) &&
		(abs(our.trans[i, 'V2'] -  recurrent.trans[j, 'V2']) <= threshold) &&
		 	(our.trans[i, 'V3'] == recurrent.trans[j, 'V3']) &&
		(abs(our.trans[i, 'V4'] -  recurrent.trans[j, 'V4']) <= threshold)) {
		
			# if they do, they are recurrent & change their value in our master list
			recurrent.totals[i] <- FALSE
			
			# exit this loop level
			break 
		}
	}
}

# split the input file path by the / and put it in a variable
extract_filename <- strsplit(sys.argv[1], '/')[[1]]

# grab the very last bit of the path, which is the previous file's name and steal it
filename <- extract_filename[length(extract_filename)]

# write the output to a table

write.table(our.trans[recurrent.totals, ], paste0(sys.argv[2], "/", ".bed-non-recurrent-translocation-from-raw-unfiltered-list.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
