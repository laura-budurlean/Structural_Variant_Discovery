# Variables
NTASKS = 1 
TIME = 24:00:00 
MEM = 16G # specify memory
PARTITION = standard # set this to the partition you are using
NUM_SAMPLES = 10 # set this to the number of samples you are running

# Targets
.PHONY: all help 1-part1 remove_recurrent 1-part2 2-compare 3-stratification filter merge

all: 1-part1 remove_recurrent 1-part2 2-compare 3-stratification filter merge

help:
	@echo "Makefile targets:"
	@echo "  all                - Run all parts of the workflow in sequence."
	@echo "  1-part1            - Run the 1-Optical-Mapping-SV-part1.sh script."
	@echo "  remove_recurrent   - Run the Bionano-remove-recurrent.sh script."
	@echo "  1-part2            - Run the 1-Optical-Mapping-SV-part2.sh script."
	@echo "  2-compare          - Run the 2-compare.sh script."
	@echo "  3-stratification   - Run the 3-stratification.sh script."
	@echo "  filter             - Run the filter-sv.sh script."
	@echo "  merge              - Run the merge-sv.sh script."

1-part1:
	sbatch --ntasks=$(NTASKS) -t $(TIME) --mem=$(MEM) --partition=$(PARTITION) --job-name=1-Optical-Mapping-part1 \
	-o /scripts/Bionano/script.pack/out_err_files/%J.%N.out -e /scripts/Bionano/script.pack/out_err_files/%J.%N.err \
	--array=1-$(NUM_SAMPLES) 1-Optical-Mapping-SV-part1.sh

remove_recurrent: 1-part1
	sbatch --ntasks=$(NTASKS) -t $(TIME) --mem=$(MEM) --partition=$(PARTITION) --job-name=remove-recurrent \
	-o /scripts/Bionano/script.pack/out_err_files/%J.%N.out -e /scripts/Bionano/script.pack/out_err_files/%J.%N.err \
	Bionano-remove-recurrent.sh

1-part2: remove_recurrent
	sbatch --ntasks=$(NTASKS) -t $(TIME) --mem=$(MEM) --partition=$(PARTITION) --job-name=1-Optical-Mapping-part2 \
	-o /scripts/Bionano/script.pack/out_err_files/%J.%N.out -e /scripts/Bionano/script.pack/out_err_files/%J.%N.err \
	1-Optical-Mapping-SV-part2.sh

2-compare: 1-part2
	sbatch --ntasks=$(NTASKS) -t $(TIME) --mem=$(MEM) --partition=$(PARTITION) --job-name=2-compare \
	-o /scripts/Bionano/script.pack/out_err_files/%J.%N.out -e /scripts/Bionano/script.pack/out_err_files/%J.%N.err \
	--array=1-$(NUM_SAMPLES) 2-compare.sh

3-stratification: 2-compare
	sbatch --ntasks=$(NTASKS) -t $(TIME) --mem=$(MEM) --partition=$(PARTITION) --job-name=3-stratifications \
	-o /scripts/Bionano/script.pack/out_err_files/%J.%N.out -e /scripts/Bionano/script.pack/out_err_files/%J.%N.err \
	--array=1-$(NUM_SAMPLES) 3-stratifications.sh

filter: 3-stratification
	sbatch --ntasks=$(NTASKS) -t $(TIME) --mem=$(MEM) --partition=$(PARTITION) --job-name=filter-sv \
	-o /scripts/Bionano/script.pack/out_err_files/%J.%N.out -e /scripts/Bionano/script.pack/out_err_files/%J.%N.err \
	scripts/WGS/SV/filter-sv.sh

merge: filter
	sbatch --ntasks=$(NTASKS) -t $(TIME) --mem=$(MEM) --partition=$(PARTITION) --job-name=merge-sv \
	-o /scripts/Bionano/script.pack/out_err_files/%J.%N.out -e /scripts/Bionano/script.pack/out_err_files/%J.%N.err \
	scripts/WGS/SV/merge-sv.sh