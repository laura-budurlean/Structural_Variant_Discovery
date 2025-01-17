# Structural Variant Detection Pipeline Combining Optical Genome Mapping and Whole Genome Sequencing
by Laura Budurlean

A combination of optical genome mapping with Bionano and whole genome sequencing short-read data. This pipeline was created to help integrate structural variant calling from these two technologies. These results were published in:

Budurlean L, Tukaramrao DB, Zhang L, Dovat S, Broach J. Integrating optical genome mapping and whole genome sequencing in somatic structural variant detection. *Journal of Personalized Medicine. 2024*;14(3):291.

# SV Calling Pipeline Summary:
![Fig1](https://github.com/user-attachments/assets/be599605-4406-4032-b91a-4fa5d652ba7a)

# Example Pipeline Output of an IKZF1 Deletion in Leukemia Patient:
![Fig5](https://github.com/user-attachments/assets/5dbdc675-29a0-4e02-a1ad-992555c61599)
Bionano optical genome maps are aligned to precise coordinates using WGS reads to confirm a heterozygous IKZF1 deletion.

# How To Run:

## Option 1: Makefile run

You can use the makefile: ```make -f Makefile.mak help``` to run the help file, and ```make -f Makefile.mak all``` to run all the steps at once. You will still need to specify some file names/directories inside the 1- 2- 3- main .sh scripts first depending on your unique sample naming conventions.

## Option 2: Manual run

You can run each of the scripts yourself in the order specified below.

### 1. Run SV Calling Pipeline on WGS Files

You can run any SV calling pipeline you wish for your project on your WGS files. You need to obtain a VCF file with your SV calls from WGS. For our data, we opted to use LUMPY and DELLY. We use the SpeedSeq pipeline for running the SVs, which already incorporates LUMPY. We run DELLY separately. A separate repository for doing SV calling with SpeedSeq is available here: [SV-calling-with-SpeedSeq](https://github.com/laura-budurlean/SV-calling-with-SpeedSeq). For copy number alterations, we use FreeCNV.

### 2. Obtain VCF and .SMAP Files from Bionano Data

From the Bionano data, obtain the VCF and .SMAP files for each sample from the pipeline you ran and give each sample its own folder with the respective VCF and .SMAP files. We used the Bionano Rare Variant Pipeline. You can either use filtered or unfiltered files. The scripts here have their own built-in filtering parameters that you can adjust, similar to the Bionano ones available online, but you may also pre-filter Bionano files before downloading them if you wish.

### 3. Run the Scripts

Now we can run the scripts. Note: If you set up your directories differently, make sure to change the file paths in the scripts.

### 4. Load Required Modules

Load any required modules including: Python, Bedtools, R.

### 5. Adjust Headers, Parameters, and Paths

Adjust headers, parameters, and paths as required for your own project and run the following scripts in order:

1. `1-Optical-Mapping-SV-part1.sh`
2. `Bionano-remove-recurrent.sh`
3. `1-Optical-Mapping-SV-part2.sh`

### 6. Run Filtering Script

Run `filter-sv.sh` in the `/.../scripts/WGS/SV/` scripts folder.

### 7. Run Merging Script

Run `merge-sv.sh` in the `/.../scripts/WGS/SV/` scripts folder.

### 8. Notes

- Ensure your FreeCNV `.bam_CNVs` files have `chr` prefixes and no trailing tabs.
- If you encounter Python pandas errors, ensure your Python version is higher than 3.7. We used Python 3.7.1.