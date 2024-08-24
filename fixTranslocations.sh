#!/bin/bash

script= /scripts/Bionano/script.pack
OutDir= /Bionano_RareVariantData/scripts_output

cd $OutDir

for sample in #list samples
do

python /scripts/Bionano/script.pack/convertTLs.py $OutDir/translocation.${sample}.bed $OutDir/fixedTLs.${sample}.bed 

done
