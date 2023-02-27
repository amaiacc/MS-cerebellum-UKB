#!/bin/bash

ukb_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/UKB/BIG40/
lat_file=${ukb_dir}IDPs_lat_summary.txt

#-----------------------------------------
# define all phenotypes
#-----------------------------------------
# define rois as cerebellum - parcellation
rois=$(grep -e cerebellum -e Cerebellum ${ukb_dir}IDPs_summary.csv | grep -v -e intensity -e pheno | awk -F "," '{print $26}' | sed 's/"//g')

## FAST
name=FAST_VOLUME
echo ${name}
p=$(echo $name | awk -F "_" '{print $1}') # parcellation
m=$(echo $name | awk -F "_" '{print $2}') # measure
m2=$(echo $m | tr '[:upper:]' '[:lower:]')
## phens2: list of language related cortical volumes (FAST parcellation)
phens2=$( grep -v cerebell ${lat_file} | grep ${p} | grep -e ${m} -e ${m2} | awk '{print $4,$5,$6}' | sed 's/ /\n/g' | grep -v NA | sort | uniq )
# subcortical vols
# phens3: FAST
phens3=$( grep -e pallidum -e putamen -e striatum -e accumb -e caudate -e thalamus -e amygdala -e hippocampus -e brain_stem ${lat_file} | grep -e ${p} | grep ${m} | awk '{print $4,$5,$6}' | sed 's/ /\n/g' | grep -v NA | sort | uniq )
# phens3: aseg
name=aseg_VOLUME
p=$(echo $name | awk -F "_" '{print $1}') # parcellation
m=$(echo $name | awk -F "_" '{print $2}') # measure
m2=$(echo $m | tr '[:upper:]' '[:lower:]')
phens4=$( grep -e Pallidum -e Putamen -e Striatum -e Accumb -e Caudate -e Thalamus -e Amygdala -e Hippocampus -e Brain-Stem ${lat_file} | \
 grep -e ${p} | grep -e ${m} -e ${m2} | awk '{print $4,$5,$6}' | sed 's/ /\n/g' | grep -v NA | sort | uniq )
#-----------------------------------------
out_dir=/export/home/acarrion/acarrion/projects/cerebellum_UKB/data/
cd ${out_dir}
# get all phenos from IDPs_summary.csv file
## quick and dirty ... but well enough
head -n 1 ${ukb_dir}IDPs_summary.csv > UKBphenos4TableS1.csv
for p in  ${rois} ${phens2} ${phens3} ${phens4}
 do
 grep '"'${p}'"' ${ukb_dir}IDPs_summary.csv >> UKBphenos4TableS1.csv
 # awk -F ","  '($26 == $p) { print $0 }' ${ukb_dir}IDPs_summary.csv
 done
