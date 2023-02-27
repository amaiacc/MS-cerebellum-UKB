#!/bin/bash
# phenoSpD: to compute number of inpedendent phenotypes based on GWAS summary statistics
#----------------------------------------------------------------------
phenospd_dir=/export/home/acarrion/acarrion/projects/resources/github/PhenoSpD/script/
#----------------------------------------------------------------------
project=cerebellum_UKB
scripts_dir=/export/home/acarrion/acarrion/projects/${project}/scripts/PhenoSpD/

out_dir=/export/home/acarrion/acarrion/projects/${project}/data/phenoSpD/
#----------------------------------------------------------------------
mkdir -p ${out_dir}/script/
cd ${out_dir}
# copy scripts here... since relative paths are given from the phenospd scripts
cp ${phenospd_dir}/* script/
#----------------------------------------------------------------------
# create input file
## will use the LDSC input (mungestats format) to combine phenotypes in a single file, as required for PhenoSpD

# Rscript --verbose ${scripts_dir}/input4phenoSpD.R
#----------------------------------------------------------------------
# clean

# select 32 phenos (excluding totalvol by Chambers2022):
awk '{print $1,$67,$68}' input4phenoSpD_${project}.table  | head # confirm these are the cols to remove

sed 's/allele_0/snp allele_0/g' input4phenoSpD_${project}.table | awk '{$68="";print}' | awk '{$68="";print}' | sed 's/snp allele_0/allele_0/g' > input4phenoSpD_${project}_32phenos.table
grep -v NA input4phenoSpD_${project}_32phenos.table > input4phenoSpD_${project}_32phenos_clean.table

# 33 phenos (but excluding SNPs with NA values)
grep -v NA input4phenoSpD_${project}.table > input4phenoSpD_${project}_33phenos_clean.table
#----------------------------------------------------------------------
# run phenoSpD

## for 33 phenos
Rscript ./script/phenospd_ed.r --sumstats input4phenoSpD_${project}_33phenos_clean.table --out phenoSpD_${project}_33phenos

## for 32 phenos (more SNPs)
Rscript ./script/phenospd_ed.r --sumstats input4phenoSpD_${project}_32phenos_clean.table --out phenoSpD_${project}_32phenos
