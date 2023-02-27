#!/bin/bash
# PRS
#$ -N ldsc_rg
#$ -cwd
#$ -q short.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m as
#--------
# Run genetic correlation between two traits using LDSC
#--------

#--------
# get arguments
#--------
p1=${1}
p2=${2}
name=${p1}_${p2}
#--------

#--------
# load modules and define tools/directories with resources
source activate ldsc
PATH=/export/home/acarrion/acarrion/projects/resources/github/ldsc/:${PATH}
resource_dir=/export/home/acarrion/acarrion/projects/resources/
## info on LD score formats/ estimation
# https://github.com/bulik/ldsc/wiki/LD-File-Formats
# https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
ldscores_dir=${resource_dir}/LDscores/
#--------
# project specific parameters, point to data location and working directory
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/UKB_BIG40/
cd ${working_dir}
#--------

if [ ! -f ${working_dir}/${name}_rg.log ]
 then
 if [ -f ${working_dir}/${p1}.sumstats.gz ] && [ -f ${working_dir}/${p2}.sumstats.gz ] 
 then
   echo 'Run LDSC to calculate genetic correlation between: '${p1}' and '${p2}
   ldsc.py \
   --rg ${p1}.sumstats.gz,${p2}.sumstats.gz \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${name}_rg
  else
   echo Please generate ${p1}.sumstats.gz and ${p2}.sumstats.gz first.
  fi
  else
  echo File ${name}_rg.log already exists so... rg has probably been computed already.
fi
