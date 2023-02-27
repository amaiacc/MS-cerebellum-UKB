#!/bin/bash
# PRS
#$ -N ldsc_rg_cognitive
#$ -o /export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/sge/
#$ -e /export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/sge/
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
#p1=${1}
#suffix=${2}
#--------
name=${p1}_${suffix}

echo ${name}
echo ${phenos2test}
#------------------------------------
# GenCor with GWAS summary statistics
#------------------------------------ 
# define phenotypes --> as variables when submitting the job using qsub
#p2=Lee2018_EA3
#p3=Lee2018_CP
#p4=Savage2018_IQ
#p4=Grove_ASD
#p5=DDyslexia_Gialluisi2020

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
# define variable to test rg --> all sumstats
var4rg=$(echo ${phenos2test} | sed "s|-AND-|.sumstats.gz,../|g" )
var4rg=$(echo ${p1}.sumstats.gz,../${var4rg}.sumstats.gz)
#var4rg=$(echo ${p1}.sumstats.gz,${var4rg}) # to edit!
#
echo ${var4rg}
#--------

if [ ! -f ${working_dir}/${name}_rg.log ]
 then
 if [ -f ${working_dir}/${p1}.sumstats.gz ]
 then
   echo 'Run LDSC to calculate genetic correlation between: '${p1}' and cognitive traits of interest.'
   ldsc.py \
   --rg ${var4rg} \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${name}_rg
  else
   echo Please generate ${p1}.sumstats.gz first.
  fi
  else
  echo File ${name}_rg.log already exists so... rg has probably been computed already.
fi
