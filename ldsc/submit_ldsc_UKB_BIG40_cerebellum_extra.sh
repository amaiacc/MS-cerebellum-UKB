#!/bin/bash
# additional global rg runs required for local genetic correlations (sample overlap estimation from gcov_int)
#---------------------------------------------------------------------- 
# From UKB BIG40 study
## https://open.win.ox.ac.uk/ukbiobank/big40/
#---------------------------------------------------------------------- 
analysis_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/
primary_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/UKB_BIG40/
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/UKB_BIG40/
working_dir2=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/
ukb_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/UKB/BIG40/
out_dir=/export/home/acarrion/acarrion/projects/cerebellum_UKB/data/ldsc/
#----------------------------------------------------------------------
mkdir -p ${working_dir}
#------------------------------------
# 
suffix=disorders_cog
for base_pheno in 0102 0104
 do
  if [ -f ${working_dir}${base_pheno}.sumstats.gz ]
  then
  if [ ! -f ${working_dir}/${base_pheno}_${suffix}_rg.log ] && [ ! -f ${working_dir}/logs/${base_pheno}_${suffix}_rg.log ]
  then
   echo Compute rg between ${base_pheno} and traits of interest using LDSC
   qsub -v phenos2test='Grove_ASD-AND-PGC3_SCZ_primary-AND-Lee2018_CP',p1=${base_pheno},suffix=${suffix} ${analysis_dir}ldsc_rg_cognitive.sh
  fi
  fi
 done

 
#------------------------------------
# create summary tables - for all phenos
#------------------------------------
# move logs to folder
mkdir -p ${working_dir}logs
mv ${working_dir}*log ${working_dir}/logs/

#------------------------------------
# double check log files that were wrong
#------------------------------------
cd ${working_dir}/logs/
grep -A11 p2 *disorders_rg.log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > summary_disorders_rg_ldsc.table
grep -A11 p2 *disorders_cog_rg.log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > summary_disorders_cog_rg_ldsc.table

