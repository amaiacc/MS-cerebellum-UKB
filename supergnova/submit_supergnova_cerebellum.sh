#!/bin/bash
# example run for supergnova
#---------------------------------------------------------------------- 
analysis_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/supergnova/
#---------------------------------------------------------------------- 
ukb_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/UKB/BIG40/
ldsc_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/
#----------------------------------------------------------------------
ldscores_dir=/export/home/acarrion/acarrion/projects/resources/LDscores/
ldscores_dir2=/export/home/acarrion/acarrion/projects/resources/genomic_annotations/LDscores/
#----------------------------------------------------------------------
working_dir=/export/home/acarrion/acarrion/projects/cerebellum_UKB/data/supergnova/
mkdir -p ${working_dir} ${analysis_dir}/sge/
#----------------------------------------------------------------------
pheno_file=${ukb_dir}IDPs_summary.csv
lat_file=${ukb_dir}IDPs_lat_summary.txt
dos2unix ${pheno_file} ${lat_file}
#----------------------------------------------------------------------
# define rois as cerebellum - parcellation
rois2=$(grep -e Cerebellum -e cerebellum ${pheno_file} | grep -v -e INTENSITY -e intensity | grep -e "Crus I " -e "VI " | grep "Right" | awk -F "," '{print $26}' | sed 's/"//g' )
# define temp_occ phenos
name=FAST_VOLUME
echo ${name}
p=$(echo $name | awk -F "_" '{print $1}') # parcellation
m=$(echo $name | awk -F "_" '{print $2}') # measure
m2=$(echo $m | tr '[:upper:]' '[:lower:]')
rois3=$( grep -v cerebell ${pheno_file} | grep ${p} | grep -e ${m} -e ${m2} | grep -e occ_fusif | grep -e "Left" | awk -F "," '{print $26}' | sed 's/"//g' )

#----------------------------------------------------------------------
cd ${working_dir}
#----------------------------------------------------------------------
# between both cereb. measures: VI and crus I (right)
p1_name=0143



#----------------------------------------------------------------------
# between cereb. measures and cortical measures
for p1_name in ${rois2}
do
 p1=${ldsc_dir}/UKB_BIG40/${p1_name}.sumstats.gz
 for p2_name in ${rois2}
 do
 if [ ${p1_name} != ${p2_name} ] 
  then
  if [ ! -f ${working_dir}${p2_name}_${p1_name}.supergnova_lavaPartition_results.txt ] && [ ! -f ${working_dir}${p1_name}_${p2_name}.supergnova_lavaPartition_results.txt ]
  then
   p2=${ldsc_dir}/UKB_BIG40/${p2_name}.sumstats.gz
    echo ${p1} ${p2}
   qsub -v p1=${p1},p2=${p2},out=${working_dir}${p2_name}_${p1_name}.supergnova_lavaPartition_results.txt ${analysis_dir}supergnova_lavaPartition_run.sh
  fi; fi
 done
done



#----------------------------------------------------------------------
# with ASD 
p1=${ldsc_dir}/Grove_ASD.sumstats.gz
p1_name=ASD
# with SCZ
p1=${ldsc_dir}/PGC3_SCZ_primary.sumstats.gz
p1_name=SCZ
# with CP
p1=${ldsc_dir}/Lee2018_CP.sumstats.gz
p1_name=CP
# run
for p2_name in ${rois2}
do
 if [ ! -f ${working_dir}${p2_name}_${p1_name}.supergnova_lavaPartition_results.txt ] && [ ! -f ${working_dir}${p1_name}_${p2_name}.supergnova_lavaPartition_results.txt ]
 then
  p2=${ldsc_dir}/UKB_BIG40/${p2_name}.sumstats.gz
  echo ${p2}
  qsub -v p1=${p1},p2=${p2},out=${working_dir}${p2_name}_${p1_name}.supergnova_lavaPartition_results.txt ${analysis_dir}supergnova_run.sh
 fi
 done

#----------------------------------------------------------------------