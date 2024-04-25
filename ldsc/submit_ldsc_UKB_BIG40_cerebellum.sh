#!/bin/bash
# Heritability and genetic correlation analyses of selected imaging derived phenotypes ( UKB (BIG40) )
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

## requirements
# Select regions to format:  select_phenos2download.R
# Run $analysis_dir/QC_sumstats/submit_publicGWAS_QC_UKB_BIG40.sh

cd ${primary_dir}

pheno_file=${ukb_dir}IDPs.csv
dos2unix ${pheno_file}
# define rois as cerebellum - parcellation
rois=$(grep -e cerebellum -e Cerebellum ${ukb_dir}IDPs_summary.csv | grep -v -e intensity -e pheno | awk -F "," '{print $27}' | sed 's/"//g')
#------------------------------------
# run
mkdir -p ${analysis_dir}sge
cd ${analysis_dir}sge
#-----------------------------------
# LDSC format: munge sumstats 
for base_pheno in ${rois}
 do
  base_file=${base_pheno}.txt.gz
  if [ ! -f ${working_dir}${base_pheno}.sumstats.gz ]
  then
   echo Submit LDSC formatting ${base_pheno}
   qsub ${analysis_dir}ldsc_format_UKB_BIG40.sh ${base_pheno} ${pheno_file}
  fi
 done

#-----------------------------------
# heritability
## add total cerebellar volume, downloaded and formatted in GWASsumstats_Chambers2022.sh

# copy sumstats to working dir, so that all are placed in the same directory
if [ -f ${working_dir2}Chambers2022_TotalCerebellarVolume.sumstats.gz ]
 then
 cp ${working_dir2}Chambers2022_TotalCerebellarVolume.sumstats.gz ${working_dir}Chambers2022_TotalCerebellarVolume.sumstats.gz
fi

rois=$(echo Chambers2022_TotalCerebellarVolume $rois)
for base_pheno in ${rois}
 do
  base_file=${base_pheno}.txt.gz
  if [ -f ${working_dir}${base_pheno}.sumstats.gz ]
  then
   if [ ! -f ${working_dir}/${base_pheno}_h2.log ] && [ ! -f ${working_dir}/logs/${base_pheno}_h2.log ]
   then
    echo Compute h2 using LDSC ${base_file}
    qsub ${analysis_dir}ldsc_h2.sh ${base_pheno}
   fi
  fi
done
# Run genetic correlations for L and R
lat_file=${ukb_dir}IDPs_lat_summary.txt
dos2unix ${lat_file}
phensLR=$(grep -e cerebellum -e Cerebellum ${lat_file} | grep -v -e intensity | awk '{print $3}' | grep -v -e "NA$")
for name in ${phensLR}
do
if [ ! -f ${working_dir}${name}_LR_rg.log ] && [ ! -f ${working_dir}/logs/${name}_LR_rg.log ]
 then
 
 name2=$(echo ${name}"\t")
 p1=$(grep -P ${name2} ${lat_file} | awk '{print $5}') # use -P to grep including tab, otherwise it may give multiple matches (e.g. Crus_I and Crus_II)
 p2=$(grep -P ${name2} ${lat_file} | awk '{print $6}')
 if [ -f ${working_dir}/${p1}.sumstats.gz ] && [ -f ${working_dir}/${p2}.sumstats.gz ] 
 then
  echo 'Run LDSC to calculate genetic correlation between left and right for ' ${name}': '${p1}' and '${p2}
  qsub ${analysis_dir}ldsc_rg_LR.sh ${name} ${p1} ${p2}
 fi
fi
done

# Run genetic correlations within atlas/parcellation type (for left and right separatelly)
parcs=$(grep -e cerebellum -e Cerebellum ${lat_file} | grep -v -e intensity -e INTENSITY | awk '{print $1,$2}' | sort | uniq | sed 's/ /_/g')
for name in ${parcs}
do
 echo ${name}
 p=$(echo $name | awk -F "_" '{print $1}') # parcellation
 m=$(echo $name | awk -F "_" '{print $2}') # measure

 for col in {4..6}
  do
   phens=($(grep -e cerebellum -e Cerebellum ${lat_file} | grep ${p} | grep ${m} | awk -v n=${col} '{print $n}' | grep -v 'NA' | sort | uniq)) # generate array, the content will be accessed by the indices $i and $j
   n=${#phens[@]}
   if [ ${n} -ge 2 ]
    then 
    # Run pairwise genetic correlations
    for i in $(seq 0 $(($n-1))); do
    for j in $(seq 0 $(($n-1))); do
     if [ "$j" -gt "$i" ]; then
      # define p1 and p2 as getting indices from phens
      p1=${phens[${i}]}
      p2=${phens[${j}]}
      if [ ! -f ${working_dir}${p1}_${p2}_rg.log ] && [ ! -f ${working_dir}/logs/${p1}_${p2}_rg.log ]
       then
       echo Run rg between: ${p1} ${p2}
       qsub ${analysis_dir}ldsc_rg.sh ${p1} ${p2}
      fi
     fi
    done; done
   else 
    echo 'Skip, cannot run genetic correlations for '${n}' trait(s)'.
   fi
  done

done

## for the cerebellum, run also vermis vs left and right
name=FAST_VOLUME
 echo ${name}
 p=$(echo $name | awk -F "_" '{print $1}') # parcellation
 m=$(echo $name | awk -F "_" '{print $2}') # measure
 phens=($(grep -e cerebellum ${lat_file} | grep ${p}  | grep ${m} | grep -e cer -e Cer | awk '{print $4,$5,$6}' | sed 's/ /\n/g' | grep -v NA | sort | uniq )) # generate array, the content will be accessed by the indices $i and $j
 n=${#phens[@]}
 if [ ${n} -ge 2 ]
    then 
    # Run pairwise genetic correlations
    for i in $(seq 0 $(($n-1))); do
    for j in $(seq 0 $(($n-1))); do
     if [ "$j" -gt "$i" ]; then
      # define p1 and p2 as getting indices from phens
      p1=${phens[${i}]}
      p2=${phens[${j}]}
      if [ ! -f ${working_dir}${p1}_${p2}_rg.log ] && [ ! -f ${working_dir}/logs/${p1}_${p2}_rg.log ]
       then
       echo Run rg between: ${p1} ${p2}
       qsub ${analysis_dir}ldsc_rg.sh ${p1} ${p2}
      fi
     fi
    done; done
   else 
    echo 'Skip, cannot run genetic correlations for '${n}' trait(s)'.
 fi

#------------------------------------
# Genetic correlations of total cerebellar volume with regional cerebellar volumes
p1=Chambers2022_TotalCerebellarVolume
for p2 in ${rois}
 do
  if test "$p1" != "$p2"
  then
    if [ -f ${working_dir}${p2}.sumstats.gz ]
    then
      if [ ! -f ${working_dir}${p1}_${p2}_rg.log ] && [ ! -f ${working_dir}/logs/${p1}_${p2}_rg.log ]
       then
       echo Run rg between: ${p1} ${p2}
       qsub ${analysis_dir}ldsc_rg.sh ${p1} ${p2}
  fi; fi; fi
 done

#------------------------------------
# Genetic correlations with disorders/other traits
suffix=disorders_cog
for base_pheno in ${rois}
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

# suffix=handedness
# for base_pheno in ${rois}
 # do
  # if [ -f ${working_dir}${base_pheno}.sumstats.gz ]
  # then
  # if [ ! -f ${working_dir}/${base_pheno}_${suffix}_rg.log ] && [ ! -f ${working_dir}/logs/${base_pheno}_${suffix}_rg.log ]
  # then
   # echo Compute rg between ${base_pheno} and traits of interest using LDSC
   # qsub -v phenos2test='LeftHandedness_CuellarPartida2020-AND-Ambidextrous_CuellarPartida2020',p1=${base_pheno},suffix=${suffix} ${analysis_dir}ldsc_rg_cognitive.sh
  # fi
  # fi
 # done

#------------------------------------
# define sets of regions of interest (subcortical and selection of cortical) for rg
#------------------------------------
name=FAST_VOLUME
echo ${name}
p=$(echo $name | awk -F "_" '{print $1}') # parcellation
m=$(echo $name | awk -F "_" '{print $2}') # measure
m2=$(echo $m | tr '[:upper:]' '[:lower:]')
## phens2: list of language related cortical volumes (FAST parcellation)
phens2=$( grep -v cerebell ${lat_file} | grep ${p} | grep -e ${m} -e ${m2} | \
 grep -v -e thal -e stem -e stria -e caudate -e hippocampus -e pallidum -e putamen | \
  awk '{print $4,$5,$6}' | sed 's/ /\n/g' | grep -v NA | sort | uniq )
phens3=$( grep -e pallidum -e putamen -e striatum -e accumb -e caudate -e thalamus -e amygdala -e hippocampus -e brain_stem ${lat_file} | grep -e ${p} | grep ${m} | awk '{print $4,$5,$6}' | sed 's/ /\n/g' | grep -v NA | sort | uniq )
# aseg
name=aseg_VOLUME
p=$(echo $name | awk -F "_" '{print $1}') # parcellation
m=$(echo $name | awk -F "_" '{print $2}') # measure
m2=$(echo $m | tr '[:upper:]' '[:lower:]')
phens4=$( grep -e Pallidum -e Putamen -e Striatum -e Accumb -e Caudate -e Thalamus -e Amygdala -e Hippocampus -e Brain-Stem ${lat_file} | \
 grep -e ${p} | grep -e ${m} -e ${m2} | awk '{print $4,$5,$6}' | sed 's/ /\n/g' | grep -v NA | sort | uniq )


#------------------------------------
# subcortical volumes
#------------------------------------
# Run pairwise genetic correlations with subcortical volumes
## for all cerebellar phenotypes
for p2 in ${phens3} ${phens4}; do
 for p1 in ${rois}; do
  if [ ! -f ${working_dir}${p1}_${p2}_rg.log ] && [ ! -f ${working_dir}/logs/${p1}_${p2}_rg.log ]
   then
   echo Run rg between: ${p1} ${p2}
   qsub ${analysis_dir}ldsc_rg.sh ${p1} ${p2}
  else
   echo File already exists: ${p1}_${p2}_rg.log
  fi
 done
done

#------------------------------------
# cortical language network / whole cortical FAST
#------------------------------------
# Run pairwise genetic correlations with language related cortical volumes
## for selected cerebellar phenotypes
rois2=$( grep -e Cerebellum -e cerebellum ${ukb_dir}IDPs_summary.csv | grep -v -e intensity -e pheno | awk -F "," '{print $27}' | sed 's/"//g' )
rois2=$(echo Chambers2022_TotalCerebellarVolume $rois2)

# genetic correlations with cortical volumes related to language
## + with rest of the HO atlas, for completeness + comply with reviewers
for p2 in ${phens2}; do
 for p1 in ${rois}; do
  if [ ! -f ${working_dir}/${p1}_${p2}_rg.log ] && [ ! -f ${working_dir}/logs/${p1}_${p2}_rg.log ]
   then # **to check: this is not working properly **
   echo Run rg between: ${p1} ${p2}
   qsub ${analysis_dir}ldsc_rg.sh ${p1} ${p2}
   else
   echo File already exists: ${p1}_${p2}_rg.log
  fi 
 done
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

#------------------------------------
## h2 estimates
# columns to get:
##Total Observed scale h2: 0.0103 (0.0027)
##Lambda GC: 1.0466
##Mean Chi^2: 1.0469
##Intercept: 1.0093 (0.006)
##Ratio: 0.1984 (0.1285)
cd ${working_dir}logs
h2_file=$(ls -1 *h2.log )
paste <(echo 'File') <(echo 'h2 (se)') <(echo 'Lambda GC') <(echo 'Mean Chi^2') <(echo 'Intercept (se)') <(echo 'Ratio (se)')> summary_h2_ldsc.table
paste <(ls -1 ${h2_file}) \
      <(grep "h2:" ${h2_file} | awk '{print $5,$6}' ) \
      <(grep "Lambda GC:" ${h2_file} | awk '{print $3}') \
      <(grep "Mean Chi^2" ${h2_file} | awk '{print $3}') \
      <(grep "Intercept" ${h2_file} | awk '{print $2,$3}') \
      <(grep "Ratio" ${h2_file} | awk '{print $2,$3}') >> summary_h2_ldsc.table

## genetic correlations between left and right
grep -A11 p2 *LR_rg.log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > summary_LRrg_ldsc.table
grep -A11 p2 *rg.log | grep -v -e LR -e cognitive -e disorders | grep -v Analysis | grep -e p2 -e sumstats | \
grep -v -e 'log-$' -e 'time elapsed' -e '.py' -e '_sumstats' -e 'ps.sumstats' -e 'NA' > summary_rg_ldsc.table
grep -A11 p2 *disorders_rg.log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > summary_disorders_rg_ldsc.table
grep -A11 p2 *disorders_cog_rg.log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > summary_disorders_cog_rg_ldsc.table

#grep -A11 p2 *cognitive_rg.log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > summary_cognitive_rg_ldsc.table
#grep -A11 p2 *handedness_rg.log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > summary_handedness_rg_ldsc.table
#grep -A11 p2 *readingGenLang_rg.log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > summary_readingGenLang_rg_ldsc.table
mkdir -p ${working_dir}/output/
mv ${working_dir}logs/*table ${working_dir}/output/
## format summary and convert tables to csv
Rscript --verbose ${analysis_dir}/ldsc_parse_output_UKB_BIG40.R
#------------------------------------