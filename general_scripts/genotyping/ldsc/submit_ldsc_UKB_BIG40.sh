#!/bin/bash
# Heritability and genetic correlation analyses of selected imaging derived phenotypes ( UKB (BIG40) )
#---------------------------------------------------------------------- 
# From UKB BIG40 study
## https://open.win.ox.ac.uk/ukbiobank/big40/
#---------------------------------------------------------------------- 
analysis_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/
primary_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/UKB_BIG40/
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/UKB_BIG40/
#----------------------------------------------------------------------
mkdir -p ${working_dir}

## requirements
# Select regions to format:  select_phenos2download.R
# Run $analysis_dir/QC_sumstats/submit_publicGWAS_QC_UKB_BIG40.sh

cd ${primary_dir}
files2format=$(ls *QC.txt.gz)

pheno_file=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/UKB/BIG40/IDPs.csv
dos2unix ${pheno_file}
#------------------------------------
# run

mkdir -p ${analysis_dir}sge
cd ${analysis_dir}sge
#-----------------------------------
# LDSC format: munge sumstats 
for base_file in ${files2format}
 do
  base_pheno=$(echo ${base_file} | sed 's/.QC.txt.gz//g')
  if [ ! -f ${working_dir}${base_pheno}.sumstats.gz ]
  then
   echo Submit LDSC formatting ${base_pheno}
   qsub ${analysis_dir}ldsc_format_UKB_BIG40.sh ${base_pheno} ${pheno_file}
  fi
 done
# heritability
for base_file in ${files2format}
 do
  base_pheno=$(echo ${base_file} | sed 's/.QC.txt.gz//g')
  if [ -f ${working_dir}${base_pheno}.sumstats.gz ]
  then
  if [ ! -f ${working_dir}/${base_pheno}_h2.log ]
  then
   echo Compute h2 using LDSC ${base_file}
   qsub ${analysis_dir}ldsc_h2.sh ${base_pheno}
  fi
  fi
 done
# Run genetic correlations for L and R
lat_file=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/UKB/BIG40/IDPs_lat_summary.txt
dos2unix ${lat_file}
phensLR=$(awk '{print $3}' ${lat_file} | grep -v -e "NA$")
for name in ${phensLR}
do
if [ ! -f ${working_dir}${name}_LR_rg.log ]
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
parcs=$(awk '{print $1,$2}' ${lat_file} | sort | uniq | sed 's/ /_/g')
for name in ${parcs}
do
 echo ${name}
 p=$(echo $name | awk -F "_" '{print $1}') # parcellation
 m=$(echo $name | awk -F "_" '{print $2}') # measure

 for col in {4..6}
  do
   phens=($(grep ${p} ${lat_file} | grep ${m} | awk -v n=${col} '{print $n}' | grep -v 'NA' | sort | uniq)) # generate array, the content will be accessed by the indices $i and $j
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
      if [ ! -f ${working_dir}${p1}_${p2}_rg.log ]
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
 phens=($(grep ${p} ${lat_file} | grep ${m} | grep -e cer -e Cer | awk '{print $4,$5,$6}' | sed 's/ /\n/g' | grep -v NA | sort | uniq )) # generate array, the content will be accessed by the indices $i and $j
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
      if [ ! -f ${working_dir}${p1}_${p2}_rg.log ]
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

#------------------------------------
# create summary tables
#------------------------------------
cd ${working_dir}
## h2 estimates
# columns to get:
##Total Observed scale h2: 0.0103 (0.0027)
##Lambda GC: 1.0466
##Mean Chi^2: 1.0469
##Intercept: 1.0093 (0.006)
##Ratio: 0.1984 (0.1285)
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
grep -A11 p2 *rg.log | grep -v -e LR -e cognitive | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > summary_rg_ldsc.table
grep -A11 p2 *cognitive_rg.log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > summary_cognitive_rg_ldsc.table


#------------------------------------
# move logs to folder
#mkdir -p logs
#mv *log ./logs/
#------------------------------------