#!/bin/bash
# Heritability and genetic correlation analyses of selected imaging derived phenotypes ( UKB (BIG40) )
#---------------------------------------------------------------------- 
# From UKB BIG40 study
## https://open.win.ox.ac.uk/ukbiobank/big40/
#---------------------------------------------------------------------- 
analysis_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/UKB_BIG40/
ukb_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/UKB/BIG40/
out_dir=/export/home/acarrion/acarrion/projects/cerebellum_UKB/data/ldsc/partitioned/
#----------------------------------------------------------------------
ldscores_dir=/export/home/acarrion/acarrion/projects/resources/LDscores/
ldscores_dir2=/export/home/acarrion/acarrion/projects/resources/genomic_annotations/LDscores/
#----------------------------------------------------------------------
mkdir -p ${analysis_dir}sge
mkdir -p ${out_dir}
## requirements
# Select regions to format:  select_phenos2download.R
# Run $analysis_dir/QC_sumstats/submit_publicGWAS_QC_UKB_BIG40.sh

pheno_file=${ukb_dir}IDPs.csv
dos2unix ${pheno_file}
# define rois as cerebellum - parcellation
rois=$(grep -e Cerebellum -e cerebellum ${ukb_dir}IDPs_summary.csv | grep -v -e INTENSITY -e intensity -e pheno | awk -F "," '{print $26}' | sed 's/"//g')
rois=$(echo Chambers2022_TotalCerebellarVolume $rois)
rois2=$( grep -e Cerebellum -e cerebellum ${ukb_dir}IDPs_summary.csv | grep -v -e INTENSITY -e intensity | grep -e crus_I -e crus_II -e _X -e '_VI"' -e Cortex -e White-Matter | awk -F "," '{print $26}' | sed 's/"//g' )

#------------------------------------
## for selected phenotypes
# run partitioned h2 analysis

# Adult_enhancers_Cerebellum Adult_promoters_Cerebellum NeanderthalDepleted_Chen2020 NeanderthalMatchSNPs_sprime_Browning2018 ncHAR SelectiveSweep Fetal_brain_GSE63648_7pcw
annots2test="Fetal_brain_GSE63648_Hu_gain Adult_enhancers_promoters_Cerebellum NeanderthalDepleted_Chen2020 NeanderthalMatchSNPs_sprime_Browning2018 ncHAR SelectiveSweep"

#qsub ${analysis_dir}ldsc_h2_partitioned.sh \
   ${p} ${s_file} \
   ${ldscores_dir}baseline_v1.2/baseline. \
   baseline \
   ${out_dir}/\
   cerebellum_baseline

# using h2 to get enrichment: two-tailed test
for p in ${rois}
 do
 s_file=${working_dir}/${p}.sumstats.gz

 for annot in ${annots2test}
  do
  name=$(grep ${annot} ${ldscores_dir2}/*ldcts |grep -v evol_annots_all.ldcts | awk -F":" '{print $2}' | awk '{print $1}' | sort | uniq)
  ref_ld=$(grep ${annot} ${ldscores_dir2}/*ldcts |grep -v evol_annots_all.ldcts | awk -F":" '{print $2}' | awk '{print $2}' | sort | uniq )
  if [ ! -f ${out_dir}/${p}_annots_cerebellum_${name}.log ] && [ ! -f ${out_dir}/logs/${p}_annots_cerebellum_${name}.log ]
  then
  echo ${p1}_cerebellum_${name}
  qsub ${analysis_dir}ldsc_h2_partitioned.sh\
    ${p} ${s_file}\
    ${ldscores_dir}/baselineLD_v2.2/baselineLD.\
    ${ref_ld}\
    ${out_dir}/\
    cerebellum_${name}
   sleep 2s
  fi
 done
done

# combine all results to process later
cd ${out_dir}
grep -e Category -e L2_1 *cerebellum*results | sort | uniq > ${out_dir}evol_annots_cerebellum.table

#------------------------------------
# to run as cell-type s-ldsc
for p in ${rois2}
 do
 s_file=${working_dir}/${p}.sumstats.gz
 # using h2-cts
 #qsub ${analysis_dir}ldsc_h2cts_partitioned.sh ${p} ${s_file} ${ldscores_dir}baselineLD_v2.2/baselineLD. ${ldscores_dir2}/evol_annots_set1.ldcts ${out_dir}/ cerebellum_evol_annots
 #qsub ${analysis_dir}ldsc_h2cts_partitioned.sh ${p} ${s_file} ${ldscores_dir}baselineLD_v2.2/baselineLD. ${ldscores_dir2}/evol_annots_adultHGEpfc.ldcts ${out_dir}/ cerebellum_evol_annots_HGEpfc
 #qsub ${analysis_dir}ldsc_h2cts_partitioned.sh ${p} ${s_file} ${ldscores_dir}baselineLD_v2.2/baselineLD. ${ldscores_dir2}/evol_annots_fetalHGE.ldcts ${out_dir}/ cerebellum_evol_annots_fetalHGE
done


#------------------------------------

#------------------------------------
# move logs to folder
mkdir -p logs results
mv *log ./logs/
mv *results* ./results/

for f in $(ls | grep log | grep -v annots)
 do
 mv ${f} logs/
 done


#------------------------------------