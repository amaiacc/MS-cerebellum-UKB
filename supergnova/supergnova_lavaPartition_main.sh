#!/bin/bash
project=cerebellum_fusiform
# Run supergnova, using genome partition as defined by LAVA
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
## create genome partiton as defined by LAVA
# bash ${analysis_dir}/supergenova_partition_lavaPart.sh


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

## create a list of pairs, or iterate over two loops (without repeating combinations of phenotypes!)
rm ${project}_pheno.pairs.txt
echo ${rois2} >> ${project}_pheno.pairs.txt
echo "Chambers2022_TotalCerebellarVolume 0143
Chambers2022_TotalCerebellarVolume 0146" >> ${project}_pheno.pairs.txt

for P1 in ${rois2} Chambers2022_TotalCerebellarVolume
 do 
 for P2 in ${rois3} Grove_ASD Lee2018_CP PGC3_SCZ_primary
 do 
  echo ${P1} ${P2} >> ${project}_pheno.pairs.txt
 done
done
# ensure that no pairs are repetead, just in case
sort ${project}_pheno.pairs.txt | uniq > ${project}_unique_pheno.pairs.txt


#----------------------------------------------------------------------
# run local bivar analysis

N_PAIRS=$(wc -l ${project}_unique_pheno.pairs.txt | awk '{print $1}')

for I in $(seq 1 $N_PAIRS); do
 # extract phenotype IDs
 P1=$(awk 'NR=='$I' {print $1}' ${project}_unique_pheno.pairs.txt)
 P2=$(awk 'NR=='$I' {print $2}' ${project}_unique_pheno.pairs.txt)
 if [ ${P1} != ${P2} ]; then
    # define sumstats files
    p1=$(find ${ldsc_dir} -name ${P1}.sumstats.gz | head -n 1 ) # sumstats p1
    p2=$(find ${ldsc_dir} -name ${P2}.sumstats.gz | head -n 1 ) # sumstats p2
    out_file=${working_dir}${P2}_${P1}.supergnova_lavaPartition_results.txt
    # check if output files for these phenotypes already exist
    if [ ! -f ${working_dir}${P2}_${P1}.supergnova_lavaPartition_results.txt ] && [ ! -f ${working_dir}${P1}_${P2}.supergnova_lavaPartition_results.txt ]
    then
        # if not, SUBMIT JOB! 
        echo "Run" ${p1} ${p2}
        qsub -pe by_node 9 -v p1=${p1},p2=${p2},out=${working_dir}${P1}_${P2}.supergnova_lavaPartition_results.txt ${analysis_dir}supergnova_lavaPartition_run.sh 
    else
        echo "*  output files already exist for: " ${P1} ${P2}
        #ls ${P2}_${P1}.supergnova_lavaPartition_results.txt

    fi
 fi
 done
