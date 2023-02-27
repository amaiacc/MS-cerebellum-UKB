#!/bin/bash
#----------------------------------------------------------------------
# workflow for running local genetic correlations using LAVA
## https://ctg.cncr.nl/software/lava
## https://github.com/josefin-werme/LAVA
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# ## preliminary steps ##
#----------------------------------------------------------------------
# define working directoy
analysis_dir=/export/home/acarrion/acarrion/projects/cerebellum_UKB/scripts/lava/
dos2unix ${analysis_dir}/*

#----------------------------------------------------------------------
project=cerebellum_fusiform
source ${analysis_dir}/${project}_settings.sh
#----------------------------------------------------------------------
# locus files -> partitions from lava and from supergnova
## copy them to working_dir
cp ${lava_dir}/support_data/*.locfile ${working_dir}
# create partition file using the bed files that are used by supergnova (to make results comparable)
#bash ${analysis_dir}/lava_partition_supergnovaPart.sh
cp /export/home/acarrion/acarrion/projects/general_scripts/genotyping/lava/*locfile ${working_dir}
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Create sample overlap file (optional input for lava --> get sample overlap from ldsc)
cd ${working_dir}
## get names of sumstats, as defined in the info file
awk '{print $4}' ${info_file} | tail -n +2 | sed 's/.sumstats.gz//g' > ${project}_sumstats2read.txt
# get sample overlap from LDSC rg --> gcov_int
Rscript ${analysis_dir}/get_sample_overlap_from_ldsc.R ${project}_sumstats2read.txt
#----------------------------------------------------------------------
# Run LAVA
# Submit separate jobs to cluster
## 1- UNIV: run h2 univ --> in parallel for each phenotype
## 2- BIVAR: then read univ results for each pair --> filter loci based on univ h2 --> run bivar only on those loci
## 3- PCOR: read bivar and univ for target and conditional phenotype --> filter loci that are non-zero for bivar + univ cond --> run conditional
#----------------------------------------------------------------------
# UNIV
# run univ analysis for each phenotype separatelly
# to avoid computing these multiple times when phenotypes are repeteated across pairs for rg analysis
#----------------------------------------------------------------------
while read P1; do
    # check if output files for these phenotypes already exist
    if [[ -f ${out_dir}/${P1}.univ ]]; then
        echo "*  output files already exist for: "${P1}
    else
        # if not, SUBMIT JOB! 
        echo "Run" ${P1}
        qsub -pe by_node 9 ${scripts_dir}lava_univ_sge.job ${analysis_dir}/${project}_settings.sh ${P1}
    fi
 done < ${project}_sumstats2read.txt

#----------------------------------------------------------------------
# BIVAR
#----------------------------------------------------------------------
## create a list of pairs, or iterate over two loops (without repeating combinations of phenotypes!)
rm ${project}_pheno.pairs.txt

echo "0143 0146
Chambers2022_TotalCerebellarVolume 0143
Chambers2022_TotalCerebellarVolume 0146" >> ${project}_pheno.pairs.txt
for P1 in 0143 0146 Chambers2022_TotalCerebellarVolume
 do 
 for P2 in 0102 0104 Grove_ASD Lee2018_CP PGC3_SCZ_primary
 do 
  echo ${P1} ${P2} >> ${project}_pheno.pairs.txt
 done
done

# ensure that no pairs are repeatad, just in case
sort ${project}_pheno.pairs.txt | uniq > ${project}_unique_pheno.pairs.txt
# iterate over all unique phenotype pairs
N_PAIRS=$(wc -l ${project}_unique_pheno.pairs.txt | awk '{print $1}')

for I in $(seq 1 $N_PAIRS); do
    # extract phenotype IDs
    P1=$(awk 'NR=='$I' {print $1}' ${project}_unique_pheno.pairs.txt)
    P2=$(awk 'NR=='$I' {print $2}' ${project}_unique_pheno.pairs.txt)

    # check if output files for these phenotypes already exist
    if [[ -f ${out_dir}/${P1}.${P2}.bivar ]]; then
        echo "*  output files already exist for: "${P1} ${P2}
    else
        # if not, SUBMIT JOB! 
        echo "Run" ${P1} ${P2}
        qsub -pe by_node 9 ${scripts_dir}lava_rg_sge.job ${analysis_dir}/${project}_settings.sh ${P1} ${P2}
        
        #(note: this is written for slurm, adapt the line below as necessary)
        ## sbatch -t $WT:00 -p $NODE -J $P1.$P2.lava -o slurm.$P1-$P2.%A.out lava_rg.job $SCRIPTS/settings.sh $P1 $P2
    fi
 done
#----------------------------------------------------------------------
# CONDITIONAL analyses

# qsub ${scripts_dir}lava_pcor_sge.job ${analysis_dir}/${project}_settings.sh 0143 0146 0102
# qsub ${scripts_dir}lava_pcor_sge.job ${analysis_dir}/${project}_settings.sh 0143 0146 0104


for P1 in 0143 0146 Chambers2022_TotalCerebellarVolume
 do
 for P2 in Grove_ASD Lee2018_CP PGC3_SCZ_primary
  do
   for P3 in 0143 0146 Chambers2022_TotalCerebellarVolume #  
    do
    if [ "${P1}" != "${P3}" ]; then
     # check if output files for these phenotypes already exist
    if [[ -f ${out_dir}/${P1}.${P2}_condBy.${P3}.pcor ]]; then
        echo "*  output files already exist for: "${P1} ${P2}
    else
        # if not, SUBMIT JOB! 
        echo "Run" ${P1} ${P2} "conditioned by" ${P3}
        qsub ${scripts_dir}lava_pcor_sge.job ${analysis_dir}/${project}_settings.sh ${P1} ${P2} ${P3}
    fi;
    fi
   done
  done
 done


#----------------------------------------------------------------------
# clean intermediate files
rm ${project}_sumstats2read.txt
