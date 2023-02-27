#!/bin/bash
# bed2LDscore
#$ -N bed2LDscore
#$ -o /export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/sge/
#$ -e /export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/sge/
#$ -t 1-22
#$ -cwd
#$ -q short.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m as
#--------------------------
# create annotation files from bed files
# following: https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial#step-1-creating-an-annot-file
# An annot file typically consists of CHR, BP, SNP, and CM columns, followed by one column per annotation, with the value of the annotation for each SNP (0/1 for binary categories, arbitrary numbers for continuous annotations). 
# Make sure you have the same SNPs in the same order as the .bim file used for the computation of LD scores
#--------------------------
cat_name=$1
bedfile=$2
chr=$SGE_TASK_ID
# cat_name=Fetal_brain_GSE63648
# bedfile=Fetal_brain_GSE63648_12Opcw_Hu_gain.bed
#--------------------------
ldsc_dir=/export/home/acarrion/acarrion/projects/resources/github/ldsc/
ldscores_dir=/export/home/acarrion/acarrion/projects/resources/LDscores/
bed_dir=/export/home/acarrion/acarrion/projects/resources/genomic_annotations/bedfiles/
# out dir and name prefix for output files
out_dir=/export/home/acarrion/acarrion/projects/resources/genomic_annotations/LDscores/${cat_name}/
f2=$(echo $bedfile | sed "s/.bed//g")
#-------------------------------------------
source activate ldsc
#-------------------------------------------
# sort bedfile
mkdir -p ${out_dir}
cd ${bed_dir}
if [ ! -f sorted.${bedfile} ]
 then
 echo 'Sort ${bedfile}'
 sed 's/ /\t/g' ${bedfile} | grep -v -e start -e stop > tmp.${bedfile}
 bedtools sort -chrThenSizeA -i tmp.${bedfile} > sorted.${bedfile}
 fi

if [ ! -f ${out_dir}/${f2}.${chr}.annot.gz ]
then
 ## Step 1: Creating an annot file
 python ${ldsc_dir}make_annot.py \
     --bed-file sorted.${bedfile} \
     --bimfile ${ldscores_dir}/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
     --annot-file ${out_dir}/${f2}.${chr}.annot.gz
 ## Step 2: Computing LD scores with an annot file.
 python ${ldsc_dir}ldsc.py\
     --l2\
     --bfile ${ldscores_dir}/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}\
     --ld-wind-cm 1\
     --annot ${out_dir}/${f2}.${chr}.annot.gz\
     --print-snps ${ldscores_dir}/listHM3.txt \
     --thin-annot\
     --out ${out_dir}/${f2}.${chr}\
     --print-snps ${ldscores_dir}/listHM3.txt # need to use this list of snps in order to match the baseline model
#     --print-snps ${ldscores_dir}/hapmap3_snps/hm.${chr}.snp 
fi
