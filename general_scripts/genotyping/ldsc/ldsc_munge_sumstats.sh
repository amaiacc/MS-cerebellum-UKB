#!/bin/bash
#$ -N ldsc_munge
#$ -o /export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/sge/
#$ -e /export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/sge/
#$ -cwd
#$ -q long.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m as
#--------
# Munge sumstats
#--------
source activate ldsc
PATH=/export/home/acarrion/acarrion/projects/resources/github/ldsc/:${PATH}
resource_dir=/export/home/acarrion/acarrion/projects/resources/
ldscores_dir=${resource_dir}/LDscores/
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/

# arguments defined as variables when submitting job
# p_file
# N1
# p
# snp_id

# check arguments
echo ${p_file}
echo ${p}
echo ${snp_id} ${N1}
# run
munge_sumstats.py \
   --sumstats ${p_file} \
   --snp ${snp_id} \
   --N ${N1} \
   --p p_value \
   --ignore direction,Zscore \
   --out ${working_dir}/${p} \
   --a1 effect_allele --a2 other_allele \
   --merge-alleles ${ldscores_dir}/eur_w_ld_chr/w_hm3.snplist
