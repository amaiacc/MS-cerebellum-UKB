#!/bin/bash
# PRS
#$ -N supergnova
#$ -o /export/home/acarrion/acarrion/projects/general_scripts/genotyping/supergnova/sge/
#$ -e /export/home/acarrion/acarrion/projects/general_scripts/genotyping/supergnova/sge/
#$ -cwd
#$ -q long.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m as
#--------
PATH=/export/home/acarrion/acarrion/projects/resources/github/SUPERGNOVA/:${PATH}
module load python/python3.6
supergnova_dir=/export/home/acarrion/acarrion/projects/resources/github/SUPERGNOVA/

# to define as variables when submitting job
# p1: sumstats for trait1 ## files are in the standard format that ldsc understands.
# n1: sample size for trait1
# p2: sumstats for trait2
# n2: sample size for trait2
# out: output name
# thr: thread to be used when running

#--N1 ${n1} \
#--N2 ${n2} \

python3 ${supergnova_dir}supergnova.py ${p1} ${p2} \
--bfile ${supergnova_dir}/data/bfiles/eur_chr@_SNPmaf5 \
--partition ${supergnova_dir}/data/partition/eur_chr@.bed \
--out ${out}

