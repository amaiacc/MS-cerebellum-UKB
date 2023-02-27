#!/bin/bash
#$ -N ldsc_h2_partitioned
#$ -o /export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/sge/
#$ -e /export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/sge/
#$ -cwd
#$ -q long.q
#$ -S /bin/bash
#$ -M acarrion@bcbl.eu
#$ -m as
#--------
# Run heritability enrichment for a given annotation using ldsc h2 command
#--------

#--------
# get arguments
#--------
p1=${1} # phenotype name, pointing at summary statistics file
sumstats_file=${2} # sumstats file, path to file
baseline=${3} # define baseline to be used, provide absolute paths e.g. ${ldscores_dir}baselineLD_v2.2/baselineLD.
annots2test=${4} # annotations to be tested, provide absolute paths e.g. ${ldscores_dir2}SelectiveSweep/SelectiveSweep_Peyregne2017.
out_dir=${5} # output directory
suffix=${6} # suffix to append to the output files
#--------

#--------
# load modules and define tools/directories with resources
source activate ldsc
PATH=/export/home/acarrion/acarrion/projects/resources/github/ldsc/:${PATH}
resource_dir=/export/home/acarrion/acarrion/projects/resources/
## info on LD score formats/ estimation
# https://github.com/bulik/ldsc/wiki/LD-File-Formats
# https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
ldscores_dir=${resource_dir}/LDscores/
#--------
# project specific parameters, point to data location and working directory
#ldscores_dir2=/export/home/acarrion/acarrion/projects/resources/genomic_annotations/LDscores/
## create variable name including the concatenation of all the possible annotaitons that may be within annots2test
#annots2test_paths=$(echo $annots2test | sed "s|,|,${ldscores_dir2}|g") # this could be used if relative paths were given
#--------


if [ ${annots2test} == 'baseline' ]; then allannots2test=${baseline}; else allannots2test=${baseline},${annots2test}; fi

if [ ! -f ${out_dir}/${p1}_annots_${suffix}.log ]
 then
 if [ -f ${sumstats_file} ]
 then
   echo "Running S-LDSC to calculate h2 enrichment for ${suffix} annotations in trait ${p1}"
   # using h2, per annot that we may want to check
    ldsc.py \
     --h2 ${sumstats_file} \
     --ref-ld-chr ${allannots2test}\
     --w-ld-chr ${ldscores_dir}/weights_hm3_no_hla/weights.\
     --overlap-annot \
     --frqfile-chr ${ldscores_dir}/1000G_Phase3_frq/1000G.EUR.QC.\
     --out ${out_dir}/${p1}_annots_${suffix}
  else
   echo Please generate ${sumstats_file} first.
  fi
  else
  echo File ${p1}_annots_${suffix}.log already exists so... this analysis has probably been ran already.
fi
