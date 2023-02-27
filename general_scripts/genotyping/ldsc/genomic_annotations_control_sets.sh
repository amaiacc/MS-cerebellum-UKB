# define control sets for the human gained enhancers/promoters S-LDSC runs
scripts_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/
working_dir=/export/home/acarrion/acarrion/projects/resources/genomic_annotations/
# load modules and define tools/directories with resources
source activate ldsc
PATH=/export/home/acarrion/acarrion/projects/resources/github/ldsc/:${PATH}
resource_dir=/export/home/acarrion/acarrion/projects/resources/
## info on LD score formats/ estimation
# https://github.com/bulik/ldsc/wiki/LD-File-Formats
# https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
ldscores_dir=${resource_dir}/LDscores/
#--------------------------------------------------------------------------------
cd ${working_dir}
# "Multi_tissue_chromatin" (includes both Roadmap and ENTEX data), from Finucane et al. 2018
## get names for brain "epigenome names"
rm *all*.bed
brain_names=$(grep -e "^Brain" -e "Fetal_Brain" ${ldscores_dir}Multi_tissue_chromatin.ldcts | awk -F"__" '{print $1}' | sort | uniq)
for b in ${brain_names}
 do
 grep -e ${b} ${ldscores_dir}Multi_tissue_chromatin.ldcts | awk '{print $2}' | awk -F"," '{print $1}' | sed 's/.$/.bed/g' | sed "s|^|${ldscores_dir}|g" > ${b}.bed2merge.txt
 while read line
  do
   cat ${line} >> ${b}.all.bed
  done < ${b}.bed2merge.txt
  # bedtools merge requires that you presort your data by chromosome and then by start position 
  sort -k1,1 -k2,2n ${b}.all.bed > ${b}.all.sorted.bed
  # merge intervals if overlapping
  bedtools merge -i ${b}.all.sorted.bed > ${working_dir}/bedfiles/${b}.all.sorted.merged.bed
 done

## EID: epigenome name
## E073: brain dorsolateral prefrontal cortex
## E081: fetal brain male
## E082: fetal brain female
