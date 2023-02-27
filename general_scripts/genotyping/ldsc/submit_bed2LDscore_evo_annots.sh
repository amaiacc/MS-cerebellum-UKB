# submit jobs to generate LDscores for evolutionary annotations
scripts_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/
bed_dir=/export/home/acarrion/acarrion/projects/resources/genomic_annotations/bedfiles/
#
cd ${bed_dir}

control_annots='Fetal_Brain Brain'
test_annots='Fetal_brain_GSE63648 Adult_enhancers Adult_promoters dhsHAR SelectiveSweep  ncHAR NeanderthalDepleted NeanderthalDenisovanDepleted NeanderthalMatchSNPs'
## to run

for cat_name in Adult_enhancers_promoters Fetal_brain_GSE63648_Hu_gain # ${control_annots} ${test_annots}
do
 cd ${bed_dir}
 bedfiles=$(ls ${cat_name}*bed| grep -v '^sorted')
 echo $bedfiles
 #--------------------------------
 # create annotations for each bed file within category
 for bed in ${bedfiles}
  do
  qsub ${scripts_dir}genomic_annotations_bed2LDscore_array.sh ${cat_name} ${bed}
  done
 done
#--------------------------------
# create ldcts files -> defining all annotations and paths
ldscores_dir=/export/home/acarrion/acarrion/projects/resources/LDscores/
ldscores_dir2=/export/home/acarrion/acarrion/projects/resources/genomic_annotations/LDscores/
rm ${ldscores_dir2}/evol_annots_all.ldcts
touch ${ldscores_dir2}/evol_annots_all.ldcts

for p in $(ls ${ldscores_dir2}/*/*.1.l2.ldscore.gz)
do 
 name=$(echo $p | sed "s|${ldscores_dir2}/||g" | sed 's/.1.l2.ldscore.gz//g' | awk -F"/" '{print $2}')
 echo "${name} ${p}" | sed 's/1.l2.ldscore.gz//g' | sed 's/ /\t/g' >> ${ldscores_dir2}/evol_annots_all.ldcts
 done
#--------------------------------
# create ldcts files with control set defined
grep -e HAR -e Neanderthal -e Selective ${ldscores_dir2}/evol_annots_all.ldcts > ${ldscores_dir2}/evol_annots_set0.ldcts
grep -e ncHAR -e Chen2020 -e Selective ${ldscores_dir2}/evol_annots_all.ldcts > ${ldscores_dir2}/evol_annots_set1.ldcts
# adult HGE, prefrontal cortex, also controlling for dorsolateral prefrontal cortex epigenetic markers
control=$(grep -e Dorsolateral_Prefrontal ${ldscores_dir2}/evol_annots_all.ldcts | awk '{print $2}')
grep -e PFC ${ldscores_dir2}/evol_annots_all.ldcts | awk -v c=${control} 'BEGIN{OFS=",";} {print $0,c}' > ${ldscores_dir2}/evol_annots_adultHGEpfc.ldcts
# adult HGE, cerebellum, also controlling for dorsolateral prefrontal cortex epigenetic markers
control=$(grep -e Dorsolateral_Prefrontal ${ldscores_dir2}/evol_annots_all.ldcts | awk '{print $2}')
grep -e Cerebellum ${ldscores_dir2}/evol_annots_all.ldcts | awk -v c=${control} 'BEGIN{OFS=",";} {print $0,c}' > ${ldscores_dir2}/evol_annots_adultHGEcerebellum.ldcts
# fetal
control=$(grep -e Fetal_Brain_Male ${ldscores_dir2}/evol_annots_all.ldcts | awk '{print $2}')
grep -e Fetal_brain ${ldscores_dir2}/evol_annots_all.ldcts | awk -v c=${control} 'BEGIN{OFS=",";} {print $0,c}' > ${ldscores_dir2}/evol_annots_fetalHGE.ldcts

#--------------------------------
# check that the number of SNPs in chr1 is the same as in the baseline model, for all the custom annotations
#zcat ${ldscores_dir}/baselineLD_v2.2/baselineLD.1.annot.gz | wc
#zcat ${ldscores_dir2}*/*.1.annot.gz | wc
#--------------------------------