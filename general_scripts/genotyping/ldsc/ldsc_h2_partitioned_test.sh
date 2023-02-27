#!/bin/bash
## test run partitioned heritability analysis, following: https://github.com/bulik/ldsc/wiki/Partitioned-Heritability
scripts_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/

# load modules and define tools/directories with resources
source activate ldsc
PATH=/export/home/acarrion/acarrion/projects/resources/github/ldsc/:${PATH}
resource_dir=/export/home/acarrion/acarrion/projects/resources/
## info on LD score formats/ estimation
# https://github.com/bulik/ldsc/wiki/LD-File-Formats
# https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
ldscores_dir=${resource_dir}/LDscores/
ldscores_dir2=/export/home/acarrion/acarrion/projects/resources/genomic_annotations/LDscores/
## info on LD score formats/ estimation
# https://github.com/bulik/ldsc/wiki/LD-File-Formats
# https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
#--------------------------
working_dir=${resource_dir}datasets/GWAS_sumstats/ldsc/partitioned/
mkdir -p ${working_dir}
cd ${scripts_dir}/sge/

#--------------------------
# partitioned heritability
#--------------------------
## cell-type specific analysis
#cts_name=Cahoy
#cts_name=Multi_tissue_chromatin
# edit file ${cts_name}.ldcts to add full paths -> ${cts_name}_ed.ldcts
# sed "s|${cts_name}_1000Gv3_ldscores|${ldscores_dir}${cts_name}_1000Gv3_ldscores|g" ${ldscores_dir}${cts_name}.ldcts > ${ldscores_dir}${cts_name}_ed.ldcts

for p in Okbay2022_EA4 Lee2018_EA3 Grove_ASD PASS_Autism_Hujoel2019 Ripke2014_SCZ # GIANT_BMI_Speliotes2010 Grove_ASD PASS_Autism_Hujoel2019 Ripke2014_SCZ # PASS_Autism_Grove2019_Hujoel2019
 do
 s_file=${working_dir}/../${p}.sumstats.gz
 # using h2-cts - e.g. cell-type Cahoy, as ref.
 #qsub ${scripts_dir}ldsc_h2cts_partitioned.sh ${p} ${s_file} ${ldscores_dir}baselineLD_v2.2/baselineLD. ${ldscores_dir}${cts_name}_ed.ldcts ${working_dir}/ ${cts_name}
 # using h2-cts
 qsub ${scripts_dir}ldsc_h2cts_partitioned.sh ${p} ${s_file} ${ldscores_dir}baselineLD_v2.2/baselineLD. ${ldscores_dir2}/evol_annots_set1.ldcts ${working_dir}/ evol_annots
 qsub ${scripts_dir}ldsc_h2cts_partitioned.sh ${p} ${s_file} ${ldscores_dir}baselineLD_v2.2/baselineLD. ${ldscores_dir2}/evol_annots_adultHGEpfc.ldcts ${working_dir}/ evol_annots_HGEpfc
 qsub ${scripts_dir}ldsc_h2cts_partitioned.sh ${p} ${s_file} ${ldscores_dir}baselineLD_v2.2/baselineLD. ${ldscores_dir2}/evol_annots_fetalHGE.ldcts ${working_dir}/ evol_annots_fetalHGE
 done

# using h2 to get enrichment: two-tailed test
for p in Okbay2022_EA4 Lee2018_EA3 Grove_ASD PASS_Autism_Hujoel2019 Ripke2014_SCZ # GIANT_BMI_Speliotes2010 Grove_ASD PASS_Autism_Hujoel2019 Ripke2014_SCZ # PASS_Autism_Grove2019_Hujoel2019
 do
 s_file=${working_dir}/../${p}.sumstats.gz
 for annot in Fetal_brain_GSE63648_7pcw # ncHAR SelectiveSweep 
  do
  name=$(grep ${annot} ${ldscores_dir2}/*ldcts |grep -v evol_annots_all.ldcts | sort | uniq | awk -F":" '{print $2}' | awk '{print $1}')
  ref_ld=$(grep ${annot} ${ldscores_dir2}/*ldcts |grep -v evol_annots_all.ldcts | sort | uniq | awk -F":" '{print $2}' | awk '{print $2}')
  qsub ${scripts_dir}ldsc_h2_partitioned.sh \
   ${p} ${s_file} \
   ${ldscores_dir}/baselineLD_v2.2/baselineLD. \
   ${ref_ld}\
   ${working_dir}/\
   ${name}
 done
done
 
#--------------------------
# Eising et al. 2021 (bioRxiv) Supp.Table 17
## but note that I'm testing word reading here, while the GenLang study reports the enrichment results for the MTAG sumstats
#--------------------------
working_dir2=/export/home/acarrion/acarrion/projects/coeduca/data/working_data/GenLang/projects/GWASMA_freeze_v1/data/ldsc/
p=WR_Z_combined_STERR_GCOFF
s_file=${working_dir2}${p}.sumstats.gz

# using h2-cts - cell-type Cahoy, as ref.
qsub ${scripts_dir}ldsc_h2cts_partitioned.sh ${p} ${s_file} ${ldscores_dir}baselineLD_v2.2/baselineLD. ${ldscores_dir}${cts_name}_ed.ldcts ${working_dir}/ ${cts_name}
# using h2-cts
qsub ${scripts_dir}ldsc_h2cts_partitioned.sh ${p} ${s_file} ${ldscores_dir}baselineLD_v2.2/baselineLD. ${ldscores_dir2}/evol_annots_set1.ldcts ${working_dir}/ evol_annots
# using h2
qsub ${scripts_dir}ldsc_h2_partitioned.sh \
  ${p} ${s_file} \
  ${ldscores_dir}/baselineLD_v2.2/baselineLD. \
  ${ldscores_dir2}/SelectiveSweep/SelectiveSweep_Peyregne2017.,${ldscores_dir2}/ncHAR/ncHAR_Capra2013. \
  ${working_dir}/\
  SelectiveSweep_nCHAR
  
qsub ${scripts_dir}ldsc_h2_partitioned.sh \
  ${p} ${s_file} \
  ${ldscores_dir}/baselineLD_v2.2/baselineLD. \
  ${ldscores_dir2}/NeanderthalMatchSNPs/NeanderthalMatchSNPs_sprime_Browning2018. \
  ${working_dir}/\
  NeanderthalMatchSNPs

qsub ${scripts_dir}ldsc_h2_partitioned.sh \
  ${p} ${s_file} \
  ${ldscores_dir}/baselineLD_v2.2/baselineLD. \
  ${ldscores_dir2}/NeanderthalDepleted/NeanderthalDepleted_Vernot2016. \
  ${working_dir}/\
  NeanderthalDepleted_Vernot2016
#--------------------------

#--------------------------
# Tilot eta  al. 2020
#--------------------------
p=Mean_Full_SurfArea_noGC
p=Mean_Full_SurfArea
p=Mean_parsorbitalis_surfavg
s_file=${working_dir}/../ENIGMA3/${p}.sumstats.gz
cd ${scripts_dir}/sge/
# using h2-cts - cell-type Cahoy, as ref.
#qsub ${scripts_dir}ldsc_h2cts_partitioned.sh ${p} ${s_file} ${ldscores_dir}baselineLD_v2.2/baselineLD. ${ldscores_dir}${cts_name}_ed.ldcts ${cts_name}
# using h2-cts
qsub ${scripts_dir}ldsc_h2cts_partitioned.sh ${p} ${s_file} ${ldscores_dir}baselineLD_v2.2/baselineLD. ${ldscores_dir2}/evol_annots_set1.ldcts ${working_dir}/ evol_annots
qsub ${scripts_dir}ldsc_h2cts_partitioned.sh ${p} ${s_file} ${ldscores_dir}baselineLD_v2.2/baselineLD. ${ldscores_dir2}/evol_annots_adultHGEpfc.ldcts ${working_dir}/ evol_annots_HGEpfc
qsub ${scripts_dir}ldsc_h2cts_partitioned.sh ${p} ${s_file} ${ldscores_dir}baselineLD_v2.2/baselineLD. ${ldscores_dir2}/evol_annots_fetalHGE.ldcts ${working_dir}/ evol_annots_fetalHGE

 # using h2 to get enrichment
qsub ${scripts_dir}ldsc_h2_partitioned.sh \
  ${p} ${s_file} \
  ${ldscores_dir}/baselineLD_v2.2/baselineLD. \
  ${ldscores_dir2}/SelectiveSweep/SelectiveSweep_Peyregne2017.,${ldscores_dir2}/ncHAR/ncHAR_Capra2013. \
  ${working_dir}/\
  SelectiveSweep_nCHAR

qsub ${scripts_dir}ldsc_h2_partitioned.sh \
  ${p} ${s_file} \
  ${ldscores_dir}/baselineLD_v2.2/baselineLD. \
  ${ldscores_dir2}/Fetal_brain_GSE63648/Fetal_brain_GSE63648_7pcw_Hu_gain.,${ldscores_dir2}/Fetal_Brain/Fetal_Brain_Male.all.sorted.merged. \
  ${working_dir}/\
  Fetal_7pcw
