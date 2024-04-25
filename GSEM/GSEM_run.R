# GSEM
#----------------------------------------------------------
# running genomic SEM analysis for cerebellar substructures
# following: https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects
#----------------------------------------------------------
rm(list=ls())

#----------------------------------------------------------
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
working_dir=paste0(dir,"resources/datasets/GWAS_sumstats/")
out_dir=paste0(dir,"cerebellum_UKB/data/GSEM/")
scripts_dir=paste0(dir,"cerebellum_UKB/scripts/GSEM/")
source(paste0(scripts_dir,"semPlotModel_GSEM.R"))
#----------------------------------------------------------
# could read directly the munged data..
setwd(paste0(working_dir,"ldsc/UKB_BIG40/"))
#----------------------------------------------------------
# define possible subset for analyses:
# - all: 33 cerebellar measures
# - VR: FAST cerebellar volumes from right cereb. and vermis
# - VL: FAST cerebellar volumes from left cereb. and vermis
for (analysis_run in c("VL","VR")){
  idps<-read.csv(paste0(working_dir,"downloaded_data/UKB/BIG40/IDPs_summary.csv")) %>% 
    filter(Cerebellum==1&(parcellation=="FAST"|parcellation=="aseg")&lobe=="SUBCORTICAL"&measure=="VOLUME")
  
  source(paste0(scripts_dir,"genomicSEM.R"))
  gc()
}
rm(analysis_run)
#----------------------------------------------------------
