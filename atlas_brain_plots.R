# clean workspace
rm(list=ls())
# libraries and custom functions
# library(ggsegHO) # to visualize regions
library(ggseg)
library(ggsegExtra)
library(ggsegHO) # to visualize regions
# fix ggsegHO labelling issue
## https://github.com/ggseg/ggesgHO/issues/1
swap_labels<-function(d,l1,l2){
  # get indices and geometry
  w1<-which(d$label==l1)
  r1<-d$region[w1] %>% unique()
  l1<-d$label[w1] %>% unique()
  h1<-d$hemi[w1] %>% unique()
  
  w2<-which(d$label==l2)
  r2<-d$region[w2] %>% unique()
  l2<-d$label[w2] %>% unique()
  h2<-d$hemi[w2] %>% unique()
  
  
  # swap
  d$region[w1]<-r2
  d$label[w1]<-l2
  d$hemi[w1]<-h2
  
  d$region[w2]<-r1
  d$label[w2]<-l1
  d$hemi[w2]<-h1
  
  return(d)
  
}
# check labels that have been switched
labels2swap<-rbind(
  cbind("lh_Juxtapositional.Lobule.Cortex..formerly.Supplementary.Motor.Cortex.","lh_Subcallosal.Cortex"),
  cbind("lh_Lingual.Gyrus","lh_Parahippocampal.Gyrus..posterior.division"),
  cbind("lh_Angular.Gyrus","lh_Supramarginal.Gyrus..posterior.division"),
  cbind("rh_Inferior.Temporal.Gyrus..anterior.division","rh_Inferior.Temporal.Gyrus..posterior.division"),
  cbind("rh_Frontal.Operculum.Cortex","rh_Occipital.Fusiform.Gyrus" )
) %>% data.frame()

for (i in 1:NROW(labels2swap)){
  hoCort$data<-swap_labels(d=hoCort$data,l1=labels2swap$X1[i],l2=labels2swap$X2[i])
}
rm(i)
# to check
# tmp2<-hoCort$data %>% filter(region %in% hoCort$data$region[w])
# tmp2 %>%
#   mutate(region2plot=region) %>% 
#   brain_join(hoCort) %>%
#   ggplot() + geom_sf(aes(fill=region2plot)) +
#   theme_void() +
#   theme(legend.position="bottom",legend.text=element_text(size=16)) +
#   NULL
# rm(tmp2)
# rm(w)
#
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(RColorBrewer)
library(MetBrewer)
library(dplyr)
library(tidyr)
library(DT)

#---------------------------------------------------
options(stringsAsFactors = FALSE)
root="UKB_BIG40"
project="cerebellum_UKB"
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
primary_dir=paste(dir,"resources/datasets/GWAS_sumstats/downloaded_data/UKB/","BIG40","/",sep="")
#---------------------------------------------------
options(stringsAsFactors = FALSE)
root="UKB_BIG40"
project="cerebellum_UKB"
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
primary_dir=paste(dir,"resources/datasets/GWAS_sumstats/downloaded_data/UKB/","BIG40","/",sep="")

out_dir=paste(dir,project,"data/atlas_figures/",sep="/")
if(!dir.exists(out_dir)){dir.create(out_dir)}
# set working dir
setwd(out_dir)
#---------------------------------------------------
# some functions
source(paste(dir,"general_scripts/helper_functions.R",sep=""))
source(paste(dir,"general_scripts/genotyping/ldsc/helper_functions_ldsc.R",sep=""))
#---------------------------------------------------
# read 
## phenotype information
idps<-read.csv(paste(primary_dir,"IDPs_summary.csv",sep=""),header=TRUE)
# %>%   filter(!is.na(Cerebellum)|!is.na(Subcortical))

# define order of regions:
cer_regions<-c("TOTAL","CORTEX","WHITE-MATTER",
               "I-IV","V","VI","CRUS I","CRUS II",
               "VIIB","VIIIA","VIIIB","IX","X")
#---------------------------------------------------
# cerebellum aseg
# associate color to region
target_labels<- ggseg(atlas="aseg")$data %>%
  dplyr::select(label,region) %>% distinct() %>%
  mutate(matchLabel=tolower(label) %>%
           gsub("\\.|-","_",.) %>%
           gsub("__","_",.) %>%
           gsub("_division|_part","",.)
  )
# Cerebellar measures from atlas 'aseg'
cer<-target_labels %>% 
  mutate(cerebellum=if_else(grepl("cereb",label),
                            as.character(region) %>% 
                              gsub("right-|left-","",.) %>% 
                              gsub("-"," ",.) %>% sapply(.,simpleCap),
                            "")) %>% distinct()
cer$cerebellum[cer$cerebellum==""]<-NA
cer_plot <- cer %>% brain_join(aseg) %>%
  ggplot() +
  geom_sf(aes(fill=cerebellum),color="grey65",size=0.5) +
  # scale_fill_brain("aseg",na.value="transparent") +
  scale_fill_manual(values = c("darkblue","lightblue"),na.value="transparent",
                    labels=c(unique(cer$cerebellum[!is.na(cer$cerebellum)]),"")) +
  theme_void() +
  theme(legend.position="right",legend.text=element_text(size=6)) +
  guides(fill=guide_legend(ncol=1,title="",
                           override.aes = list(color="white"))) +
  # theme_darkbrain() +
  NULL


ggsave(cer_plot,file="cerebellum_aseg_regions.png",width=4,height=2)
rm(cer,cer_plot)

#---------------------------------------------------
# associate color to region
target_labels<- ggseg(atlas="hoCort")$data %>% 
  dplyr::select(label,region) %>% distinct() %>%
  mutate(matchLabel=tolower(label) %>%
           gsub("\\.|-","_",.) %>%
           gsub("__","_",.) %>%
           gsub("_division|_part","",.)
  )
#' Cortical regions included in analysis (language-related):
cort_network<-target_labels %>% 
  mutate(langNet=if_else(grepl("cereb",label),
                            as.character(region) %>% 
                              gsub("right-|left-","",.) %>% 
                              gsub("-"," ",.) %>% sapply(.,simpleCap),
                            "")) %>% distinct()
cort_network$langNet[cort_network$langNet==""]<-NA
# define order of cortical measures, by lobe...
order_cortical<-sort(unique(cort_network$region))
order_cortical<-order_cortical[c(1:2, # IFG
                                 10:11, # SUPRAMARGINAL
                                 7, # PRECUNEUS
                                 8:9, # STG
                                 3:5, # MTG
                                 12:14,6 # Fusif. Gyrus, Temp and Occ
)]



cort_network_plot <- cort_network %>%
  # select(label,region.p2,region.cortical) %>% distinct() %>% 
  brain_join(hoCort) %>%
  mutate(region.cortical=factor(region,levels=order_cortical)) %>%
  ggplot() +
  geom_sf(aes(fill=region),color="lightgrey",size=0.5) +
  # scale_fill_brain() +
  # scale_fill_viridis_d(option = "inferno",na.value="transparent") +
  scale_fill_met_d(name="Signac") +
  # scale_fill_manual(na.value="transparent") +
  theme_void() +
  theme(legend.position="bottom",legend.text=element_text(size=12)) +
  guides(fill=guide_legend(ncol=3,title="")) +
  # labs(title="Cortical volumes") +
  # theme_darkbrain() +
  NULL

