#' ---
#' title: 'Genetic correlation of cerebellar structures and subcortical volumes'
#' author: Amaia Carrion Castillo
#' date: "`r Sys.Date()`"
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 5
#'     toc_float: true
#'     highlight: "textmate"
#'   pdf_document:
#'     keep_tex: true
#' ---
#' <style type="text/css">
#' .main-container {
#' max-width: 1800px;
#' margin-left: auto;
#' margin-right: auto;
#' }
#' </style>
#' 
#' ## Data
#' 
#' Summary statistics for the cerebellar subregions were obtained from the UK Biobank BIG40 server:
#' 
#'    - https://open.win.ox.ac.uk/ukbiobank/big40/
#'    - GWAS summary stats computed from N~31,000
#'
#' Summary statistics for other traits of interest were obtained from publicly available GWASes for:
#'
#'
#' * subcortical volumes
#' 
#' ## Analysis
#' 
#' **Genetic correlation** is the proportion of variance explained by common genetic factors. 
#' 
#' It was computed using GWAS summary statistics by running the LDSC tool.
#' 
#' * This measure ranges from -1 and 1. 
#' * I tested whether the genetic correlation estimate is significantly higher than 0.
#' 
#---------------------------------------------------
# clean workspace
rm(list=ls())
# libraries and custom functions
# library(ggsegHO) # to visualize regions
library(ggseg)
library(ggsegExtra)
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
ldsc_dir=paste(dir,"resources/datasets/GWAS_sumstats/ldsc/",root,"/output/",sep="")
out_dir=paste(dir,project,"data/ldsc/",sep="/")
if(!dir.exists(out_dir)){dir.create(out_dir)}
# set working dir
setwd(ldsc_dir)
#----------------------------------------------------------------------
# some functions
source(paste(dir,"general_scripts/helper_functions.R",sep=""))
source(paste(dir,"general_scripts/genotyping/ldsc/helper_functions_ldsc.R",sep=""))
#----------------------------------------------------------------------
# defined using phenoSpD
# Effective Number of Independent Variables [VeffLi] (using Equation 5 of Li and Ji 2005): 19
ntests<-19
#----------------------------------------------------------------------
# read 
## phenotype information
idps<-read.csv(paste(primary_dir,"IDPs_summary.csv",sep=""),header=TRUE)
# %>%   filter(!is.na(Cerebellum)|!is.na(Subcortical))

# define order of regions:
cer_regions<-c("TOTAL","CORTEX","WHITE-MATTER",
               "I-IV","V","VI","CRUS I","CRUS II",
               "VIIB","VIIIA","VIIIB","IX","X")
# read estimates: rg across FAST atlas measures
# and subset only cerebellum measures
#
rg1_subcort<-read.csv(paste0(ldsc_dir,"ldsc_rg_atlas_UKB_BIG40.csv")) 
rg2_subcort<-read.csv(paste0(ldsc_dir,"ldsc_rg_not_atlas_UKB_BIG40.csv"))
# fill in info for total cerebellum volume, from Chambers et al. 2022
rg3_subcort<-read.csv(paste0(ldsc_dir,"ldsc_rg_other_UKB_BIG40.csv"))
rg3_subcort <- rg3_subcort %>%
  mutate(measure.p1=if_else(p1=="Chambers2022_TotalCerebellarVolume","VOLUME",measure.p1),
         region.p1=if_else(p1=="Chambers2022_TotalCerebellarVolume","TotalCerebellarVolume",region.p1)
         ) %>%
  mutate(measure=if_else(measure.p1==measure.p2,measure.p1,""))
# filter(measure.p1=="VOLUME" & region.p1=="TotalCerebellarVolume")

# combine all
rg_subcort<-merge(rg1_subcort,rg2_subcort,all=TRUE)
rg_subcort<-merge(rg_subcort,rg3_subcort,all=TRUE)
rm(rg1_subcort,rg2_subcort,rg3_subcort)

# filter to keep only rg's of interest
rg_subcort <- rg_subcort %>% 
  filter(measure=="VOLUME") %>%
  # keep only rgs where first p is cerebellar and second is not
  filter(region.p1 %in% c("TotalCerebellarVolume",idps$region[!is.na(idps$Cerebellum)])) %>%
  filter(region.p2 %in% idps$region[!is.na(idps$Subcortical)])
# edit some columns
rg_subcort <- rg_subcort %>% 
  mutate(hemisphere=if_else(grepl("^V_",region.p1),"V",hemisphere.p1),
         region=gsub("^V_","",region.p1) %>% 
           gsub("CEREBELLUM_|CEREBELLUM-|CerebellarVolume","",.) %>% 
           gsub("_"," ",.) %>% toupper ) %>%
  mutate(hemisphere=if_else(is.na(hemisphere),"BiLat",hemisphere)) %>%
  mutate(hemisphere=factor(hemisphere,levels=c("BiLat","V","L","R")),
         type=if_else((region=="CORTEX"|region=="WHITE-MATTER"|region=="TOTAL"),"Global","Region"),
         region=factor(region,levels=cer_regions),
         region.subcortical=gsub("_|-"," ",region.p2) %>% tolower() %>% sapply(simpleCap) %>%
           gsub(" Area","",.)) %>%
  mutate(laterality=if_else(hemisphere=="V","Vermis",
                            if_else(hemisphere=="BiLat","Bilateral",
                                    if_else(hemisphere.p1==hemisphere.p2,"Ipsilateral","Contralateral")
                                    )),
         hemisphere.p2=if_else(is.na(hemisphere.p2),"NotLat",hemisphere.p2)) %>%
  mutate(hemisphere.p2=factor(hemisphere.p2,levels=c("L","R","NotLat")))

rg_subcort <- rg_subcort %>% 
  mutate(parcellation.subcortical=if_else(is.na(parcellation.p2),parcellation,parcellation.p2)) %>%
  mutate(parcellation.subcortical=factor(parcellation.subcortical))

rg_subcort <- rg_subcort %>%
  mutate(sig=if_else(parcellation.subcortical=="FAST",if_else(p<(0.05/(13*ntests)),"*",""),
                     if_else(parcellation.subcortical=="aseg",if_else(p<(0.05/(15*ntests)),"*",""),
                             "")))

# define subcortical structures as factor, define order
rg_subcort$region.subcortical<-
  factor(rg_subcort$region.subcortical,
            levels=c("Brain Stem","Amygdala","Hippocampus",
                     "Thalamus",
                     "Pallidum",
                     # Striatum
                     ## dorsal
                     "Putamen","Caudate",
                     ## ventral
                     "Accumbens",
                     "Ventral Striatum"))

# define only columns of interest
rg_subcort <- rg_subcort %>% dplyr::select(contains("parcellation"),contains("measure"),
                                          contains("hemisphere"),contains("region"),
                                          contains("rg"),p,contains("gcov"),
                                          sig,
                                          type,laterality)

# wide format to compare both segmentations (FAST and aseg)
rg_subcort_wide <- rg_subcort %>%
  mutate(region.subcortical=as.character(region.subcortical) %>%
           gsub("Accumbens|Ventral Striatum","Accumbens/Ventral Striatum",.)) %>%
  dplyr::select(parcellation.subcortical,hemisphere.p1,hemisphere.p2,hemisphere,
                region,region.subcortical,type,laterality,
                contains("rg"),contains("gcov"),p) %>%
  pivot_wider(names_from=parcellation.subcortical,values_from=c(contains("rg"),p,contains("gcov")))
rg_subcort_wide<-rg_subcort_wide %>% filter(!(is.na(rg_aseg)|is.na(rg_FAST)))

# save files
write.csv(rg_subcort,file=paste(out_dir,"/ldsc_rg_subcortical","_",root,"_",project,".csv",sep=""),row.names = FALSE)

#----------------------------------------------------------------------
# generate plots to visualize genetic correlations
rg_subcort_plot_fast <-rg_subcort %>% 
  filter(parcellation.subcortical=="FAST") %>%
  mutate(region=factor(region,levels=rev(cer_regions))) %>%
  mutate(sig_val=if_else(rg<0,rg.CI95lower*1.5,rg.CI95upper*1.5)) %>%
  ggplot(aes(x=region,y=rg,width=0.7,color=hemisphere,fill=hemisphere,shape=hemisphere.p2,alpha=sig)) +
  geom_hline(yintercept = 0,color="black",alpha=0.7,linetype="dashed") +
  geom_errorbar(aes(ymin=rg.CI95upper , ymax=rg.CI95lower,width=.5), position=position_dodge(.7)) +
  geom_point(position=position_dodge(.7), stat="identity",size=3,color="white") +
  geom_text(aes(x=region,y=sig_val,label=sig,color=hemisphere),alpha=0.8,
            fontface=2,size=6,position=position_dodge(.8),stat="identity") +
  labs(color="Cerebellar hemisphere",
       fill="Cerebellar hemisphere",
       shape="Subcortical hemisphere",
       x=bquote('Estimate ('*rho*')'),
       y="Cerebellar volume",
       title="Genetic correlation with subcortical volumes (FAST)") +
  facet_grid(type~region.subcortical, scales = "free",space="free") +
  coord_flip() +
  # scale_color_met_d(name="Morgenstern") +
  scale_fill_viridis_d() + scale_color_viridis_d() +
  scale_alpha_manual(values=c(0.3,0.9)) +
  scale_shape_manual(values=c(21,24,22)) +
  guides(shape = guide_legend(override.aes = list(color="black") )
         ) +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="grey85",size=1),
        legend.position = "bottom") +
  panel_border(size = 1) +
  NULL
rg_subcort_plot_aseg <- rg_subcort %>% 
  filter(parcellation.subcortical=="aseg") %>%
  mutate(region=factor(region,levels=rev(cer_regions))) %>%
  mutate(sig_val=if_else(rg<0,rg.CI95lower*1.5,rg.CI95upper*1.5)) %>%
  ggplot(aes(x=region,y=rg,width=0.7,color=hemisphere,fill=hemisphere,shape=hemisphere.p2,alpha=sig)) +
  geom_hline(yintercept = 0,color="black",alpha=0.7,linetype="dashed") +
  geom_errorbar(aes(ymin=rg.CI95upper , ymax=rg.CI95lower,width=.5), position=position_dodge(.7)) +
  geom_point(position=position_dodge(.7), stat="identity",size=3,color="white") +
  geom_text(aes(x=region,y=sig_val,label=sig,color=hemisphere),alpha=0.8,
            fontface=2,size=6,position=position_dodge(.8),stat="identity") +
  labs(color="Cerebellar hemisphere",
       fill="Cerebellar hemisphere",
       shape="Subcortical hemisphere",
       x=bquote('Estimate ('*rho*')'),
       y="Cerebellar volume",
       title="Genetic correlation with subcortical volumes (aseg)") +
  facet_grid(type~region.subcortical, scales = "free",space="free") +
  coord_flip() +
  # scale_color_met_d(name="Morgenstern") +
  scale_fill_viridis_d() + scale_color_viridis_d() +
  scale_alpha_manual(values=c(0.3,0.9)) +
  scale_shape_manual(values=c(21,24,22)) +
  guides(shape = guide_legend(override.aes = list(color="black") )
  ) +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="grey85",size=1),
        legend.position = "bottom") +
  panel_border(size = 1) +
  NULL

plot_grid(rg_subcort_plot_fast,rg_subcort_plot_aseg,ncol=1)

## compare rg values for aseg and FAST subcortical parcellations
rg_subcort_comparison <- rg_subcort_wide %>% 
  ggplot(aes(color=hemisphere,fill=hemisphere,shape=hemisphere.p2)) +
  geom_abline(slope=1,intercept=0,linetype="dashed") +
  geom_point(aes(x=rg_aseg,y=rg_FAST),
             alpha=0.7,size=4,color="white") +
  labs(fill="Cerebellar hemisphere",
       color="Cerebellar hemisphere",
       shape="Subcortical hemisphere",
       x=bquote('Estimate ('*rho*') Subcortical parcellation FAST'),
       y=bquote('Estimate ('*rho*') Subcortical parcellation aseg'),
       title="Genetic correlation with subcortical volumes",
       subtitle="Subcortical volume parcellation comparison: aseg and FAST") +
  facet_wrap(.~region.subcortical) +
  # scale_color_met_d(name="Morgenstern") +
  scale_fill_viridis_d() + scale_color_viridis_d() +
  scale_shape_manual(values=c(21,24,22)) +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="grey85",size=1),
        legend.position = "bottom") +
  panel_border(size = 1) +
  NULL


# save plots
ggsave(rg_subcort_plot_fast, 
       file=paste(out_dir,"/ldsc_rg_subcorticalFAST","_",root,"_",project,".png",sep=""),
       width=12,height=7)
ggsave(rg_subcort_plot_aseg, 
       file=paste(out_dir,"/ldsc_rg_subcorticalAseg","_",root,"_",project,".png",sep=""),
       width=12,height=7)
ggsave(rg_subcort_comparison , 
       file=paste(out_dir,"/ldsc_rg_subcortical_comparison","_",root,"_",project,".png",sep=""),
       width=7,height=7)

leg<-get_legend(
  rg_subcort_comparison + geom_point(aes(x=rg_aseg,y=rg_FAST,color=hemisphere)) +
    guides(fill = "none")
  )

plot_grid(rg_subcort_plot_aseg + theme(legend.position = "none"),
          rg_subcort_comparison + theme(legend.position = "none"),
          leg,
          ncol=1,rel_heights = c(1,1,0.1),
          labels=c("A","B")) %>% 
  ggsave(.,file=paste(out_dir,"/Figure_ldsc_rg_subcortical_Aseg_comparison","_",root,"_",project,".png",sep=""),
         width=14,height=12)


rg_subcort_plot_fast %>% 
  ggsave(.,file=paste(out_dir,"/Figure_ldsc_rg_subcortical_FAST","_",root,"_",project,".png",sep=""),
         width=12,height=8)

#--------------------------------------------
# clustering of the cerebellum-subcortical genetic correalations
library(corrplot)
## FAST
fast_rg<-rg_subcort_wide %>% 
  mutate(subcortical=paste0(region.subcortical," (",hemisphere.p2,")")%>% gsub(" (NA)","",.),
         cerebellum=paste0(region," (",hemisphere,")") %>% gsub(" (NA)","",.)) %>%
  dplyr::select(cerebellum,subcortical,rg_FAST) %>%
  pivot_wider(names_from=cerebellum,values_from=rg_FAST)

fast_rg_matrix<-as.matrix(fast_rg[,-1])
rownames(fast_rg_matrix)<-fast_rg$subcortical
# correlation between patterns --> to get how similar each of them are
# get order for subcortical structures
hc_sub<-hclust(dist(fast_rg_matrix))
order_sub<-rownames(fast_rg_matrix)[hc_sub$order]
# get order for cerebellar measurs
hc_cer<-hclust(dist(t(fast_rg_matrix)))
order_cer<-colnames(fast_rg_matrix)[hc_cer$order]

png(paste0(out_dir,"corrplot_fast_cerebellum_rg.png"),width=1000,height=700)
corrplot(fast_rg_matrix[order_sub,order_cer],method="square",
         cl.ratio = 0.2, tl.col = "black", tl.srt = 90, 
         col = COL2('PuOr', 10)
)
dev.off()

## aseg
aseg_rg<-rg_subcort_wide %>% 
  mutate(subcortical=paste0(region.subcortical," (",hemisphere.p2,")")%>% gsub(" (NA)","",.),
         cerebellum=paste0(region," (",hemisphere,")") %>% gsub(" (NA)","",.)) %>%
  dplyr::select(cerebellum,subcortical,rg_aseg) %>%
  pivot_wider(names_from=cerebellum,values_from=rg_aseg)

aseg_rg_matrix<-as.matrix(aseg_rg[,-1])
rownames(aseg_rg_matrix)<-aseg_rg$subcortical
# correlation between patterns --> to get how similar each of them are
# get order for subcortical structures
hc_sub<-hclust(dist(aseg_rg_matrix))
order_sub<-rownames(aseg_rg_matrix)[hc_sub$order]
# get order for cerebellar measurs
hc_cer<-hclust(dist(t(aseg_rg_matrix)))
order_cer<-colnames(aseg_rg_matrix)[hc_cer$order]

png(paste0(out_dir,"corrplot_aseg_cerebellum_rg.png"),width=1000,height=700)
corrplot(aseg_rg_matrix[order_sub,order_cer],method="square",
         cl.ratio = 0.2, tl.col = "black", tl.srt = 90, col = COL2('PuOr', 10),
        
         )
dev.off()
