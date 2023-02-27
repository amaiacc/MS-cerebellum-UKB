#' ---
#' title: 'Genetic correlation of cerebellar structures and psychiatric disorders'
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
#' Summary statistics for other traits of interest were obtained from publicly available GWASes.
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
idps$pheno<-numeric_nchar(idps$Pheno)
# define lateralized measures
idps_lat<-read.table(paste(primary_dir,"IDPs_lat_summary.txt",sep=""),header=TRUE) %>% filter(is.na(NotLat)) %>% dplyr::select(-NotLat)
idps_lat$L[!is.na(idps_lat$L)]<-numeric_nchar(idps_lat$L[!is.na(idps_lat$L)])
idps_lat$R[!is.na(idps_lat$R)]<-numeric_nchar(idps_lat$R[!is.na(idps_lat$R)])
idps_lat<-merge(idps_lat,idps[,c("pheno","region","lobe")],by.x="L",by.y="pheno")


# define order of regions:
cer_regions<-c("TOTAL","CORTEX","WHITE-MATTER",
               "I-IV","V","VI","CRUS I","CRUS II",
               "VIIB","VIIIA","VIIIB","IX","X")

# define order for traits -> in plot
traits_order<-c("Grove_ASD","PGC3_SCZ_primary","Lee2018_CP")

# read data
rg_dis<-read.csv(paste0(ldsc_dir,"ldsc_rg_disorders_cog_UKB_BIG40.csv"))
rg_dis<-rg_dis %>% mutate(measure=if_else(p1=="Chambers2022_TotalCerebellarVolume","VOLUME",measure),
               region=if_else(p1=="Chambers2022_TotalCerebellarVolume","TotalCerebellarVolume",region)
)
# filter only cerebellum
rg_dis<-rg_dis[grep("CEREBELL",toupper(rg_dis$region)),]
# edit some columns
rg_dis <- rg_dis %>% 
  filter(measure!="INTENSITY") %>%
  mutate(hemisphere=if_else(grepl("^V_",region),"V",hemisphere),
         region=gsub("^V_","",region) %>% 
           gsub("CEREBELLUM_|CEREBELLUM-|CerebellarVolume","",.) %>% 
           gsub("_"," ",.) %>% toupper ) %>%
  mutate(hemisphere=if_else(is.na(hemisphere),"BiLat",hemisphere)) %>%
  mutate(hemisphere=factor(hemisphere,levels=c("BiLat","V","L","R")),
         type=if_else((region=="CORTEX"|region=="WHITE-MATTER"|region=="TOTAL"),"Global","Region"),
         region=factor(region,levels=cer_regions),
         p2=factor(p2,levels=traits_order))

rg_dis <- rg_dis %>%
  mutate(sig=if_else(p<(0.05/(3*ntests)),"*",""))


# save files
write.csv(rg_dis,file=paste(out_dir,"/ldsc_rg_disorders","_",root,"_",project,".csv",sep=""),row.names = FALSE)
#----------------------------------------------------------------------
# generate plots to visualize genetic correlations
traits_order2<-c("ASD","SCZ","CP")
levels(rg_dis$p2)<-traits_order2
rg_plot_list<-lapply(c("rg_dis"),function(rg){
  get(rg) %>%  mutate(region=factor(region,levels=rev(cer_regions))) %>%
    ggplot(aes(x=region,y=rg,width=0.7,color=hemisphere,alpha=sig)) +
    geom_hline(yintercept = 0,color="black",alpha=0.7,linetype="dashed") +
    geom_point(position=position_dodge(.7), stat="identity",size=3) +
    # geom_bar(aes(y=rg),position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=rg.CI95upper , ymax=rg.CI95lower,width=.5), position=position_dodge(.7)) +
    labs(y=bquote('Estimate ('*rho*')'),
         x="",
         color="Cerebellar hemisphere") +
    coord_flip(ylim=c(-0.25,0.25)) +
    facet_grid(type ~ p2, scales = "free",space ="free") +
    # scale_color_met_d(name="Morgenstern") +
    scale_color_viridis_d() +
    scale_alpha_manual(values=c(0.3,0.9)) +
    guides(alpha="none") +
    theme(panel.spacing=unit(0, "lines"),
          strip.background =element_rect(fill="white",color="grey85",size=1),
          legend.position = "none") +
    panel_border(size = 1)+
    NULL
})

rg_legend<-get_legend(rg_plot_list[[1]] + theme(legend.position = "bottom"))

prow<-plot_grid(plotlist=rg_plot_list,nrow=1,rel_widths = c(2,1,1))

combined_plot<-plot_grid(prow,rg_legend,rel_heights = c(1,0.1),ncol=1)
combined_plot

ggsave(combined_plot,file=paste0(out_dir,"Figure_ldsc_rg_disorders_cog.png"),width=7,height=5)

# to do's for plot:
# - flag nom. sig.
# - legend center, make legend name smaller.
# - fix label overlap in x-axis.
#--------------------------------------------


#----------------------------------------------------------------------
#' ### Cerebellar measures
#----------------------------------------------------------------------
#' These cerebellar volumes were selected for the genetic correlation analysis, and there were three measures per region (left, right and vermis):
#' 
#' Global measures:
#' 
#' * Cortex
#' * White-matter
#' 
#' Lobules:
#' 
#' * Crus I
#' * Crus II
#' * X
#' * VI

#----------------------------------------------------------------------
#' ## Disorders
#----------------------------------------------------------------------
#' The following GWAS summary statistics were used:
#' 
#' * SCZ
#' * ASD
#' * ADHD
#' 
# The only nominally significant genetic correlation with disorders (p= `r subset(rg_dis,p<0.05)$p`).
rg_dis %>% rename(Behaviour=p2) %>% 
  select(-IDP_short_name,-measure,-gcov_int,-gcov_int_se) %>%
  arrange(p) %>%
  mutate(across(where(is.numeric),round,digits=2)) %>%
  datatable(rownames=FALSE,caption="Genetic correlation between handedness and  cerebellar structures.")

rg_plot_list[[1]] + theme(legend.position = "bottom")

