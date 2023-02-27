#' ---
#' title: 'Partitioned heritability analysis: evolutionary annotations'
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
#' 
#' ## Analysis
#' 
#' **Partitioned SNP heritability** is the proportion of variance explained by common genetic factors, and was computed using GWAS summary statistics by running the LDSC tool.
#' 
#' * This measure ranges from 0-1.
#' * We typically test whether the heritability estimate is significantly higher than 0.
#' 

#---------------------------------------------------
# clean workspace
rm(list=ls())
# libraries and custom functions
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
# library(corrplot)
library(ggcorrplot)
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
# ldsc_dir=paste(dir,"resources/datasets/GWAS_sumstats/ldsc/",root,"/",sep="")
out_dir=paste(dir,project,"data/ldsc/partitioned/",sep="/")
if(!dir.exists(out_dir)){dir.create(out_dir)}
# set working dir
setwd(out_dir)
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
idps<-read.csv(paste(primary_dir,"IDPs_summary.csv",sep=""),header=TRUE) %>% filter(Cerebellum==1)
idps$pheno<-numeric_nchar(idps$Pheno)
idps_short<-idps[,c("pheno","region","IDP_short_name","IDP_description","hemisphere")]
 
# define lateralized measures
idps_lat<-read.table(paste(primary_dir,"IDPs_lat_summary.txt",sep=""),header=TRUE) %>% filter(is.na(NotLat)) %>% dplyr::select(-NotLat)
idps_lat$L[!is.na(idps_lat$L)]<-numeric_nchar(idps_lat$L[!is.na(idps_lat$L)])
idps_lat$R[!is.na(idps_lat$R)]<-numeric_nchar(idps_lat$R[!is.na(idps_lat$R)])
idps_lat<-merge(idps_lat,idps[,c("pheno","region","lobe")],by.x="L",by.y="pheno")
# define order of regions:
cer_regions<-c("TOTAL","CORTEX","WHITE-MATTER",
               "I-IV","V","VI","CRUS I","CRUS II",
               "VIIB","VIIIA","VIIIB","IX","X")
#----------------------------------------------------------------------
# read enrichment estimates from the S-LDSC run
sh2<-read.table(paste0(out_dir,"/evol_annots_cerebellum.table"),header=TRUE)
# clean file
colnames(sh2)<-c("Category",colnames(sh2)[2:NCOL(sh2)])
sh2<- sh2 %>% filter(Enrichment!="Enrichment") %>%
  mutate(file=gsub(".Category","",Category),
         pheno=strsplit(Category,"_") %>% sapply("[[",1),
         annot=strsplit(Category,"annots_cerebellum_") %>% sapply("[[",2) %>% gsub(".results:L2_1","",.)
         ) %>%
  dplyr::select(-Category)

# combine with idp info
sh2<-merge(sh2,idps_short,by="pheno",all.x=TRUE)
# fix cols for Chambers et al. 2022
sh2<-sh2 %>% mutate(region=if_else(pheno=="Chambers2022","TotalCerebellarVolume",region) )

#
sh2 <- sh2 %>% 
  mutate(hemisphere=if_else(grepl("^V_",region),"V",hemisphere),
         region=gsub("^V_","",region) %>% 
           gsub("CEREBELLUM_|CEREBELLUM-|CerebellarVolume","",.) %>% 
           gsub("_"," ",.) %>% toupper ) %>%
  mutate(hemisphere=if_else(is.na(hemisphere),"BiLat",hemisphere)) %>%
  mutate(hemisphere=factor(hemisphere,levels=c("BiLat","V","L","R")),
         type=if_else((region=="CORTEX"|region=="WHITE-MATTER"|region=="TOTAL"),"Global","Region"),
         region=factor(region,levels=cer_regions) ) %>%
  # define adjusted p-val
  mutate(p.BNFadj=as.numeric(Enrichment_p)*ntests*6,
         p.FDRadj=p.adjust(Enrichment_p,method = "fdr"))

# make columns numeric
sh2<-sh2 %>%
  mutate(stat="h2-stratified", program="LDSC") %>%
  mutate(across(contains("Prop"), as.numeric)) %>%
  mutate(across(contains("Enrichment"), as.numeric)) %>%
  # compute 95% CIs
  mutate(Enrichment.CI95upper=Enrichment+(1.96*Enrichment_std_error),
         Enrichment.CI95lower=Enrichment-(1.96*Enrichment_std_error)
  )
#----------------------------------------------------------------------
# save files
write.csv(sh2,file=paste0(out_dir,"/sldsc_","evol_annots","_",root,"_",project,".csv"),row.names = FALSE)
#----------------------------------------------------------------------
## generate plots
# a="Fetal_brain_GSE63648_7pcw_Hu_gain"
annots<-unique(sh2$annot) %>% sort()
# sort order of annotations
annots<-annots[c(2,1,3,6,5,4)]
sh2<-sh2 %>% mutate(
  annot=factor(annot,levels=annots),
  clean_annot=gsub("_hg19.merged|.all_sorted.merged|_Cerebellum","",annot) %>% 
    gsub("_Chen2020|_Browning2018|_Peyregne2017|_Capra2013|_GSE63648","",.) %>%
    gsub("_"," ",.) %>%
    gsub("Adult enhancers promoters","Adult brain\n(Cerebellum)\nHuman gained\nenhancers/promoters ",.) %>%
    gsub("Fetal brain Hu gain","Fetal brain\nHuman gained\nenhancers",.) %>%
    gsub("NeanderthalDepleted","Archaic Depleted\nRegions",.) %>%
    gsub("NeanderthalMatchSNPs sprime", "Neanderthal\nIntrogressed SNPs\n(Sprime)",.) %>%
    gsub("SelectiveSweep","Ancient\nSelective Sweeps",.) %>%
    gsub("ncHAR","Human Accelerated\nRegions",.)
  ) %>% 
  arrange(annot) %>%
  mutate(clean_annot=factor(clean_annot,unique(clean_annot))) %>%
  mutate(sig=if_else(p.FDRadj<0.05,"*",""))

sh2_cer_plotlist <- lapply(annots, function(a){
    tmp<-sh2 %>% filter(annot==a) %>% mutate(region=factor(region,levels=rev(cer_regions)))
    clean_annot<-unique(tmp$clean_annot)
    ggplot(tmp,aes(x=region,y=Enrichment,width=0.7,
                   color=hemisphere)) +
      geom_hline(yintercept = 1,color="black",alpha=0.7,linetype="dashed") +
      geom_point(position=position_dodge(.7), stat="identity",size=3) +
      geom_errorbar(aes(ymin=Enrichment.CI95upper, 
                        ymax=Enrichment.CI95lower,width=.5), position=position_dodge(.7)) +
      labs(subtitle=unique(clean_annot),
           x="",
           y="Enrichment", #bquote('Estimate ('*h^2*')'),
           color="",shape="") +
      # ylim(c(0,0.5)) +
      facet_grid(type~.,scales = "free",space ="free") +
      coord_flip() +
      # scale_color_met_d(name="Morgenstern") +
      scale_color_viridis_d() +
      theme(panel.spacing=unit(0, "lines"),
            strip.background =element_rect(fill="white",color="grey85",size=1),
            legend.position = "none") +
      panel_border(size = 1)+
      NULL
  })
leg<-get_legend(sh2_cer_plotlist[[1]]+ theme(legend.position="bottom") + 
                  labs(color="Hemisphere",shape="Hemisphere") )
sh2_cer_plot<-plot_grid(
  plot_grid(plotlist=sh2_cer_plotlist,nrow=2),leg,ncol=1,rel_heights = c(1,0.1))
  
ggsave(sh2_cer_plot, 
       file=paste0(out_dir,"/sldsc_","evol_annots","_",root,"_",project,"_v1.png"),
       width=12,height=10)
# different figure for same plot
sh2_cer_plot2<-sh2 %>%
  mutate(region=factor(region,levels=rev(cer_regions))) %>%
  mutate(sig_val=if_else(Enrichment<1,Enrichment.CI95lower*0.8,Enrichment.CI95upper*1.5)) %>%
  ggplot(.,aes(x=region,y=Enrichment,width=0.7,
               # color=measure,
               # shape=hemisphere,
               color=hemisphere, alpha=sig)) +
  geom_hline(yintercept = 1,color="black",alpha=0.7,linetype="dashed") +
  geom_point(position=position_dodge(.7), stat="identity",size=3) +
  geom_errorbar(aes(ymin=Enrichment.CI95upper, 
                    ymax=Enrichment.CI95lower,width=.5), position=position_dodge(.7)) +
  geom_text(aes(x=region,y=-0.1,label=sig,color=hemisphere),
            fontface=2,size=6,position=position_dodge(.8),stat="identity") +
  labs(x="",
       y="Enrichment", #bquote('Estimate ('*h^2*')'),
       color="",shape="") +
  # ylim(c(0,0.5)) +
  facet_grid(type~clean_annot,scales = "free",space="free_y") +
  coord_flip() +
  # scale_color_met_d(name="Morgenstern") +
  scale_color_viridis_d() +
  scale_alpha_manual(values=c(0.3,0.9)) +
  guides(alpha="none") +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="grey85",size=1),
        legend.position = "bottom") +
  panel_border(size = 1) +
  NULL

ggsave(sh2_cer_plot2, 
       file=paste0(out_dir,"/Figure_sldsc_","evol_annots","_",root,"_",project,".png"),
       width=12,height=6)

# compare Enrichment estimate for:
## archaic deserts and fetal brain HGE

sh2_wide<- sh2 %>% 
  dplyr::select(Enrichment,Enrichment.CI95upper,Enrichment.CI95lower,Enrichment_p,p.BNFadj,annot,region,hemisphere,type) %>%
  pivot_wider(id_cols = c(region,hemisphere,type),
              names_from = annot,values_from = c(Enrichment,Enrichment.CI95upper,Enrichment.CI95lower,Enrichment_p,p.BNFadj))

annots<-unique(sh2$annot) %>% as.character()

enrichment_cors <- sh2_wide %>% group_by(hemisphere) %>% 
  select(starts_with("Enrichment_")) %>% psych::corr.test()


sh2_wide2<-sh2_wide %>%
  mutate(sig=if_else(
    (Enrichment_p_Fetal_brain_GSE63648_Hu_gain.all_sorted.merged>0.05&Enrichment_p_NeanderthalDepleted_Chen2020),
    "n.s.", 
    if_else(
      (p.BNFadj_Fetal_brain_GSE63648_Hu_gain.all_sorted.merged<0.05|p.BNFadj_NeanderthalDepleted_Chen2020<0.05),
      "BFadj p<0.05",
      "unadj p<0.05" 
    )
  ) ) %>%
  mutate(sig=factor(sig,levels=c("n.s.","unadj p<0.05","BFadj p<0.05")))

sh2_cer_ArchaicDeserts_FetalHGE <- sh2_wide2 %>%
 ggplot(.,aes(x=Enrichment_Fetal_brain_GSE63648_Hu_gain.all_sorted.merged,
              y=Enrichment_NeanderthalDepleted_Chen2020,
              alpha=sig,
              color=hemisphere)) +
  geom_hline(yintercept = 1,color="black",alpha=0.7,linetype="dashed") +
  geom_vline(xintercept = 1,color="black",alpha=0.7,linetype="dashed") +
  geom_point(position=position_dodge(.7), stat="identity",size=3) +
  geom_errorbarh(aes(xmin=Enrichment.CI95upper_Fetal_brain_GSE63648_Hu_gain.all_sorted.merged, 
                    xmax=Enrichment.CI95lower_Fetal_brain_GSE63648_Hu_gain.all_sorted.merged,width=.5), position=position_dodge(.7)) +
  geom_errorbar(aes(ymin=Enrichment.CI95upper_NeanderthalDepleted_Chen2020, 
                     ymax=Enrichment.CI95lower_NeanderthalDepleted_Chen2020,width=.5), position=position_dodge(.7)) +
  facet_grid(type~hemisphere,scales = "free",space="free_y") +
  # geom_text(aes(x=region,y=-0.1,label=sig,color=hemisphere),
  #           fontface=2,size=6,position=position_dodge(.8),stat="identity") +
  labs(x="Enrichment (fetal brain HGE)",
       y="Enrichment (Archaic Depleted Regions)", #bquote('Estimate ('*h^2*')'),
       color="",shape="") +
  coord_cartesian(xlim=c(-2,8),ylim=c(-1,2)) +
  scale_color_viridis_d() +
  scale_alpha_manual(values=c(0.1,0.6,1)) +
  guides() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="grey85",size=1),
        legend.position = "bottom") +
  panel_border(size = 1) +
  NULL


ggsave(sh2_cer_ArchaicDeserts_FetalHGE,
       file=paste0("sldsc_fetalHGE_ArchaicDepleted_cerebellum_UKB.png"),
       width=8,height=6)
