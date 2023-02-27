#' ---
#' title: 'Heritability and genetic correlation of cerebellar structures'
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
#' **SNP heritability** is the proportion of variance explained by common genetic factors, and was computed using GWAS summary statistics by running the LDSC tool.
#' 
#' * This measure ranges from 0-1.
#' * We typically test whether the heritability estimate is significantly higher than 0.
#' 
#' **Genetic correlation** is the proportion of variance explained by common genetic factors. 
#' 
#' It was computed using GWAS summary statistics by running the LDSC tool.
#' 
#' * This measure ranges from 0-1. 
#' * When testing the genetic correlation of the left and right of the same structure, I tested whether the genetic correlation estimate is significantly lower than 1.
#'     * 1-tailed test compared with a normal distribution (z = (1-rg)/SE): `dnorm(z)/2`
#' * For the rest of the correlations, I tested whether the genetic correlation estimate is significantly higher than 0.
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
idps_lat<-read.table(paste(primary_dir,"IDPs_lat_summary.txt",sep=""),header=TRUE) %>% filter(is.na(NotLat)) %>% select(-NotLat)
idps_lat$L[!is.na(idps_lat$L)]<-numeric_nchar(idps_lat$L[!is.na(idps_lat$L)])
idps_lat$R[!is.na(idps_lat$R)]<-numeric_nchar(idps_lat$R[!is.na(idps_lat$R)])
idps_lat<-merge(idps_lat,idps[,c("pheno","region","lobe")],by.x="L",by.y="pheno")
# define order of regions:
cer_regions<-c("TOTAL","CORTEX","WHITE-MATTER",
               "I-IV","V","VI","CRUS I","CRUS II",
               "VIIB","VIIIA","VIIIB","IX","X")
# read estimates: h2 and rg
# and subset only cerebellum measures
h2 <- read.csv(paste(ldsc_dir,"ldsc_h2_UKB_BIG40.csv",sep=""))
h2 <- h2 %>%
  mutate(measure=if_else(pheno=="Chambers2022_TotalCerebellarVolume","VOLUME",measure),
         region=if_else(pheno=="Chambers2022_TotalCerebellarVolume","TotalCerebellarVolume",region)
  )
#
h2 <- h2[grep("cereb",tolower(h2$IDP_short_name)),] %>% filter(measure=="VOLUME") %>% arrange(parcellation,IDP_short_name)
h2 <- h2 %>% mutate(hemisphere=if_else(grepl("^V_",region),"V",hemisphere),
                    region=gsub("^V_","",region) %>% 
                      gsub("CEREBELLUM_|CEREBELLUM-|CerebellarVolume","",.) %>% 
                      gsub("_"," ",.) %>% toupper ) %>%
  mutate(hemisphere=if_else(is.na(hemisphere),"BiLat",hemisphere)) %>%
  mutate(hemisphere=factor(hemisphere,levels=c("BiLat","V","L","R")),
         type=if_else((region=="CORTEX"|region=="WHITE-MATTER"|region=="TOTAL"),"Global","Region"),
         region=factor(region,levels=cer_regions)) %>%
  # define adjusted p-val
  mutate(p.BNFadj=p.adjust(p,method = "bonferroni"))

h2 <- h2 %>%  mutate(sig=if_else(p<(0.05/ntests),"*",""))

# genetic correlations
# total cerebellum with regions
rg_total <-read.csv(paste0(ldsc_dir,"ldsc_rg_other_UKB_BIG40.csv"))
rg_total <- rg_total %>% 
  filter(grepl("TotalCerebellarVolume", p1)) %>%
  filter(grepl("cereb",tolower(IDP_short_name.p2)))
  # filter(measure.p2=="VOLUME") %>% arrange(parcellation.p2,IDP_short_name_lat.p2)

# left and right
rg_lr <- read.csv(paste(ldsc_dir,"ldsc_rgLR_UKB_BIG40.csv",sep=""))
rg_lr <- rg_lr[grep("cereb",tolower(rg_lr$IDP_short_name_lat)),] %>% filter(measure=="VOLUME") %>% arrange(parcellation,IDP_short_name_lat)
rg_lr <- rg_lr %>% mutate(region=gsub("^V_","",region) %>% gsub("CEREBELLUM_|CEREBELLUM-","",.) %>% gsub("_"," ",.) ) %>%
  mutate(
    hemisphere="LR",
    region=factor(region,levels=cer_regions),
    type=if_else((region=="CORTEX"|region=="WHITE-MATTER"),"Global","Region")
  ) %>%
  # define adjusted p-val
  mutate(p.BNFadj=p.adjust(p2,method = "bonferroni"))
#
rg_atlas <- read.csv(paste(ldsc_dir,"ldsc_rg_atlas_UKB_BIG40.csv",sep="")) %>% filter(measure=="VOLUME") %>%
            filter(IDP_short_name.p1!=IDP_short_name.p2)
rg_atlas <- rg_atlas[intersect(grep("cer",tolower(rg_atlas$IDP_short_name.p2)),grep("cer",tolower(rg_atlas$IDP_short_name.p1))),]
rg_atlas <- rg_atlas %>% 
  mutate(hemisphere.p1=if_else(grepl("^V_",region.p1),"V",hemisphere.p1),
         hemisphere.p2=if_else(grepl("^V_",region.p2),"V",hemisphere.p2),
         region.p1=gsub("^V_","",region.p1) %>% gsub("CEREBELLUM_|CEREBELLUM-","",.) %>% gsub("_"," ",.),
         region.p2=gsub("^V_","",region.p2) %>% gsub("CEREBELLUM_|CEREBELLUM-","",.) %>% gsub("_"," ",.) 
                ) %>%
  mutate(hemisphere=if_else(hemisphere.p1==hemisphere.p2,hemisphere.p1,paste0(hemisphere.p1,hemisphere.p2))) %>%
  mutate(hemisphere=gsub("RL","LR",hemisphere) %>% gsub("LV","VL",.) %>% gsub("RV","VR",.)) %>%
  mutate(hemisphere=factor(hemisphere,levels=c("V","L","R","LR","VL","VR")),
         region.p1=factor(region.p1,levels=cer_regions),
         region.p2=factor(region.p2,levels=cer_regions)
         )
# filter out vermis/hemisphere (if not the same substructure)
rg_atlas2 <- rg_atlas %>% filter(region.p1!=region.p2&(hemisphere %in% c("LR","VL","VR")) )
rg_atlas <- rg_atlas %>% filter(!(region.p1!=region.p2&(hemisphere %in% c("LR","VL","VR")) ) )
# define adjusted p-val
rg_atlas <- rg_atlas %>% 
  mutate(p.BNFadj=p.adjust(p,method = "bonferroni"))

# define left-vermis and right-vermis rg for each measure
rg_lr2<-subset(rg_atlas,region.p1==region.p2&hemisphere %in% c("LR","VL","VR")) %>%
  mutate(region=region.p1)
rg_lr2<-merge(rg_lr2,rg_lr,by=c("parcellation","measure","region","hemisphere",
                                "rg","rg.SE","rg.CI95upper","rg.CI95lower","gcov_int","gcov_int_se"),
              suffixes=c("",".diff1"),all=TRUE)
rg_lr2 <- rg_lr2 %>% 
  mutate(region=factor(region,levels=cer_regions),
         type=if_else((region=="CORTEX"|region=="WHITE-MATTER"),"Global","Region")
         )



# save files
write.csv(h2,file=paste(out_dir,"/ldsc_h2","_",root,"_",project,".csv",sep=""),row.names = FALSE)
write.csv(rg_lr2,file=paste(out_dir,"/ldsc_rg_LR_VL_VR","_",root,"_",project,".csv",sep=""),row.names = FALSE)
write.csv(rg_atlas,file=paste(out_dir,"/ldsc_rg_atlas","_",root,"_",project,".csv",sep=""),row.names = FALSE)
write.csv(rg_total,file=paste(out_dir,"/ldsc_rg_total","_",root,"_",project,".csv",sep=""),row.names = FALSE)

#----------------------------------------------------------------------

## generate plots
h2_cer_plot <- h2 %>% mutate(region=factor(region,levels=rev(cer_regions))) %>%
  ggplot(.,aes(x=region,y=h2,width=0.7,
                                # color=measure,
                                # shape=hemisphere,
                                color=hemisphere, alpha=sig)) +
  geom_hline(yintercept = 0,color="black",alpha=0.7,linetype="dashed") +
  geom_point(position=position_dodge(.7), stat="identity",size=3)+
  geom_errorbar(aes(ymin=h2.CI95upper , ymax=h2.CI95lower,width=.5), position=position_dodge(.7)) +
  geom_text(aes(x=region,y=h2.CI95lower*0.8,label=sig,color=hemisphere),alpha=0.8,
            fontface=2,size=6,position=position_dodge(.8),stat="identity") +
  labs(title="",x="",y=bquote('Estimate ('*h^2*')'),color="",shape="") +
  ylim(c(0,0.5)) +
  facet_grid(type~.,scales = "free",space ="free") +
  coord_flip() +
  # scale_color_met_d(name="Morgenstern") +
  scale_color_viridis_d() +
  scale_alpha_manual(values=c(0.9)) +
  guides(alpha="none") +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="grey85",size=1),
        legend.position = "bottom") +
  panel_border(size = 1)+
  NULL

# generate plots to visualize genetic correlations between left and right
w<-which(is.na(rg_lr2$p.diff1))
rg_lr2$p.diff1[w]<-rg_lr2$p2.diff1[w]
rg_lr2 <- rg_lr2 %>% 
  mutate( pval=if_else(!is.na(p2.diff1),p2.diff1,p),
          sig=if_else(pval<(0.05/NROW(rg_lr2)),"*","") 
          )

rg_cer_plot <- rg_lr2 %>% mutate(region=factor(region,levels=rev(cer_regions))) %>%
  mutate(sig_val=if_else(rg<0.9,rg.CI95lower*0.8,rg.CI95upper*1.1)) %>%
  ggplot(data=.,aes(x=region,y=rg,width=0.7,
         # color=measure,
         # shape=hemisphere,
         color=hemisphere,alpha=sig) ) +
  geom_hline(yintercept = c(0,1),color="black",alpha=0.7,linetype="dashed") +
  geom_point(position=position_dodge(.7), stat="identity",size=3,shape="diamond") +
  geom_errorbar(aes(ymin=rg.CI95upper , ymax=rg.CI95lower,width=.5), position=position_dodge(.7)) +
  geom_text(aes(x=region,y=sig_val,label=sig,color=hemisphere),alpha=0.8,
            fontface=2,size=6,position=position_dodge(.8),stat="identity") +
  # labs(title="",x="",y=bquote('Estimate ('*rho[LR]*')')) +
  labs(title="",x="",y="Estimate(rg)",color="",shape="") +
  ylim(c(-0.5,1.1)) +
  facet_grid(type~.,scales = "free",space ="free") +
  coord_flip() +
  scale_color_met_d(name="Demuth") +
  scale_alpha_manual(values=c(0.3,0.9)) +
  guides(alpha="none") +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="grey85",size=1),
        legend.position = "bottom") +
  panel_border(size = 1) +
  NULL
h2_rgLR_combined <- plot_grid(h2_cer_plot,rg_cer_plot,labels=c("A","B"))
# leg<-get_legend(h2_cer_plot + theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))
# comb<-plot_grid(h2_cer_plot + theme(legend.position = "none"),
#           rg_cer_plot, #+ theme( axis.text.y=element_blank() ),
#           rel_widths = c(1,1),nrow=1,
#           labels=c("A","B"))
# h2_rgLR_combined <- plot_grid(comb,
#           plot_grid(leg,NULL,nrow=1,rel_widths=c(1,0.4)),
#           ncol=1,rel_heights =c(1,0.1))
# rm(leg,comb)
#----------------------------------------------------------------------
#' ## Heritability
#----------------------------------------------------------------------
#' All cerebellar measures are heritable (estimates significantly greater than 0). 
#' Most estimates are in the range of 0.20-0.35, with the exception of: _Crus I (vermis)_ and the total white matter intenstiy measures.
#' 
#' For the lateral measures, the estimates of the left and the right are in the same range.
#' 
#' 
#' 
h2_cer_plot

h2 %>% select(-IDP_short_name,-Lambda.GC,-h2.CI95upper,-h2.CI95lower,-lobe) %>% datatable(rownames=FALSE,caption="Heritability estimates of cerebellar structures.")


#----------------------------------------------------------------------
#' ## Genetic correlations
#----------------------------------------------------------------------
#' ### Left and right
#' 
#' The genetic correlation between left and right substructures is close to 1.
#' 
#' However, it is significantly lower than 1 for several substructures, most notably: _X_, _IX_ and _I-IV_ (where p-values < 0.05 / `r NROW(rg_lr)`).
rg_cer_plot

rg_lr %>%  
  select(-IDP_short_name_lat,-rg.CI95upper,-rg.CI95lower,-lobe,-L,-R,-p,-gcov_int,-gcov_int_se) %>%  
  arrange(p2) %>% rename(p=p2) %>%
  datatable(rownames=FALSE,caption="Genetic correlation estimates for L and R sides of cerebellar structures.")

#--------------------------------------------
#' ### Within atlas
parc="FAST"
#' #### `r parc`
m="VOLUME"
#par(mfrow=c(1,2))
#' Genetic correlation for `r parc` atlas substructures (within hemisphere). 

# corrplot_rg_atlas(atlas=rg_atlas,lr=rg_lr,parc=parc,m=m) # previous version using corrplot
rg_atlas_lat<-rg_atlas %>%
  filter(hemisphere=="L"|hemisphere=="R") %>%
  # define corrected p-val to be used in plot
  mutate(p.unadj=p,
         p=p.BNFadj)


corrplot_lr<-ggcorrplot_rg_atlas(atlas=rg_atlas %>% mutate(p=p.BNFadj),
                                 lr=rg_lr %>% mutate(p2=p.BNFadj),
                                 parc=parc,m=m)
corrplot_lr

#' The diagonal shows the genetic correlations between left and right.
#' The upper triangle shows the correlations across left hemisphere substructures
#' and the lower triangle shows the correlations across right hemisphere substructures.
#' Non-significant correlations are in white.
#' 
#' 
#' Genetic correlation within `r parc` atlas substructures (non-lateral measures). 
rg_atlas_notlat<-rg_atlas %>% filter(hemisphere=="V") %>%
  # define corrected p-val to be used in plot
  mutate(p.unadj=p,
         p=p.BNFadj)
corrplot_notlat<-ggcorrplot_rg_notlat_atlas(atlas=rg_atlas_notlat,lr=rg_lr,parc=parc,m=m)
corrplot_notlat
#par(mfrow=c(1,1))
rg_atlas %>%  
  select(-IDP_short_name.p1,-IDP_short_name.p2,-rg.CI95upper,-rg.CI95lower,-p1,-p2, -gcov_int,-gcov_int_se) %>%  
  arrange(p) %>%
  datatable(rownames=FALSE,caption="Genetic correlation estimates for cerebellar structures within atlas.")
#--------------------------------------------
leg<-get_legend(corrplot_lr + theme(legend.position="bottom",legend.justification="center"))
corrplots_combined <- plot_grid(plot_grid(
  corrplot_lr + theme(legend.position="none"),
  # NULL,
  corrplot_notlat + theme(legend.position="none"),
  nrow=1,rel_widths = c(1,1),labels=c("C","D")),
  leg, ncol=1,rel_heights = c(1,0.2)
)

fig4ms<-plot_grid(h2_rgLR_combined,corrplots_combined,ncol=1,rel_heights = c(1,1))
ggsave(fig4ms,file=paste0(out_dir,"/Figure_cerebellar_h2_rg.png"),width=12,height=12)
