#' ---
#' title: 'Heritability analysis of the *reading network*'
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
#' ## Description and aim
#' Heritability analysis based on the Imaging Derived Phenotypes (IDPs) from the UK Biobank BIG40 server:
#' 
#'    - https://open.win.ox.ac.uk/ukbiobank/big40/
#'    - GWAS summary stats computed from N~31,000
#'    
#' SNP heritability (proportion of variance explained by common genetic factors) was computed using GWAS summary statistics by running the LDSC tool.
#' This measure ranges from 0-1, and we typically test whether the heritability estimate is significantly higher than 0.
#' 
#' The aim of the current analysis is to prioritize regions from the reading network based on their heritability estimates for further downstream analyses (e.g. for polygenic risk scores).
#' 
#' 
#---------------------------------------------------
# clean workspace
rm(list=ls())
# libraries and custom functions
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(DT)

#---------------------------------------------------
options(stringsAsFactors = FALSE)
root="UKB_BIG40"
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
primary_dir=paste(dir,"resources/datasets/GWAS_sumstats/downloaded_data/UKB/","BIG40","/",sep="")
ldsc_dir=paste(dir,"resources/datasets/GWAS_sumstats/ldsc/",root,"/",sep="")
# set working dir
setwd(ldsc_dir)
#----------------------------------------------------------------------
# some functions
source(paste(dir,"general_scripts/helper_functions.R",sep=""))
source(paste(dir,"general_scripts/genotyping/ldsc/helper_functions_ldsc.R",sep=""))
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
# read estimates: h2 and rg
h2<-read.csv(paste(ldsc_dir,"ldsc_h2_UKB_BIG40.csv",sep=""))
rg<-read.csv(paste(ldsc_dir,"ldsc_rgLR_UKB_BIG40.csv",sep=""))
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#' ### Selected Imaging Derived Phenotypes (IDPs)
#' `r NROW(h2)` IDPs were selected broadly based on their relevance for the reading network (see below and attached for the full table).
#' 
#' Whenever a lateral measure was selected, the contralateral measure was also included in the analysis, 
#' and the genetic correlation between the left and right measures was computed  (i.e. `r NROW(idps_lat)` measures).
#' 
#' __Genetic correlation__ how much of the correlation between two traits can be explained by genetic variance.
#' 
#' 
#' Note that:
#' 
#' * If two traits have low phenotypic correlation it is still possible that they have a high genetic correlation. It would mean that a large proportion of that phenotypic correlation will be explained by genetic factors (even if the total amount of genetics is small).
#' * If two traits have high phenotypic correlation (e.g. left and right measures of a brain region), and similar heritability measures: the genetic correlation between the two may still be different to 1, because the genetic contribution to each of the traits may be different.
#' 
#--------------------------------------------

#----------------------------------------------------------------------
# generate plots
#----------------------------------------------------------------------
# define common palette to be used across plots
set12_pal<-c(brewer.pal(n = ceiling(length(levels(h2$measure))/2), "Set1"),brewer.pal(n = ceiling(length(levels(h2$measure))/2), "Set2"))
cols_subc<-set12_pal[which(levels(rg$measure) %in% unique(subset(h2,lobe=="SUBCORTICAL")$measure))]
cols_GM<-set12_pal[which(levels(rg$measure) %in% unique(subset(h2,lobe!="dMRI"&lobe!="SUBCORTICAL")$measure))]
cols_WM<-set12_pal[which(levels(h2$measure) %in% unique(subset(h2,lobe=="dMRI")$measure))]


# just to get the legend
h2_all1_plot<-ggplot(data=h2,aes(x=region,y=h2,width=0.7,color=measure)) +
  geom_bar(aes(y=h2),position=position_dodge(), stat="identity") +
  scale_color_manual(values=set12_pal) +
  theme(legend.position = "right",legend.direction = "horizontal")  +
  NULL
h2_all2_plot<-ggplot(data=h2,aes(x=region,y=h2,width=0.7,alpha=hemisphere)) +
  geom_bar(aes(y=h2),position=position_dodge(), stat="identity") +
  scale_color_manual(values=set12_pal) +
  theme(legend.position = "right",legend.direction = "horizontal")  +
  NULL

h2_plot_fill_legend <-  get_legend(h2_all1_plot)
h2_plot_alpha_legend <-  get_legend(h2_all2_plot)

h2_subc_plot<-ggplot(data=subset(h2,lobe=="SUBCORTICAL"),aes(x=region,y=h2,width=0.7,color=measure,alpha=hemisphere)) +
  geom_hline(yintercept = 0,color="black",alpha=0.7,linetype="dashed") +
  geom_point(position=position_dodge(.7), stat="identity",size=3)+
  # geom_bar(aes(y=h2),position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=h2.CI95upper , ymax=h2.CI95lower,width=.5), position=position_dodge(.7)) +
  ylab(bquote('Estimate ('*h^2*')')) +
  labs(title="Subcortical structures") +
  facet_grid(.~parcellation,scales = "free_y",space="free") +
  scale_color_manual(values=cols_subc) +
  scale_alpha_discrete(range=c(0.6,1)) +
  coord_flip() +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 0, hjust = 1, vjust=0.5), axis.title.y=element_blank()) +
  NULL

h2_GM_plot<- ggplot(data=subset(h2,lobe!="SUBCORTICAL"&lobe!="dMRI"),aes(x=region,y=h2,width=0.7,color=measure,alpha=hemisphere)) +
  geom_hline(yintercept = 0,color="black",alpha=0.7,linetype="dashed") +
  # geom_bar(aes(y=h2),position=position_dodge(), stat="identity") +
  geom_point(position=position_dodge(.7), stat="identity",size=3)+
  geom_errorbar(aes(ymin=h2.CI95upper , ymax=h2.CI95lower,width=.5), position=position_dodge(.7)) +
  ylab(bquote('Estimate ('*h^2*')')) +
  labs(title="Grey matter regions") +
  facet_grid(lobe~parcellation,scales = "free_y",space="free") +
  scale_color_manual(values=cols_GM) +
  scale_alpha_discrete(range=c(0.6,1)) +
  coord_flip() +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 0, hjust = 1, vjust=0.5), axis.title.y=element_blank()) +
  NULL

h2_cc_plot<-ggplot(data=h2[grep("CORPUS_CALLOSUM",h2$region),],aes(x=region,y=h2,width=0.7,color=measure)) +
  geom_hline(yintercept = 0,color="black",alpha=0.7,linetype="dashed") +
  # geom_bar(aes(y=h2),position=position_dodge(), stat="identity") +
  geom_point(position=position_dodge(.7), stat="identity",size=3)+
  geom_errorbar(aes(ymin=h2.CI95upper , ymax=h2.CI95lower,width=.5), position=position_dodge(.7)) +
  ylab(bquote('Estimate ('*h^2*')')) +
  # labs(title="White matter tracts") +
  facet_wrap(.~parcellation,scales = "free_y") +
  scale_color_manual(values=cols_WM) +
  scale_alpha_discrete(range=c(0.6,1)) +
  # coord_cartesian(ylim = c(0,max(h2$h2.CI95upper))) +
  coord_flip() +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 0, hjust = 1, vjust=0.5), axis.title.y=element_blank()) +
  NULL

h2_WMrest_plot<-ggplot(data=h2[grep("CORPUS_CALLOSUM",h2$region,invert=TRUE),] %>% filter(lobe=="dMRI"),
                       aes(x=region,y=h2,width=0.7,color=measure,alpha=hemisphere)) +
  geom_hline(yintercept = 0,color="black",alpha=0.7,linetype="dashed") +
  # geom_bar(aes(y=h2),position=position_dodge(), stat="identity") +
  geom_point(position=position_dodge(.7), stat="identity",size=3)+
  geom_errorbar(aes(ymin=h2.CI95upper , ymax=h2.CI95lower,width=.5), position=position_dodge(.7)) +
  ylab(bquote('Estimate ('*h^2*')')) +
  # labs(title="White matter tracts") +
  facet_wrap(parcellation~.,scales = "free_y",ncol=1) +
  scale_color_manual(values=cols_WM) +
  scale_alpha_discrete(range=c(0.6,1)) +
  # coord_cartesian(ylim = c(0,max(h2$h2.CI95upper))) +
  coord_flip() +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 0, hjust = 1, vjust=0.5), axis.title.y=element_blank()) +
  NULL

WMtitle <- ggdraw() + 
  draw_label(
    "White matter tracts",
    fontface = 'bold'
    # x = 0,
    # hjust = 1
  ) +
  # theme(
  #   # add margin on the left of the drawing canvas,
  #   # so title is aligned with left edge of first plot
  #   plot.margin = margin(0, 0, 20, 0)
  # ) +
  NULL

h2_WM_plot<-plot_grid(WMtitle,
                      plot_grid(h2_WMrest_plot,h2_cc_plot,align="hv",ncol=1,rel_heights = c(0.5,1.2)),
                      ncol=1,
                      # align="v",
                      rel_heights = c(0.1, 1)
                      )

h2_WM_plot2<-plot_grid(WMtitle,
                      plot_grid(h2_cc_plot,
                                h2_WMrest_plot + theme(legend.position = "bottom"),
                                align="hv",ncol=1,rel_heights = c(0.5,1.2)),
                      ncol=1,
                      # align="v",
                      rel_heights = c(0.1, 1)
                      )

h2_all_plot<-plot_grid(h2_GM_plot, h2_WM_plot,
          h2_subc_plot,
          # plot_grid(h2_plot_fill_legend,h2_plot_alpha_legend,ncol=2,rel_widths = c(1,0.4)),
          align="hv",axis="tblr",
          ncol=1,rel_heights = c(0.74,0.26,0.15,0.05))
#--------------------------------------------
# generate plots to visualize genetic correlations between left and right
rg_subc_plot<-ggplot(data=subset(rg,lobe=="SUBCORTICAL"),aes(x=region,y=rg,width=0.7,color=measure)) +
  geom_hline(yintercept = 1,color="black",alpha=0.7,linetype="dashed") +
  geom_point(position=position_dodge(.7), stat="identity",size=3,shape="triangle")+
  # geom_bar(aes(y=rg),position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=rg.CI95upper , ymax=rg.CI95lower,width=.5), position=position_dodge(.7)) +
  ylab(bquote('Estimate ('*rho*')')) +
  labs(title="Subcortical structures") +
  facet_grid(.~parcellation,scales = "free_y",space="free") +
  scale_color_manual(values=cols_subc) +
  # coord_cartesian(ylim = c(0,max(rg$rg.CI95upper))) +
  coord_flip() +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 0, hjust = 1, vjust=0.5), axis.title.y=element_blank()) +
  NULL

rg_GM_plot<- ggplot(data=subset(rg,lobe!="SUBCORTICAL"&lobe!="dMRI"),aes(x=region,y=rg,width=0.7,color=measure)) +
  geom_hline(yintercept = 1,color="black",alpha=0.7,linetype="dashed") +
  geom_point(position=position_dodge(.7), stat="identity",size=3,shape="triangle")+
  # geom_bar(aes(y=rg),position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=rg.CI95upper , ymax=rg.CI95lower,width=.5), position=position_dodge(.7)) +
  ylab(bquote('Estimate ('*rho*')')) +
  labs(title="Grey matter regions") +
  facet_grid(lobe~parcellation,scales = "free_y",space="free") +
  scale_color_manual(values=cols_GM) +
  # coord_cartesian(ylim = c(0,max(rg$rg.CI95upper))) +
  coord_flip() +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 0, hjust = 1, vjust=0.5), axis.title.y=element_blank()) +
  NULL


rg_WM_plot<-ggplot(data=rg %>% filter(lobe=="dMRI"),
                       aes(x=region,y=rg,width=0.7,color=measure)) +
  geom_hline(yintercept = 1,color="black",alpha=0.7,linetype="dashed") +
  geom_point(position=position_dodge(.7), stat="identity",size=3,shape="triangle")+
  # geom_bar(aes(y=rg),position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=rg.CI95upper , ymax=rg.CI95lower,width=.5), position=position_dodge(.7)) +
  ylab(bquote('Estimate ('*rho*')')) +
  # labs(title="White matter tracts") +
  facet_wrap(parcellation~.,scales = "free_y",ncol=1) +
  scale_color_manual(values=cols_WM) +
  # coord_cartesian(ylim = c(0,max(rg$rg.CI95upper))) +
  coord_flip() +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 0, hjust = 1, vjust=0.5), axis.title.y=element_blank()) +
  NULL

rg_all_plot<-plot_grid(rg_GM_plot, rg_WM_plot,
                       # plot_grid(rg_plot_fill_legend,rg_plot_alpha_legend,ncol=2,rel_widths = c(1,0.4)),
                       align="hv",axis="tblr",
                       ncol=1,rel_heights = c(0.74,0.35,0.05))


# #--------------------------------------------
# # # save plots and summary table
# ggsave(h2_subc_plot + theme(legend.position="right"),file=paste("h2_subcortical_",root,".png",sep=""))
# ggsave(h2_WM_plot + theme(legend.position="right"),file=paste("h2_WM_",root,".png",sep=""),width=15,height=8)
# ggsave(h2_GM_plot + theme(legend.position="right"),file=paste("h2_GM_",root,".png",sep=""),width=15,height=15)
# # ggsave(h2_all_plot ,file=paste("h2_all_",root,".png",sep=""),width=15,height=15)



# describe results

#' ## Subcortical structures
#' ### Heritability
opts_chunk$set(fig.width=7, fig.height=5 )
#' `r subset(h2,lobe=="SUBCORTICAL") %>% NROW` 
#' measures of the hippocampus and the thalamus were selected:
#' 
#'    - regional and tissue volume (left and right)
#'    - regional and tissue intensity (left and right)
#'    - regional T2* (left and right)
h2_subc_plot + theme(legend.position="bottom")
subset(h2,lobe=="SUBCORTICAL")  %>% select(-Lambda.GC,-h2.CI95upper,-h2.CI95lower,-lobe) %>% datatable(rownames=FALSE,caption="Heritability estimates of subcortical structures.")

min<- subset(h2,lobe=="SUBCORTICAL") %>% pull(h2) %>% min
max<- subset(h2,lobe=="SUBCORTICAL") %>% pull(h2) %>% max

#' `r (subset(h2,lobe=="SUBCORTICAL") %>% NROW) - (subset(h2,lobe=="SUBCORTICAL") %>% filter(p>0.05) %>% NROW)` measures were significantly heritable, with heritabily (h2) values ranging between 
#'
#'    - min: `r min`
#' (`r paste("_",paste(subset(h2,lobe=="SUBCORTICAL"&h2==min)$IDP_short_name, collapse="_ and _"),"_",sep="")`)
#'
#'    - max: `r max`
#' (`r paste("_",paste(subset(h2,lobe=="SUBCORTICAL"&h2==max)$IDP_short_name, collapse="_ and _"),"_",sep="")`)
#'
#' Just `r subset(h2,lobe=="SUBCORTICAL"&p>(0.05/NROW(h2))) %>% NROW()` measure does not survive Bonferroni multiple comparison correction (i.e. p>0.05/`r NROW(h2)`):
#' 
#' `r subset(h2,lobe=="SUBCORTICAL"&p>(0.05/NROW(h2)))  %>% select(-Lambda.GC,-h2.CI95upper,-h2.CI95lower,-lobe) %>% kable(row.names=FALSE)`

rm(min,max)


#' ### Genetic correlation: L/R
rg_subc_plot

#' Most measures have genetic correlations of 1 between the left and right. 
#' There is just one for which the estimate is significantly different to 1 (nominal, p2<0.05)
subset(rg,lobe=="SUBCORTICAL") %>% filter(p2<0.05) %>%  
  select(-rg.CI95upper,-rg.CI95lower,-lobe,-L,-R,-p) %>%  
  datatable(rownames=FALSE,caption="Genetic correlation estimates for L and R sides of subcortical structures.")


#' ## Grey matter regions
#' ### Heritability
opts_chunk$set(fig.width=15, fig.height=20 )
#' `r subset(h2,lobe!="SUBCORTICAL"&lobe!="dMRI") %>% NROW` 
#' measures of regional grey matter cortical measures were included.
#' 
subset(h2,lobe!="SUBCORTICAL"&lobe!="dMRI") %>% pull(measure) %>% table() %>% data.frame() %>% filter(Freq!=0) %>% kable(caption="Summary of cortical regional measure types.")

h2_GM_plot + theme(legend.position="bottom")
subset(h2,lobe!="SUBCORTICAL"&lobe!="dMRI")  %>% select(-Lambda.GC,-h2.CI95upper,-h2.CI95lower,-lobe) %>% datatable(rownames=FALSE,caption="Heritability estimates of grey matter regions.")


min<- subset(h2,lobe!="SUBCORTICAL"&lobe!="dMRI"&p<0.05) %>% pull(h2) %>% min
max<- subset(h2,lobe!="SUBCORTICAL"&lobe!="dMRI"&p<0.05) %>% pull(h2) %>% max

#' `r (subset(h2,lobe!="SUBCORTICAL"&lobe!="dMRI") %>% NROW) - (subset(h2,lobe!="SUBCORTICAL"&lobe!="dMRI") %>% filter(p>0.05) %>% NROW)` 
#' measures were (nominally, p<0.05) significantly heritable, with heritabily (h2) values ranging between 
#'
#'    - min: `r min`
#' (`r paste("_",paste(subset(h2,lobe!="SUBCORTICAL"&lobe!="dMRI"&h2==min)$IDP_short_name, collapse="_ and _"),"_",sep="")`)
#'
#'    - max: `r max`
#' (`r paste("_",paste(subset(h2,lobe!="SUBCORTICAL"&lobe!="dMRI"&h2==max)$IDP_short_name, collapse="_ and _"),"_",sep="")`)
#'
#' `r subset(h2,lobe!="SUBCORTICAL"&lobe!="dMRI"&p>(0.05/NROW(h2))) %>% NROW()` measures do not survive Bonferroni multiple comparison correction (i.e. p>0.05/`r NROW(h2)`):
#' 
#' `r subset(h2,lobe!="SUBCORTICAL"&lobe!="dMRI"&p>(0.05/NROW(h2)))  %>% select(-Lambda.GC,-h2.CI95upper,-h2.CI95lower,-lobe) %>% kable(row.names=FALSE)`
#' 
#' 
#' ### Genetic correlation: L/R
rg_GM_plot
#' Most measures have genetic correlations of 1 between the left and right. 
#' There are a few for which the estimate is significantly different to 1 (nominal, p2<0.05)
subset(rg,lobe!="SUBCORTICAL"&lobe!="dMRI") %>% filter(p2<0.05) %>% 
  select(-rg.CI95upper,-rg.CI95lower,-lobe,-L,-R,-p) %>%  
  datatable(rownames=FALSE,caption="Genetic correlation estimates for L and R sides of grey matter regions.")

#' 
#' ## White matter tracts
#' ### Heritability
opts_chunk$set(fig.width=15, fig.height=20 )
#' `r subset(h2,lobe=="dMRI")  %>% NROW` 
#' measures of white matter tract measures were included.
subset(h2,lobe=="dMRI") %>% pull(measure) %>% table() %>% data.frame() %>% filter(Freq!=0) %>% kable(caption="Summary of dMRI measure types.")

h2_WM_plot2
subset(h2,lobe=="dMRI")  %>% select(-Lambda.GC,-h2.CI95upper,-h2.CI95lower,-lobe) %>% datatable(rownames=FALSE,caption="Heritability estimates of white matter tract measures.")


min<- subset(h2,lobe=="dMRI"&p<0.05) %>% pull(h2) %>% min
max<- subset(h2,lobe=="dMRI"&p<0.05) %>% pull(h2) %>% max

#' `r (subset(h2,lobe=="dMRI") %>% NROW) - (subset(h2,lobe=="dMRI") %>% filter(p>0.05) %>% NROW)` 
#' measures were (nominally, p<0.05) significantly heritable, with heritabily (h2) values ranging between 
#'
#'    - min: `r min`
#' (`r paste("_",paste(subset(h2,lobe=="dMRI"&h2==min)$IDP_short_name, collapse="_ and _"),"_",sep="")`)
#'
#'    - max: `r max`
#' (`r paste("_",paste(subset(h2,lobe=="dMRI"&h2==max)$IDP_short_name, collapse="_ and _"),"_",sep="")`)
#'
#' `r subset(h2,lobe=="dMRI"&p>(0.05/NROW(h2))) %>% NROW()` measure does not survive Bonferroni multiple comparison correction (i.e. p>0.05/`r NROW(h2)`):
#' 
#' `r subset(h2,lobe=="dMRI"&p>(0.05/NROW(h2)))  %>% select(-Lambda.GC,-h2.CI95upper,-h2.CI95lower,-lobe) %>% kable(row.names=FALSE)`
#' 
#' 
#' ### Genetic correlation: L/R
rg_WM_plot
#' Most measures have genetic correlations of 1 between the left and right. 
#' There are a few for which the estimate is significantly different to 1 (nominal, p2<0.05)
subset(rg,lobe=="dMRI") %>% filter(p2<0.05) %>% 
  select(-rg.CI95upper,-rg.CI95lower,-lobe,-L,-R,-p) %>%  
  datatable(rownames=FALSE,caption="Genetic correlation estimates for L and R sides of white matter tract measures.")

#' ## Summary
#' 
#' Overall, most measure are heritable, and there seems to be no a clear difference between the different atlas definitions for a given region.
#' 
#' In order to select the regions, we should therefore focus on a relevant parcellation (given arguments more related to the neurobiological relevance rather than the power derived from this heritability analysis).
#' 


