#' ---
#' title: 'Genetic correlation of cerebellar structures and cortical volumes'
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
#' * cortical measures related to language
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
library(ggseg)
#
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(gghighlight)
library(RColorBrewer)
library(MetBrewer)
library(dplyr)
library(tidyr)
library(DT)
#---------------------------------------------------
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
#---------------------------------------------------

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
idps<-read.csv(paste(primary_dir,"IDPs_summary.csv",sep=""),header=TRUE) %>% filter(!is.na(Cerebellum))

# define order of regions:
cer_regions<-c("TOTAL","CORTEX","WHITE-MATTER",
               "I-IV","V","VI","CRUS I","CRUS II",
               "VIIB","VIIIA","VIIIB","IX","X")
# read estimates: rg across FAST atlas measures
# and subset only cerebellum measures
#
rg_cort1<-read.csv(paste0(ldsc_dir,"ldsc_rg_atlas_UKB_BIG40.csv")) %>% 
  filter(region.p1 %in% idps$region) %>%
  filter(region.p2 %in% idps$region) %>%
  filter(parcellation=="FAST")

# fill in info for total cerebellum volume, from Chambers et al. 2022
rg_cort3<-read.csv(paste0(ldsc_dir,"ldsc_rg_other_UKB_BIG40.csv"))
rg_cort3 <- rg_cort3 %>%
  mutate(measure.p1=if_else(p1=="Chambers2022_TotalCerebellarVolume","VOLUME",measure.p1),
         region.p1=if_else(p1=="Chambers2022_TotalCerebellarVolume","TotalCerebellarVolume",region.p1)
  ) %>%
  mutate(measure=if_else(measure.p1==measure.p2,measure.p1,""))
# combine all
rg_cort<-merge(rg_cort1,rg_cort3,all=TRUE)
rm(rg_cort1,rg_cort3)

# keep only rgs where first p is cerebellar and second is not
w1<-grep("CEREBELL",toupper(rg_cort$region.p1))
w2<-grep("CEREBELL",toupper(rg_cort$region.p2),invert=TRUE)
w0<-intersect(w1,w2) # both
w3<-grep("ACCUM|HIPPO|CAUDATE|STEM|AMYG|PALLID|THALAMUS|STRIATUM|PUTAMEN|PLANUM_TEMP|ANGULAR_G",toupper(rg_cort$region.p2),invert=TRUE)
w0<-intersect(w0,w3)
rg_cort<-rg_cort[w0,]
rm(w1,w2,w0)

# edit some columns
rg_cort <- rg_cort %>% 
  mutate(hemisphere=if_else(grepl("^V_",region.p1),"V",hemisphere.p1),
         region=gsub("^V_","",region.p1) %>% 
           gsub("CEREBELLUM_|CEREBELLUM-|CerebellarVolume","",.) %>% 
           gsub("_"," ",.) %>% toupper ) %>%
  mutate(hemisphere=if_else(is.na(hemisphere),"BiLat",hemisphere)) %>%
  mutate(hemisphere=factor(hemisphere,levels=c("BiLat","V","L","R")),
         type=if_else((region=="CORTEX"|region=="WHITE-MATTER"|region=="TOTAL"),"Global","Region"),
         region=factor(region,levels=cer_regions),
         region.cortical=gsub("_|-"," ",region.p2) %>% tolower() %>% sapply(simpleCap)
         ) 
rg_cort<-rg_cort %>% mutate(
  region.cortical=gsub("Gyrus","\nGyrus",region.cortical) %>% gsub("Cortex","\nCortex",.) 
  )


# define order of cortical measures, by lobe...
order_cortical<-sort(unique(rg_cort$region.cortical))
order_cortical<-order_cortical[c(1:2, # IFG
                                 10:11, # SUPRAMARGINAL
                                 7, # PRECUNEUS
                                 8:9, # STG
                                 3:5, # MTG
                                 12:14,6 # Fusif. Gyrus, Temp and Occ
)]

# define sig after multiple testing correction
rg_cort <- rg_cort %>%
  mutate(sig=if_else(p<(0.05/(ntests*length(order_cortical))),"*",""))


rg_cort<-rg_cort %>% mutate(region.cortical=factor(region.cortical,levels=order_cortical))
# save files
write.csv(rg_cort,file=paste(out_dir,"/ldsc_rg_cortical","_",root,"_",project,".csv",sep=""),row.names = FALSE)



#----------------------------------------------------------------------
#' ### Cerebellar volumes
#----------------------------------------------------------------------
#' These cerebellar volumes were selected for the genetic correlation analysis, and there were three measures per region (left, right and vermis):
#' 
#' * VI
#     - Associated with reading (fMRI): RH
#' * Crus I
#     - Associated with reading (fMRI): RH


#----------------------------------------------------------------------
#' ### Cortical measures
#----------------------------------------------------------------------
# plot(hoCort)

# associate color to region
target_labels<- ggseg(atlas="hoCort")$data %>% select(label) %>% distinct() %>%
  mutate(matchLabel=tolower(label) %>%
           gsub("\\.|-","_",.) %>%
           gsub("__","_",.) %>%
           gsub("_division|_part","",.)
         )

rg_cort<-rg_cort %>% mutate(hemi.p2=gsub("L","lh",hemisphere.p2) %>% gsub("R","rh",.)) %>%
  mutate(matchLabel=paste(hemi.p2,region.p2,sep="_") %>% 
           gsub("area_|thick_|smri_|cort.","",.) %>% 
           gsub(paste0("FAST","_"),"",.) %>%
           gsub("\\.|-","_",.) %>% tolower() %>%
           gsub("occ_","occipital_",. ) %>%
           gsub("temp_","temporal_",.) %>%
           gsub("front_","frontal_",.) %>%
           gsub("inf_","inferior_",.) %>%
           gsub("sup_","superior_",.) %>%
           gsub("mid_","middle_",.) %>%
           gsub("_ant","_anterior",.) %>%
           gsub("_post","_posterior",.) %>%
           gsub("pars","pars_",.) %>%
           gsub("tri$","triangularis",.) %>%
           gsub("op$","opercularis",.) %>%
           gsub("orb$","orbitalis",.) %>%
           gsub("fusif_","fusiform_",.) %>%
           gsub("supramarg_","supramarginal_",.) %>%
           gsub("precun_","precuneous_",.) %>%
           gsub("tempocc","temporooccipital",.)
        )
# add color label
rg_cort_colors<-rg_cort %>% dplyr::select(region.cortical) %>% distinct() %>%
  mutate(color=
           met.brewer(name="Signac",n=length(unique(rg_cort_colors$region.cortical)))
  )

rg_cort<-merge(rg_cort,rg_cort_colors,all=TRUE)

# w<-which(rg_cort$matchLabel %in% target_labels$matchLabel)
# rg_cort$matchLabel[-w]

rg_cort<-merge(rg_cort,target_labels)
#' Cortical regions included in analysis (language-related):
rg_cort_network<-rg_cort %>% select(label,region.p2,region.cortical,color) %>% distinct() %>% brain_join(hoCort) %>%
  mutate(region.cortical=factor(region.cortical,levels=order_cortical)) %>%
  ggplot() +
  geom_sf(aes(fill=region.cortical),color="lightgrey",size=0.5) +
  # scale_fill_brain() +
  # scale_fill_viridis_d(option = "inferno",na.value="transparent") +
  # scale_fill_manual(na.value="transparent") +
  # scale_fill_met_d(name="Signac") +
  scale_fill_manual(values=unique(rg_cort$color),na.value="transparent",
                    labels=c(gsub("\n","",levels(rg_cort$region.cortical)),"")) +

  theme_void() +
  theme(legend.position="bottom",legend.text=element_text(size=12)) +
  guides(fill=guide_legend(ncol=3,title="",override.aes = list(color="white"))) +
  # labs(title="Cortical volumes") +
  # theme_darkbrain() +
  NULL

rg_cort_network

rg_cort_network2<-rg_cort_network + theme(legend.position="left") +
  facet_grid(hemi~.) +
  guides(fill=guide_legend(ncol=1,title="",override.aes = list(color="white")))

rg_cort_network %>% ggsave(file=paste0(out_dir,"../atlas_figures/HarvardOxford_vols_langNetwork.png"),height=6,width=8)

rg_cort_network2 %>%  ggsave(file=paste0(out_dir,"../atlas_figures/HarvardOxford_vols_langNetwork_2.png"),height=4,width=10)

#----------------------------------------------------------------------
#' ## Results
#----------------------------------------------------------------------
#' ### All results
#+ fig.width=12, fig.height=10, fig.align="center"
rg_cort_plot 

rg_cort %>%
  arrange(p) %>%
  rename(Cerebellar.Region=region,
         Cerebellar.Hemis=hemisphere,
         Cortical.Region=region.cortical,
         Corticall.Hemis=hemisphere.p2) %>%
  select(contains("Cerebellar"),contains("Cortical"),rg,rg.SE,p) %>%
  mutate(p=format(p,scientific=TRUE,digits=2)) %>%
  mutate(across(where(is.numeric),round,digits=3)) %>%
  datatable(rownames=FALSE,caption="Table of all genetic correlations with cortical volumes (FAST).")

#----------------------------------------------------------------------
# generate plots to visualize genetic correlations
rg_cort_plot <- rg_cort %>% 
  mutate(region=factor(region,levels=rev(cer_regions))) %>%
  mutate(sig_val=if_else(rg<0,rg.CI95lower*1.5,rg.CI95upper*1.5)) %>%
  ggplot(aes(x=region,y=rg,width=0.7,color=hemisphere,fill=hemisphere,shape=hemisphere.p2,alpha=sig)) +
  geom_hline(yintercept = 0,color="black",alpha=0.7,linetype="dashed") +
  geom_errorbar(aes(ymin=rg.CI95upper, ymax=rg.CI95lower,width=.5), position=position_dodge(.7)) +
  geom_point(position=position_dodge(.7), stat="identity",size=3,color="white") +
  geom_text(aes(x=region,y=sig_val,label=sig,color=hemisphere),alpha=0.8,
            fontface=2,size=6,position=position_dodge(.8),stat="identity") +
  labs(color="Cerebellar hemisphere",
       fill="Cerebellar hemisphere",
       shape="Cortical hemisphere",
       y=bquote('Estimate ('*rho*')'),
       x=""
       # title="Genetic correlation with cortical language network"
  ) +
  facet_wrap(region.cortical ~.,ncol=7) +
  coord_flip() +
  # scale_color_met_d(name="Morgenstern") +
  scale_fill_viridis_d() + scale_color_viridis_d() +
  scale_alpha_manual(values=c(0.3,0.9)) +
  scale_shape_manual(values=c(21,24,22)) +
  guides(shape = guide_legend(override.aes = list(color="black") ),
         alpha="none"
  ) +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="transparent",color="grey85",size=1),
        legend.position = "bottom") +
  panel_border(size = 1) +
  NULL
rg_cort_plot %>% 
  ggsave(.,file=paste0(out_dir,"/Figure_ldsc_rg_cortical.png"),width=9.5,height=7)

# now add color to facets to reflect which cortical region they reflect
## following: https://stackoverflow.com/questions/41631806/change-facet-label-text-and-background-colour/60046113#60046113

g <- ggplot_gtable(ggplot_build(rg_cort_plot))
stript <- which(grepl('strip-t', g$layout$name))
n<-length(unique(rg_cort_colors$color))
fills <- c(unique(rg_cort_colors$color)[((n/2)+1):n],
           unique(rg_cort_colors$color)[1:(n/2)]
)
           
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid::grid.newpage()
grid::grid.draw(g)
rm(stript,fills,k)

ggsave(g,file=paste0(out_dir,"/Figure_ldsc_rg_cortical_2.png"),width=9.5,height=7)

# rg_cort_plot_combined<-plot_grid(rg_cort_network2,grid::grid.draw(g),nrow=1,labels=c("A","B"))
# ggsave(rg_cort_plot_combined,file=paste0(out_dir,"/Figure_ldsc_rg_cortical_combined.png"),width=15,height=7)


#--------------------------------------------




#' ### Significant after multiple comparison correction

#' After correcting for multiple testing comparison (i.e. p-value<0.05/`r (ntests*length(order_cortical)))` = `r 0.05/(ntests*length(order_cortical)))`)
#' the positive genetic correlations between 
#' *Crus I/lobule VI* and 
#' temporo-occipital volumes (*Temp Occ Fusif Cortex/Occ Fusif Gyrus*) are significant.
rg_cort %>%
  arrange(p) %>% filter(p<0.05/(ntests*length(order_cortical))) %>%
  rename(Cerebellar.Region=region,
         Cerebellar.Hemis=hemisphere,
         Cortical.Region=region.cortical,
         Corticall.Hemis=hemisphere.p2) %>%
  select(contains("Cerebellar"),contains("Cortical"),rg,rg.SE,p) %>%
  mutate(p=format(p,scientific=TRUE,digits=2)) %>%
  mutate(across(where(is.numeric),round,digits=3)) %>%
  datatable(rownames=FALSE,caption="Table of significant (p<0.05/336) genetic correlations with cortical volumes (FAST).")


#' ### Crus I overview
roi="CRUS I"

p1<-rg_cort %>% filter(hemisphere=="R"&region==roi) %>% select(rg,p,label) %>% 
  brain_join(hoCort) %>% 
  # reposition_brain(hemi ~ side) %>%
  ggplot() +
  geom_sf(aes(fill=rg),colour="black") +
  scale_fill_distiller(type="div",limit=c(-0.5,0.5)) +
  # geom_sf_label(aes(label=ifelse(p<0.05,region,NA)),alpha=0.8) +
  labs(subtitle="Right") +
  theme_void() +
  theme(legend.position="none") +
  NULL
p2<-rg_cort %>% filter(hemisphere=="L"&region==roi) %>% select(rg,p,label) %>% 
  brain_join(hoCort) %>% 
  # reposition_brain(hemi ~ side) %>%
  ggplot() +
  geom_sf(aes(fill=rg),colour="black") +
  scale_fill_distiller(type="div",limit=c(-0.5,0.5)) +
  # geom_sf_label(aes(label=ifelse(p<0.05,region,NA)),alpha=0.8) +
  labs(subtitle="Left") +
  theme_void() +
  theme(legend.position="none") +
  NULL
p3<-rg_cort %>% filter(hemisphere=="V"&region==roi) %>% select(rg,p,label) %>% 
  brain_join(hoCort) %>% 
  # reposition_brain(hemi ~ side) %>%
  ggplot() +
  geom_sf(aes(fill=rg),colour="black") +
  scale_fill_distiller(type="div",limit=c(-0.5,0.5)) +
  # geom_sf_label(aes(label=ifelse(p<0.05,region,NA)),alpha=0.8) +
  labs(subtitle="Vermis") +
  theme_void() +
  theme(legend.position="none") +
  NULL
leg<-get_legend(p1 +  theme(legend.position="right"))
region_plot<-plot_grid(plot_grid(p1,p2,p3,ncol=1),leg,ncol=2,labels=roi,rel_widths=c(1,0.1))
ggsave(region_plot,file=paste0(out_dir,"/cereb_",gsub(" ","_",roi),"_cortical_rg.pdf"),height=6,width=8)

#+ fig.width=14, fig.height=10, fig.align="center"
region_plot
# clean intermediate
rm(p1,p2,p3,roi)


#' ### Lobule VI overview
roi="VI"

p1<-rg_cort %>% filter(hemisphere=="R"&region==roi) %>% select(rg,p,label) %>% 
  brain_join(hoCort) %>% 
  # reposition_brain(hemi ~ side) %>%
  ggplot() +
  geom_sf(aes(fill=rg),colour="black") +
  scale_fill_distiller(type="div",limit=c(-0.5,0.5)) +
  # geom_sf_label(aes(label=ifelse(p<0.05,region,NA)),alpha=0.8) +
  labs(subtitle="Right") +
  theme_void() +
  theme(legend.position="none") +
  NULL
p2<-rg_cort %>% filter(hemisphere=="L"&region==roi) %>% select(rg,p,label) %>% 
  brain_join(hoCort) %>% 
  # reposition_brain(hemi ~ side) %>%
  ggplot() +
  geom_sf(aes(fill=rg),colour="black") +
  scale_fill_distiller(type="div",limit=c(-0.5,0.5)) +
  # geom_sf_label(aes(label=ifelse(p<0.05,region,NA)),alpha=0.8) +
  labs(subtitle="Left") +
  theme_void() +
  theme(legend.position="none") +
  NULL
p3<-rg_cort %>% filter(hemisphere=="V"&region==roi) %>% select(rg,p,label) %>% 
  brain_join(hoCort) %>% 
  # reposition_brain(hemi ~ side) %>%
  ggplot() +
  geom_sf(aes(fill=rg),colour="black") +
  scale_fill_distiller(type="div",limit=c(-0.5,0.5)) +
  # geom_sf_label(aes(label=ifelse(p<0.05,region,NA)),alpha=0.8) +
  labs(subtitle="Vermis") +
  theme_void() +
  theme(legend.position="none") +
  NULL
leg<-get_legend(p1 +  theme(legend.position="right"))
region_plot<-plot_grid(plot_grid(p1,p2,p3,ncol=1),leg,ncol=2,labels=roi,rel_widths=c(1,0.1))
ggsave(region_plot,file=paste0(out_dir,"/cereb_",gsub(" ","_",roi),"_cortical_rg.pdf"),height=6,width=8)

#+ fig.width=14, fig.height=10, fig.align="center"
region_plot
# clean intermediate
rm(p1,p2,p3,roi)

