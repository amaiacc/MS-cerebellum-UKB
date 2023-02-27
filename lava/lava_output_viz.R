library(dplyr)
library(tidyr)
# library(circlize) 
# library(chorddiag)  #devtools::install_github("mattflor/chorddiag")
library(ggplot2)
library(cowplot)
library(viridis); library(MetBrewer)
# clean environment
rm(list=ls())
gc()
#------------------------------------
# Set directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
#------------------------------------
# some functions
source(paste(dir,"general_scripts/helper_functions.R",sep=""))
# working dir
working_dir=paste0(dir,"cerebellum_UKB/data/lava/output/")
setwd(working_dir)
#------------------------------------
# read bivar file with all results
bivar<-read.csv("lava_bivar_all_traits.csv") %>% mutate(end=stop)
bivar_wide<-read.csv("lava_bivar_wide_all_traits.csv") #%>%
  # mutate(chr=strsplit())

# read locus file
loci<-read.table("../blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile",header=TRUE) %>%
  rename(chr=CHR,start=START,end=STOP,loc=LOC) %>%
  dplyr::select(chr,start,end,loc) %>%
  mutate(locus_id=paste0(chr,":",start,"-",end))

# define genome
genome<-data.frame(chr=sort(unique(loci$chr)) ) %>%
  mutate( 
    start=0,
    min= sapply(chr,function(x){subset(loci,chr==x) %>% pull(start) %>% min() } ),
    max= sapply(chr,function(x){subset(loci,chr==x) %>% pull(end) %>% max() } )
)
#------------------------------------
# visualize!
n<-NROW(bivar)
# edit names, and make new variables to be used in the plot
bivar2<-bivar %>% mutate(trait1=if_else(trait1=="CEREBELLUM_CRUS_I_R" & trait2=="CEREBELLUM_VI_R",
                        "CRUS_I_R_VI_R",trait1) %>% 
                          gsub("CEREBELLUM_","",.) %>%
                          gsub("Chambers2022_","",.) %>%
                          gsub("TotalCerebellarVolume","Total Cerebellar Volume",.) %>%
                          gsub("_"," ",.) %>%
                          gsub(" L"," (L)",.) %>% gsub(" R"," (R)",.) %>%
                          gsub("CRUS","Crus",.) %>% gsub("VI","Lobule VI",.),
                        trait2=if_else( (trait1=="Crus I (R) Lobule VI (R)"),
                                        "Crus I (R) Lobule VI (R)",trait2) %>%
                          gsub("CEREBELLUM_","",.) %>%
                          gsub("Chambers2022_","",.) %>%
                          gsub("TotalCerebellarVolume","Total Cerebellar Volume",.) %>%
                          gsub("Grove_|Lee2018_|PGC3_|_primary","",.) %>%
                          gsub("_"," ",.) %>%
                          gsub("GYRUS","\nGYRUS",.) %>% gsub("CORTEX","\nCORTEX",.) %>%
                          tolower() %>% sapply(.,simpleCap) %>%
                          gsub(" L$"," (L)",.) %>% 
                          gsub(" R$"," (R)",.) %>% gsub(" \\(r\\)"," (R)",.) %>%
                          gsub("Vi","VI",.) %>% #gsub("VI","Lobule VI",.) %>%
                          gsub("Asd","ASD",.) %>%
                          gsub("Cp","CP",.) %>% gsub("Scz","SCZ",.)
                        ) %>%
  mutate(sig1=if_else(p.adjFDR<0.05,"FDR adjusted p-value < 0.05","FDR adjusted p-value > 0.05"),
         sig2=if_else(p.adjBF<0.05,"Bonferroni adjusted p-value < 0.05","Bonferroni adjusted p-value > 0.05"),
         sig3=if_else(p.adjBF<0.05,paste0("p-value < 0.05/",n),if_else(p<0.05,"p-value < 0.05","p-value > 0.05"))
         ) %>%
  mutate(sig3=factor(sig3,levels=c(paste0("p-value < 0.05/",n),"p-value < 0.05","p-value > 0.05")) )

# to ensure same x-axis for all plots **try to ensure** --> all loci are present
loci$trait1<-NA
bivar2 <- merge(bivar2,loci,all=TRUE,
      by.x=c("locus","chr","start","stop",
             "locus_id","trait1"),
      by.y=c("loc","chr","start","end",
             "locus_id","trait1"))

# levels(bivar2$locus)<-loci$loc
#------------------------------------
## plots
#------------------------------------
a<-bivar2 %>% filter(trait1=="Crus I (R) Lobule VI (R)" | 
                       (trait1=="Total Cerebellar Volume" & (trait2=="Crus I (R)"|trait2=="VI (R)")) |
                       is.na(trait1) ) %>%
  # placeholders, to keep loci
    mutate(trait1=if_else(is.na(trait1),"",trait1),
           trait2=if_else(is.na(trait2),"",trait2),
           rho=if_else(is.na(rho),-100,rho),
           rho.lower=if_else(is.na(rho.lower),-99,rho.lower),
           rho.upper=if_else(is.na(rho.upper),-101,rho.upper),
           sig3=if_else(is.na(sig3),"",as.character(sig3))
         ) %>%
  mutate(sig3=factor(sig3,levels=c(paste0("p-value < 0.05/",n),"p-value < 0.05","p-value > 0.05","")),
         trait1=factor(trait1,levels=c("Total Cerebellar Volume","Crus I (R) Lobule VI (R)","") ),
         trait2=factor(trait2,levels=c("Crus I (R)","VI (R)","Crus I (R) Lobule VI (R)","") )
  )
  # mutate(sig3=factor(sig3,levels=c(paste0("p-value < 0.05/",n),"p-value < 0.05","p-value > 0.05","")),
  #        locus=factor(locus))
pA <- a %>%
  ggplot(., aes(x=start,y=rho,alpha=sig3),color="black") +
  geom_hline(yintercept = c(-1,0,1),color="grey85") +
  geom_errorbar(aes(ymin=rho.lower,ymax=rho.upper,color=trait1)) +
  geom_point(size=3,color="white",
             aes(fill=trait1,
                 shape=trait2)
  ) +
  facet_grid(.~chr,scales="free_x",space = "free_x") +
  labs(x="loci") +
  # scale_alpha_manual(values=c(1,0.2,0.05,0)) +
  scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_manual(values=c(viridis(3)[c(1,3)],"#ffffff")) +
  scale_color_manual(values=c(viridis(3)[c(1,3)],"#ffffff")) +
  scale_shape_manual(values=c(21,24,23,22),
                     labels=c("Total Cerebellar Volume  ~ Crus I (R)","Total Cerebellar Volume ~ Lobule VI (R)","Crus I (R) ~ Lobule VI (R)","")
  ) +
  scale_x_discrete(
    breaks = unique(bivar2$locus),
    labels = element_blank(),
    expand = c(0,0)
  ) +
  coord_cartesian(ylim=c(-1,1)) +
  guides(fill="none",color="none",
         alpha=guide_legend(override.aes = list(color="grey85",fill="black", shape=18),title="",order=1),
         shape=guide_legend(override.aes = list(color="white",
                                                fill=c(viridis(3)[1],viridis(3)[1],viridis(3)[3],"#ffffff")
         ), title="",order=2)
  ) +
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="grey85",size=1),
        legend.position = "none") +
  panel_border() +
  NULL
# rm(a)
#------------------------------------
# define colors to match HO atlas of language network (Figure Sx)
colors_fusif<-met.brewer(name="Signac",n=14)[c(14,13)]
b <- bivar2 %>% filter(grepl("Fusif",trait2) | is.na(trait1) ) %>%
  # placeholders, to keep loci
  mutate(trait1=if_else(is.na(trait1),"",trait1),
         trait2=if_else(is.na(trait2),"",trait2),
         rho=if_else(is.na(rho),-100,rho),
         rho.lower=if_else(is.na(rho.lower),-99,rho.lower),
         rho.upper=if_else(is.na(rho.upper),-101,rho.upper),
         sig3=if_else(is.na(sig3),"",as.character(sig3))
  ) %>%
  mutate(sig3=factor(sig3,levels=c(paste0("p-value < 0.05/",n),"p-value < 0.05","p-value > 0.05","")),
         locus=factor(locus)) %>%
  mutate(
    trait1=factor(trait1,levels=c("Crus I (R)","Lobule VI (R)","Total Cerebellar Volume","")),
    trait2=factor(trait2,levels=c("Occ Fusif \ngyrus (L)","Temp Occ Fusif \ncortex (L)",""))
  )
pB <- b %>%
  ggplot(aes(x=start,y=rho,
             color=trait2,
             fill=trait2,
             shape=trait1,
             alpha=sig3)) +
  geom_hline(yintercept = c(-1,0,1),color="grey85") +
  geom_errorbar(aes(ymin=rho.lower,ymax=rho.upper)) +
  geom_point(size=3,color="white") +
  facet_grid(.~chr,scales="free_x",space = "free_x") +
  labs(x="loci") +
  # scale_alpha_manual(values=c(1,0.2,0.05,0)) +
  scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_manual(values=c(colors_fusif,"#ffffff")) +
  scale_color_manual(values=c(colors_fusif,"#ffffff")) +
  scale_shape_manual(values=c(21,24,22,23)) +
  scale_x_discrete(
    breaks = unique(bivar2$locus),
    labels = element_blank(),
    expand = c(0,0)
  ) +
  coord_cartesian(ylim=c(-1,1)) +
  guides(alpha="none", color="none",
         fill=guide_legend(override.aes = list(color=c(colors_fusif,"#ffffff"), shape=18),order=2, title="Cortical volume"),
         shape=guide_legend(override.aes = list(color=c("black","black","black","white")),order=1, title="Cerebellar volume"))+
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="grey85",size=1),
        legend.position = "none") +
  panel_border() +
  NULL
# rm(b)
#------------------------------------
c <- bivar2 %>% filter(grepl("CP|SCZ|ASD",trait2) | is.na(trait1) ) %>%
  # placeholders, to keep loci
  mutate(trait1=if_else(is.na(trait1),"",trait1),
         trait2=if_else(is.na(trait2),"",trait2),
         rho=if_else(is.na(rho),-100,rho),
         rho.lower=if_else(is.na(rho.lower),-99,rho.lower),
         rho.upper=if_else(is.na(rho.upper),-101,rho.upper),
         sig3=if_else(is.na(sig3),"",as.character(sig3))
  ) %>%
  mutate(sig3=factor(sig3,levels=c(paste0("p-value < 0.05/",n),"p-value < 0.05","p-value > 0.05","")),
         locus=factor(locus)) %>%
  mutate(
    trait1=factor(trait1,levels=c("Crus I (R)","Lobule VI (R)","Total Cerebellar Volume","")),
    trait2=factor(trait2,levels=c("ASD","SCZ","CP",""))
    )

pC <- c %>%
  ggplot(aes(x=start,y=rho,
             color=trait2,
             fill=trait2,
             shape=trait1,
             alpha=sig3)) +
  geom_hline(yintercept = c(-1,0,1),color="grey85") +
  geom_errorbar(aes(ymin=rho.lower,ymax=rho.upper)) +
  geom_point(size=3,color="white") +
  facet_grid(.~chr,scales="free_x",space = "free_x") +
  labs(x="loci") +
  # scale_alpha_manual(values=c(1,0.2,0.05,0)) +
  scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_met_d(name="Austria") +
  scale_color_met_d(name="Austria") +
  scale_shape_manual(values=c(21,24,22,23)) +
  scale_x_discrete(
    breaks = unique(bivar2$locus),
    labels = element_blank(),
    expand = c(0,0)
  ) +
  coord_cartesian(ylim=c(-1,1)) +
  guides(alpha="none", color="none",
         fill=guide_legend(override.aes = list(color=c(met.brewer(name="Austria",n=3),"#ffffff"), shape=18), title="Cognitive trait/\nPyschiatric disorder"),
         shape="none"
         # shape=guide_legend(override.aes = list(color="black"),order=1, title="Cerebellar volume")
         ) +
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="grey85",size=1),
        legend.position = "none") +
  panel_border() +
  NULL
#------------------------------------
# combine plots and get legends

p<-plot_grid(pA,
             pB,
             pC,#NULL,
             ncol=1,labels=c("A","B","C"))
l<-plot_grid(get_legend(pA+theme(legend.position="right")+guides(shape="none")),
             plot_grid(
               get_legend(pB+theme(legend.position="right")+guides(fill="none")),
               get_legend(pA+theme(legend.position="right")+guides(alpha="none")),
               ncol=1
             ),
             get_legend(pB+theme(legend.position="right")+guides(shape="none")),
             get_legend(pC+theme(legend.position="right")),
             nrow=1,
             rel_widths = c(1,1.4,1,1),
             align="v"
)
p2<-plot_grid(p,NULL,l,ncol=1,rel_heights = c(2,0.1,0.5))

ggsave(p2,file="Figure_rg_local_lava_bivar_all.png",width=12,height=10)
write.csv(bivar2,file="lava_bivar_all_traits_clean.csv",row.names=FALSE)

#------------------------------------
loci_lava<-bivar2 %>% filter(p.adjBF<0.05) %>% pull(locus_id) %>% unique()
a1<-a %>% filter(locus_id %in% loci_lava) %>% filter(rho!=(-100)) %>%
  mutate(locus_id=factor(locus_id,levels=loci_lava))
b1<-b %>% filter(locus_id %in% loci_lava) %>% filter(rho!=(-100)) %>%
  mutate(locus_id=factor(locus_id,levels=loci_lava))
c1<-c %>% filter(locus_id %in% loci_lava) %>% filter(rho!=(-100)) %>%
  mutate(locus_id=factor(locus_id,levels=loci_lava))

pA1 <- a1 %>%
  ggplot(., aes(x=locus_id,y=rho,alpha=sig3),color="black") +
  geom_hline(yintercept = c(-1,0,1),color="grey85") +
  geom_errorbar(aes(ymin=rho.lower,ymax=rho.upper,color=trait1)) +
  geom_point(size=3,color="white",
             aes(fill=trait1,
                 shape=trait2)
  ) +
  # facet_grid(.~chr,scales="free_x",space = "free_x") +
  labs(x="") +
  # scale_alpha_manual(values=c(1,0.2,0.05,0)) +
  scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_manual(values=c(viridis(3)[c(1,3)],"#ffffff")) +
  scale_color_manual(values=c(viridis(3)[c(1,3)],"#ffffff")) +
  scale_shape_manual(values=c(21,24,23,22),
                     labels=c("Total Cerebellar Volume  ~ Crus I (R)","Total Cerebellar Volume ~ Lobule VI (R)","Crus I (R) ~ Lobule VI (R)","")
  ) +
  scale_x_discrete(
    limits = rev(loci_lava),
    # labels = element_blank(),
    # expand = c(0,0),
    drop=FALSE
  ) +
  coord_flip() +
  guides(fill="none",color="none",
         alpha=guide_legend(override.aes = list(color="black",fill="black", shape=18),title="",order=1),
         shape=guide_legend(override.aes = list(color="white",
                                                fill=c(viridis(3)[1],viridis(3)[1],viridis(3)[3])
         ), title="",order=2)
  ) +
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="grey85",size=1),
        legend.position = "none") +
  panel_border(color="black") +
  NULL

pB1 <- b1 %>%
  ggplot(aes(x=locus_id,y=rho,
             color=trait2,
             fill=trait2,
             shape=trait1,
             alpha=sig3)) +
  geom_hline(yintercept = c(-1,0,1),color="grey85") +
  geom_errorbar(aes(ymin=rho.lower,ymax=rho.upper)) +
  geom_point(size=3,color="white") +
  # facet_grid(.~chr,scales="free_x",space = "free_x") +
  labs(x="") +
  # scale_alpha_manual(values=c(1,0.2,0.05,0)) +
  scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_manual(values=c(colors_fusif,"#ffffff")) +
  scale_color_manual(values=c(colors_fusif,"#ffffff")) +
  scale_shape_manual(values=c(21,24,22,23)) +
  scale_x_discrete(
    limits = rev(loci_lava),
    labels = element_blank(),
    # expand = c(0,0),
    drop=FALSE
  ) +
  coord_flip() +
  guides(alpha="none", color="none",
         fill=guide_legend(override.aes = list(color=c(colors_fusif), shape=18),order=2, title="Cortical volume"),
         shape=guide_legend(override.aes = list(color="black"),order=1, title="Cerebellar volume"))+
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="grey85",size=1),
        legend.position = "none") +
  panel_border(color="black") +
  NULL

pC1 <- c1 %>%
  ggplot(aes(x=locus_id,y=rho,
             color=trait2,
             fill=trait2,
             shape=trait1,
             alpha=sig3)) +
  geom_hline(yintercept = c(-1,0,1),color="grey85") +
  geom_errorbar(aes(ymin=rho.lower,ymax=rho.upper)) +
  geom_point(size=3,color="white") +
  # facet_grid(.~chr,scales="free_x",space = "free_x") +
  labs(x="") +
  # scale_alpha_manual(values=c(1,0.2,0.05,0)) +
  scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_met_d(name="Austria") +
  scale_color_met_d(name="Austria") +
  scale_shape_manual(values=c(21,24,22,23)) +
  scale_x_discrete(
    limits = rev(loci_lava),
    # labels = element_blank(),
    # expand = c(0,0),
    position="top",
    drop=FALSE
  ) +
  coord_flip() +
  guides(alpha="none", color="none",
         fill=guide_legend(override.aes = list(color=c(met.brewer(name="Austria",n=3)), shape=18), title="Cognitive trait/\nPyschiatric disorder"),
         shape="none"
         # shape=guide_legend(override.aes = list(color="black"),order=1, title="Cerebellar volume")
  ) +
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="grey85",size=1),
        legend.position = "none") +
  panel_border(color="black") +
  NULL

# combine plots and get legends

p1<-plot_grid(pA1,
             pB1,
             pC1,#NULL,
             nrow=1,
             rel_widths = c(2,1.2,2),
             labels=c("A","B","C"))
l1<-plot_grid(get_legend(pA1+theme(legend.position="right")+guides(shape="none")),
             # plot_grid(
               get_legend(pB1+theme(legend.position="right")+guides(fill="none")),
               get_legend(pA1+theme(legend.position="right")+guides(alpha="none")),
               # ncol=1
             # ),
             get_legend(pB1+theme(legend.position="right")+guides(shape="none")),
             get_legend(pC1+theme(legend.position="right")),
             nrow=1,
             rel_widths = c(1,1,1.4,1,1)
             # align="v"
)
p11<-plot_grid(p1,l1,ncol=1,rel_heights  = c(3,0.5))

ggsave(p11,file="Figure_rg_local_lava_bivar_BF.png",width=14,height=10)

# #------------------------------------
# # attempt to visualize using circlize --> not worth it I think...
# ## following:
# ## - https://www.royfrancis.com/beautiful-circos-plots-in-r/
# ## - https://www.data-to-viz.com/graph/chord.html
# circos.clear()
# col_text <- "grey40"
# circos.par("track.height"=0.8, gap.degree=5, cell.padding=c(0, 0, 0, 0))
# 
# circos.initialize(factors=genome$chr, 
#                   xlim=matrix(c(rep(0, NROW(genome)), genome$max), ncol=2))
# # chr backbone, customized
# circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
#   chr=CELL_META$sector.index
#   xlim=CELL_META$xlim
#   ylim=CELL_META$ylim
#   circos.text(mean(xlim), mean(ylim), chr, cex=0.5, col=col_text, 
#               facing="bending.inside", niceFacing=TRUE)
# }, bg.col="grey90", bg.border=F, track.height=0.06)
# 
# # # add x axis --> maybe should be each locus?
# # brk <-seq(from=0,to=max(genome$max),by=1000*10000) # in Mb
# # circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
# #   circos.axis(h="top", major.at=brk, labels=round(brk/10000, 1), labels.cex=0.4, 
# #               col=col_text, labels.col=col_text, lwd=0.7, labels.facing="clockwise")
# # }, bg.border=F)
# 
# # label loci --> only the ones included in the bivar analysis
# w<-which(loci$locus_id %in% bivar_wide$locus_id)
# circos.genomicLabels(loci[w,], 
#                      labels.column = 4, 
#                      cex=0.25, col=col_text, line_lwd=0.5, line_col="grey80", 
#                      side="inside", connection_height=0.05, labels_height=0.04)
# 
# #
# 
# # rg
# circos.track(factors=gc$chr, x=gc$x, y=gc$value, panel.fun=function(x, y) {
#   circos.lines(x, y, col="grey50", lwd=0.6)
#   circos.segments(x0=0, x1=max(ref$V2), y0=0.3, y1=0.3, lwd=0.6, lty="11", col="grey90")
#   circos.segments(x0=0, x1=max(ref$V2), y0=0.5, y1=0.5, lwd=0.6, lty="11", col="grey90")
#   circos.segments(x0=0, x1=max(ref$V2), y0=0.7, y1=0.7, lwd=0.6, lty="11", col="grey90")
# }, ylim=range(gc$value), track.height=0.08, bg.border=F)
# # gc y axis
# circos.yaxis(at=c(0.3, 0.5, 0.7), labels.cex=0.25, lwd=0, tick.length=0, labels.col=col_text, col="#FFFFFF")
# 
