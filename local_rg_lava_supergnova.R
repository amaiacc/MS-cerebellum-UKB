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
working_dir=paste0(dir,"cerebellum_UKB/data/")
setwd(working_dir)
if(!dir.exists(paste0(working_dir,"local_rg_combined"))){
  dir.create(paste0(working_dir,"local_rg_combined"))
}


#------------------------------------
lava<-read.csv("lava/output/lava_bivar_all_traits_clean.csv") %>% filter(!is.na(trait1)) %>%
  dplyr::select(locus_id,trait1,trait2,
         n.snps,n.pcs,
         rho,rho.lower,rho.upper,
         r2,r2.lower,r2.upper,
         p,p.adjBF,p.adjFDR)
supergnova<-read.csv("supergnova/output/supergnova_bivar_all_traits_clean.csv") %>% filter(!is.na(trait1)) %>%
  dplyr::select(locus_id,trait1,trait2,
                m,
                h2_1,h2_2,
                rho,corr,var,
                p,p.adjBF,p.adjFDR)
# combine both
bivar<-merge(lava,supergnova,by=c("locus_id","trait1","trait2"),
             suffixes=c(".lava",".supergnova"),all=TRUE) %>%
  arrange(locus_id,trait1,trait2) %>%
  dplyr::select(locus_id,trait1,trait2,
                rho.lava,p.lava,p.adjBF.lava,p.adjFDR.lava,
                rho.supergnova,p.supergnova,p.adjBF.supergnova,p.adjFDR.supergnova) %>%
  mutate(chr=strsplit(locus_id,":") %>% sapply(.,"[[",1) %>% as.numeric(),
         start=strsplit(locus_id,":") %>% sapply(.,"[[",2) %>%
           strsplit(.,"-") %>% sapply(.,"[[",1) %>% as.numeric(),
         end=strsplit(locus_id,":") %>% sapply(.,"[[",2) %>%
           strsplit(.,"-") %>% sapply(.,"[[",2) %>% as.numeric()
  ) %>% 
  mutate(Analysis=if_else((trait2=="Crus I (R) Lobule VI (R)"|trait2=="Crus I (R)"|trait2=="VI (R)"),"Cerebellum",
                          if_else( (trait2=="Temp Occ Fusif \ncortex (L)"|trait2=="Occ Fusif \ngyrus (L)"),
                                   "Fusiform",
                                   if_else((trait2=="SCZ"|trait2=="ASD"|trait2=="CP"),"Psych. Disorder/ Cognitive",""
                                   ))) ) %>%
  mutate(chr=as.numeric(chr)) %>%
  arrange(chr,start) %>%
  mutate(cov.supergnova=round(rho.supergnova,digits=4),
         trait1=gsub("\n"," ",trait1) %>% 
           gsub("Crus I (R) Lobule VI (R)","Crus I (R)",.),
         trait2=gsub("\n"," ",trait2) %>%
           gsub("Crus I (R) Lobule VI (R)","Lobule VI (R)",.)
         ) %>%
  dplyr::select(-rho.supergnova) %>%
  rename(Trait1=trait1,Trait2=trait2)
# clean columns etc
bivar_clean <- bivar %>%
  dplyr::select(Analysis,locus_id,Trait1,Trait2,
                rho.lava,p.lava,p.adjBF.lava,
                cov.supergnova,p.supergnova,p.adjBF.supergnova)


# filter based on significance
bivar_sig <- bivar %>% filter(p.adjBF.lava<0.05|p.adjBF.supergnova<0.05)
bivar_sig1 <-bivar_clean %>% 
  filter((p.adjBF.lava<0.05&p.supergnova<0.05)|
           (p.adjBF.supergnova<0.05&p.lava<0.05)) %>%
  mutate(rho.lava=round(rho.lava,digits=2),
         cov.supergnova=round(cov.supergnova,digits=4),
         across(contains("p.adjBF"), ~ format(.,digits=2)),
         across(c("p.lava","p.supergnova"), ~ format(.,scientific=TRUE,digits=2)),) %>%
  dplyr::select(Analysis,locus_id,Trait1,Trait2,
                rho.lava,p.lava,p.adjBF.lava,
                cov.supergnova,p.supergnova,p.adjBF.supergnova)

bivar_sig2<-bivar %>% 
  filter((p.adjFDR.lava<0.05&p.supergnova<0.05)|
           (p.adjFDR.supergnova<0.05&p.lava<0.05))

bivar_sig3<-bivar %>% 
  filter((p.adjFDR.lava<0.05&p.adjFDR.supergnova<0.05))

# for whichever loci sig1 or sig2 --> get all local rgs
bivar_sig_all<-bivar %>% filter(locus_id %in% bivar_sig$locus_id)
bivar_sig1_all<-bivar_clean %>% filter(locus_id %in% bivar_sig1$locus_id)
bivar_sig2_all<-bivar %>% filter(locus_id %in% bivar_sig2$locus_id)
bivar_sig3_all<-bivar %>% filter(locus_id %in% bivar_sig3$locus_id)

# save tables
if(!dir.exists(paste0(working_dir,"local_rg_combined"))){
  dir.create(paste0(working_dir,"local_rg_combined"))
}
setwd(paste0(working_dir,"local_rg_combined"))
# save tables
write.csv(bivar,"local_rg_combined_all.csv",row.names = FALSE)
write.csv(bivar_sig,"local_rg_combined_sig.csv",row.names = FALSE)
write.csv(bivar_sig1,"local_rg_combined_sigBFconsistent.csv",row.names = FALSE)
write.csv(bivar_sig2,"local_rg_combined_sigFDRconsistent.csv",row.names = FALSE)
write.csv(bivar_sig1_all,"local_rg_combined_sigBFconsistent_all.csv",row.names = FALSE)
write.csv(bivar_sig2_all,"local_rg_combined_sigFDRconsistent_all.csv",row.names = FALSE)
write.csv(bivar_sig3_all,"local_rg_combined_sigFDRboth_all.csv",row.names = FALSE)

# in latex format
library(xtable)
sink("local_rg_combined_sigBFconsistent.tex")
cat("\\begin{landscape}\n")
print(
  xtable(bivar_sig1,digits=4,
  caption="Local genetic correlations. Consistent signals across analysis tools (LAVA and SUPERGNOVA) are shown, i.e. estimates that were significant 
after Bonferroni multiple comparison correction in at least one of the tools and were nominally significant in the other one."
  ),
  include.rownames=FALSE)
cat("\\end{landscape}\n")
sink()

sink("local_rg_combined_sigFDRconsistent.tex")
print(xtable(bivar_sig2,digits=4,
       caption="Local genetic correlations.
Consistent signals across analysis tools (LAVA and SUPERGNOVA) are shown, i.e. estimates that were significant 
after FDR multiple comparison correction in at least one of the tools and were nominally significant in the other one.
       "),include.rownames=FALSE)
sink()

#------------------------------------
setwd(paste0(working_dir,"local_rg_combined"))
# visualize...!?
thr<-data.frame(program=toupper(c("lava","supergnova")),n=c(1387,40536)) %>%
  mutate(BFthr=0.05/n,
         FDRthr=c(subset(bivar,p.adjFDR.lava<0.05)$p.lava %>% max(),
                  subset(bivar,p.adjFDR.supergnova<0.05)$p.supergnova %>% max()
         ),
         minlog10BFthr=(-log10(BFthr)),
         minlog10FDRthr=(-log10(FDRthr))
  )
#------------------------------------
# add 
bivar2 <-bivar %>%
  mutate(minLog10p.lava=(-log10(p.lava)),
         minLOG10p.supergnova=(-log10(p.supergnova))
         ) %>%
  mutate(minLog10p.lava=if_else(is.na(p.lava),-1,minLog10p.lava),
         minLOG10p.supergnova=if_else(is.na(p.supergnova),-1,minLOG10p.supergnova)
         )  %>%
  mutate(col=if_else( (is.na(p.lava)|is.na(p.supergnova)),"blue","black"))
comparison_plot <- bivar2 %>% # bivar_sig_all %>%
  ggplot(., aes(
    x=minLog10p.lava,
    y=minLOG10p.supergnova,
    color=col,
    # fill=trait2,
    # shape=trait1
    ) ) +
  geom_vline(data=thr %>% filter(program=="LAVA"),aes(xintercept = minlog10FDRthr),linetype="dashed",color="grey55") +
  geom_vline(data=thr %>% filter(program=="LAVA"),aes(xintercept = minlog10BFthr),linetype="solid",color="grey55") +
  geom_hline(data=thr %>% filter(program=="SUPERGNOVA"),aes(yintercept = minlog10FDRthr),linetype="dashed",color="grey55") +
  geom_hline(data=thr %>% filter(program=="SUPERGNOVA"),aes(yintercept = minlog10BFthr),linetype="solid",color="grey55") +
  geom_point(size=2,alpha=0.5) +
  labs(x=bquote(-log[10] ~ "(p-value LAVA)"),
       y=bquote(-log[10] ~ "(p-value SUPERGNOVA)")) +
  scale_color_manual(values=c("black","blue")) +
  guides(color="none") +
  NULL
  
comparison_plot %>% ggsave(.,file="comparison_minLog10Pval_lava_supergnova.png",width = 5,height = 5)
#------------------------------------
bivar_long<-bivar_sig_all %>% 
  pivot_longer(cols=c(-locus_id,-trait1,-trait2,-chr,-start,-end),names_to="tmp") %>%
  mutate(program=gsub("p.adj","pAdj",tmp) %>%
           strsplit(.,"\\.") %>% sapply(.,"[[",2),
         stat=gsub("p.adj","pAdj",tmp) %>%
           strsplit(.,"\\.") %>% sapply(.,"[[",1)) %>%
  dplyr::select(-tmp) %>%
  pivot_wider(names_from="stat",values_from="value") %>%
  mutate(minlog10P=-log10(p),
         traitPair=paste0(trait1,"~",trait2)) %>%
  arrange(chr,start) %>%
  mutate(locus_id=factor(locus_id,levels=unique(locus_id)),
         program=toupper(program))
#

# ggplot(data=bivar_long) +
#   geom_hline(yintercept=0,linetype="solid",color="black") +
#   geom_hline(yintercept = (-log10(0.05)),linetype="dashed") +
#   geom_hline(data=thr,aes(yintercept = minlog10FDRthr),linetype="dashed",color="blue") +
#   geom_hline(data=thr,aes(yintercept = minlog10BFthr),linetype="dashed",color="red") +
#   geom_point(aes(x=locus_id,y=minlog10P,shape=trait1,color=trait2)) +
#   facet_grid(program~.) +
#   coord_flip() +
#   NULL
# 

# define subsets to present
a<-bivar_long %>% filter(trait1=="Crus I (R) Lobule VI (R)" | 
                       (trait1=="Total Cerebellar Volume" & (trait2=="Crus I (R)"|trait2=="VI (R)")) |
                       is.na(trait1) ) %>%
  mutate(
    # sig3=factor(sig3,levels=c(paste0("p-value < 0.05/",n),"p-value < 0.05","p-value > 0.05","")),
         trait1=factor(trait1,levels=c("Total Cerebellar Volume","Crus I (R) Lobule VI (R)","") ),
         trait2=factor(trait2,levels=c("Crus I (R)","VI (R)","Crus I (R) Lobule VI (R)","") )
  )
colors_fusif<-met.brewer(name="Signac",n=14)[c(14,13)]
b <- bivar_long %>% filter(grepl("Fusif",trait2) | is.na(trait1) ) %>%
  mutate(
    trait1=factor(trait1,levels=c("Crus I (R)","Lobule VI (R)","Total Cerebellar Volume","")),
    trait2=factor(trait2,levels=c("Occ Fusif \ngyrus (L)","Temp Occ Fusif \ncortex (L)",""))
  )
c <- bivar_long %>% filter(grepl("CP|SCZ|ASD",trait2) | is.na(trait1) ) %>%
  mutate(
    trait1=factor(trait1,levels=c("Crus I (R)","Lobule VI (R)","Total Cerebellar Volume","")),
    trait2=factor(trait2,levels=c("ASD","SCZ","CP",""))
  )

#------------------------------------
pA <- a %>%
  ggplot(., aes(x=locus_id,y=minlog10P,
             fill=trait1,
             shape=trait2)) +
  # geom_hline(yintercept=0,linetype="solid",color="black") +
  # geom_hline(yintercept = (-log10(0.05)),linetype="dashed") +
  geom_hline(data=thr,aes(yintercept = minlog10FDRthr),linetype="dashed",color="grey55") +
  geom_hline(data=thr,aes(yintercept = minlog10BFthr),linetype="solid",color="grey55") +
  geom_point(size=3,color="white",alpha=0.8, position="dodge") +
  facet_grid(.~program,scales="free_x",space = "free_x") +
  labs(
    x="",
    y=bquote(-log[10] ~ "(p-value)")) +
  coord_flip() +
  # scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_manual(values=c(viridis(3)[c(1,3)],"#ffffff")) + # ,
  # scale_color_manual(values=c(viridis(3)[c(1,3)])) + # ,"#ffffff"
  scale_shape_manual(values=c(21,24,23), #,22
                     labels=c("Total Cerebellar Volume  ~ Crus I (R)",
                              "Total Cerebellar Volume ~ Lobule VI (R)",
                              "Crus I (R) ~ Lobule VI (R)","")
  ) +
  guides(fill="none",color="none",
         shape=guide_legend(override.aes = list(color="white",
                                                fill=c(viridis(3)[1],
                                                       viridis(3)[1],
                                                       viridis(3)[3])
         ), title="",order=2)
  ) +
  scale_x_discrete(
    # breaks = ,
    limits = rev(levels(c$locus_id))#,
    # expand = c(0,0)
  ) +
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="black",size=1),
        legend.position = "none") +
    panel_border(color="black") +
  NULL


pB <- b %>%
  ggplot(aes(x=locus_id,y=minlog10P,
             color=trait2,
             fill=trait2,
             shape=trait1)) +
  # geom_hline(yintercept=0,linetype="solid",color="black") +
  # geom_hline(yintercept = (-log10(0.05)),linetype="dashed") +
  geom_hline(data=thr,aes(yintercept = minlog10FDRthr),linetype="dashed",color="grey55") +
  geom_hline(data=thr,aes(yintercept = minlog10BFthr),linetype="solid",color="grey55") +
  geom_point(size=3,color="white",alpha=0.8) +
  facet_grid(.~program,scales="free_x",space = "free_x") +
  labs(
    x="",
    y=bquote(-log[10] ~ "(p-value)")) +
  coord_flip() +
  # scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_manual(values=c(colors_fusif)) + #,"#ffffff"
  scale_color_manual(values=c(colors_fusif)) + # ,"#ffffff"
  scale_shape_manual(values=c(21,24,22,23)) +
  scale_x_discrete(
    # breaks = ,
    limits = rev(levels(c$locus_id)),
    labels = element_blank()#,
    # expand = c(0,0)
  ) +
  guides(alpha="none", color="none",
         fill=guide_legend(override.aes = list(color=c(colors_fusif), shape=18),order=2, title="Cortical volume"),
         shape=guide_legend(override.aes = list(color=c("black","black","black")),order=1, title="Cerebellar volume")
         )+
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="black",size=1),
        legend.position = "none") +
    panel_border(color="black") +
  NULL


pC <- c %>%
  ggplot(aes(x=locus_id,y=minlog10P,
             color=trait2,
             fill=trait2,
             shape=trait1)) +
  # geom_hline(yintercept=0,linetype="solid",color="black") +
  # geom_hline(yintercept = (-log10(0.05)),linetype="dashed") +
  geom_hline(data=thr,aes(yintercept = minlog10FDRthr),linetype="dashed",color="grey55") +
  geom_hline(data=thr,aes(yintercept = minlog10BFthr),linetype="solid",color="grey55") +
  geom_point(size=3,color="white",alpha=0.8) +
  facet_grid(.~program,scales="free_x",space = "free_x") +
  labs(
    x="",
    y=bquote(-log[10] ~ "(p-value)")) +
  coord_flip() +
  scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_met_d(name="Austria") +
  scale_color_met_d(name="Austria") +
  scale_shape_manual(values=c(21,24,22,23)) +
  scale_x_discrete(
    # breaks = ,
    limits = rev(levels(c$locus_id)),
    # labels = element_blank(),
    # expand = c(0,0),
    position="top"
  ) +
  guides(alpha="none", 
         # color="none",
         fill=guide_legend(override.aes = list(color=c(met.brewer(name="Austria",n=3)), 
                                               shape=18), title="Cognitive trait/\nPyschiatric disorder"),
         shape="none"
         # shape=guide_legend(override.aes = list(color="black"),order=1, title="Cerebellar volume")
  ) +
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="black",size=1),
        legend.position = "none") +
  panel_border(color="black") +
  NULL
# combine plots
p<-plot_grid(pA,
             pB,
             pC,#NULL,
             nrow=1,
             rel_widths = c(2,1.4,2),
             labels=c("A","B","C"))
l<-plot_grid(
  # get_legend(pA+theme(legend.position="right")+guides(shape="none")),
             # plot_grid(
               get_legend(pB+theme(legend.position="right")+guides(fill="none")),
               get_legend(pA+theme(legend.position="right")+guides(alpha="none")),
               # ncol=1
             # ),
             get_legend(pB+theme(legend.position="right")+guides(shape="none")),
             get_legend(pC+theme(legend.position="right")),
             ncol=1
             # align="v"
)
p2<-plot_grid(p,l,nrow=1,rel_widths = c(2,0.5))
ggsave(p2,file="Figure_rg_local_lava_supergnova_bivar_BFone.png",width=20,height=10)
#------------------------------------
# show only lava
loci_lava<-bivar %>% filter(p.adjBF.lava<0.05) %>% pull(locus_id) %>% unique()
bivar_long1<-subset(bivar_long,program=="LAVA") %>%
  filter(locus_id %in% loci_lava)
bivar_long1$locus_id<-droplevels(bivar_long1$locus_id)
thr1<-subset(thr,program=="LAVA")
# define subsets to present
a<-bivar_long1 %>% filter(trait1=="Crus I (R) Lobule VI (R)" | 
                           (trait1=="Total Cerebellar Volume" & (trait2=="Crus I (R)"|trait2=="VI (R)")) |
                           is.na(trait1) ) %>%
  mutate(
    # sig3=factor(sig3,levels=c(paste0("p-value < 0.05/",n),"p-value < 0.05","p-value > 0.05","")),
    trait1=factor(trait1,levels=c("Total Cerebellar Volume","Crus I (R) Lobule VI (R)","") ),
    trait2=factor(trait2,levels=c("Crus I (R)","VI (R)","Crus I (R) Lobule VI (R)","") )
  )
colors_fusif<-met.brewer(name="Signac",n=14)[c(14,13)]
b <- bivar_long1 %>% filter(grepl("Fusif",trait2) | is.na(trait1) ) %>%
  mutate(
    trait1=factor(trait1,levels=c("Crus I (R)","Lobule VI (R)","Total Cerebellar Volume","")),
    trait2=factor(trait2,levels=c("Occ Fusif \ngyrus (L)","Temp Occ Fusif \ncortex (L)",""))
  )
c <- bivar_long1 %>% filter(grepl("CP|SCZ|ASD",trait2) | is.na(trait1) ) %>%
  mutate(
    trait1=factor(trait1,levels=c("Crus I (R)","Lobule VI (R)","Total Cerebellar Volume","")),
    trait2=factor(trait2,levels=c("ASD","SCZ","CP",""))
  )

#------------------------------------
pA <- a %>%
  ggplot(., aes(x=locus_id,y=minlog10P,
                fill=trait1,
                shape=trait2)) +
  # geom_hline(yintercept=0,linetype="solid",color="black") +
  # geom_hline(yintercept = (-log10(0.05)),linetype="dashed") +
  geom_hline(data=thr1,aes(yintercept = minlog10FDRthr),linetype="dashed",color="grey55") +
  geom_hline(data=thr1,aes(yintercept = minlog10BFthr),linetype="solid",color="grey55") +
  geom_point(size=3,color="white",alpha=0.8) +
  # facet_grid(.~program,scales="free_x",space = "free_x") +
  labs(
    x="",
    y=bquote(-log[10] ~ "(p-value)")) +
  coord_flip() +
  # scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_manual(values=c(viridis(3)[c(1,3)],"#ffffff")) + # ,
  # scale_color_manual(values=c(viridis(3)[c(1,3)])) + # ,"#ffffff"
  scale_shape_manual(values=c(21,24,23), #,22
                     labels=c("Total Cerebellar Volume  ~ Crus I (R)",
                              "Total Cerebellar Volume ~ Lobule VI (R)",
                              "Crus I (R) ~ Lobule VI (R)","")
  ) +
  scale_x_discrete(
    # breaks = ,
    limits = rev(levels(c$locus_id))#,
    # expand = c(0,0)
  ) +
  guides(fill="none",color="none",
         shape=guide_legend(override.aes = list(color="white",
                                                fill=c(viridis(3)[1],
                                                       viridis(3)[1],
                                                       viridis(3)[3])
         ), title="",order=2)
  ) +
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="black",size=1),
        legend.position = "none") +
    panel_border(color="black") +
  NULL


pB <- b %>%
  ggplot(aes(x=locus_id,y=minlog10P,
             color=trait2,
             fill=trait2,
             shape=trait1)) +
  # geom_hline(yintercept=0,linetype="solid",color="black") +
  # geom_hline(yintercept = (-log10(0.05)),linetype="dashed") +
  geom_hline(data=thr1,aes(yintercept = minlog10FDRthr),linetype="dashed",color="grey55") +
  geom_hline(data=thr1,aes(yintercept = minlog10BFthr),linetype="solid",color="grey55") +
  geom_point(size=3,color="white",alpha=0.8) +
  # facet_grid(.~program,scales="free_x",space = "free_x") +
  labs(
    x="",
    y=bquote(-log[10] ~ "(p-value)")) +
  coord_flip() +
  # scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_manual(values=c(colors_fusif)) + #,"#ffffff"
  scale_color_manual(values=c(colors_fusif)) + # ,"#ffffff"
  scale_shape_manual(values=c(21,24,22,23)) +
  scale_x_discrete(
    # breaks = ,
    labels = element_blank(),
    limits = rev(levels(c$locus_id))#,
    # expand = c(0,0)
  ) +
  guides(alpha="none", color="none",
         fill=guide_legend(override.aes = list(color=c(colors_fusif), shape=18),order=2, title="Cortical volume"),
         shape=guide_legend(override.aes = list(color=c("black","black","black")),order=1, title="Cerebellar volume")
  )+
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="black",size=1),
        legend.position = "none") +
    panel_border(color="black") +
  NULL


pC <- c %>%
  ggplot(aes(x=locus_id,y=minlog10P,
             color=trait2,
             fill=trait2,
             shape=trait1)) +
  # geom_hline(yintercept=0,linetype="solid",color="black") +
  # geom_hline(yintercept = (-log10(0.05)),linetype="dashed") +
  geom_hline(data=thr1,aes(yintercept = minlog10FDRthr),linetype="dashed",color="grey55") +
  geom_hline(data=thr1,aes(yintercept = minlog10BFthr),linetype="solid",color="grey55") +
  geom_point(size=3,color="white",alpha=0.8) +
  # facet_grid(.~program,scales="free_x",space = "free_x") +
  labs(
    x="",
    y=bquote(-log[10] ~ "(p-value)")) +
  coord_flip() +
  scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_met_d(name="Austria") +
  scale_color_met_d(name="Austria") +
  scale_shape_manual(values=c(21,24,22,23)) +
  scale_x_discrete(
    # breaks = ,
    limits = rev(levels(c$locus_id)),
    # labels = element_blank(),
    # expand = c(0,0),
    position="top"
  ) +
  guides(alpha="none", 
         # color="none",
         fill=guide_legend(override.aes = list(
           color=c(met.brewer(name="Austria",n=3)),
           shape=18), title="Cognitive trait/\nPyschiatric disorder"),
         shape="none"
         # shape=guide_legend(override.aes = list(color="black"),order=1, title="Cerebellar volume")
  ) +
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="black",size=1),
        legend.position = "none") +
    panel_border(color="black") +
  NULL
# combine plots
# combine plots
p<-plot_grid(pA,
             pB,
             pC,#NULL,
             nrow=1,
             rel_widths = c(2,1.2,2),
             labels=c("A","B","C"))
l<-plot_grid(
  # get_legend(pA+theme(legend.position="right")+guides(shape="none")),
  # plot_grid(
  get_legend(pB+theme(legend.position="right")+guides(fill="none")),
  get_legend(pA+theme(legend.position="right")+guides(alpha="none")),
  # ncol=1
  # ),
  get_legend(pB+theme(legend.position="right")+guides(shape="none")),
  get_legend(pC+theme(legend.position="right")),
  ncol=1
  # align="v"
)
p2<-plot_grid(p,l,nrow=1,rel_widths = c(1,0.25))
if(!dir.exists(paste0(working_dir,"local_rg_combined"))){
  dir.create(paste0(working_dir,"local_rg_combined"))
}
setwd(paste0(working_dir,"local_rg_combined"))
ggsave(p2,file="Figure_rg_local_lavaBF_bivar.png",width=12,height=8)

#------------------------------------
# show only lava
loci_supergnova<-bivar %>% filter(p.adjBF.supergnova<0.05) %>% pull(locus_id) %>% unique()
bivar_long1<-subset(bivar_long,program=="SUPERGNOVA") %>%
  filter(locus_id %in% loci_supergnova)
bivar_long1$locus_id<-droplevels(bivar_long1$locus_id)
thr1<-subset(thr,program=="SUPERGNOVA")
# define subsets to present
a<-bivar_long1 %>% filter(trait1=="Crus I (R) Lobule VI (R)" | 
                            (trait1=="Total Cerebellar Volume" & (trait2=="Crus I (R)"|trait2=="VI (R)")) |
                            is.na(trait1) ) %>%
  mutate(
    # sig3=factor(sig3,levels=c(paste0("p-value < 0.05/",n),"p-value < 0.05","p-value > 0.05","")),
    trait1=factor(trait1,levels=c("Total Cerebellar Volume","Crus I (R) Lobule VI (R)","") ),
    trait2=factor(trait2,levels=c("Crus I (R)","VI (R)","Crus I (R) Lobule VI (R)","") )
  )
colors_fusif<-met.brewer(name="Signac",n=14)[c(14,13)]
b <- bivar_long1 %>% filter(grepl("Fusif",trait2) | is.na(trait1) ) %>%
  mutate(
    trait1=factor(trait1,levels=c("Crus I (R)","Lobule VI (R)","Total Cerebellar Volume","")),
    trait2=factor(trait2,levels=c("Occ Fusif \ngyrus (L)","Temp Occ Fusif \ncortex (L)",""))
  )
c <- bivar_long1 %>% filter(grepl("CP|SCZ|ASD",trait2) | is.na(trait1) ) %>%
  mutate(
    trait1=factor(trait1,levels=c("Crus I (R)","Lobule VI (R)","Total Cerebellar Volume","")),
    trait2=factor(trait2,levels=c("ASD","SCZ","CP",""))
  )

#------------------------------------
pA <- a %>%
  ggplot(., aes(x=locus_id,y=minlog10P,
                fill=trait1,
                shape=trait2)) +
  # geom_hline(yintercept=0,linetype="solid",color="black") +
  # geom_hline(yintercept = (-log10(0.05)),linetype="dashed") +
  geom_hline(data=thr1,aes(yintercept = minlog10FDRthr),linetype="dashed",color="grey55") +
  geom_hline(data=thr1,aes(yintercept = minlog10BFthr),linetype="solid",color="grey55") +
  geom_point(size=3,color="white",alpha=0.8) +
  # facet_grid(.~program,scales="free_x",space = "free_x") +
  labs(
    x="",
    y=bquote(-log[10] ~ "(p-value)")) +
  coord_flip() +
  # scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_manual(values=c(viridis(3)[c(1,3)],"#ffffff")) + # ,
  # scale_color_manual(values=c(viridis(3)[c(1,3)])) + # ,"#ffffff"
  scale_shape_manual(values=c(21,24,23), #,22
                     labels=c("Total Cerebellar Volume  ~ Crus I (R)",
                              "Total Cerebellar Volume ~ Lobule VI (R)",
                              "Crus I (R) ~ Lobule VI (R)","")
  ) +
  scale_x_discrete(
    # breaks = ,
    limits = rev(levels(c$locus_id))#,
    # expand = c(0,0)
  ) +
  guides(fill="none",color="none",
         shape=guide_legend(override.aes = list(color="white",
                                                fill=c(viridis(3)[1],
                                                       viridis(3)[1],
                                                       viridis(3)[3])
         ), title="",order=2)
  ) +
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="black",size=1),
        legend.position = "none") +
  panel_border(color="black") +
  NULL


pB <- b %>%
  ggplot(aes(x=locus_id,y=minlog10P,
             color=trait2,
             fill=trait2,
             shape=trait1)) +
  # geom_hline(yintercept=0,linetype="solid",color="black") +
  # geom_hline(yintercept = (-log10(0.05)),linetype="dashed") +
  geom_hline(data=thr1,aes(yintercept = minlog10FDRthr),linetype="dashed",color="grey55") +
  geom_hline(data=thr1,aes(yintercept = minlog10BFthr),linetype="solid",color="grey55") +
  geom_point(size=3,color="white",alpha=0.8) +
  # facet_grid(.~program,scales="free_x",space = "free_x") +
  labs(
    x="",
    y=bquote(-log[10] ~ "(p-value)")) +
  coord_flip() +
  # scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_manual(values=c(colors_fusif)) + #,"#ffffff"
  scale_color_manual(values=c(colors_fusif)) + # ,"#ffffff"
  scale_shape_manual(values=c(21,24,22,23)) +
  scale_x_discrete(
    # breaks = ,
    labels = element_blank(),
    limits = rev(levels(c$locus_id))#,
    # expand = c(0,0)
  ) +
  guides(alpha="none", color="none",
         fill=guide_legend(override.aes = list(color=c(colors_fusif), shape=18),order=2, title="Cortical volume"),
         shape=guide_legend(override.aes = list(color=c("black","black","black")),order=1, title="Cerebellar volume")
  )+
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="black",size=1),
        legend.position = "none") +
  panel_border(color="black") +
  NULL


pC <- c %>%
  ggplot(aes(x=locus_id,y=minlog10P,
             color=trait2,
             fill=trait2,
             shape=trait1)) +
  # geom_hline(yintercept=0,linetype="solid",color="black") +
  # geom_hline(yintercept = (-log10(0.05)),linetype="dashed") +
  geom_hline(data=thr1,aes(yintercept = minlog10FDRthr),linetype="dashed",color="grey55") +
  geom_hline(data=thr1,aes(yintercept = minlog10BFthr),linetype="solid",color="grey55") +
  geom_point(size=3,color="white",alpha=0.8) +
  # facet_grid(.~program,scales="free_x",space = "free_x") +
  labs(
    x="",
    y=bquote(-log[10] ~ "(p-value)")) +
  coord_flip() +
  scale_alpha_manual(values=c(1,0.2,0.01,0)) +
  scale_fill_met_d(name="Austria") +
  scale_color_met_d(name="Austria") +
  scale_shape_manual(values=c(21,24,22,23)) +
  scale_x_discrete(
    # breaks = ,
    limits = rev(levels(c$locus_id)),
    # labels = element_blank(),
    # expand = c(0,0),
    position="top"
  ) +
  guides(alpha="none", 
         # color="none",
         fill=guide_legend(override.aes = list(
           color=c(met.brewer(name="Austria",n=3)),
           shape=18), title="Cognitive trait/\nPyschiatric disorder"),
         shape="none"
         # shape=guide_legend(override.aes = list(color="black"),order=1, title="Cerebellar volume")
  ) +
  theme_cowplot() +
  theme(panel.spacing=unit(0, "lines"),
        strip.background =element_rect(fill="white",color="black",size=1),
        legend.position = "none") +
  panel_border(color="black") +
  NULL
# combine plots
# combine plots
p<-plot_grid(pA,
             pB,
             pC,#NULL,
             nrow=1,
             rel_widths = c(2,1.4,2),
             labels=c("A","B","C"))
l<-plot_grid(
  # get_legend(pA+theme(legend.position="right")+guides(shape="none")),
  # plot_grid(
  get_legend(pB+theme(legend.position="right")+guides(fill="none")),
  get_legend(pA+theme(legend.position="right")+guides(alpha="none")),
  # ncol=1
  # ),
  get_legend(pB+theme(legend.position="right")+guides(shape="none")),
  get_legend(pC+theme(legend.position="right")),
  ncol=1
  # align="v"
)
p2<-plot_grid(p,l,nrow=1,rel_widths = c(2,0.5))
if(!dir.exists(paste0(working_dir,"local_rg_combined"))){
  dir.create(paste0(working_dir,"local_rg_combined"))
}
setwd(paste0(working_dir,"local_rg_combined"))
ggsave(p2,file="Figure_rg_local_supergnovaBF_bivar.png",width=15,height=10)


#------------------------------------


