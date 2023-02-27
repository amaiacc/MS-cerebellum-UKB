library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
# clean environment
rm(list=ls())
gc()
#------------------------------------
# Set directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
#------------------------------------
# some functions
source(paste(dir,"general_scripts/helper_functions.R",sep=""))
## phenotype information
primary_dir=paste(dir,"resources/datasets/GWAS_sumstats/downloaded_data/UKB/","BIG40","/",sep="")
idps<-read.csv(paste(primary_dir,"IDPs_summary.csv",sep=""),header=TRUE)
idps$pheno<-numeric_nchar(idps$Pheno)
#------------------------------------
arg = commandArgs(T)
arg<-c(paste0(dir,"cerebellum_UKB/data/supergnova"),
       paste0(dir,"cerebellum_UKB/data/supergnova/output/"),
       "local_rg.csv"
       )
input_dir=arg[1];out_dir=arg[2];name=arg[3]
#------------------------------------
setwd(input_dir)
# read all files in input directory
for (f in list.files(input_dir,pattern="lavaPartition")){
  tmp<-read.table(f,sep=" ",header=TRUE)
  tmp$file<-f
  if(exists("d")){d<-rbind(d,tmp)} else {d<-tmp}
}
rm(tmp)
setwd(out_dir)
#------------------------------------
# get phenotypes 1 and 2
d <- d %>% mutate(
  tmp=strsplit(file,"\\.") %>% sapply("[[",1) %>% gsub("Chambers2022_","",.),
  p1=strsplit(tmp,"_") %>% sapply("[[",1),
  p2=strsplit(tmp,"_") %>% sapply("[[",2),
  locus_id=paste0(chr,":",start,"-",end)
) %>% dplyr::select(-tmp,-file)

# merge with idps to get names
d<-merge(d,idps[,c("pheno","region","hemisphere")],by.x="p1",by.y="pheno",all.x=TRUE) %>% 
  merge(.,idps[,c("pheno","region","hemisphere")],by.x="p2",by.y="pheno",all.x=TRUE,suffixes=c(".p1",".p2")) %>%
  mutate(p2_name=if_else(!is.na(region.p2),paste0(region.p2,"_",hemisphere.p2),p2))
#
d <- d %>% 
  mutate(trait1=paste(region.p1,hemisphere.p1,sep="_"),
         trait2=paste(region.p2,hemisphere.p2,sep="_")) %>%
  mutate(
    trait1=if_else(is.na(region.p1),p1,trait1),
    trait2=if_else(is.na(region.p2),p2,trait2)
  )

#------------------------------------
# correct for multiple comparisons
# BF corrected rg's
d <- d %>% mutate(p.adjBF=p.adjust(p,method="bonferroni"),p.adjFDR=p.adjust(p,method="fdr"))
# only if non-zero heritability is present for one of the traits --> i.e. corr is not NA
d2<- d %>% filter(!is.na(corr)) %>% 
  mutate(p.adjBF2=p.adjust(p,method="bonferroni"),p.adjFDR2=p.adjust(p,method="fdr"))
#------------------------------------
d_clean <- d2 %>%
  mutate(traits=paste(trait1,"-",trait2)  ) %>%
  dplyr::select(#locus,
                locus_id,traits,
                # p1,p2,
                h2_1,h2_2,
                m,
                # n.pcs,
                rho,
                # rho_lower_upper,
                # r2,r2_lower_upper,
                corr,var,
                p,p.adjBF,p.adjFDR)

d_wide<- d_clean %>% 
  # dplyr::select(-p1,-p2) %>%
  pivot_wider(names_from=traits,values_from = c(h2_1,h2_2,m,rho,corr,var,p,p.adjBF,p.adjFDR)) %>%
  select(locus_id,contains(unique(d_clean$traits))) %>% arrange(locus_id)


# subset only sig
d2_f<-subset(d_clean,p.adjBF<0.05) %>% arrange(locus_id)
d2_f_wide<- d2_f %>% 
  # dplyr::select(-p1,-p2) %>%
  pivot_wider(names_from=traits,values_from = c(h2_1,h2_2,m,rho,corr,var,p,p.adjBF,p.adjFDR)) %>%
  select(locus_id,contains(unique(d_clean$traits))) %>% arrange(locus_id)

d2_f2<-subset(d_clean,p.adjFDR<0.05) %>% arrange(locus_id)
d2_f2_wide<- d2_f2 %>% 
  # dplyr::select(-p1,-p2) %>%
  pivot_wider(names_from=traits,values_from = c(h2_1,h2_2,m,rho,corr,var,p,p.adjBF,p.adjFDR)) %>%
  select(locus_id,contains(unique(d_clean$traits))) %>% arrange(locus_id)

# save
write.csv(d,"supergnova_bivar_all_traits.csv",row.names = FALSE)
write.csv(d2,"supergnova_bivar_all_traits_nonzeroh2.csv",row.names = FALSE)
write.csv(d_wide,"supergnova_bivar_wide_all_traits_nonzeroh2.csv",row.names = FALSE)

write.csv(d2_f,"supergnova_bivar_all_traits_nonzeroh2_BFsig.csv",row.names = FALSE)
write.csv(d2_f_wide,"supergnova_bivar_wide_all_traits_nonzeroh2_BFsig.csv",row.names = FALSE)

write.csv(d2_f2,"supergnova_bivar_all_traits_nonzeroh2_FDRsig.csv",row.names = FALSE)
write.csv(d2_f2_wide,"supergnova_bivar_wide_all_traits_nonzeroh2_FDRsig.csv",row.names = FALSE)

#------------------------------------
# check correlation of p-vals for 0143 and 0146
dw143_146 <- d2 %>% filter(p1 %in% c("0143","0146")) %>%
  dplyr::select(p1,p2,p2_name,locus_id,p, rho,corr) %>%
  pivot_wider(id_cols=c(p2,p2_name,locus_id),names_from=p1,values_from = c(p, rho,corr))

tmp<-dw143_146 %>% filter(! p2 %in% c("0143","0146"))

comp_VI_CrusI_pvals <- tmp %>% ggplot() +
  geom_point(aes(x=log10(p_0143),y=log10(p_0146),color=p2_name),alpha=0.5) +
  scale_color_brewer(palette=1,type="qual",name="") +
  theme_cowplot()  + theme(legend.position="none")
comp_VI_CrusI_rho <- tmp %>% ggplot() +
  geom_point(aes(x=rho_0143,y=rho_0146,color=p2_name),alpha=0.5) +
  scale_color_brewer(palette=1,type="qual",name="") +
  theme_cowplot() + theme(legend.position="none")
leg<-get_legend(comp_VI_CrusI_rho + theme(legend.position="bottom"))
plot_grid(plot_grid(comp_VI_CrusI_pvals,comp_VI_CrusI_rho),leg,ncol=1,rel_heights = c(1,0.1))

cor(tmp$rho_0143,tmp$rho_0146,use ="pairwise.complete.obs")
#------------------------------------
# subset: shared loci between cereb. crus I and cereb. lobule VI (right)
## loci sig at FDR 0.05
d2 <- d %>% 
  filter(p1=="0146" & p2=="0143") %>%
  mutate(p.adjFDR=p.adjust(p,method="fdr"))
d2_f<-d2 %>% filter(p.adjFDR<0.05)
#------------------------------------
# rg for those loci: for cereb and fusif phenos
d3 <- d %>% filter(grepl("FUSIF",p2_name)) %>%
  filter(locus_id %in% d2_f$locus_id) %>%
  mutate(p.adjFDR=p.adjust(p,method="fdr")) 
d3_f<-d3 %>% filter(p.adjFDR<0.05)

# check if those loci are associated with ASD or CP
d4 <- d %>% filter(p2 %in% c("ASD","SCZ","CP")) %>%
  filter(locus_id %in% d3_f$locus_id) %>%
  mutate(p.adjFDR=p.adjust(p,method="fdr"))
d4_f <- d4 %>% filter(p.adjFDR<0.05)
# check all values for these loci
lapply(unique(d4$locus_id),function(l){
  subset(d,locus_id==l) %>% arrange(p) %>% dplyr::select(locus_id,p1,p2,rho,corr,p)
} )
