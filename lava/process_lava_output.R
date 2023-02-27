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
arg<-c(paste0(dir,"cerebellum_UKB/data/lava/output/"),
       paste0("cerebellum_UKB/data/lava/output/"),
       "local_rg.csv"
       )
input_dir=arg[1];out_dir=arg[2];name=arg[3]
#------------------------------------
setwd(input_dir)
# read all files in input directory, combine per type
#------------------------------------
## UNIV
univ_files<-list.files(input_dir,pattern="univ$")
univ_files<-univ_files[which(lapply(strsplit(univ_files,"\\."),length)==2)]
for (f in univ_files){
  tmp<-read.table(f,sep=" ",header=TRUE)
  tmp$file<-f
  if(exists("univ")){univ<-rbind(univ,tmp)} else {univ<-tmp}
}
rm(tmp,univ_files)

# get phenotypes 1 and 2
univ <- univ %>% mutate(
  p1=phen,
  locus_id=paste0(chr,":",start,"-",stop)
)
# merge with idps to get names
univ<-merge(univ,idps[,c("Pheno","region","hemisphere")],by.x="p1",by.y="Pheno",all.x=TRUE)

univ<- univ %>%
  mutate(trait=paste(region,hemisphere,sep="_")) %>%
  mutate(trait=if_else(is.na(region),phen,trait)) %>% 
  arrange(locus)
# format to wide
univ_wide<- univ %>%
  dplyr::select(locus,locus_id,trait,n.snps,n.pcs,h2.obs,p) %>%
  pivot_wider(names_from=trait,values_from = c(n.snps,n.pcs,h2.obs,p)) %>%
  select(locus,locus_id,contains(unique(univ$trait))) %>% arrange(locus)
# save univ files
write.csv(univ,file="lava_univ_all_traits.csv",row.names = FALSE)
write.csv(univ_wide,file="lava_univ_wide_all_traits.csv",row.names = FALSE)


#------------------------------------
## bivar
for (f in list.files(input_dir,pattern="bivar$")){
  if(length(readLines(f))>1){
    tmp<-read.table(f,sep=" ",header=TRUE)
    tmp$file<-f
    if(exists("bivar")){bivar<-rbind(bivar,tmp)} else {bivar<-tmp}
  }
}
rm(tmp)

bivar <- bivar %>% mutate(
  tmp=strsplit(file,"\\.bivar") %>% sapply("[[",1),
  p1=strsplit(tmp,"\\.") %>% sapply("[[",1),
  p2=strsplit(tmp,"\\.") %>% sapply("[[",2),
  locus_id=paste0(chr,":",start,"-",stop)
) %>% dplyr::select(-tmp)

# merge with idps to get names
bivar<-merge(bivar,idps[,c("pheno","region","hemisphere")],by.x="p1",by.y="pheno",all.x=TRUE) %>% 
  merge(.,idps[,c("pheno","region","hemisphere")],by.x="p2",by.y="pheno",all.x=TRUE,suffixes=c(".p1",".p2")) %>%
  mutate(p2_name=if_else(!is.na(region.p2),paste0(region.p2,"_",hemisphere.p2),p2))

# BF corrected rg's
bivar <- bivar %>% 
  mutate(p.adjBF=p.adjust(p,method="bonferroni"),p.adjFDR=p.adjust(p,method="fdr")) %>%
  mutate(trait1=paste(region.p1,hemisphere.p1,sep="_"),
         trait2=paste(region.p2,hemisphere.p2,sep="_")) %>%
  mutate(
    trait1=if_else(is.na(region.p1),phen1,trait1),
    trait2=if_else(is.na(region.p2),phen2,trait2)
         )
bivar_clean <- bivar %>%
  mutate(traits=paste(trait1,"-",trait2),
         rho_lower_upper=paste0(round(rho.lower,digits=2),";",round(rho.upper,digits=2)),
         r2_lower_upper=paste0(round(r2.lower,digits=2),";",round(r2.upper,digits=2))
         ) %>%
  dplyr::select(locus,locus_id,traits,n.snps,n.pcs,rho,rho_lower_upper,r2,r2_lower_upper,p,p.adjBF,p.adjFDR) %>%
  distinct()

bivar_wide<- bivar_clean %>%
  pivot_wider(names_from=traits,values_from = c(n.snps,n.pcs,rho,rho_lower_upper,r2,r2_lower_upper,p,p.adjBF,p.adjFDR)) %>%
  select(locus,locus_id,contains(unique(bivar_clean$traits))) %>% arrange(locus)

corrected_loci<-bivar$locus[bivar$p.adjBF<0.05]
corrected_lociFDR<-bivar$locus[bivar$p.adjFDR<0.05]
bivar_wide_corrected <- bivar_clean %>% filter(locus %in% corrected_loci) %>%
  pivot_wider(names_from=traits,values_from = c(n.snps,n.pcs,rho,rho_lower_upper,r2,r2_lower_upper,p,p.adjBF,p.adjFDR)) %>%
  select(locus,locus_id,contains(unique(bivar_clean$traits))) %>% arrange(locus)

bivar_wide_correctedFDR <- bivar_clean %>% filter(locus %in% corrected_lociFDR) %>%
  pivot_wider(names_from=traits,values_from = c(n.snps,n.pcs,rho,rho_lower_upper,r2,r2_lower_upper,p,p.adjBF,p.adjFDR)) %>%
  select(locus,locus_id,contains(unique(bivar_clean$traits))) %>% arrange(locus)


bivar_clean %>% filter(locus %in% corrected_lociFDR) %>% arrange(locus_id) %>% View()

# save bivar files
write.csv(bivar,file="lava_bivar_all_traits.csv",row.names = FALSE)
write.csv(bivar_wide,file="lava_bivar_wide_all_traits.csv",row.names = FALSE)
write.csv(bivar_wide_corrected,file="lava_bivar_wide_all_traits_BFsig.csv",row.names = FALSE)
write.csv(bivar_wide_correctedFDR,file="lava_bivar_wide_all_traits_FDRsig.csv",row.names = FALSE)

#------------------------------------
## PCOR
## partial correlations
for (f in list.files(input_dir,pattern="pcor$")){
  tmp<-read.table(f,sep=" ",header=TRUE)
  tmp$file<-f
  if(exists("pcor")){pcor<-rbind(pcor,tmp)} else {pcor<-tmp}
}
rm(tmp)

pcor <- pcor %>% mutate(
  tmp=strsplit(file,"\\.pcor") %>% 
    sapply("[[",1) %>% 
    gsub("_condBy","_condBy\\.",.) %>%
    gsub("\\.\\.",".",.),
  p1=strsplit(tmp,"\\.") %>% sapply("[[",1),
  p2=strsplit(tmp,"\\.") %>% sapply("[[",2) %>% gsub("_condBy","",.),
  p3=strsplit(tmp,"\\.") %>% sapply("[[",3),
  locus_id=paste0(chr,":",start,"-",stop)
) %>% dplyr::select(-tmp)
# merge with idps to get names
pcor<-merge(pcor,idps[,c("pheno","region","hemisphere")],by.x="p1",by.y="pheno",all.x=TRUE) %>% 
  merge(.,idps[,c("pheno","region","hemisphere")],by.x="p3",by.y="pheno",all.x=TRUE,suffixes=c(".p1",".p3")) %>%
  mutate(p3_name=if_else(!is.na(region.p3),paste0(region.p3,"_",hemisphere.p3),p3))


pcor <- pcor %>% 
  mutate(trait1=paste(region.p1,hemisphere.p1,sep="_"),
         # trait2=paste(region.p2,hemisphere.p2,sep="_"),
         trait3=paste(region.p3,hemisphere.p3,sep="_")) %>%
  mutate(
    trait1=if_else(is.na(region.p1),p1,trait1),
    # trait2=if_else(is.na(region.p2),phen2,trait2),
    trait3=if_else(is.na(region.p3),p3,trait3)
  )

pcor_clean <- pcor %>%
  mutate(traits=paste(trait1,"-",p2,trait3),
         pcor_lower_upper=paste0(round(ci.lower,digits=2),";",round(ci.upper,digits=2))
  ) %>%
  dplyr::select(locus,locus_id,
                traits,
                trait1,p2,trait3,
                n.snps,n.pcs,
                r2.phen1_z,r2.phen2_z,
                pcor, pcor_lower_upper,
                p) %>%
  distinct()
write.csv(pcor_clean,file="lava_pcor_cerebellum-cognitive_conditioned.csv",row.names = FALSE)
