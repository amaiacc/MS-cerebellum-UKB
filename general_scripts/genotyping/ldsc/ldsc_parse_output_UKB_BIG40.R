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
ldsc_dir=paste(dir,"resources/datasets/GWAS_sumstats/ldsc/",root,"/output/",sep="")
# set working dir
setwd(ldsc_dir)
#----------------------------------------------------------------------
# some functions
source(paste(dir,"general_scripts/helper_functions.R",sep=""))
source(paste(dir,"general_scripts/genotyping/ldsc/helper_functions_ldsc.R",sep=""))
#----------------------------------------------------------------------
# read phenotype information
idps<-read.csv(paste(primary_dir,"IDPs_summary.csv",sep=""),header=TRUE)
idps$pheno<-numeric_nchar(idps$Pheno)
idps<-idps %>% mutate(measure=if_else(is.na(measure)&Units=="mm3","VOLUME",measure))
# lateralized measures
idps_lat<-read.table(paste(primary_dir,"IDPs_lat_summary.txt",sep=""),header=TRUE) %>% filter(is.na(NotLat)) %>% dplyr::select(-NotLat)
idps_lat$L[!is.na(idps_lat$L)]<-numeric_nchar(idps_lat$L[!is.na(idps_lat$L)])
idps_lat$R[!is.na(idps_lat$R)]<-numeric_nchar(idps_lat$R[!is.na(idps_lat$R)])
idps_lat<-merge(idps_lat,idps[,c("pheno","region","lobe")],by.x="L",by.y="pheno")

#----------------------------------------------------------------------
# read summary result form LDSC run and combine with IDP info
#----------------------------------------------------------------------
## heritability
h2<-read_ldsc_h2(h2file=paste(ldsc_dir,"summary_h2_ldsc.table",sep=""))
# combine with IDP info and select columns to keep
h2<-merge(h2,idps,by=c("pheno"),all.x=TRUE)
h2<-h2 %>% dplyr::select(pheno,IDP_short_name,parcellation,measure,lobe,region,hemisphere,N,h2,h2.SE,Intercept,Intercept.SE,Lambda.GC,h2.CI95upper,h2.CI95lower) %>%
  arrange(lobe,region,measure,parcellation) %>% filter(parcellation!="wg"|is.na(parcellation)) %>%
  mutate(IDP_short_name=if_else(is.na(IDP_short_name),pheno,IDP_short_name))
h2$p<-pchisq((h2$h2/h2$h2.SE)^2,1,lower.tail = FALSE) %>% format(scientific=TRUE,digits=2) %>% as.numeric()
h2$h2<-h2$h2 %>% format(digits=1)  %>% as.numeric()
h2$h2.SE<-h2$h2.SE %>% format(digits=2)  %>% as.numeric()
h2$measure<-as.factor(h2$measure)
#----------------------------------------------------------------------
## genetic correlation between left and right (for lateral measures)
## left and right
rg<-read_ldsc_rg(rgfile=paste(ldsc_dir,"summary_LRrg_ldsc.table",sep=""))
rg<-rg%>% rename(R=p1,L=p2)
# combine with IDP info and select columns to keep
rg<-merge(idps_lat,rg,by.x=c("L","R","IDP_short_name_lat"),by.y=c("L","R","file"),all=TRUE)
rg<-rg %>% dplyr::select(IDP_short_name_lat,parcellation,measure,lobe,region,L,R,rg,rg.SE,rg.CI95upper,rg.CI95lower,p,gcov_int,gcov_int_se) %>%
  arrange(lobe,region,measure,parcellation) %>% filter(parcellation!="wg")
rg$p2<-pchisq(((1-rg$rg)/rg$rg.SE)^2,1,lower.tail = FALSE) %>% format(scientific=TRUE,digits=2) %>% as.numeric() # statistical significance for different to 1
# rg$rg<-rg$rg %>% format(digits=3)  %>% as.numeric()
# rg$rg.SE<-rg$rg.SE %>% format(digits=2)  %>% as.numeric()
rg$measure<-as.factor(rg$measure)
rgLR<-rg
rm(rg)
#----------------------------------------------------------------------
## rg within parcellation and measure type (left and right)
rg_all<-read_ldsc_rg(rgfile=paste(ldsc_dir,"summary_rg_ldsc.table",sep=""))
# combine with IDP info and select columns to keep
rg_all<-merge(rg_all,idps,by.x=c("p1"),by.y=c("pheno"),all.x=TRUE) %>% filter(p1!=p2)
# combine to get info for second, given that the category/measure type is the same
rg_all2<-merge(rg_all,idps,
                by.x=c("p2"),
                by.y=c("pheno"),
                suffixes=c(".p1",".p2"),all.x=TRUE) %>% 
  mutate(parcellation=if_else(parcellation.p1==parcellation.p2,parcellation.p1,"Different parcellation")) %>%
  filter(!is.na(rg)) %>%
  # keep hemisphere as such, if it's the same for both regions
  mutate(hemisphere=if_else(hemisphere.p1==hemisphere.p2,
                            hemisphere.p1,"contralateral")) %>%
  # for same measure L and R, keep only one instance, so p1=L, p2=R
  filter(!( (region.p1==region.p2)&(hemisphere.p1=="R"&hemisphere.p2=="L")) |
           (is.na(region.p1)|is.na(region.p2)))

# define the measures within the same atlas
rg_atlas<- rg_all2 %>%
  filter(parcellation!="Different parcellation") %>%
  dplyr::select(-parcellation.p1,-parcellation.p2) %>%
  # keep only if the measure type is the same
  filter(measure.p1==measure.p2) %>%
  mutate(measure=measure.p1) %>% dplyr::select(-measure.p2,-measure.p1)
#
rg_atlas<-rg_atlas %>% dplyr::select(parcellation,
                              measure,
                              contains("IDP_short_name"),contains("hemisphere"),contains("region"),
                              p1,p2,
                              rg,rg.SE,rg.CI95upper,rg.CI95lower,p,gcov_int,gcov_int_se)  %>%
  arrange(measure,parcellation) %>% filter(parcellation!="wg")
# rg_atlas$rg<-rg_atlas$rg %>% format(digits=1)  %>% as.numeric()
# rg_atlas$rg.SE<-rg_atlas$rg.SE %>% format(digits=2)  %>% as.numeric()
rg_atlas$measure<-as.factor(rg_atlas$measure)

# rg between same measure (e.g. vol) but different parcellations
rg_not_atlas <- rg_all2 %>% 
  filter(parcellation=="Different parcellation") %>%
  # keep only if the measure type is the same
  filter(measure.p1==measure.p2) %>%
  mutate(measure=measure.p1) %>% dplyr::select(-measure.p2,-measure.p1)

# rg of one UKB measure with some other type of measure
rg_other <- rg_all2 %>% 
  filter(is.na(parcellation))


#----------------------------------------------------------------------

#----------------------------------------------------------------------
## rg with other measures
rg_files<-list.files(ldsc_dir,pattern="rg.*.table")
rg_files<-rg_files[grep("LRrg|summary_rg_ldsc.table",invert=TRUE,rg_files)]

for (rg_f in rg_files){
  # f<-strsplit(rg_f,"_") %>% sapply("[[",2)
  f<-gsub("summary_|_rg_ldsc.table","",rg_f)
  rg<-read_ldsc_rg(rgfile=paste(ldsc_dir,rg_f,sep=""))
  # combine with IDP info and select columns to keep
  rg<-merge(rg,idps,by.x=c("p1"),by.y=c("pheno"),all.x=TRUE)
  
  rg<-rg %>% dplyr::select(p1,IDP_short_name,parcellation,hemisphere,measure,region,p2,rg,rg.SE,rg.CI95upper,rg.CI95lower,p,gcov_int,gcov_int_se)  %>%
    arrange(measure,parcellation)
  # rg$rg<-rg$rg %>% format(digits=1)  %>% as.numeric()
  # rg$rg.SE<-rg$rg.SE %>% format(digits=2)  %>% as.numeric()
  rg$measure<-as.factor(rg$measure)
  write.csv(rg,paste(ldsc_dir,"ldsc_rg_",f,"_UKB_BIG40.csv",sep=""),row.names=FALSE)
  rm(rg)
}


#----------------------------------------------------------------------
write.csv(h2,paste(ldsc_dir,"ldsc_h2_UKB_BIG40.csv",sep=""),row.names=FALSE)
write.csv(rgLR,paste(ldsc_dir,"ldsc_rgLR_UKB_BIG40.csv",sep=""),row.names=FALSE)
write.csv(rg_atlas,paste(ldsc_dir,"ldsc_rg_atlas_UKB_BIG40.csv",sep=""),row.names=FALSE)
write.csv(rg_not_atlas,paste(ldsc_dir,"ldsc_rg_not_atlas_UKB_BIG40.csv",sep=""),row.names=FALSE)
write.csv(rg_other,paste(ldsc_dir,"ldsc_rg_other_UKB_BIG40.csv",sep=""),row.names=FALSE)

# write.csv(rg_all,paste(ldsc_dir,"ldsc_rg_all_UKB_BIG40.csv",sep=""),row.names=FALSE)
#--------------------------------------------