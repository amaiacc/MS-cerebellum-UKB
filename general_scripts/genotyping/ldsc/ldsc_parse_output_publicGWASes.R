# clean workspace
rm(list=ls())
# libraries and custom functions
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
library(ggrepel)
library(dplyr)
library(tidyr)

#---------------------------------------------------
options(stringsAsFactors = FALSE)
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
ldsc_dir=paste(dir,"resources/datasets/GWAS_sumstats/ldsc/",sep="")
# set working dir
setwd(ldsc_dir)
#----------------------------------------------------------------------
# some functions
source(paste(dir,"general_scripts/helper_functions.R",sep=""))
source(paste(dir,"general_scripts/genotyping/ldsc/helper_functions_ldsc.R",sep=""))
## Convert observed h2 to liability scale
## We need the inverse complementary error function to do the conversion
## Here's the function
erfcinv <- function(x) { qnorm(x/2, lower = FALSE)/sqrt(2) }

## https://github.com/saralpulit/Afib-Stroke-Overlap/blob/master/h2.obs2liab.R
## prevalence =  estimated trait prevalence
## cases = total number of cases, controls = total number of controls
## h2obs = h2 on the observed scale
h2liab <- function(prevalence, cases, controls, h2obs) {
  
  P <- cases/(cases + controls)
  Z <- sqrt(2)*erfcinv(2*prevalence)
  f <- 1/sqrt(2*pi) * exp(-Z^2/2)
  obs_to_liab <- (prevalence*(1-prevalence))^2/(P*(1-P))/f^2
  h2.liab <- h2obs * obs_to_liab
  
  return(h2.liab)
  
}
#----------------------------------------------------------------------
# Read data
gwas_info<-read.csv(paste(ldsc_dir,"/../","PublicSummary_stats_overview_info.csv",sep=""))
gwas_info<-gwas_info %>% 
  mutate(FirstAuthor=strsplit(First.Author," ") %>% sapply("[[",1) %>% gsub("-","",.)) %>%
  mutate(ref=paste0(FirstAuthor,Year)) %>%
  mutate(Trait_name_ref=paste(Trait_name,ref,sep="\n")) %>%
  mutate(pheno1=paste(Trait_name,"_",FirstAuthor,Year,sep=""),
         pheno2=paste(FirstAuthor,Year,"_",Trait_name,sep=""),
         pheno3=paste(FirstAuthor,"_",Trait_name,sep="")
         )

#----------------------------------------------------------------------
# read summary result form LDSC run and combine with GWAS info
## heritability
h2<-read_ldsc_h2(h2file=paste(ldsc_dir,"summary_h2_ldsc.table",sep=""))
# rename some files to match the reference..
h2<- h2 %>% mutate(
  pheno2= gsub("PASS_|_Hujoel2019","",pheno) %>%
    gsub("Autism","ASD_",.) %>%
    gsub("ASD_","ASD_PGCCrossTrait2013",.) %>%
    gsub("__","_",.)
)

#----------------------------------------------------------------------
# combine with info
h2_1<-merge(h2,gwas_info,by.x=c("pheno2"),by.y="pheno1")
h2_2<-merge(h2,gwas_info,by.x=c("pheno2"),by.y="pheno2")
h2_3<-merge(h2,gwas_info,by.x=c("pheno2"),by.y="pheno3")
# combine all
h2_all<-merge(h2_1,h2_2,all=TRUE)
h2_all<-merge(h2_all,h2_3,all=TRUE)
# clean intermediate
rm(h2_1,h2_2,h2_3)
rm(h2)
#----------------------------------------------------------------------
# select columns to keep
h2_all<-h2_all %>% select(ref,pheno,Trait_name,Trait_name_ref,Trait,Reported.trait,
                          PopPrevalence,
                          ref_h2,ref_se,
                          N,Ncases,Ncontrols,
                          h2,h2.SE,Intercept,Intercept.SE,Lambda.GC,h2.CI95upper,h2.CI95lower)

h2_all$p<-pchisq((h2_all$h2/h2_all$h2.SE)^2,1,lower.tail = FALSE) %>% format(scientific=TRUE,digits=2) %>% as.numeric()
h2_all$h2<-h2_all$h2 %>% format(digits=1)  %>% as.numeric()
h2_all$h2.SE<-h2_all$h2.SE %>% format(digits=2)  %>% as.numeric()
h2_all$PopPrevalence<-as.numeric(h2_all$PopPrevalence)
h2_all$ref_h2<-as.numeric(h2_all$ref_h2)
# convert binary estimates to liability scale
h2_all$h2liab<-NA
h2_all$h2liab[which(h2_all$Trait=="Binary")]<- sapply(which(h2_all$Trait=="Binary"),function(w){
  h2liab(prevalence = h2_all$PopPrevalence[w],
         cases=h2_all$Ncases[w],
         controls=h2_all$Ncontrols[w],
         h2obs = h2_all$h2[w]
         )
} )

h2_all$estimate<-h2_all$h2
h2_all$estimate[which(h2_all$Trait=="Binary")]<-h2_all$h2liab[which(h2_all$Trait=="Binary")]
#----------------------------------------------------------------------
# plot
h2_plot<-ggplot(h2_all %>% filter(!is.na(ref_h2)), aes(x=estimate,y=ref_h2,color=)) + 
  geom_point() +
  geom_text_repel(aes(label=Trait_name_ref)) +
  ylab("Heritability estimate in reference publication") +
  xlab("Heritability estimate (LDSC) from public sumstats") +
  xlim(c(0,0.3)) + ylim(c(0,0.3)) +
  geom_abline()


#----------------------------------------------------------------------
# save output
write.csv(h2_all,paste(ldsc_dir,"ldsc_h2_publicGWASes.csv",sep=""),row.names=FALSE)
ggsave(h2_plot,file=paste(ldsc_dir,"ldsc_h2_comparison_publicGWASes.png",sep=""))
#--------------------------------------------