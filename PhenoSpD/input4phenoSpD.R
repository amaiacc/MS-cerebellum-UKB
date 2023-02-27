# Create input files for phenoSpD by combining munge sumstats files
#---------------------------------------------------
library(dplyr); library(tidyr)
library(data.table)

#---------------------------------------------------
# clean workspace
rm(list=ls())
#---------------------------------------------------
options(stringsAsFactors = FALSE)
root="UKB_BIG40"
project="cerebellum_UKB"
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
primary_dir=paste(dir,"resources/datasets/GWAS_sumstats/downloaded_data/UKB/","BIG40","/",sep="")
ldsc_dir=paste(dir,"resources/datasets/GWAS_sumstats/ldsc/",root,"/",sep="")
out_dir=paste(dir,project,"data/phenoSpD/",sep="/")
if(!dir.exists(out_dir)){dir.create(out_dir)}
# set working dir
setwd(ldsc_dir)
#----------------------------------------------------------------------
# define IDPs to include
idps<-read.csv(paste(primary_dir,"IDPs_summary.csv",sep=""),header=TRUE) %>%
  filter(!is.na(Cerebellum)) %>% filter(measure=="VOLUME") %>%
  filter(grepl("cerebell",tolower(IDP_short_name))) %>%
  mutate(idp=if_else(nchar(pheno)==3,paste0("0",pheno),as.character(pheno)))
phenos2check<-c(idps$idp,"Chambers2022_TotalCerebellarVolume")
#----------------------------------------------------------------------
sumstats_files<-list.files(paste0(primary_dir,"/sumstats/"),pattern="txt.gz")
sumstats_files2<-list.files(gsub("BIG40","Chambers2022",primary_dir),pattern="gz")
# cannot use only munged sumstats because they have Z but not beta/se which are required by phenoSpD

# get SNPs from munged sumstats to subset a number of more managable SNPs
sumstats<-paste0(ldsc_dir,phenos2check[1],".sumstats.gz") %>% fread(.) %>% dplyr::select(SNP,A1,A2)

# read sequentially for all phenos, and combine
for (p in phenos2check){
  w<-grep(p,sumstats_files)
  if(length(w)>0){
    f<-paste0(primary_dir,"/sumstats/",sumstats_files[w])
    cat('Reading data', grep(p,phenos2check),'out of',length(phenos2check),'.\n')
    cat('File:',f,'\n')
    cat('-----------------------------------------------------------\n')
    tmp<-fread(f)
    tmp<-tmp %>% dplyr::select(rsid,a1,a2,beta,se) %>%
      filter(rsid %in% sumstats$SNP)
    
    # GWAS-es from the UKB-BIG40 were run using BGENIE
    # BGENIE: 'the regression model we code the first and second alleles as 0 and 1 respectively, so the beta coefficient refers to the effect of having an extra copy of the second allele.'
    # a_2 The character code for the first allele (a string). --> effect allele (because BIG40 released betas as direction of a_2!)
    # a_1 The character code for the second allele (a string). --> not-effect allele
    ## https://jmarchini.org/bgenie-usage/
    
    colnames(tmp)<- colnames(tmp) %>% 
      gsub("rsid","SNP",.) %>% gsub("a2","A1",.) %>% gsub("a1","A2",.) %>%
      gsub("^beta$",paste0(p,"_b"),.) %>%
      gsub("^se$",paste0(p,"_se"),.)
    sumstats<-merge(sumstats,tmp,all.x=TRUE)
    rm(f,tmp)
  } else {
    f<-paste0(gsub("BIG40","Chambers2022",primary_dir),sumstats_files2)
    cat('Reading data', grep(p,phenos2check),'out of',length(phenos2check),'.\n')
    cat('File:',f,'\n')
    cat('-----------------------------------------------------------\n')
    tmp<-fread(f)
    tmp<-tmp %>% dplyr::select(variant_id,effect_allele,other_allele,beta,standard_error) %>%
      filter(variant_id %in% sumstats$SNP)
    colnames(tmp)<- colnames(tmp) %>% 
      gsub("variant_id","SNP",.) %>% gsub("effect_allele","A1",.) %>% gsub("other_allele","A2",.) %>%
      gsub("^beta$",paste0(p,"_b"),.) %>%
      gsub("^standard_error$",paste0(p,"_se"),.)
    tmp<-tmp %>%
      mutate(A1=toupper(A1),
             A2=toupper(A2))
    sumstats<-merge(sumstats,tmp,all.x=TRUE)
    rm(f,tmp)
  }
  rm(w)
  gc()
}
rm(p)

# remove extra "_" from the colnames
colnames(sumstats)<-colnames(sumstats) %>% gsub("Chambers2022_","",.)
apply(sumstats,2,function(x) sum(is.na(x)))


# input for phenoSpD
sumstats4phenoSpD<-sumstats %>%
  rename(allele_0=A1,
         allele_1=A2) %>%
  mutate(across(contains("beta"),as.numeric)) %>%
  dplyr::select(-SNP)
rownames(sumstats4phenoSpD)<-sumstats$SNP
# convert all to numeric

tmp<-sumstats %>% head(100) %>%
  rename(allele_0=A1,
         allele_1=A2) %>%
  mutate(across(contains("beta"),as.numeric)) %>%
  dplyr::select(-SNP)
rownames(tmp)<-sumstats$SNP
# save
write.table(sumstats4phenoSpD,
            file=paste0(out_dir,"/","input4phenoSpD_",project,".table"),
            col.names = TRUE,row.names = TRUE,quote=FALSE)
