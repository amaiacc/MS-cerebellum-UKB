library("LDlinkR")

# preparation prior to looking at association profiles where signif. local genetic correlations are
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
primary_dir=paste0(dir,"resources/datasets/GWAS_sumstats/QC/")
ldsc_dir=paste(dir,"resources/datasets/GWAS_sumstats/ldsc/",sep="")
lava_dir=paste(dir,project,"data/lava/output/",sep="/")
out_dir=paste(dir,project,"data/viz_loci/",sep="/")
if(!dir.exists(out_dir)){dir.create(out_dir)}

setwd(primary_dir)

#----------------------------------------------------------------------
# check significant genetic correlation signals (i.e. as in Table S8, lava_bivar_wide_all_traits_BFsig.csv)
#----------------------------------------------------------------------
# define IDPs to include
idps<-read.csv(paste(primary_dir,"../downloaded_data/UKB/BIG40/IDPs_summary.csv",sep=""),header=TRUE) %>%
  filter(!is.na(Cerebellum)) %>% filter(measure=="VOLUME") %>%
  filter(grepl("cerebell",tolower(IDP_short_name))) %>%
  mutate(idp=if_else(nchar(pheno)==3,paste0("0",pheno),as.character(pheno)))
idps<-subset(idps,(region=="CEREBELLUM_VI"&hemisphere=="R")|
               (region=="CEREBELLUM_CRUS_I"&hemisphere=="R")
               )
phenos2check<-c(idps$idp,
                "0102","0104",
                "Chambers2022_TotalCerebellarVolume","SCZ","CP")
#----------------------------------------------------------------------
sumstats_files1<-list.files(paste0(primary_dir,"/UKB_BIG40/"),pattern="QC.txt.gz")
sumstats_files1<-paste0(paste0(primary_dir,"/UKB_BIG40/"),
                        sumstats_files1[grep(paste0(phenos2check,collapse="|"),sumstats_files1)]
                        )
sumstats_files2<-list.files(primary_dir,pattern="QC.gz")
sumstats_files2<-paste0(primary_dir,
                        sumstats_files2[grep(paste0(phenos2check,collapse="|"),sumstats_files2)]
                        )
#
# sumstats_files1<-list.files(paste0(ldsc_dir,root),pattern="sumstats.gz")
# sumstats_files1<-sumstats_files1[grep(paste0(phenos2check,collapse="|"),sumstats_files1)]
# sumstats_files2<-list.files(ldsc_dir,pattern="sumstats.gz")
# sumstats_files2<-sumstats_files2[grep(paste0(phenos2check,collapse="|"),sumstats_files2)]

sumstats_files<-c(sumstats_files1,sumstats_files2) %>% unique()
# rm(sumstats_files1,sumstats_files2)
#----------------------------------------------------------------------
# read sumstats
# get SNPs from munged sumstats to subset a number of more managable SNPs
# sumstats<-paste0(ldsc_dir,root,"/",phenos2check[1],".sumstats.gz") %>% fread(.) %>% dplyr::select(SNP,A1,A2)

# read sequentially for all phenos, and combine
for (f in sumstats_files1){
  # if(file.exists(paste0(ldsc_dir,root,"/",p))){f<-paste0(ldsc_dir,root,"/",p)
  # } else {f<-paste0(ldsc_dir,"/",p)}
  cat('Reading data', grep(f,sumstats_files1),'out of',length(sumstats_files1),'.\n')
  cat('File:',f,'\n')
  cat('-----------------------------------------------------------\n')
  # GWAS-es from the UKB-BIG40 were run using BGENIE
  # BGENIE: 'the regression model we code the first and second alleles as 0 and 1 respectively, so the beta coefficient refers to the effect of having an extra copy of the second allele.'
  # a_2 The character code for the first allele (a string). --> effect allele (because BIG40 released betas as direction of a_2!)
  # a_1 The character code for the second allele (a string). --> not-effect allele
  ## https://jmarchini.org/bgenie-usage/
  p2<- gsub(primary_dir,"",f) %>% gsub(root,"",.) %>%
    gsub(".sumstats.gz|.QC|.txt.gz","",.) %>% gsub("/","",.)
  
  tmp<-fread(f)
  
  tmp<-tmp %>% dplyr::select(rsid,chr,pos,a_1,a_2,beta,se) 
  # %>%  filter(rsid %in% sumstats$rsid)
  
  colnames(tmp)<- colnames(tmp) %>% 
    gsub("rsid","SNP",.) %>% gsub("a_2","A1",.) %>% gsub("a_1","A2",.) %>%
    gsub("^Z$",paste0(p2,"_Z"),.) %>%
    gsub("^N$",paste0(p2,"_N"),.) %>%
    gsub("^beta$",paste0(p2,"_b"),.) %>%
    gsub("^se$",paste0(p2,"_se"),.)
  

  # combine
  if(exists("sumstats")){
    sumstats<-merge(sumstats,tmp,all.x=TRUE) 
  } else { 
    sumstats<-tmp
  }
  # check dimensions of sumstats
  dim(sumstats)
  # clean intermediate
  rm(p2,tmp)
  gc()
}

sumstats<-sumstats %>% 
  dplyr::select(SNP,chr,pos,A1,A2,contains("_b"),contains("_se")) %>%
  mutate(chr=as.numeric(chr)) %>% 
  filter(!is.na(chr))

#----------------------------------------------------------------------
# load separately/manually, and combine
#-------------------
## Total Cerebellar Volume
p2="TotalCerebellarVolume"
f<-paste0(primary_dir,"../downloaded_data/UKB/Chambers2022/GCST90020190_buildGRCh37.tsv.gz")
tmp<-fread(f)
tmp<-tmp %>% dplyr::select(variant_id,chromosome,base_pair_location,effect_allele,other_allele,beta,standard_error)
colnames(tmp)<- colnames(tmp) %>% 
  gsub("variant_id","SNP",.) %>% 
  gsub("chromosome","chr",.) %>%
  gsub("base_pair_location","pos",.) %>%
  gsub("effect_allele","A1",.) %>% gsub("other_allele","A2",.) %>%
  gsub("^Z$",paste0(p2,"_Z"),.) %>%
  gsub("^N$",paste0(p2,"_N"),.) %>%
  gsub("^beta$",paste0(p2,"_b"),.) %>%
  gsub("^standard_error$",paste0(p2,"_se"),.)

tmp<- tmp %>% mutate(A1=toupper(A1),A2=toupper(A2))
# combine
sumstats<-merge(sumstats,tmp,all=TRUE)
rm(p2,f,tmp)
#-------------------
## CP
p2="CP"
f<-list.files(primary_dir,pattern=p2)
f<-paste0(primary_dir,f[grep("QC.gz",f)])
tmp<-fread(f)
tmp<-tmp %>% dplyr::select(rsid,chr,pos,a_1,a_0,beta,se)

colnames(tmp)<- colnames(tmp) %>% 
  gsub("rsid","SNP",.) %>% gsub("a_1","A1",.) %>% gsub("a_0","A2",.) %>%
  gsub("^Z$",paste0(p2,"_Z"),.) %>%
  gsub("^N$",paste0(p2,"_N"),.) %>%
  gsub("^beta$",paste0(p2,"_b"),.) %>%
  gsub("^se$",paste0(p2,"_se"),.)
tmp<- tmp %>% mutate(A1=toupper(A1),A2=toupper(A2))

# combine
sumstats<-merge(sumstats,tmp,all=TRUE)
rm(p2,f,tmp)
#-------------------
## SCZ **hemen nago**
p2="SCZ"
f<-paste0(primary_dir,"../downloaded_data/PGC/SCZ/SCZ2022/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz")
tmp<-fread(f)
tmp<-tmp %>% dplyr::select(ID,CHROM,POS,A1,A2,BETA,SE)

colnames(tmp)<- colnames(tmp) %>% 
  gsub("ID","SNP",.) %>% 
  gsub("CHROM","chr",.) %>%
  gsub("POS","pos",.) %>%
  gsub("^Z$",paste0(p2,"_Z"),.) %>%
  gsub("^N$",paste0(p2,"_N"),.) %>%
  gsub("^BETA$",paste0(p2,"_b"),.) %>%
  gsub("^SE$",paste0(p2,"_se"),.)

tmp <- tmp %>% mutate(chr=as.numeric(chr)) %>% 
  filter(!is.na(chr))

# combine
sumstats<-merge(sumstats,tmp,all=TRUE,by=c("SNP","chr","pos"),
                suffixes = c("",".scz"))
rm(p2,f,tmp)
#-------------------
#----------------------------------------------------------------------
sumstats<-as.data.frame(sumstats) %>% distinct()
# compute z-scores
phenos<-colnames(sumstats)[grep("_b",colnames(sumstats))] %>% gsub("_b","",.)
for(p in phenos){
  sumstats[,paste0(p,"_z")]<-sumstats[,paste0(p,"_b")]/sumstats[,paste0(p,"_se")]
}
# get all markers
markers_all<-sumstats %>% dplyr::select(SNP,chr,pos) %>% arrange(chr,pos) %>%
  mutate(marker=SNP)
markers_all<-as.data.frame(markers_all) 

# save combined file for future!
fwrite(sumstats,file=paste0(out_dir,"sumstats_all_phenos.csv"),row.names=FALSE)
fwrite(markers_all,file=paste0(out_dir,"markers_all_phenos.csv"),row.names=FALSE)
# subset only hm3 markers
w_hm3_snps<-read.table(paste0(dir,"resources/LDscores/w_hm3.snplist"),header=TRUE)
sumstats_hm3<-subset(sumstats, SNP %in% w_hm3_snps$SNP)
markers_hm3<-subset(markers_all, SNP %in% w_hm3_snps$SNP)
fwrite(sumstats_hm3,file=paste0(out_dir,"sumstats_hm3_all_phenos.csv"),row.names=FALSE)
fwrite(markers_hm3,file=paste0(out_dir,"markers_hm3_all_phenos.csv"),row.names=FALSE)


rm(sumstats); gc()
#----------------------------------------------------------------------
# get correlation using LDmatrix() function, for each locus of interest
LDlink_token="660aece5746d"
# read loci of interest
loci2check<-read.csv(paste0(lava_dir,"lava_bivar_wide_all_traits_BFsig.csv"),stringsAsFactors = FALSE) %>% 
  dplyr::select(locus_id) %>%
  mutate(chr=strsplit(locus_id,":") %>% sapply("[[",1) %>% as.numeric(),
         tmp=strsplit(locus_id,":") %>% sapply("[[",2),
         start=strsplit(tmp,"-") %>% sapply("[[",1) %>% as.numeric(),
         end=strsplit(tmp,"-") %>% sapply("[[",2) %>% as.numeric(),
  ) %>%
  dplyr::select(-tmp)

# get correlation matrix for each locus
for (l in loci2check$locus_id){
  locus_file=paste0(out_dir,"LDmatrix_corr_",l,".txt") %>% gsub("-|:","_",.)
  if(!file.exists(locus_file)){
    cat("Get correlation matrix for:",l,"\n")
    locus=subset(loci2check,locus_id==l)
    markers_locus<-subset(markers_hm3,
                          chr==as.numeric(locus$chr) &
                            pos>=as.numeric(locus$start) &
                            pos<=as.numeric(locus$end)) %>%
      dplyr::select(marker,chr,pos)
    markers_locus<-markers_locus[grep("rs",markers_locus$marker),]
    if(NROW(markers_locus)<=1000){
      corr_locus<-LDmatrix(snps=markers_locus$marker,token=LDlink_token)
      
      # remove markers from list if not available in corr
      markers_locus<-markers_locus %>% filter(marker %in% corr_locus$RS_number)
      
      # save as output
      write.table(corr_locus,file=locus_file,sep="\t",row.names=FALSE)
      write.table(markers_locus,file= gsub("LDmatrix_","Markers_",locus_file),sep="\t",row.names=FALSE)
      # clean intermediate
      rm(corr_locus)
    }
    rm(locus,markers_locus)
    
  } else {
    cat("File already exists:",locus_file,"\n")
  }
  rm(locus_file)
}

rm(LDlink_token)
