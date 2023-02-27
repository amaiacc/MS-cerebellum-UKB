# look at association profiles where signif. local genetic correlations are
# to run after: local_rg_regions_prepare.R
#---------------------------------------------------
library(dplyr); library(tidyr)
library(data.table)
library("gassocplot")
library(ggplot2)

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
lava_dir=paste(dir,project,"data/lava/output/",sep="/")
out_dir=paste(dir,project,"data/viz_loci/",sep="/")
if(!dir.exists(out_dir)){dir.create(out_dir)}

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
# read loci of interest
loci<-read.csv(paste0(lava_dir,"lava_bivar_all_traits_clean.csv"),stringsAsFactors = FALSE) %>% 
  mutate(chr=strsplit(locus_id,":") %>% sapply("[[",1) %>% as.numeric(),
         tmp=strsplit(locus_id,":") %>% sapply("[[",2),
         start=strsplit(tmp,"-") %>% sapply("[[",1) %>% as.numeric(),
         end=strsplit(tmp,"-") %>% sapply("[[",2) %>% as.numeric(),
  ) %>%
  dplyr::select(-tmp)
loci2check<-loci %>% filter(p.adjBF<0.05) %>% pull(locus_id) %>% unique()

# read combined sumstats for all traits
sumstats<-fread(paste0(out_dir,"sumstats_hm3_all_phenos.csv")) %>% distinct() %>%
  mutate(chr=)
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# make plots!
#----------------------------------------------------------------------
setwd(out_dir)
l= "4:102544804-104384534" # "4:102544804-104384534"
l="14:57460782-58447798"
for (l in loci2check){
  loc=gsub("-|:","_",l)
  # get z-scores for phenos of interest
  phenos2check <- loci %>% filter(locus_id==l & p.adjBF<0.05) %>% 
    dplyr::select(phen1,phen2) %>% unlist() %>% unique() %>%
    gsub("146","0146",.) %>%
    gsub("143","0143",.) %>%
    gsub("102","0102",.) %>%
    gsub("104","0104",.) %>%
    gsub("Chambers2022_TotalCerebellarVolume","TotalCerebellarVolume",.) %>%
    gsub("PGC3_SCZ_primary","SCZ",.) %>%
    gsub("Lee2018_CP","CP",.) %>%
    factor(.,levels=c("TotalCerebellarVolume","0143","0146",
              "0102","0104",
              "CP","SCZ")) %>% sort()
  names_phenos<-phenos2check %>%
    gsub("0146","Crus I (R)",.) %>%
    gsub("0143","Lobule VI (R)",.) %>%
    gsub("0102","Temporal Occipital Fusiform Cortex (L)",.) %>%
    gsub("0104","Occipital Fusiform Gyrus (L)",.) %>%
    gsub("TotalCerebellarVolume","Total Cerebellar Volume",.)
  
  if(file.exists(paste0("LDmatrix_corr_",loc,".txt"))){
    # read corr
    corr_locus<-read.table(paste0("LDmatrix_corr_",loc,".txt"),header=TRUE)
    corr<-corr_locus[,-1]%>% as.matrix()
    corr<-(corr)^(1/2)
    # remove markers from list if not available in corr
    markers<-read.table(paste0("Markers_corr_",loc,".txt"),header=TRUE) %>% 
      filter(marker %in% corr_locus$RS_number) %>% distinct()
    z<-sumstats %>% filter(SNP %in% markers$marker) %>%
      arrange(chr,pos) %>%
      dplyr::select(SNP,paste0(phenos2check,"_z"))
    # ensure that the marker order is the same
    z<-z[match(markers$marker,z$SNP),] %>% distinct()
    rownames(z)<-z$SNP
    # keep only z values
    z<-as.matrix(z[,-1])
    p<-stack_assoc_plot(markers=markers,
                        z = z,
                        corr=corr,
                        traits=names_phenos,
    )
    ggsave(p,file=paste0(loc,"_all_phenos_BF0.05.png"),width=7,height=length(phenos2check)*3)
    rm(corr_locus,corr)
    } else {
    markers<-sumstats %>% dplyr::select(SNP,chr,pos) %>%
      filter(chr==as.numeric(strsplit(loc,"_") %>% sapply("[[",1)) &
               pos>=as.numeric(strsplit(loc,"_") %>% sapply("[[",2)) &
               pos<=as.numeric(strsplit(loc,"_") %>% sapply("[[",3))) %>%
      dplyr::select(SNP,chr,pos) %>%
      rename(marker=SNP) %>% distinct()
    z<-sumstats %>% filter(SNP %in% markers$marker) %>%
        arrange(chr,pos) %>%
        dplyr::select(SNP,matches(paste0(phenos2check,collapse="|"),perl=TRUE)) %>%
        dplyr::select(SNP,contains("_z"))
    # ensure that the marker order is the same
    z<-z[match(markers$marker,z$SNP),] %>% distinct()
    rownames(z)<-as.character(z$SNP)
    # keep only z values
    z<-as.matrix(z[,-1])
    p<-stack_assoc_plot(markers=markers,
                        z = z,
                        traits=names_phenos,
    )
    ggsave(p,file=paste0(loc,"_all_phenos_BF0.05_noCorr.png"),width=7,height=length(phenos2check)*3)
    }
    rm(markers,z)
    rm(p,loc,phenos2check,names_phenos)
}

#----------------------------------------------------------------------
# create latex file to concatenate all those plots
loci2check_cer<-loci %>% filter(p.adjBF<0.05) %>% 
  filter(phen2=="143"|phen2=="146") %>% pull(locus_id) %>% unique()
loci2check_fusif<-loci %>% filter(p.adjBF<0.05) %>% 
  filter(phen2=="102"|phen2=="104") %>% pull(locus_id) %>% unique()
loci2check_cp_scz<-loci %>% filter(p.adjBF<0.05) %>% 
  filter(phen2=="Lee2018_CP"|phen2=="PGC3_SCZ_primary") %>% pull(locus_id) %>% unique()

sink("local_rg_regions_BF_viz.tex")
# cerebellum rg's first
cat("% cerebellum\n")
for (l in loci2check_cer){
  loc=gsub("-|:","_",l)
  file=list.files(pattern=loc)
  file=file[grep(".png",file)]
  cat("\\begin{figure}\n")
  cat("\\centering\n")
  cat(paste0("\\includegraphics[height=0.9\\textheight]{figures/local_rg_viz_loci/",file,"}\n"))
  cat("\\caption{Regional association plot for traits with a significant local genomic correlation in locus",l,", (LAVA). HapMap3 SNPs are shown.}\n")
  cat(paste0("\\label{fig:",loc,"}\n"))
  cat("\\end{figure}\n")
  cat("\n")
}

#
cat("% fusiform\n")
for (l in loci2check_fusif[!loci2check_fusif %in% loci2check_cer]){
  loc=gsub("-|:","_",l)
  file=list.files(pattern=loc)
  file=file[grep(".png",file)]
  cat("\\begin{figure}\n")
  cat("\\centering\n")
  cat(paste0("\\includegraphics[height=0.9\\textheight]{figures/local_rg_viz_loci/",file,"}\n"))
  cat("\\caption{Regional association plot for traits with a significant local genomic correlation in locus",l,", (LAVA). HapMap3 SNPs are shown.}\n")
  cat(paste0("\\label{fig:",loc,"}\n"))
  cat("\\end{figure}\n")
  cat("\n")
}

#
cat("% CP and SCZ\n")
for (l in loci2check_cp_scz[!loci2check_cp_scz %in% loci2check_cer]){
  loc=gsub("-|:","_",l)
  file=list.files(pattern=loc)
  file=file[grep(".png",file)]
  cat("\\begin{figure}\n")
  cat("\\centering\n")
  cat(paste0("\\includegraphics[height=0.9\\textheight]{figures/local_rg_viz_loci/",file,"}\n"))
  cat("\\caption{Regional association plot for traits with a significant local genomic correlation in locus",l,", (LAVA). HapMap3 SNPs are shown.}\n")
  cat(paste0("\\label{fig:",loc,"}\n"))
  cat("\\end{figure}\n")
  cat("\n")
}
sink()

