# running genomic SEM analysis for cerebellar substructures
# following: https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects
#----------------------------------------------------------

#----------------------------------------------------------
library(dplyr)
library(tidyr)
# library(devtools)
# install_github("GenomicSEM/GenomicSEM")
require(GenomicSEM)
library(lavaan)
require(Matrix)
require(stats)
# for plotting
library(ggplot2); library(ggrepel)
library(cowplot)
#----------------------------------------------------------
# function to run CFA models
run_cfa_model<-function(nFactors,thr,ldscoutput=LDSCoutput,EFA){
  cat("---------------------------------------------\n")
  cat(paste0("Running CFA model after EFA ",EFA,"\nN factors=",nFactors,"; loading threshold=",thr,"\n"))
  cat("---------------------------------------------\n")
  
  loadings<-get(EFA)[[nFactors]]
  f_loadings<-sapply(1:NCOL(loadings),function(i){
    w<-which(abs(loadings[,paste0("Factor",i)])>thr) %>% names()
    return(w)
  })
  names(f_loadings)<-paste0("F",1:NCOL(loadings))
  
  # check which manifest variables are not assigned to any factor
  traits2assign<-trait.names2[!trait.names2 %in% (unlist(f_loadings) %>% unique())]
  traits2assignF<-sapply(traits2assign,function(x) {
    l<-loadings[x,]
    maxF<-which(l==max(l)) %>% names() %>% gsub("Factor","F",.)
    return(maxF)
    })
  
  # assign traits to corresponding F
  for (f in names(f_loadings)){
    w<-which(traits2assignF==f) %>% names()
    if(length(w)>0){
      f_loadings[[f]]<-c(f_loadings[[f]],w)
    }
    rm(w)
  }
  rm(f)
  rm(traits2assignF,traits2assign)
  
  # remove factor if there are no measures with loadings that are higher than the threshold
  w0<-which(lapply(f_loadings,length)==0)
  if(length(w0)>0){
    f_loadings[[w0]]<-NULL
  }
  
  nFactors2<-length(f_loadings)
  
  # generate formula (generic!)
  nFactorModel <- lapply(1:length(f_loadings),function(f){
    paste0("F",f, " =~ NA*", paste(f_loadings[[f]],collapse=" + "),"\n")
  }) %>% do.call("paste0",.)
  
  # model_part1<-paste0("F",1:nFactors) %>% combn(.,m=2) %>% apply(.,2,function(t) paste0(t,collapse=" ~~ 0*")) %>% paste(.,collapse="\n")
  # 
  # nFactorModel<-paste0(nFactorModel,"\n",model_part1)
  
  cat("Run model:\n")
  print(nFactorModel)
  
  #run the model
  CFA_nFactors <-usermodel(LDSCoutput, estimation = "DWLS", model = nFactorModel, 
                           CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
  
  # save output if model converged
  if(!is.null( CFA_nFactors )){
    #
    output<-CFA_nFactors$modelfit
    
    # identify heywood cases
    vars2constrain<-CFA_nFactors$results %>% filter(STD_All<(-1)|(STD_All>1)) %>% pull(rhs) %>% unique() %>% paste(.,collapse=",")
    
    if(length(vars2constrain)>0){
      output$vars2constrain_Heywood<-vars2constrain
      } else {
      output$vars2constrain_Heywood<-""
      }

  } else{
    output<-data.frame(matrix(nrow=1,ncol=7))
    colnames(output)<-c("chisq","df","p_chisq","AIC","CFI","SRMR","vars2constrain_Heywood")
  }
  # add additional parameters
  output$nFactors<-nFactors2
  output$FactorLoadingThreshold<-thr
  output$model<-nFactorModel

  # return output summary table
  return(output)
}
#----------------------------------------------------------


#----------------------------------------------------------
if(analysis_run=="all"){
  idps<-idps %>% filter(parcellation=="FAST")
  traits <- paste0(paste0("0",idps$Pheno),".sumstats.gz")
  trait.names<-idps$IDP_short_name
  # traits <- paste0(c(paste0("0",idps$Pheno),"Chambers2022_TotalCerebellarVolume"),".sumstats.gz")
  # trait.names<-c(idps$IDP_short_name,"TotalCerebellarVolume")
  } else {
    if(analysis_run=="VR"){
      idps<-idps %>% filter((hemisphere=="R"|is.na(hemisphere))&parcellation=="FAST")
    } else if(analysis_run=="VL"){
      idps<-idps %>% filter((hemisphere=="L"|is.na(hemisphere))&parcellation=="FAST")
    }
    traits <- paste0(paste0("0",idps$Pheno),".sumstats.gz")
    trait.names<-idps$IDP_short_name
  }
# make sure illegal characters are not included
trait.names2 <- trait.names %>% gsub("-","_",.) %>%
  gsub("IDP_T1_FAST_ROIs_|cerebellum_|aseg_|volume_","",.) %>%
  gsub("^L_|lh","left",.) %>%
  gsub("^R_|rh","right",.) %>%
  gsub("^V_","vermis",.) %>% 
  gsub("crus_I","CrusI",.) %>%
  gsub("I_IV","I.IV",.) %>%
  gsub("_","",.)
#----------------------------------------------------------
## Step 2: Run multivariable LDSC
# define sumstats files from munging for LDSC
# traits <- paste0("0",idps$Pheno,".sumstats.gz")
sample.prev <- rep(NA,times=length(traits))
population.prev <- rep(NA,times=length(traits))
ld<-paste0(dir,"resources/LDscores/eur_w_ld_chr/")
wld <- paste0(dir,"resources/LDscores/eur_w_ld_chr/")

#----------------------------------------------------------
# read combined object for analysis run
filename_run=paste0(out_dir,paste0("LDSC_cerebellum_",analysis_run,".RData"))
if(!file.exists(filename_run)){
  LDSCoutput <- ldsc(traits, 
                     sample.prev, 
                     population.prev, 
                     ld, 
                     wld, 
                     trait.names2)
  
  
  save(LDSCoutput, file=filename_run)
  gc()
} else {
 load(filename_run)
}
rm(filename_run)
#----------------------------------------------------------
# The output (named LDSCoutput here) is a list with 5 named variables in it:
#   
# LDSCoutput$S is the covariance matrix (on the liability scale for case/control designs).
LDSCoutput$S
# LDSCoutput$V which is the sampling covariance matrix in the format expected by lavaan.
LDSCoutput$V
# LDSCoutput$I is the matrix of LDSC intercepts and cross-trait (i.e., bivariate) intercepts.
LDSCoutput$I
# LDSCoutput$N contains the sample sizes (N) for the heritabilities and sqrt(N1N2) for the co-heritabilities. These are the sample sizes provided in the munging process.
LDSCoutput$N
# LDSCoutput$m is the number of SNPs used to construct the LD score.
LDSCoutput$m

#----------------------------------------------------------
Ssmooth<-as.matrix((nearPD(LDSCoutput$S, corr = FALSE))$mat)
# nScree to get number of factors for EFA, according to different criteria
# - kaiser: lamda > 1
# - parallel analysis: 
# - acceleration factor: elbow of the scree plot
#-  optimal coordinates: corresponds to an extrapolation of the preceeding eigenvalue by a regression line between the eigenvalue coordinates and the last eigenvalue coordinate
n_scree_crit<-nFactors::nScree(Ssmooth)
summary(n_scree_crit)
n<-n_scree_crit$Components %>% max() # get max n for EFA
#----------------------------------------------------------

#----------------------------------------------------------
## Run the SEM model without SNPs
#----------------------------------------------------------

#----------------------------------------------------------
# Common factor model
#----------------------------------------------------------
#To run using DWLS estimation#
CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")
CommonFactor_DWLS$modelfit
#----------------------------------------------------------

#----------------------------------------------------------
# Test common Factor Model
#----------------------------------------------------------
# generate formula for common factor model (all regions loading)
OneFactorModel <- paste0("F1 =~ NA*", paste(colnames(LDSCoutput$S),collapse=" + "),"\n",
                         "F1 ~~ 1*F1")

CFA_1Factors <-usermodel(LDSCoutput, estimation = "DWLS", model = OneFactorModel, 
                         CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFA_1Factors$modelfit
#----------------------------------------------------------

#----------------------------------------------------------
run_efa<-function(Ssmooth,i){
  tmp<-factanal(covmat = Ssmooth, factors = i, rotation = "promax")
  return(tmp$loadings) #print the loadings
}
EFA_loadings<-lapply(1:n, function(i){run_efa(Ssmooth,i)})
#----------------------------------------------------------

#----------------------------------------------------------
cfa<-lapply(seq(from=0.4,to=0.6,by=0.1),function(t) 
  run_cfa_model(nFactors = n,thr=t,ldscoutput=LDSCoutput,EFA="EFA_loadings")
) %>% do.call("rbind",.)

# for each model, check how many regions are included
cfa<-cfa %>% mutate(measures_in_model=
                              model %>% gsub(paste0("F",1:10,collapse="|"),"",.) %>% 
                              gsub(" =\\~ NA\\*","",.) %>% gsub("\\\n"," + ",.) %>% 
                              gsub(" \\+ ",",",.) %>% gsub(",$","",.)
) %>%
  mutate(measures_n_model=measures_in_model %>% strsplit(.,",") %>% sapply(.,function(x) unique(x) %>% length()),
         measures_in_model=measures_in_model %>% strsplit(.,",") %>% sapply(.,function(x) unique(x) %>% paste0(.,collapse=",")))

#----------------------------------------------------------
# **hemen nago**
# to do: edit run_cfa to include also variables that in whichever factor they have the largest loading for...
p_a<-cfa %>%
  ggplot(data=.,aes(x=CFI,y=SRMR,color=as.factor(measures_n_model))) + 
  geom_vline(xintercept = 0.90) + 
  geom_hline(yintercept = 0.05) +
  geom_point() + geom_text_repel(aes(label=FactorLoadingThreshold)) +
  NULL
rm(p_a)
#----------------------------------------------------------
# to re-run using defined threshold
# thr=cfa %>% arrange(-measures_n_model,SRMR,CFI) %>% head(n=1) %>% pull(FactorLoadingThreshold)
thr=0.5
# 
m<-"cfa"
subset(get(m),nFactors==n & FactorLoadingThreshold==thr) %>% print()

loadings<-get(paste0("EFA_loadings"))[[n]]
f_loadings<-sapply(1:NCOL(loadings),function(i){
  w<-which(abs(loadings[,paste0("Factor",i)])>thr) %>% names()
  return(w)
})
names(f_loadings)<-paste0("F",1:NCOL(loadings))

# check which manifest variables are not assigned to any factor
traits2assign<-trait.names2[!trait.names2 %in% (unlist(f_loadings) %>% unique())]
traits2assignF<-sapply(traits2assign,function(x) {
  l<-loadings[x,]
  maxF<-which(l==max(l)) %>% names() %>% gsub("Factor","F",.)
  return(maxF)
})

# assign traits to corresponding F
for (f in names(f_loadings)){
  w<-which(traits2assignF==f) %>% names()
  if(length(w)>0){
    f_loadings[[f]]<-c(f_loadings[[f]],w)
  }
  rm(w)
}
rm(f)
rm(traits2assignF,traits2assign)

# remove factor if there are no measures with loadings that are higher than the threshold
w0=which(lapply(f_loadings,length)==0)
if(length(w0)>0){
  f_loadings[[w0]]<-NULL
}
rm(w0)
# generate formula (generic!)
nFactorModel <- lapply(1:length(f_loadings),function(f){
  paste0("F",f, " =~ NA*", paste(f_loadings[[f]],collapse=" + "),"\n")
}) %>% do.call("paste0",.)

# model_part1<-paste0("F",1:n) %>% combn(.,m=2) %>% apply(.,2,function(t) paste0(t,collapse=" ~~ 0*")) %>% paste(.,collapse="\n")
# nFactorModel<-paste0(nFactorModel,"\n",model_part1)

#run the model
CFA_nFactors <-usermodel(LDSCoutput, estimation = "DWLS", model = nFactorModel, 
                         CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
# identify Heywood cases
vars2constrain1<-CFA_nFactors$results %>% filter(STD_All<(-1)|(STD_All>1)) %>% pull(rhs) %>% unique()
vars2constrain2<-CFA_nFactors$results %>% filter(!grepl("^F",lhs)) %>% 
  # group_by(rhs) %>% summarise(total_std_all=sum(STD_All)) %>% filter(abs(total_std_all)>1.3) %>% 
  filter(STD_All<(-0.05)) %>%
  pull(rhs) %>% unique()
vars2constrain<-c(vars2constrain1,vars2constrain2) %>% unique()
rm(vars2constrain1,vars2constrain2)

# re-run if necessary
if(length(vars2constrain)>0){
  model_part2<-paste(vars2constrain,"~~",letters[1:length(vars2constrain)],"*",vars2constrain,"\n") %>% paste(.,collapse="")
  model_part3<-paste(letters[1:length(vars2constrain)],"> .001\n") %>% paste(.,collapse="")
  
  
  nFactorModel2<- paste0(nFactorModel,"\n",
                         model_part2,model_part3)
  
  #run the model
  CFA_nFactors <-usermodel(LDSCoutput, estimation = "DWLS", model = nFactorModel2, 
                           CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
  
  rm(nFactorModel2)
  rm(model_part2,model_part3)
}

assign(paste0("CFA_",n,"Factors","_thr",thr),CFA_nFactors)
# clean intermediate
rm(m,loadings,f_loadings)
rm(vars2constrain)
rm(nFactorModel,CFA_nFactors)
#
rm(thr)

#----------------------------------------------------------
# save
setwd(out_dir)

cols <- wesanderson::wes_palette(name = "Moonrise2", n = 4, type = "discrete")
colorlist <- list(man = cols[1], lat = cols[2])

for(cfa_name in ls(pattern="CFA_")){
  # check output
  tmp<-get(cfa_name)
  n<-unique(tmp$results$lhs) %>% grep("F",.) %>% length()
  for (f in 1:n) {
    tmp$results %>% filter(lhs==paste0("F",f) & !grepl("^F",rhs)) %>% 
      # filter(p_value<0.05) %>% 
      print()
  }
  rm(n)
  # save file
  cfa_file<-paste0(cfa_name,"_",analysis_run,".xlsx")
  openxlsx::write.xlsx(get(cfa_name),file=cfa_file)
  rm(cfa_file)
  
  # figure, need to tweak
  png(paste0(cfa_name,"_",analysis_run,".png"),width=1500,height=800)
  semPaths(semPlotModel_GSEM(tmp),
           layout="tree2",whatLabels = "est",
           sizeMan=5,
           nCharNodes=10,
           label.cex=1.5,
           edge.label.cex=1,
           color = colorlist,
           style="OpenMx",
           nDigits=3,
           reorder=TRUE)
  dev.off()
  #
  rm(tmp)
}
rm(cfa_name)

save(list=ls(pattern="CFA_|FactorModel|scree|EFA|cfa"),
     file=paste0("GSEM_cerebellum_",analysis_run,"_CFA.Rdata"))

#----------------------------------------------------------
rm(list=ls(pattern="CFA_|FactorModel|scree|EFA|cfa|Factor|LDSC|Ssmooth"))
rm(trait.names,trait.names2,f,ld,population.prev,sample.prev)
# rm(filename_run)
gc()
