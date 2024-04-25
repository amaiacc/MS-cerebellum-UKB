rm(list=ls())

#----------------------------------------------------------
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
working_dir=paste0(dir,"resources/datasets/GWAS_sumstats/")
out_dir=paste0(dir,"cerebellum_UKB/data/GSEM/")
scripts_dir=paste0(dir,"cerebellum_UKB/scripts/GSEM/")
source(paste0(scripts_dir,"semPlotModel_GSEM.R"))

# organize output from GSEM analysis
cerebellum_levels<-c("I-IV","V","VI","CrusI","CrusII","VIIb","VIIIa","VIIIb","IX","X")
#-----------------------
# read results from EFA
for (f in list.files(out_dir,pattern="Rdata")){
  efa=strsplit(f,"cerebellum_") %>% sapply("[[",2) %>% gsub("_CFA.Rdata","_EFA",.)
  load(f);
  assign(efa,EFA_loadings)
  rm(efa,n_scree_crit,EFA_loadings,cfa,OneFactorModel)
}
rm(f)
#-----------------------
# model fit stats
cfa_stats<-lapply(list.files(out_dir,pattern="xlsx"),function(f){
  t<-openxlsx::read.xlsx(f,sheet=1) %>% mutate(model=gsub(".xlsx","",f))
}) %>% do.call("rbind",.) %>%
  mutate( model_nFactors=strsplit(model,"_") %>% sapply("[[",2) %>% gsub("Factors","",.),
          model_Fthreshold=if_else(grepl("thr",model),
                                   strsplit(model,"_") %>% sapply("[[",3) %>% gsub("thr","",.),
                                   "-"  )) %>%
  mutate(model_measures=if_else(grepl("thr",model),
                                strsplit(model,"Factors_",1) %>% sapply("[[",2) %>% gsub("thr","",.),
                                strsplit(model,"_") %>% sapply("[[",3)
  ) ) %>%
  mutate(model_measures=gsub(unique(model_Fthreshold) %>% paste(.,collapse="|"),"",model_measures) %>% gsub("_","",.)) %>%
  dplyr::select(model,model_measures,model_nFactors,model_Fthreshold,
                CFI,SRMR,AIC,chisq,df) %>%
  arrange(model_measures,model_nFactors)

#-----------------------
# model loadings and structure
cfa_results<-lapply(list.files(out_dir,pattern="xlsx"),function(f){
  t<-openxlsx::read.xlsx(f,sheet=2) %>% 
    mutate(model=gsub(".xlsx","",f))
}) %>% do.call("rbind",.) %>%
  mutate(across("p_value", ~format(.x,scientific=TRUE,digits=2))) %>%
  mutate(across(c(contains("Est"),contains("SE")),  ~as.numeric(.x) )) %>%
  mutate_if(is.numeric,round,2)

# rename and organize output
cfa_results<-cfa_results %>%
  mutate(#model_type=strsplit(model,"_") %>% sapply("[[",1),
    model_nFactors=strsplit(model,"_") %>% sapply("[[",2) %>% gsub("Factors","",.)
  ) %>%
  mutate(lhs=gsub("left|right","lat.",lhs)  %>% gsub("vermis","vermis.",.),
         rhs=gsub("left|right","lat.",rhs)  %>% gsub("vermis","vermis.",.)
         ) %>%
  rename(Variable1=lhs,Variable2=rhs,operator=op) %>%
  mutate(Variable1=gsub("I.IV","I-IV",Variable1),
         Variable2=gsub("I.IV","I-IV",Variable2)) %>%
  mutate(Variable1=factor(Variable1,levels=c(paste0("F",1:4),paste(rep(c("lat","vermis"),each=length(cerebellum_levels)),cerebellum_levels,sep="."))),
         Variable2=factor(Variable2,levels=c(paste0("F",1:4),paste(rep(c("lat","vermis"),each=length(cerebellum_levels)),cerebellum_levels,sep=".")))
         ) %>%
  mutate(model_Fthreshold=if_else(grepl("thr",model),
                                   strsplit(model,"_") %>% sapply("[[",3) %>% gsub("thr","",.),
                                   "-"  )) %>%
  mutate(model_measures=if_else(grepl("thr",model),
                                strsplit(model,"Factors_",1) %>% sapply("[[",2) %>% gsub("thr","",.),
                                strsplit(model,"_") %>% sapply("[[",3) 
  ) ) %>%
  mutate(model_measures=gsub(unique(model_Fthreshold) %>% paste(.,collapse="|"),"",model_measures) %>% gsub("_","",.))


# create sub-tables for the CommonFactor model and CFA's derived from EFA (nEFA)
cols_values<-colnames(cfa_results)[grep("STD|Unstand|p_value",colnames(cfa_results))]

cfa_results_OneF_wide<-subset(cfa_results,model_nFactors==1) %>% 
  dplyr::select(-model_Fthreshold,-model_nFactors,-model) %>%
  pivot_wider(id_cols=c(Variable1,operator,Variable2), 
              names_from = model_measures,
              values_from = cols_values,
              names_glue="{model_measures}.{.value}",
              ) %>%
  dplyr::select(Variable1,operator,Variable2,contains("VL"),contains("VR"))


cfa_results_nF_wide<-subset(cfa_results,model_nFactors!=1) %>%
  dplyr::select(-model) %>%
  # rename factors to make them match
  # VL F2 or VR F3 --> F2/F3
  mutate(Variable1=if_else(model_measures=="VL"&Variable1=="F2","F2/F3",
                           if_else(model_measures=="VR"&Variable1=="F3","F2/F3",Variable1)),
         Variable2=if_else(model_measures=="VL"&Variable2=="F2","F2/F3",
                            if_else(model_measures=="VR"&Variable2=="F3","F2/F3",Variable2))) %>%
  # rename factors to make them match
  # VL F3 or VR F2 --> F3/F2
  mutate(Variable1=if_else(model_measures=="VL"&Variable1=="F3","F3/F2",
                           if_else(model_measures=="VR"&Variable1=="F2","F3/F2",Variable1)),
         Variable2=if_else(model_measures=="VL"&Variable2=="F3","F3/F2",
                            if_else(model_measures=="VR"&Variable2=="F2","F3/F2",Variable2))) %>%
  # format wider
  pivot_wider(id_cols=c(Variable1,operator,Variable2,model_Fthreshold,model_nFactors), 
              names_from = model_measures,
              values_from = cols_values,
              names_glue="{model_measures}.{.value}",
  ) %>%
  dplyr::select(Variable1,operator,Variable2,contains("VL"),contains("VR"))

# combine EFA models
efa_df<-lapply(ls(pattern="_EFA"),function(f){
  efa<-get(f)
  n<-length(efa)
  efa_n<-efa[[n]]
  efa_df<-matrix(as.numeric(efa_n),attributes(efa_n)$dim,dimnames=attributes(efa_n)$dimnames) %>%
    data.frame() %>% mutate(model=f %>% gsub("_EFA","",.), Variable=rownames(.))
  return(efa_df)
}) %>% do.call("rbind",.) %>%
  mutate(Variable=gsub("left|right","lat.",Variable) %>% gsub("vermis","vermis.",.)) %>%
  pivot_wider(id_cols =Variable,names_from=model,values_from=contains("Factor"),
              names_glue="{model}.{.value}") %>%
  mutate_if(is.numeric,round,2) %>%
  mutate(Variable=gsub("I.IV","I-IV",Variable)) %>%
  mutate(measure_type=strsplit(Variable,"\\.") %>% sapply("[[",1),
         measure_cerebellum=strsplit(Variable,"\\.") %>% sapply("[[",2) %>%
           factor(.,levels=cerebellum_levels)) %>%
  dplyr::select(Variable,contains("measure"),contains("VL"),contains("VR")) %>%
  arrange(measure_type,measure_cerebellum)
  
#
all_tables<-list("Model fit indices"=cfa_stats,
                 "Common Factor Models"=cfa_results_OneF_wide,
                 "EFA loadings, n=4"=efa_df,
                 "CFA n=4 Models"=cfa_results_nF_wide)
openxlsx::write.xlsx(all_tables,file=paste0(out_dir,"GSEM_cerebellum_summary_tables.xlsx"))



# # test plot
load( paste0(out_dir,"GSEM_cerebellum_VL_CFA.Rdata" ))

png(paste0("CFA_4factors_thr0.5_VL_d.png"),width=3500,height=1500)
semPaths(semPlotModel_GSEM(CFA_4Factors_thr0.5),
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

# this fails because it's not in the right format...
# sem::pathDiagram(CFA_4Factors_thr0.5$results,
            # file="pathDiagram",
            # output.type="dot")

basic_dot<-function(output,est.label="STD_All",filename="test"){
  all_vars<-c(output$results$rhs,output$results$lhs) %>% unique()
  manifest_vars<-all_vars[grep("^F",all_vars,invert=TRUE)]
  latent_vars<-all_vars[grep("^F",all_vars)]
  
  output$results$label<-output$results[,est.label] %>% round(.,digits=2) %>% as.character()
  
  t<-output$results %>% # filter(p_value<0.05/NROW(output$results)) %>%
    mutate(path=paste(lhs,op,rhs) %>% gsub("=~","->",.) %>% gsub("~~","<->",.)) %>%
    mutate(tmp=paste(path,label,sep=","))
  m<-paste0("\n",paste(t$tmp,collapse="\n"))
  sem.model<-specifyModel(text=m)
  pathDiagram(sem.model,obs.variables = manifest_vars,
              same.rank=latent_vars,
              style="ram",
              ignore.double=FALSE,
              error.nodes=TRUE,
              output="dot",
              file=filename
              )
  
}

CFA_4Factors_thr0.5$results<- CFA_4Factors_thr0.5$results %>% mutate(lhs=gsub("left","L-",lhs) %>% gsub("vermis","V-",.),
       rhs=gsub("left","L-",rhs) %>% gsub("vermis","V-",.)
) %>%
  mutate(tmpR1=rhs %>% strsplit(.,"-") %>% sapply("[[",1),
         tmpR2=rhs %>% gsub("L-|V-","",.) %>% gsub(paste0("F",1:4,collapse="|"),"",.),
         tmpL1=lhs %>% strsplit(.,"-") %>% sapply("[[",1),
         tmpL2=lhs %>% gsub("L-|V-","",.) %>% gsub(paste0("F",1:4,collapse="|"),"",.)
  ) %>%
  mutate(lhs=paste(tmpL2,tmpL1,sep="") %>% gsub("L$"," (L)",.) %>% gsub("V$"," (V)",.),
         rhs=paste(tmpR2,tmpR1,sep="") %>% gsub("L$"," (L)",.) %>% gsub("V$"," (V)",.)
  ) %>%
  dplyr::select(-contains("tmp"))

basic_dot(output = CFA_4Factors_thr0.5,filename="CFA_4Factors_thr0.5")
