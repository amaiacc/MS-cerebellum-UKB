# attempt to visualize rg cerebellum results 
set.seed(10)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggseg)
## using circus plots
library(circlize)
library(units)
# library(wesanderson)
library(png)
#---------------------------------------------------
library(ggsegHO) # to visualize regions
# fix ggsegHO labelling issue
## https://github.com/ggseg/ggesgHO/issues/1
swap_labels<-function(d,l1,l2){
  # get indices and geometry
  w1<-which(d$label==l1)
  r1<-d$region[w1] %>% unique()
  l1<-d$label[w1] %>% unique()
  h1<-d$hemi[w1] %>% unique()
  
  w2<-which(d$label==l2)
  r2<-d$region[w2] %>% unique()
  l2<-d$label[w2] %>% unique()
  h2<-d$hemi[w2] %>% unique()
  
  
  # swap
  d$region[w1]<-r2
  d$label[w1]<-l2
  d$hemi[w1]<-h2
  
  d$region[w2]<-r1
  d$label[w2]<-l1
  d$hemi[w2]<-h1
  
  return(d)
  
}
# check labels that have been switched
labels2swap<-rbind(
  cbind("lh_Juxtapositional.Lobule.Cortex..formerly.Supplementary.Motor.Cortex.","lh_Subcallosal.Cortex"),
  cbind("lh_Lingual.Gyrus","lh_Parahippocampal.Gyrus..posterior.division"),
  cbind("lh_Angular.Gyrus","lh_Supramarginal.Gyrus..posterior.division"),
  cbind("rh_Inferior.Temporal.Gyrus..anterior.division","rh_Inferior.Temporal.Gyrus..posterior.division"),
  cbind("rh_Frontal.Operculum.Cortex","rh_Occipital.Fusiform.Gyrus" )
) %>% data.frame()

for (i in 1:NROW(labels2swap)){
  hoCort$data<-swap_labels(d=hoCort$data,l1=labels2swap$X1[i],l2=labels2swap$X2[i])
}
rm(i)
# to check
# tmp2<-hoCort$data %>% filter(region %in% hoCort$data$region[w])
# tmp2 %>%
#   mutate(region2plot=region) %>% 
#   brain_join(hoCort) %>%
#   ggplot() + geom_sf(aes(fill=region2plot)) +
#   theme_void() +
#   theme(legend.position="bottom",legend.text=element_text(size=16)) +
#   NULL
# rm(tmp2)
# rm(w)
#---------------------------------------------------

#---------------------------------------------------
options(stringsAsFactors = FALSE)
root="UKB_BIG40"
project="cerebellum_UKB"
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
primary_dir=paste(dir,"resources/datasets/GWAS_sumstats/downloaded_data/UKB/","BIG40","/",sep="")
ldsc_dir=paste(dir,"resources/datasets/GWAS_sumstats/ldsc/",root,"/output/",sep="")
out_dir=paste(dir,project,"data/ldsc/",sep="/")
if(!dir.exists(out_dir)){dir.create(out_dir)}
# set working dir
setwd(out_dir)
#----------------------------------------------------------------------
# define order of regions:
cer_regions<-c("TOTAL","CORTEX","WHITE-MATTER",
               "I-IV","V","VI","CRUS I","CRUS II",
               "VIIB","VIIIA","VIIIB","IX","X")
#----------------------------------------------------------------------
# Read
rg_cort<-read.csv(file=paste(out_dir,"/ldsc_rg_cortical","_",root,"_",project,".csv",sep=""))
# edit some columns

col_fun1=colorRamp2(c(-0.5,0,0.5), c("blue", "white","red"))
col_funPos=colorRamp2(c(0, 0.5), c("white", "blue"))

rg_cort <- rg_cort %>%
  # mutate(from=paste0(region," (",hemisphere,")") %>% gsub("\n","",.),
  #        to=paste0(region.cortical," (",hemisphere.p2,")") %>% gsub("\n","",.)
  # ) %>%
  mutate(to=paste0(region," (",hemisphere,")") %>% gsub("\n","",.),
         from=paste0(region.cortical," (",hemisphere.p2,")") %>% gsub("\n","",.)
  ) %>%
  mutate(region=factor(region,levels=cer_regions),
         color_direction=if_else(sig=="","grey",
                                 col_fun1(rg))) %>%
  arrange(region)
# simplify dataset for testing purposes
d <- rg_cort %>% filter(sig=="*")
d3<-d %>% dplyr::select(from,to,rg) %>% mutate(from=gsub("\\.","\n",from),
                                               to=gsub("\\.","\n",to)
) %>%
  mutate(tmp1=strsplit(to," \\(") %>% sapply("[[",1) %>% 
           factor(.,levels=cer_regions[cer_regions %in% unique(d$region)]),
         tmp2=strsplit(from," \\(") %>% sapply("[[",1)  %>% 
           factor(.,levels=c(
             "Lingual Gyrus",
             "Occ Fusif Gyrus",
             "Temp Occ Fusif Cortex",
             "Temp Fusif Cortex Post",
             "Mid Temp Gyrus Tempocc",
             "Inf Temp Gyrus Tempocc",
             "Mid Temp Gyrus Post",
             "Inf Temp Gyrus Post",
             "Mid Front Gyrus",
             "Paracing Gyrus",
             "Cing Gyrus Ant",
             "Parahipp Gyrus Post",
             "Insular Cortex"  ))
  ) %>%
  arrange(tmp2,tmp1)
  # arrange(tmp1,tmp2)


#----------------------------------------------------------------------
## phenotype information
idps<-read.csv(paste(primary_dir,"IDPs_summary.csv",sep=""),header=TRUE) %>% 
  mutate(region=gsub("V_CEREBELLUM","CEREBELLUM",region)) %>%
  mutate(hemisphere2=if_else(is.na(hemisphere),"BiLat",hemisphere) %>% factor(.,levels=c("BiLat","L","R","V"))) %>%
  dplyr::select(-Units,-Type,-contains("Cat"),
                -contains("N."),-contains("Npar."),-contains("Nnonpar."),
                -downloaded,-Global,-Reading_selection,-ABCD_BrainBehaviour,
                -Cerebellum,Cerebellum2,
                -Subcortical,
                -parcellation,-measure)

# subset to only keep idps present in rg_cort
idps<-idps %>% filter(pheno %in% c(rg_cort$p1,rg_cort$p2))

# idps$gColor = wes_palette("Moonrise2",4) %>% as.character %>% 
#   .[as.numeric(idps$hemisphere2)]




# associate color to region
target_labels<- ggseg(atlas="hoCort")$data %>% select(region,hemi,label) %>% distinct() %>%
  mutate(matchLabel=tolower(label) %>%
           gsub("\\.|-","_",.) %>%
           gsub("__","_",.) %>%
           gsub("_division|_part","",.)
  )
idps <- idps %>%
  mutate(hemi=gsub("L","lh",hemisphere) %>% gsub("R","rh",.)) %>%
  # filter(!grepl("CEREBELLUM",region)) %>%
  mutate(matchLabel=paste(hemi,region,sep="_") %>%
           gsub("area_|thick_|smri_|cort.","",.) %>%
           gsub(paste0("FAST","_"),"",.) %>%
           gsub("\\.|-","_",.) %>% tolower() %>%
           gsub("occ_","occipital_",. ) %>%
           gsub("temp_","temporal_",.) %>%
           gsub("front_","frontal_",.) %>%
           gsub("inf_","inferior_",.) %>%
           gsub("_inf$","_inferior",.) %>%
           gsub("sup_","superior_",.) %>%
           gsub("_sup$","_superior",.) %>%
           gsub("mid_","middle_",.) %>%
           gsub("_med_","_medial_",.) %>%
           gsub("_ant","_anterior",.) %>%
           gsub("_post","_posterior",.) %>%
           gsub("cent_","central_",.) %>%
           gsub("pars","pars_",.) %>%
           gsub("tri$","triangularis",.) %>%
           gsub("op$","opercularis",.) %>%
           gsub("orb$","orbitalis",.) %>%
           # gsub("operc","opercular",.) %>%
           gsub("parietal_operc","parietal_operculum",.) %>%
           gsub("frontal_operc","frontal_operculum",.) %>%
           gsub("central_operc","central_opercular",.) %>%
           gsub("heschl_gyrus","heschl_s_gyrus_includes_h1_and_h2_",.) %>%
           gsub("fusif_","fusiform_",.) %>%
           gsub("cing","cingulate",.) %>%
           gsub("orb","orbital",.) %>%
           gsub("intracalc","intracalcarine",.) %>%
           gsub("supramarg_","supramarginal_",.) %>%
           gsub("supracalc","supracalcarine",.) %>%
           gsub("precun_","precuneous_",.) %>%
           gsub("tempocc","temporooccipital",.) %>%
           gsub("latoccipital","lateral_occipital",.) %>%
           gsub("juxtapos_lobule_cortex","juxtapositional_lobule_cortex_formerly_supplementary_motor_cortex_",.) %>%
           gsub("parahipp","parahippocampal",.) %>%
           gsub("posteriorcentral","postcentral",.)
  )

idps2<-merge(idps,rg_cort[,c("p1","region")] %>% rename(pheno=p1) %>% distinct(),by=c("pheno"),all.y=TRUE)
idps3<-merge(idps,rg_cort[,c("p2","region.cortical")] %>% rename(pheno=p2,region=region.cortical) %>% distinct(),by=c("pheno"))

idps1<-rbind(idps2,idps3) %>% rename(region=region.y)
rm(idps2,idps3)

idps1<- idps1 %>% mutate(val=paste0(region," (",hemisphere,")") %>% gsub("\n","",.) %>%
                           gsub("\\(NA\\)","(V)",.) %>% gsub("TOTAL \\(V\\)","TOTAL (BiLat)",.)
                         )


# subset to only keep idps present in rg_cort
idps_subset<-idps1 %>% filter(pheno %in% c(d$p1,d$p2))

idps_subset$color<-factor(idps_subset$region)
levels(idps_subset$color)<-rand_color(n=length((idps_subset$color)))
#----------------------------------------------------------------------
# create a figure for each cortical region (i.e. from)
if(!dir.exists("figures_cortical_regions")){dir.create("figures_cortical_regions")}
idps_cort<-subset(idps_subset,lobe!="SUBCORTICAL")[,c("val","color","matchLabel")] %>% merge(.,target_labels,by="matchLabel")

# for (l in 1:NROW(idps_cort)){
#   f=paste0("figures_cortical_regions/",idps_cort$val[l],".png")%>% gsub(" ","_",.)
#   tmp<-idps_cort[l,] %>%  brain_join(hoCort)
#   p<-  ggplot(tmp) +
#     geom_brain(aes(fill=color),atlas=hoCort,
#                # side="lateral",
#                hemi=idps_cort$hemi[l]) +
#     scale_fill_manual(values=as.character(unique(idps_cort$color[l])),na.value="transparent",
#                       labels=idps_cort$region[l]) +
#     theme_void() +
#     theme(legend.position="bottom",legend.text=element_text(size=12)) +
#     guides(fill=guide_legend(ncol=3,title="",override.aes = list(color="white"))) +
#     NULL
#   
#   png(file=f)
#   print(p)
#   dev.off()
#   rm(f,tmp,p)
#   }


img_list2<-lapply(1:NROW(idps_cort),function(l){
  f <- paste0("figures_cortical_regions/edited/",idps_cort$val[l],"_ed.png")%>% gsub(" ","_",.)
  if(file.exists(f)){
    o <- as.raster(readPNG(f))
  } else {
    o <-NULL
  }
  return(o)
})

names(img_list2)<-idps_cort$val
subset(idps_subset,lobe=="SUBCORTICAL") # somehow add these to the list as well, so that the raster works?

img_list3<-append(img_list2,
                  vector(mode='list', length=NROW(subset(idps_subset,lobe=="SUBCORTICAL"))
                         )
)
names(img_list3)<-c(names(img_list2),subset(idps_subset,lobe=="SUBCORTICAL")$val)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# define grid colors
d2 <- merge(data.frame(val=c(d$from,d$to)),idps_subset %>% dplyr::select(val,region,color,lobe),all.x=TRUE) 

d_colors<- d2$color[match(c(d$from,d$to),d2$val)] %>% as.character()
names(d_colors)<-d2$val[match(c(d$from,d$to),d2$val)] %>% as.character()
# check
d_colors[c(d$from,d$to)]
#----------------------------------------------------------------------
# create a chord diagram that contains the rg values
pdf("rg_cerebellar_cortical_circleplot.pdf",width=6,height=12)
chordDiagram(d3,
             # numeric.column=rg,
             col=d$color_direction,
             directional=TRUE,
             grid.col=d_colors, #,
             annotationTrack = c("grid"),
             preAllocateTracks = 2,
             # preAllocateTracks = list(track.height = 0.3),
             # order=idps$hemisphere2,grid.col=setNames(idps$gColor,idps$hemisphere2)
             small.gap=3,
             big.gap = 10
             )

circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", 
              niceFacing = TRUE, adj = c(0, 0.5),cex=0.7)
  # circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, track.height = 2, bg.border = NA)

dev.off()
# add brain figures
circos.track(track.index=1, ylim = c(0, 1), panel.fun = function(x, y) {
  img = img_list3[[ CELL_META$sector.index ]]
  if(!is.null(img)){
    circos.raster(img, CELL_META$xcenter, CELL_META$ycenter, 
                  width = "1cm",
                  # width = CELL_META$xrange, height = CELL_META$yrange,
                  niceFacing = TRUE,
                  facing = "downward")
  }
}, track.height = 1, bg.border = NA)

# u=0
# for(si in get.all.sector.index()){
#   xplot=get.cell.meta.data("xplot",si)
#   u=u+1
# 
#   # Small workaround because coordinate 0 should be 360
#   if(xplot[1] == 0) xplot[1] = 360
# 
#   # x=.86*cos(as_radians((xplot[2]+xplot[1])/2))
#   # y=.86*sin(as_radians((xplot[2]+xplot[1])/2))
#   deg=(xplot[2]+xplot[1])/2
#   rad=(deg*pi)/180
# 
#   x=.86*cos(rad)
#   y=.86*sin(rad)
#   # add plot here! **to edit**
#   f=paste0("figures_cortical_regions/",si,"_ed.png") %>% gsub(" ","_",.)
#   if(file.exists(f)){
#     png::readPNG(f) %>%
#       rasterImage(.,x-05, y-0.09, x+0.09, y+0.09)
#   }
# 
# }



# # as a matrix
# d2 <- d %>% dplyr::select(from,to,rg) %>%
#   pivot_wider(names_from=to,values_from=rg)
# row.names(d2)<-d2$from
# d2<-d2 %>% dplyr::select(-from) %>% as.matrix
# 
# 
# 
# d_m<-d %>% dplyr::select(from,to,rg) %>% pivot_wider(values_from=rg,names_from=to)
# tmp<-d_m$from
# d_m<- d_m %>% dplyr::select(-from) %>% as.matrix()
# rownames(d_m)<-tmp
# cereb.region<-strsplit(tmp,"\\.") %>% sapply("[[",1)
# cereb.hemis<-strsplit(tmp,"\\.") %>% sapply("[[",2)
# rm(tmp)
# 
# chordDiagram(d_m, col=cereb.region)
# 

# library(wesanderson)
# 
# dplyr::select(p1,p2,rg,region,region.cortical,hemisphere,hemisphere.p2) %>%
#   mutate(region=factor(region,levels=cer_regions),
#          region.cortical=factor(region.cortical,
#                                 levels=c("Insular \nCortex",
#                                          "Paracing \nGyrus", 
#                                          "Cing \nGyrus Ant",
#                                          "Parahipp \nGyrus Post", 
#                                          "Mid Front \nGyrus",
#                                          # mid temporal
#                                          "Mid Temp \nGyrus Tempocc",
#                                          "Mid Temp \nGyrus Post",
#                                          # inferior temporal
#                                          "Inf Temp \nGyrus Post",
#                                          "Inf Temp \nGyrus Tempocc",
#                                          # fusiform
#                                          "Temp Occ Fusif \nCortex", 
#                                          "Temp Fusif \nCortex Post", 
#                                          "Occ Fusif \nGyrus",
#                                          "Lingual \nGyrus"
#                                 ))
#   ) 