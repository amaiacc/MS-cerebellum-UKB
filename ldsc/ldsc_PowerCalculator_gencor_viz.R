#----------------------------------------------------------------------
# GCTA power calculator - for genetic correlation
## between a behavioural trait/ brain measure
## given h2 estimates and N
#----------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
theme_set(theme_cowplot())
#----------------------------------------------------------------------
rm(list=ls())
#----------------------------------------------------------------------
options(stringsAsFactors = FALSE)
# define pattern for this run
root="UKB_BIG40"

# define working_dir
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
working_dir=paste0(dir,"cerebellum_UKB/data/ldsc/")

#-----------------------------------------------
results<-read.csv(file=paste0(working_dir,"power_calc_results_",root,".csv",sep=""))
#
results$side<-results$hemisphere
# results$side[which(results$hemisphere=="NotLat")]<-"Vermis"
# results$side<-factor(results$side,levels=c("L","R","Vermis"))
table(results$side,useNA="ifany")

# clean region name
results$region_name <- results$region %>% gsub("smri_thick_cort.|smri_area_cort.|destrieux_","",.) %>% 
  gsub("\\.","-",.) %>% gsub("-and-","+",.) %>% gsub("AREA_|area_|THICK_|thick_","",.) %>% toupper()
# 
results$region_2<-results$region_name %>% gsub("V_","",.) %>% gsub("CEREBELLUM_","",.)
results$alpha<-factor(results$alpha,levels=sort(unique(results$alpha),decreasing = TRUE))
levels(results$alpha)<-c("0.05","0.05/(19*3)")
#-----------------------------------------------
# Make plots
#-----------------------------------------------
results<-subset(results,!is.na(pwr))
# results<-results %>% filter(region_2=="X"|
#                               region_2=="CRUS_I"|
#                               region_2=="CRUS_II")
n<-results %>% pull(region_2) %>% unique() %>% length()

# colors<-brewer.pal(n+1, "Set1")[-6]

plot_power<-function(x,data){ 
  d<-subset(data,ref1==x)
  x2<-d$Reported.trait
  max_rg_pwr80<-d %>% group_by(side,alpha) %>% filter(pwr<=0.8) %>% summarise(max_rg=max(rg))
  
  p0<-ggplot(data=d,aes(x=rg,y=pwr,
                        color=alpha)) + 
    geom_hline(yintercept =0.8,color="black",linetype="dashed") + 
    geom_vline(data=max_rg_pwr80,aes(xintercept=max_rg,color=alpha),linetype="dashed") +
    geom_line(size=1,alpha=0.8) +
    ylim(c(0,1)) +
    scale_x_continuous(breaks=c(0,0.05,0.1,0.15,0.2),label=c(0,0.05,0.1,0.15,0.2),limits=c(0,0.25)) +
    # facet_wrap(region~ side,scales = "free_x",ncol=3) +
    facet_grid(. ~side) +
    labs(x = bquote(''*rho*'(G)'), y = "Power", title = x2) + 
    theme(legend.position = "bottom") +
    # scale_colour_grey() +
    ggsci::scale_color_jco(name=bquote(''*alpha*'')) +
    theme(panel.spacing=unit(0, "lines"),
          strip.background =element_rect(fill="white",color="grey85",size=1),
          legend.position = "none") +
    panel_border(size = 1)+
    NULL
    NULL
  
  return(p0)
}

#-----------------------------------------------
refs2plot<-c("PGCSCZ,2022","Grove2019","Lee2018")
plot_list<-lapply(refs2plot, function(x){
  plot_power(x = x,data=results)+  theme(legend.position = "none") 
})
# legend<-get_legend(plot_list[[1]] +  theme(legend.position = "bottom") )
# power_combined_plot<-plot_grid(plot_grid(plotlist=plot_list,row=1),
#           legend,
#           rel_heights = c(1,0.1),
#           ncol=1
# )
leg<-get_legend(plot_list[[1]]+ theme(legend.position="bottom",legend.justification = "center"))

power_combined_plot<-plot_grid(plot_grid(plotlist=plot_list,ncol=1),
                               leg,
                               ncol=1,rel_heights = c(1,0.1)
                               # align="hv"
                               )
power_combined_plot

ggsave(power_combined_plot,file=paste0(working_dir,"power_rg_GWASes_cerebellum.png"),width=12,height=12)

for(x in refs2plot){
  p<-plot_power(x = x,data=results)
  ggsave(p,file=paste0(working_dir,"power_rg_",x,"_GWASes_cerebellum.png"))
}
rm(x,p)
#-----------------------------------------------
