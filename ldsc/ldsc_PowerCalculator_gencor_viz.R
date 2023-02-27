#----------------------------------------------------------------------
# GCTA power calculator - for genetic correlation
## between a behavioural trait/ brain measure
## given h2 estimates and N
#----------------------------------------------------------------------
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
#----------------------------------------------------------------------
rm(list=ls())
#----------------------------------------------------------------------
options(stringsAsFactors = FALSE)
# define pattern for this run
root="UKB_BIG40"

# define working_dir
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
working_dir=paste(dir,"resources/datasets/GWAS_sumstats/ldsc/",root,"/cerebellum/",sep="")

#-----------------------------------------------
results<-read.csv(file=paste0(working_dir,"power_calc_results_",root,".csv",sep=""))
#
results$side<-results$hemisphere
results$side[which(results$hemisphere=="NotLat")]<-"Vermis"
results$side<-factor(results$side,levels=c("L","R","Vermis"))
table(results$side)

# clean region name
results$region_name <- results$region %>% gsub("smri_thick_cort.|smri_area_cort.|destrieux_","",.) %>% 
  gsub("\\.","-",.) %>% gsub("-and-","+",.) %>% gsub("AREA_|area_|THICK_|thick_","",.) %>% toupper()
# 
results$region_2<-results$region_name %>% gsub("V_","",.) %>% gsub("CEREBELLUM_","",.)

#-----------------------------------------------
# Make plots
#-----------------------------------------------
results<-results %>% filter(region_2=="X"|
                              region_2=="CRUS_I"|
                              region_2=="CRUS_II")
n<-results %>% pull(region_2) %>% unique() %>% length()

colors<-brewer.pal(n+1, "Set1")[-6]

plot_power<-function(x,data){ 
  d<-subset(data,trait1==x)
  x2<-d$Reported.trait
  
  p0<-ggplot(data=d,aes(x=rg,y=pwr,color=region_2,alpha=0.8)) + 
    geom_hline(yintercept =0.8,color="black",linetype="dashed") + 
    geom_line(size=1) +
    scale_color_manual(values=colors) + 
    ylim(c(0,1)) +
    xlim(c(0,0.2)) +
    facet_wrap(.~ side,scales = "free_x",ncol=3) +
    labs(x = bquote(''*rho*'(G)'), y = "Power", title = x2) + 
    theme(legend.position = "bottom") +
    NULL
  
  return(p0)
}

#-----------------------------------------------
plot_list<-lapply(c("LeftHandedness","EA3","CP","IQ"), function(x){
  plot_power(x = x,data=results)+  theme(legend.position = "none") 
})
legend<-get_legend(plot_list[[1]] +  theme(legend.position = "bottom") )
power_combined_plot<-plot_grid(plot_grid(plotlist=plot_list,ncol=2),
          legend,
          rel_heights = c(1,0.1),
          ncol=1
)
power_combined_plot
ggsave(power_combined_plot,file=paste0(working_dir,"power_rg_behGWASes_reading_ROIs.png"),width=12,height=12)

for(x in results$trait1 %>% unique()){
  p<-plot_power(x = x,data=results)
  ggsave(p,file=paste0(working_dir,"power_rg_",x,"_reading_ROIs.png"))
}
rm(x,p)
#-----------------------------------------------
