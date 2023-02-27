# check relationship between the mean volume of each cerebellar measure and the h2
#---------------------------------------------------
options(stringsAsFactors = FALSE)
root="UKB_BIG40"
project="cerebellum_UKB"
# set directories, dependent on  system:
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
working_dir=paste(dir,project,"data/",sep="/")
setwd(working_dir)
#---------------------------------------------------
vols<-read.csv(file = "UKB_cerebellum_phenotypes_descriptives.csv",skip=1) %>%
  mutate(pheno=paste0("0",as.character(UKB_BIG40_Pheno)))
h2<-read.csv(file="ldsc/ldsc_h2_UKB_BIG40_cerebellum_UKB.csv")
# combine
d<-merge(vols,h2) %>%
  mutate(se.2=sd.2/N.2,
         vol.2.CI95lower=mean.2-1.96*se.2,
         vol.2.CI95upper=mean.2+1.96*se.2,
         ) %>%
  filter(parcellation=="FAST")


# correlation
r1<-cor.test(d$h2,d$mean.2) #%>% round(.,digits=3)
p1<-d %>% ggplot(aes(x=mean.2,y=h2)) +
  ggrepel::geom_text_repel(aes(label=Trait_name),alpha=0.5) +
  geom_errorbar(aes(ymin = h2.CI95lower, ymax = h2.CI95upper),alpha=0.6) +
  geom_errorbarh(aes(xmin = vol.2.CI95lower, xmax = vol.2.CI95upper),alpha=0.6) + # this is tiny so cannot see it
  geom_point() +
  geom_label(x=8000,y=0.07,label=
              paste0("r=",round(r1$estimate,digits=3),
                     ", p-value=",round(r1$p.value,digits=3))) +
  labs(x="Mean volume",y=bquote('Estimate ('*h^2*')')) +
  ylim(c(0,0.5)) +
  NULL

# use ggMarginal function to create marginal histogram
p1a<-ggExtra::ggMarginal(p1, type="histogram")
ggsave(p1a,file="h2_volume_plot.png",width=6,height=4)


#
d2<-subset(d,!(hemisphere=="V"&region=="CRUS I"))
r2<-cor.test(d2$h2,d2$mean.2) #%>% round(.,digits=3)
p2<-d2 %>% ggplot(aes(x=mean.2,y=h2)) +
  ggrepel::geom_text_repel(aes(label=Trait_name),alpha=0.5) +
  geom_errorbar(aes(ymin = h2.CI95lower, ymax = h2.CI95upper),alpha=0.6) +
  geom_errorbarh(aes(xmin = vol.2.CI95lower, xmax = vol.2.CI95upper),alpha=0.6) + # this is tiny so cannot see it
  geom_point() +
  geom_label(x=8000,y=0.07,label=
               paste0("r=",round(r2$estimate,digits=3),
                      ", p-value=",round(r2$p.value,digits=3))) +
  labs(x="Mean volume",y=bquote('Estimate ('*h^2*')')) +
  ylim(c(0,0.5)) +
  NULL

# use ggMarginal function to create marginal histogram
p2a<-ggExtra::ggMarginal(p2, type="histogram")
ggsave(p2a,file="h2_volume_plot2.png",width=6,height=4)


