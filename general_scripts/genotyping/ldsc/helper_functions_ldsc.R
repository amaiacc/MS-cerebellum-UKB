options(stringsAsFactors = FALSE)

# funtions to read and format ldsc files
read_ldsc_h2<-function(h2file){
  tmp<-read.table(file=h2file,sep="\t",header=TRUE)
  tmp$pheno<-gsub("_h2.log","",tmp$File)
  tmp<-subset(tmp,h2..se.!="")
  tmp$h2<-sapply(strsplit(tmp$h2..se.," "),"[[",1) %>% as.numeric()
  tmp$h2.SE<-sapply(strsplit(tmp$h2..se.," "),"[[",2) %>% gsub("\\(|\\)","",.) %>% as.numeric()
  tmp$Intercept<-sapply(strsplit(tmp$Intercept..se.," "),"[[",1) %>% gsub("\\(|\\)","",.) %>% as.numeric()
  tmp$Intercept.SE<-sapply(strsplit(tmp$Intercept..se.," "),"[[",2) %>% gsub("\\(|\\)","",.) %>% as.numeric()
  tmp$Ratio<-sapply(strsplit(tmp$Ratio..se.," "),"[[",1)
  tmp$Ratio.SE<-gsub("\\(|\\)","",sapply(strsplit(tmp$Ratio..se.," "),"[[",2))
  # compute 95CI: = (mean + (1.96 x SE)) to (mean - (1.96 x SE))
  tmp$h2.CI95upper<-tmp$h2+(1.96*tmp$h2.SE)
  tmp$h2.CI95lower<-tmp$h2-(1.96*tmp$h2.SE)
  tmp$stat<-"h2"
  tmp$program<-"LDSC"
  # select columns to keep
  tmp<- tmp %>% dplyr::select(pheno,
                       h2,h2.SE,h2.CI95upper,h2.CI95lower,
                       Intercept,Intercept.SE,Ratio,Ratio.SE,
                       Lambda.GC,
                       stat,program)
  return(tmp)
}

read_ldsc_rg<-function(rgfile){
  ## genetic correlation
  tmp<-read.table(rgfile)
  colnames(tmp)<-c("p1",tmp[1,-1])
  tmp<-subset(tmp,p2!="p2")
  tmp$file<-sapply(strsplit(tmp$p1,".log-"),"[[",1) %>% gsub("_LR|_rg$","",.)
  tmp$p1<-sapply(strsplit(tmp$p1,".log-"),"[[",2)
  tmp$p1<-tmp$p1 %>% gsub(".sumstats.gz","",.) %>% gsub("..\\/","",.)
  tmp$p2<-tmp$p2 %>% gsub(".sumstats.gz","",.) %>% gsub("..\\/","",.)
  tmp$p1<-sapply(strsplit(tmp$p1,":|-"),"[[",1)
  #
  tmp$rg<-as.numeric(tmp$rg)
  tmp$rg.SE<-as.numeric(tmp$se)
  tmp$p<-as.numeric(tmp$p)
  #
  # compute 95CI: = (mean + (1.96 x SE)) to (mean - (1.96 x SE))
  tmp$rg.CI95upper<-tmp$rg+(1.96*tmp$rg.SE)
  tmp$rg.CI95lower<-tmp$rg-(1.96*tmp$rg.SE)
  tmp$stat<-"rg"
  tmp$program<-"LDSC"
  # select columns to keep
  tmp<- tmp %>% dplyr::select(file,p1,p2,
                       rg,rg.SE,rg.CI95upper,rg.CI95lower,p,
                       gcov_int,gcov_int_se,
                       stat,program)
  return(tmp)
}

# functions to plot genetic correlations within atlas
corrplot_rg_atlas<-function(atlas=rg_atlas,lr=rg_lr,parc="FAST",m="VOLUME"){
  tmp<-subset(atlas,parcellation==parc & measure==m)
  l<- tmp %>% filter(hemisphere=="L") %>% dplyr::select(region.p1,region.p2,rg,p,hemisphere) # upper triangle: L
  r<- tmp %>% filter(hemisphere=="R") %>% dplyr::select(region.p2,region.p1,rg,p,hemisphere) # lower triangle: R
  notlat<- tmp %>% filter(is.na(hemisphere)) %>% dplyr::select(region.p2,region.p1,rg,p,hemisphere)
  lr<-subset(lr,parcellation==parc& measure==m) %>% mutate(region2=region,hemisphere="LR") %>%
    dplyr::select(region,region2,rg,p2,hemisphere) %>% rename(p=p2) # get pval for different to 1
  colnames(l)<-colnames(r)<-colnames(lr)<-c("p1","p2","rg","p","hemisphere")
  # t<-rbind(l,r,lr) #;rm(l,r,lr)
  # t<-subset(t,!is.na(rg))
  # if(sum(t$rg>1,na.rm=TRUE)>0){
  #   t$rg[t$rg>1]<-1 # keep range of rg between -1 and 1 (otherwise corrplot complains)
  # }
  # create empty matrices
  t_rg<-p_rg<-matrix(nrow = NROW(lr),ncol=NROW(lr))
  # upper triangle: right
  t_r<- r %>% dplyr::select(p1,p2,rg) %>% spread(p1,rg) %>% dplyr::select(-1) %>% as.matrix()
  p_r<- r %>% dplyr::select(p1,p2,p) %>% spread(p1,p) %>% dplyr::select(-1) %>% as.matrix()
  # lower triangle: left
  t_l<- l %>% dplyr::select(p1,p2,rg) %>% spread(p1,rg) %>% dplyr::select(-1) %>% as.matrix()
  p_l<- l %>% dplyr::select(p1,p2,p) %>% spread(p1,p) %>% dplyr::select(-1) %>% as.matrix()
  # diagonal: left-right
  t_lr<-lr %>% dplyr::select(p1,p2,rg) %>% spread(p1,rg) %>% dplyr::select(-1) %>% as.matrix()
  p_lr<-lr %>% dplyr::select(p1,p2,p) %>% spread(p1,p) %>% dplyr::select(-1) %>% as.matrix()
  
  # fill in matrices
  t_rg[upper.tri(t_rg)]<-t_r[upper.tri(t_r)]
  t_rg[lower.tri(t_rg)]<-t_l[lower.tri(t_l)]
  diag(t_rg)<-diag(t_lr)
  # if rg is larger than 1, fix to 1
  t_rg[which(t_rg>1)]<-1
  #
  p_rg[upper.tri(p_rg)]<-p_r[upper.tri(p_r)]
  p_rg[lower.tri(p_rg)]<-p_l[lower.tri(p_l)]
  diag(p_rg)<-diag(p_lr)
  #
  colnames(t_rg)<-rownames(t_rg)<-colnames(p_rg)<-rownames(p_rg)<-colnames(t_l)
  # plot
  corrplot(corr=t_rg, p.mat=p_rg,
           method="circle",
           # order="hclust",
           # style
           # col=brewer.pal(n=8, name="PuOr"),
           col = COL2('PuOr',n = 10),
           insig="blank",
           is.corr=TRUE,
           addCoef.col=TRUE,
           addCoefasPercent=TRUE,
           type="full",
           diag=TRUE,
           tl.col="black",
           tl.srt=45,
           mar=c(0,0,1,0),
           # title=paste(m, " (",parc,")", sep="")),
           sig.level=0.05
  )
  # return(c)
}

# for not lateralized measures
corrplot_rg_notlat_atlas<-function(atlas=rg_atlas,lr=rg_lr,parc="FAST",m="VOLUME"){
  tmp<-subset(atlas,parcellation==parc & measure==m)
  t1<- tmp %>% 
    filter(  (toupper(hemisphere)!="R"&toupper(hemisphere)!="L") & 
             (toupper(hemisphere)!="RH"&toupper(hemisphere)!="LH") &
             (toupper(hemisphere)!="RIGHT"&toupper(hemisphere)!="LEFT")
             ) %>%
    # filter(!is.na(hemisphere)) %>%
    dplyr::select(region.p2,region.p1,rg,p)
  colnames(t1)<-c("p1","p2","rg","p")
  #
  t_rownames<- t1 %>% dplyr::select(p1,p2,rg) %>% spread(p1,rg) %>% dplyr::select(1) %>% pull() %>% as.character()
  t_rg<-t1 %>% dplyr::select(p1,p2,rg) %>% spread(p1,rg) %>% dplyr::select(-1) %>% as.matrix()
  p_rg<-t1 %>% dplyr::select(p1,p2,p) %>% spread(p1,p) %>% dplyr::select(-1) %>% as.matrix()
  rownames(t_rg) <- rownames(p_rg) <- colnames(t_rg) <- colnames(p_rg)<- t_rownames
  # sort by level
  t_levs<-levels(t1$p1)[levels(t1$p1) %in% t_rownames]
  t_rg<-t_rg[t_levs,t_levs]
  p_rg<-p_rg[t_levs,t_levs]
  #
  corrplot(corr=t_rg, p.mat=p_rg,
           method="circle",
           # order="hclust",
           # style
           # col=brewer.pal(n=8, name="PuOr"),
           col = COL2('PuOr',n = 10),
           insig="blank",
           is.corr=TRUE,
           addCoef.col=TRUE,
           addCoefasPercent=TRUE,
           type="upper",
           diag=FALSE,
           tl.col="black",
           tl.srt=45,
           mar=c(0,0,1,0),
          # title=paste(m, " (",parc,")", sep="")),
           sig.level=0.05
  )
  # return(c)
}


## using ggcorrplot instead of corrplot (to be able to combine plots with cowplot)
ggcorrplot_rg_atlas<-function(atlas=rg_atlas,lr=rg_lr,parc="FAST",m="VOLUME"){
  tmp<-subset(atlas,parcellation==parc & measure==m)
  l<- tmp %>% filter(hemisphere=="L"&(hemisphere.p1==hemisphere.p2)) %>% dplyr::select(region.p1,region.p2,rg,p,hemisphere) # upper triangle: L
  r<- tmp %>% filter(hemisphere=="R"&(hemisphere.p1==hemisphere.p2)) %>% dplyr::select(region.p2,region.p1,rg,p,hemisphere) # lower triangle: R
  notlat<- tmp %>% filter(is.na(hemisphere)) %>% dplyr::select(region.p2,region.p1,rg,p,hemisphere)
  lr<-subset(lr,parcellation==parc& measure==m) %>% mutate(region2=region,hemisphere="LR") %>%
    dplyr::select(region,region2,rg,p2,hemisphere) %>% rename(p=p2) # get pval for different to 1
  colnames(l)<-colnames(r)<-colnames(lr)<-c("p1","p2","rg","p","hemisphere")
  # create empty matrices
  t_rg<-p_rg<-matrix(nrow = NROW(lr),ncol=NROW(lr))
  # upper triangle: right
  t_r<- r %>% dplyr::select(p1,p2,rg) %>% spread(p1,rg) %>% dplyr::select(-1) %>% as.matrix()
  p_r<- r %>% dplyr::select(p1,p2,p) %>% spread(p1,p) %>% dplyr::select(-1) %>% as.matrix()
  # lower triangle: left
  t_l<- l %>% dplyr::select(p1,p2,rg) %>% spread(p1,rg) %>% dplyr::select(-1) %>% as.matrix()
  p_l<- l %>% dplyr::select(p1,p2,p) %>% spread(p1,p) %>% dplyr::select(-1) %>% as.matrix()
  # diagonal: left-right
  t_lr<-lr %>% dplyr::select(p1,p2,rg) %>% spread(p1,rg) %>% dplyr::select(-1) %>% as.matrix()
  p_lr<-lr %>% dplyr::select(p1,p2,p) %>% spread(p1,p) %>% dplyr::select(-1) %>% as.matrix()
  
  # fill in matrices
  t_rg[upper.tri(t_rg)]<-t_r[upper.tri(t_r)]
  t_rg[lower.tri(t_rg)]<-t_l[lower.tri(t_l)]
  diag(t_rg)<-diag(t_lr)
  # if rg is larger than 1, fix to 1
  t_rg[which(t_rg>1)]<-1
  #
  p_rg[upper.tri(p_rg)]<-p_r[upper.tri(p_r)]
  p_rg[lower.tri(p_rg)]<-p_l[lower.tri(p_l)]
  diag(p_rg)<-diag(p_lr)
  #
  colnames(t_rg)<-rownames(t_rg)<-colnames(p_rg)<-rownames(p_rg)<-colnames(t_l)
  # plot
  w<-rev(colnames(t_rg)) # to reverse order
  ggcorrplot(corr=t_rg[w,w], p.mat=p_rg[w,w],
             # cluster="hclust",
             # method="circle",
             outline.col = "white",
             insig="blank",
             show.diag=TRUE,
             type="full",
             lab=TRUE,
             legend.title = "rg",
             sig.level=0.05
  ) +
    ggplot2::scale_fill_gradient2(limit = c(-1,1), 
                                  low=brewer.pal(8,name="PuOr")[1],
                                  mid="white",
                                  high=brewer.pal(8,name="PuOr")[8],
                                  midpoint = 0) +
    cowplot::theme_minimal_grid(font_size = 10) + panel_border() +
    labs(x="Right",y="Left",fill="rg") +
    NULL
  # return(c)
}
# for not lateralized measures
ggcorrplot_rg_notlat_atlas<-function(atlas=rg_atlas,lr=rg_lr,parc="FAST",m="VOLUME"){
  tmp<-subset(atlas,parcellation==parc & measure==m)
  t1<- tmp %>% 
    filter(  (toupper(hemisphere)!="R"&toupper(hemisphere)!="L") & 
               (toupper(hemisphere)!="RH"&toupper(hemisphere)!="LH") &
               (toupper(hemisphere)!="RIGHT"&toupper(hemisphere)!="LEFT")
    ) %>%
    # filter(!is.na(hemisphere)) %>%
    dplyr::select(region.p2,region.p1,rg,p)
  colnames(t1)<-c("p1","p2","rg","p")
  #
  t_rownames<- t1 %>% dplyr::select(p1,p2,rg) %>% spread(p1,rg) %>% dplyr::select(1) %>% pull() %>% as.character()
  t_rg<-t1 %>% dplyr::select(p1,p2,rg) %>% spread(p1,rg) %>% dplyr::select(-1) %>% as.matrix()
  p_rg<-t1 %>% dplyr::select(p1,p2,p) %>% spread(p1,p) %>% dplyr::select(-1) %>% as.matrix()
  # names
  rownames(t_rg) <- rownames(p_rg) <- colnames(t_rg) <- colnames(p_rg)<- t_rownames
  # sort by level
  t_levs<-levels(t1$p1)[levels(t1$p1) %in% t_rownames]
  t_rg<-t_rg[t_levs,t_levs]
  p_rg<-p_rg[t_levs,t_levs]
  #
  w<-rev(colnames(t_rg)) # to reverse order
  ggcorrplot(corr=t_rg[w,w], p.mat=p_rg[w,w],
             # method="circle",
             outline.col = "white",
             insig="blank",
             show.diag=TRUE,
             type="lower",
             lab=TRUE,
             sig.level=0.05
  ) +
    ggplot2::scale_fill_gradient2(limit = c(-1,1), 
                         low=brewer.pal(8,name="PuOr")[1],
                         mid="white",
                         high=brewer.pal(8,name="PuOr")[8],
                         midpoint = 0) +
    cowplot::theme_minimal_grid() + panel_border() +
    labs(x="Vermis",y="Vermis",fill="rg") +
    NULL
  # return(c)

}

