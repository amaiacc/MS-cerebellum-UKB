# get sample overlap from ldsc output and generate matrix to input lava
## adapted from sample_overlap.md, LAVA TUTORIAL (Josefin Werme, CTG Lab, VU Amsterdam)
#------------------------------------
# Set directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="F:/projects/"} else {dir="/export/home/acarrion/acarrion/projects/"}
# working dir
ldsc_dir=paste0(dir,"resources/datasets/GWAS_sumstats/ldsc/UKB_BIG40/")
lava_dir=paste0(dir,"cerebellum_UKB/data/lava/")
# define phenos to include
phenos2check<-c("0141","0143","0146","Grove_ASD","Lee2018_CP")

# read data for all rg runs
# from UKB IDPs
input_files1<-list.files(ldsc_dir,pattern="table")
input_files1<-input_files1[grep("rg",input_files1)]
input_files1<-paste0(ldsc_dir,input_files1)
# from other GWAS sumstats
input_files2<-list.files(paste0(ldsc_dir,"../"),pattern="table")
input_files2<-input_files2[grep("rg",input_files2)]
input_files2<-paste0(ldsc_dir,"../",input_files2)

for (f in c(input_files1,input_files2)){
  tmp<-read.table(f,header=FALSE)
  if (exists("scor")){scor<-rbind(scor,tmp)} else {scor<-tmp}
  rm(tmp)
}
rm(f)

# clean first col
tmp<-strsplit(scor$V1,".log")
scor$V1<-gsub("^-|^:","",sapply(tmp,"[[",2))
rm(tmp)
# get header
colnames(scor)<-scor[1,]
# clean dataset
scor<-subset(scor,p1!="p1")
#### Extract the intercepts and create a sampling correlation matrix (in R)
scor = scor[,c("p1","p2","gcov_int")]             # retain key headers
# clean pheno cols
# if path is avail, get after pathonly
scor$p1[grep("ldsc/",scor$p1)]<-sapply(strsplit(scor$p1[grep("ldsc/",scor$p1)],"ldsc/"),"[[",2)
scor$p2[grep("ldsc/",scor$p2)]<-sapply(strsplit(scor$p2[grep("ldsc/",scor$p2)],"ldsc/"),"[[",2)

# clean names
scor$p1 = gsub(".sumstats.gz|../|^/","",scor$p1)   # assuming the munged files have format [phenotype]_munge.sumstats.gz
scor$p2 = gsub(".sumstats.gz|../|^/","",scor$p2)   # (adapt as necessary)

# select only phenotypes of interest
w1<-which(scor$p1 %in% phenos2check)
w2<-which(scor$p2 %in% phenos2check)
w<-intersect(w1,w2)
scor<-scor[w,]
rm(w1,w2,w)

# duplicate to get full matrix..
scor2<-scor[,c("p2","p1","gcov_int")]
colnames(scor2)<-c("p1","p2","gcov_int")
scor<-rbind(scor,scor2)
rm(scor2)

# remove duplicates
w<-which(duplicated(scor))
if(length(w)>0){
  scor<-scor[-w,]
}
rm(w)
#
phen = unique(c(scor$p1,scor$p2))                  # obtain list of all phenotypes (assuming all combinations have been analysed)
n = length(phen)
mat = matrix(NA,n,n)                    # create matrix
rownames(mat) = colnames(mat) = phen    # set col/rownames

for (i in phen) {
  for (j in phen) {
    if (NROW(subset(scor, p1==i & p2==j))>0){
      mat[i,j] = subset(scor, p1==i & p2==j)$gcov_int
    } else if (i==j){
      mat[i,j]<-1 # if it's the same sumstats, assume 1 even if not estimated
    }
  }
}
rm(i,j)
# make numeric matrix
mat<-matrix(as.numeric(mat),ncol=n)
rownames(mat) = colnames(mat) = phen    # set col/rownames

if (!all(t(mat)==mat,na.rm = TRUE)) { mat[lower.tri(mat)] = t(mat)[lower.tri(mat)] }  # sometimes there might be small differences in gcov_int depending on which phenotype was analysed as the outcome / predictor
mat = round(cov2cor(mat),5)                   # standardise
write.table(mat, paste0(lava_dir,"sample.overlap.txt"), quote=F)   # save
