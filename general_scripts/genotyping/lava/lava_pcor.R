# This script will perform partial correlations:
## in loci that have significant bivar correlations for two traits, conditioning on trait 3

# call from the command line as: Rscript lava_pcor.R $REF_PREFIX $LOCUS_FILE $INPUT_INFO_FILE $SAMPLE_OVERLAP_FILE "$PHENO1;$PHENO2" $OUTPUT_FILENAME
# (e.g. Rscript lava_rg.R g1000_eur.maf01 blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile input.info.txt sample.overlap.txt 'mdd;scz' 'mdd.scz')
arg = commandArgs(T); refdat = arg[1]; locfile = arg[2]; infofile = arg[3]; overlapfile = arg[4]; phenos = unlist(strsplit(arg[5],";")); pheno_cond=arg[6] ; outname = arg[7]
# # to test runs
# setwd("F:/projects/cerebellum_UKB/data/lava/TMP~1/")
# refdat="g1000_eur"; locfile="blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
# infofile="cerebellum_fusiform_input.info.txt"; overlapfile="cerebellum_fusiform_sample.overlap.txt"
# phenos=c("0146","PGC3_SCZ_primary"); pheno_cond="0102"
# outname=paste0(paste(phenos,collapse="."),"_condBy",pheno_cond)

library(LAVA); library(parallel)
library(dplyr); library(tidyr)

# bivar.thresh
bivar.thresh = 0.05/1387 # .05

### READ IN DATA ###

#### read univariate analysis for all phenos
u<-lapply(c(phenos,pheno_cond),function(p){
  t<-read.table(paste0(p,".univ"),header=TRUE)
  t$pheno<-p
  return(t)
})
univ<-do.call("rbind",u)
rm(u)

univ_cond<-subset(univ,pheno==pheno_cond)
#### read bivariate analysis for pairs of phenos, to define the loci to be analysed
bivar<-read.table(paste0(phenos[1],".",phenos[2],".bivar"),header=TRUE)
# define loci to include in analysis, from bivar results
bivar2<-subset(bivar,p<bivar.thresh)
# loci to test
loci2test=univ_cond %>% filter(locus  %in% bivar2$locus) %>% filter(p<0.05) %>% pull(locus)

if(length(loci2test)>0){
  
#### read input files for analysis
input = process.input(input.info.file=infofile, 
                      # sample.overlap.file=NULL,
                      sample.overlap.file=overlapfile,
                      ref.prefix=refdat, phenos=c(phenos,pheno_cond))

# define loci to be included in the bivar analysis
loci = read.loci(locfile)
# subset loci to analyze
loci=loci %>% filter(LOC %in% loci2test)
n.loc=nrow(loci)

#### ANALYSE ####
print(paste("Starting LAVA analysis for",n.loc,"loci",
            "for phenotypes", paste(phenos,collapse=","),
            "conditioned by: ",pheno_cond)
      )

progress = ceiling(quantile(1:n.loc, seq(.05,1,.05)))   # (used for printing the progress)
pcor=NULL

# # Using mclapply to analyse the loci in parallel
# # adjust the mc.cores argument according to the number of available cores (e.g. N avail cores - 1)
# ncores=detectCores()-1
# out = mclapply(1:n.loc, mc.cores=ncores, function(i) {
# 
# # if you dont want to parallelise you can use lapply or a for loop instead
out = lapply(1:n.loc, function(i) {
	if (i %in% progress) print(paste("..",names(progress[which(progress==i)])))  # print progress

	# process locus
	locus = process.locus(loci[i,], input)

	# in some cases the locus object cannot be created due to e.g too few SNPs or negative variances in all analysed phenotypes, hence this check
	if (!is.null(locus)) {
		# extract general locus info for output
		loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)

		# # run univ & bivar analysis functions & store results
		# ub = run.univ.bivar(locus)
		# u = cbind(loc.info, ub$univ)
		# b = cbind(loc.info, ub$bivar)
		
		partcor=run.pcor(locus, target=phenos, phenos=pheno_cond,
		                 # CIs=FALSE,
		                 # p.values=FALSE,
		                 param.lim = Inf)
		if (!is.null(partcor)) {
			p = cbind(loc.info, partcor)
		}
	}
	return(list(pcor=p))
})

# head(out)
t<-do.call(rbind, lapply(out, "[[", "pcor"))
### WRITE OUTPUT ###
write.table(t, paste0(outname,".pcor"), row.names=F, quote=F)
} else {
  print(paste("There are no loci that have a non-zero local genetic correlation between phenotypes:" , phenos))
}



