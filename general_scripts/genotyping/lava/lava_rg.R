# This script will analyse the pariwise bivariate local genetic correlations between any two phenotypes of interest
# call from the command line as: Rscript lava_rg.R $REF_PREFIX $LOCUS_FILE $INPUT_INFO_FILE $SAMPLE_OVERLAP_FILE "$PHENO1;$PHENO2" $OUTPUT_FILENAME
# (e.g. Rscript lava_rg.R g1000_eur.maf01 blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile input.info.txt sample.overlap.txt 'mdd;scz' 'mdd.scz')
arg = commandArgs(T); refdat = arg[1]; locfile = arg[2]; infofile = arg[3]; overlapfile = arg[4]; phenos = unlist(strsplit(arg[5],";")); outname = arg[6]
library(LAVA); library(parallel)
library(dplyr); library(tidyr)

# univ threshold for bivar analysis
univ.thresh =1e-4 # .05 / 2495 # 0.05/n.loc(loci)

### READ IN DATA ###

#### read univariate results for each phenotype, to define the loci to be analysed
u<-lapply(phenos,function(p){
  t<-read.table(paste0(p,".univ"),header=TRUE)
  t$pheno<-p
  return(t)
})
univ<-do.call("rbind",u)
rm(u)

# filter univ base don univ.thresh
univ<-subset(univ, p<univ.thresh)

# wide format for univ
univ_wide<-univ %>% pivot_wider(id_cols = c(locus,chr,start,stop),
                                names_from=pheno,
                                values_from=c(h2.obs,p,n.snps,n.pcs))
univ_wide$pvals_NA<-apply((univ_wide[,paste0("p_",phenos)]),1,function(x) sum(is.na(x)))
# remove any locus with NAs for one of the phenos
univ_wide<-subset(univ_wide,pvals_NA==0)

#### read input files for analysis
input = process.input(infofile, overlapfile, refdat, phenos)
# define loci to be included in the bivar analysis
loci = read.loci(locfile)
# subset loci to analyze
loci= loci %>% filter(LOC %in% univ_wide$locus)
n.loc=nrow(loci)

#### ANALYSE ####
print(paste("Starting LAVA analysis for",n.loc,"loci",
            "with non-zero heritabilities (p<",univ.thresh,") for phenotypes:", paste(phenos,collapse=","))
      )
progress = ceiling(quantile(1:n.loc, seq(.05,1,.05)))   # (used for printing the progress)
u = b = NULL

# Using mclapply to analyse the loci in parallel
# adjust the mc.cores argument according to the number of available cores (e.g. N avail cores - 1)
ncores=detectCores()-1
# if you dont want to parallelise you can use lapply or a for loop instead
out = mclapply(1:n.loc, mc.cores=ncores, function(i) {
# out = lapply(1:n.loc, function(i) {
	if (i %in% progress) print(paste("..",names(progress[which(progress==i)])))  # print progress

	# process locus
	locus = process.locus(loci[i,], input)

	# in some cases the locus object cannot be created due to e.g too few SNPs or negative variances in all analysed phenotypes, hence this check
	if (!is.null(locus)) {
		# extract general locus info for output
		loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)

		# run univ & bivar analysis functions & store results
		# ub = run.univ.bivar(locus, univ.thresh=univ.thresh)
		# u = cbind(loc.info, ub$univ)
		ub=run.bivar(locus)
		if (!is.null(ub)) {
			b = cbind(loc.info, ub)
		}
	}
	return(list(bivar=b))
})

head(out)

### WRITE OUTPUT ###
write.table(do.call(rbind, lapply(out, "[[", "bivar")), paste0(outname,".bivar"), row.names=F, quote=F)
write.table(univ, paste0(outname,".univ"), row.names=F, quote=F)


