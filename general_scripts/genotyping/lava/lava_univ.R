# This script will analyse the univariate local heriatability

# call from the command line as: 
# Rscript lava_h2.R $REF_PREFIX $LOCUS_FILE $INPUT_INFO_FILE $SAMPLE_OVERLAP_FILE "$PHENO1;$PHENO2" $OUTPUT_FILENAME
# (e.g. Rscript lava_h2.R g1000_eur.maf01 blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile input.info.txt sample.overlap.txt 'mdd;scz' 'mdd.scz')
arg = commandArgs(T); refdat = arg[1]; locfile = arg[2]; infofile = arg[3]; overlapfile = arg[4]; pheno = arg[5]; outname = arg[6]
library(LAVA); library(parallel)

### check that output file does not exist
if (!file.exists(paste0(outname,".univ"))) {

### READ IN DATA ###
loci = read.loci(locfile);  n.loc = nrow(loci)
input = process.input(infofile, overlapfile, refdat, pheno)


# univ threshold for bivar analysis
# univ.thresh =1e-4 # .05 / 2495 # 0.05/n.loc(loci)


#### ANALYSE ####
print(paste("Starting LAVA analysis"))

#### univariate analysis, for each phenotype ####
print(paste("Univariate analysis for",n.loc,"loci"))

progress = ceiling(quantile(1:n.loc, seq(.05,1,.05)))   # (used for printing the progress)
univ=NULL

# Using mclapply to analyse the loci in parallel
# adjust the mc.cores argument according to the number of available cores (e.g. N avail cores - 1)
ncores=detectCores()-1

# # to test
# ncores=1
# n.loc=10

# if you dont want to parallelise you can use lapply or a for loop instead
out = mclapply(1:n.loc, mc.cores=ncores, function(i) {
	if (i %in% progress) print(paste("..",names(progress[which(progress==i)])))  # print progress

	# process locus
	locus = process.locus(loci[i,], input)

	# in some cases the locus object cannot be created due to e.g too few SNPs or negative variances in all analysed phenotypes, hence this check
	if (!is.null(locus)) {
		# extract general locus info for output
		loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)

		# run univ analysis functions
		univ = run.univ(locus)
		univ = cbind(loc.info, univ)
	}
	return(list(univ=univ))
})

# head(out)

### WRITE OUTPUT for univ###
write.table(do.call(rbind, lapply(out,"[[","univ")), paste0(outname,".univ"), row.names=F, quote=F)
}