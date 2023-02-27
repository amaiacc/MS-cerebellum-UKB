# download files required to run LAVA
## from: https://ctg.cncr.nl/software/lava

##  Reference genotype data in plink format (.bim, .bed, .fam), used for the estimation of LD
## 1000 genomes (pre-processed input files can be found here), hg19
working_dir=/export/home/acarrion/acarrion/projects/resources/reference_data/magma/
mkdir -p ${working_dir}
cd ${working_dir}
# wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip
## manual download, wget gave error
unzip g1000_eur.zip
rm g1000_eur.zip


# clone from github to get example files etc
cd /export/home/acarrion/acarrion/projects/resources/github
git clone https://github.com/josefin-werme/LAVA.git
# Locus definition file
## File that defines loci either based on genomic coordinates or a list of SNPs
## (the locus file that we used in the LAVA preprint can be found in the support_data folder; 
#### this file was obtained via https://github.com/cadeleeuw/lava-partitioning using the g1000 data phase 3, build GRCh37/hg19). 
## The locus file requires the following headers: LOC CHR, START, STOP SNPS