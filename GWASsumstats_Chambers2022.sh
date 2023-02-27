# Download GWAS sumstats from Chambers et al. 2022 --> total cerebellar volume
## PMID: 35079123 PMCID: PMC9126806 DOI: 10.1038/s41380-022-01443-8 
#--------------------------------------------------------------------------------------
## GWAS catalog
## https://www.ebi.ac.uk/gwas/downloads/summary-statistics
# Study accession: Reported trait
# GCST90020190: Total cerebellar volume (excluding Crus I vermis)
# GCST90020191: Anterior lobe of cerebellar volume (including I to V lobules)
# GCST90020192: Superior Posterior lobe of cerebellar volume (including VI to Crus I hemispheric lobules)
# GCST90020193: Superior Posterior lobe of cerebellar volume (including VI to Crus I vermal lobules, excluding Crus I vermis)
# GCST90020194: Inferior Posterior lobe of cerebellar volume (including Crus II to IX hemispheric lobules)
# GCST90020195: Inferior Posterior lobe of cerebellar volume (including Crus II to IX vermal lobules)
# GCST90020196: Flocculonodular lobe of cerebellar volume (including X hemispheric lobules)
# GCST90020197: Flocculonodular lobe of cerebellar volume (including X vermal lobules)
gwas_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/UKB/Chambers2022/
mkdir -p ${gwas_dir}
cd ${gwas_dir}
# download total cerebellar volume only
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90020001-GCST90021000/GCST90020190/GCST90020190_buildGRCh37.tsv.gz
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90020001-GCST90021000/GCST90020190/md5sum.txt
# check md5sum
md5sum GCST90020190_buildGRCh37.tsv.gz
cat md5sum.txt
#--------------------------------------------------------------------------------------
# Format for LDSC:
#source activate ldsc
PATH=/export/home/acarrion/acarrion/projects/resources/github/ldsc/:${PATH}
scripts_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/ldsc/
resource_dir=/export/home/acarrion/acarrion/projects/resources/
ldscores_dir=${resource_dir}/LDscores/
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/
# 
cd ${gwas_dir}
p_file=GCST90020190_buildGRCh37.tsv.gz
#zless ${p_file}
# variant_id      chromosome      base_pair_location      effect_allele   other_allele    beta    standard_error  p_value direction
# rs2326918       6       130840091       a       g       0.0017  0.0055  0.7516  -+
# rs112634005     3       104998275       a       c       -0.0045 0.0055  0.4167  --
p=Chambers2022_TotalCerebellarVolume
N1=33265
if [ ! -f ${working_dir}/${p}_h2.log ]
then
 qsub -v p_file=${p_file},N1=${N1},p=${p},snp_id=variant_id ${scripts_dir}/ldsc_munge_sumstats.sh 
 fi

 