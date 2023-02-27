# this 'settings' file contains relevant directories and file names used for the analyses
## adapted from settings.sh

project=cerebellum_fusiform

# key directories
lava_dir=/export/home/acarrion/acarrion/projects/resources/github/LAVA # lava files
scripts_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/lava/ # custom general scripts to run lava
sumstats_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/ # general sumstats dir

analysis_dir=/export/home/acarrion/acarrion/projects/cerebellum_UKB/scripts/lava/ # project scripts scripts to run lava
working_dir=/export/home/acarrion/acarrion/projects/cerebellum_UKB/data/lava/ # 
out_dir=${working_dir}/output/

# names of key input files
loc_file=blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile
info_file=${project}_input.info.txt
overlap_file=${project}_sample.overlap.txt

# reference data dir & prefix
ref_dir=/export/home/acarrion/acarrion/projects/resources/reference_data/magma/ # (usually this is not in the main data dir since I only keep a single copy of these files)
ref_dat=g1000_eur

# define TMPDIR
# if not available in local environment 
TMPDIR=${working_dir}TMP
mkdir -p ${TMPDIR}