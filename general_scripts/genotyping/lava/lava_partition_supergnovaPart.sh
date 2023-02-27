# create partition file using the bed files that are used by supergnova (to make results comparable)
working_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/lava/
locfile=/export/home/acarrion/acarrion/projects/resources/github/LAVA/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile
#
cd /export/home/acarrion/acarrion/projects/resources/github/SUPERGNOVA/data/partition/

head ${locfile} -n 1 > ${working_dir}/supergnova_partition.hg19.locfile
cat *bed | sort -nk1 -nk3 | grep -v chr | awk '{print NR,$0}' >> ${working_dir}/supergnova_partition.hg19.locfile
