# create partition file using the bed files that are used by supergnova (to make results comparable)
working_dir=/export/home/acarrion/acarrion/projects/general_scripts/genotyping/supergnova/partition/
lava_locfile=/export/home/acarrion/acarrion/projects/resources/github/LAVA/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile
#
mkdir -p ${working_dir}
cd ${working_dir}

awk '{print $2,$3,$4}' ${lava_locfile} > lava_partition_tmp.txt

for chr in {1..22}; do
    echo "chr start stop\n" >> lava_chr${chr}.bed
    grep -e "^${chr} " lava_partition_tmp.txt >> lava_chr${chr}.bed
done

rm lava_partition_tmp.txt
