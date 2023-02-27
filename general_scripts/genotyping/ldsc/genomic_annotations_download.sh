# donwload evolutionary annotations from different studies
# goal: to partition heritability using these annotations
# an overview of the annotations is in: annotations2download.xlsx
#-------------------------------------------
source activate ldsc

working_dir=/export/home/acarrion/acarrion/projects/resources/genomic_annotations/
mkdir -p ${working_dir}
mkdir -p ${working_dir}/bedfiles/
#-------------------------------------------
# Fetal human gained enhancers
## ref: Reilly et al. 2015, doi: 10.1126/science.1260943
#-------------------------------------------
cd ${working_dir}
mkdir -p Reilly2015_GSE63648
cd Reilly2015_GSE63648
ftp_dir=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63648/suppl/
 
wget ${ftp_dir}GSE63648_12Fpcw_ac_Hu_gain.bed.gz
wget ${ftp_dir}GSE63648_12Fpcw_me2_Hu_gain.bed.gz
wget ${ftp_dir}GSE63648_12Opcw_ac_Hu_gain.bed.gz
wget ${ftp_dir}GSE63648_12Opcw_me2_Hu_gain.bed.gz
wget ${ftp_dir}GSE63648_7pcw_ac_Hu_gain.bed.gz
wget ${ftp_dir}GSE63648_7pcw_me2_Hu_gain.bed.gz
wget ${ftp_dir}GSE63648_8_5pcw_ac_Hu_gain.bed.gz
wget ${ftp_dir}GSE63648_8_5pcw_me2_Hu_gain.bed.gz
wget ${ftp_dir}GSE63648_ReadMe_Bed_files_contents.txt

# extract gz files
for f in $(ls *gz)
 do
 gunzip $f
 done

# autosomic regions only
for t in 7pcw 8_5pcw 12Fpcw 12Opcw
 do
 grep -v chrX GSE63648_*${t}*_Hu_gain.bed | awk 'BEGIN {FS=":";OFS="\t"} {print $2,$1}' | sed 's/ /\t/g' > ${working_dir}/bedfiles/Fetal_brain_GSE63648_${t}_Hu_gain.bed
 done

# merge across all fetal tiempoints
# bedtools merge requires that you presort your data by chromosome and then by start position 
cat Fetal_brain_GSE63648*bed | sort -k1,1 -k2,2n > Fetal_brain_GSE63648_Hu_gain.all.sorted.bed
bedtools merge -i Fetal_brain_GSE63648_Hu_gain.all.sorted.bed > Fetal_brain_GSE63648_Hu_gain.all_sorted.merged.bed
rm Fetal_brain_GSE63648_Hu_gain.all.sorted.bed
#-------------------------------------------
# Adult human gained enhancers
## ref: Vermunt et al. 2016
#-------------------------------------------
# to extract bed files and liftover to hg19
PATH=/export/home/acarrion/acarrion/projects/resources/bin/liftOver/:${PATH}
chain_file=/export/home/acarrion/acarrion/projects/resources/bin/liftOver/hg38ToHg19.over.chain
mkdir -p ${working_dir}Vermunt2016_STables
cd ${working_dir}Vermunt2016_STables

# small script to reorganize STables, per region
for f in $(ls *xlsx)
 do
 Rscript S_Tables_xls2bed.R ${f}
 done
# lift over coordinates to hg19
for f in $(ls *hg38.bed )
 do
 f2=$(echo $f | sed 's/hg38.bed/hg19.bed/g' )
 liftOver ${f} ${chain_file} ${f2} ${f2}_unmapped
 cp ${f2} ${working_dir}/bedfiles/Adult_${f2}
 done

# merge enhancers and promoters
for f1 in $(ls Adult_enhancers_*.bed | grep -v merged.bed)
 do
 f2=$(echo ${f1} | sed 's/_enhancers_/_promoters_/g')
 f=$(echo ${f1} | sed 's/_enhancers_/_enhancers_promoters_/g' | sed 's/.bed//g')
 cat ${f1} ${f2} | sort -k1,1 -k2,2n > ${f}.bed
 bedtools merge -i ${f}.bed > ${f}.merged.bed
 rm ${f}.bed
 done

#-------------------------------------------
# Human accelerated regions
## ref: Capra2013
#-------------------------------------------
# supplementary Table S5 -> copied to txt file
mkdir -p ${working_dir}Capra2013_STables
cd ${working_dir}Capra2013_STables
grep -v chrX Capra2013_TableS5.txt | sed 's/ /\t/g' > ${working_dir}/bedfiles/ncHAR_Capra2013.bed


#-------------------------------------------
# Human accelerated regions
## ref: Gittelman2015
#-------------------------------------------
# supplementary Table S2 -> copied to txt file
mkdir -p ${working_dir}Gittelman2015_STables
cd ${working_dir}Gittelman2015_STables
grep -v chrX Gittelman2015_TableS2.txt | sed 's/ /\t/g' > ${working_dir}/bedfiles/dhsHAR_Gittelman2015.bed

#-------------------------------------------
# Depleted regions (neanderthal and neanderthal + denisovan)
## ref: Vernot2016
#-------------------------------------------
mkdir -p ${working_dir}Vernot2016_STables
cd ${working_dir}Vernot2016_STables
# downloaded supplementary files to this directory from: https://www.science.org/doi/10.1126/science.aad9416

# generate bed files from txt files
grep -v Chr Table_S8.txt | sed 's/^/chr/g' | sed 's/ /\t/g' > ${working_dir}/bedfiles/NeanderthalDepleted_Vernot2016.bed
grep -v Chr Table_S9.txt | sed 's/^/chr/g' | sed 's/ /\t/g' > ${working_dir}/bedfiles/NeanderthalDenisovanDepleted_Vernot2016.bed

#-------------------------------------------
# Depleted regions (neanderthal, IBDmix)
## ref: Chen et al. 2020
#-------------------------------------------
mkdir -p ${working_dir}Chen2020_IBDmix
cd ${working_dir}Chen2020_IBDmix
# download file from 

# generate bed file from Supplementary Table S8:
## also reported in Buisan et al. 2022
echo "chr1 105400000 120600000
chr3 74100000 89300000
chr7 106200000 123200000
chr8 49400000 66500000" | sed 's/ /\t/g' > ${working_dir}/bedfiles/NeanderthalDepleted_Chen2020.bed

#-------------------------------------------
# Selective sweeps
## ref: Peyregne2017
#-------------------------------------------
mkdir -p ${working_dir}Peyregne2017_STables
cd ${working_dir}Peyregne2017_STables
# download supplementary Tables S1 and S2 here
# use Table S2: extended set
grep -v "#" Supplemental_File_S2.tsv | sed 's/ /\t/g' > ${working_dir}/bedfiles/SelectiveSweep_Peyregne2017.bed


#-------------------------------------------
# Neanderthal introgressed SNPs
#-------------------------------------------
## Browning et al. 2018, doi: 10.1016/j.cell.2018.02.031
## Browning, Sharon (2018), “Sprime results for 1000 Genomes non-African populations and SGDP Papuans”, Mendeley Data, V1, doi: 10.17632/y7hyt83vxr.1
mkdir -p ${working_dir}Browning2018_Sprime
cd ${working_dir}Browning2018_Sprime
## downloaded manually from: https://data.mendeley.com/datasets/y7hyt83vxr/1
#unzip Sprime\ results\ for\ 1000\ Genomes\ non-African\ populations\ and\ SGDP\ Papuans.zip

for f in $(ls *gz)
 do
 f2=$(echo $f | sed 's/.tar.gz/.txt/g')
 zcat ${f} > ${f2}
 done
# get all Neanderthal match, across populations
awk '{ if ($9 == "match") { print $1,$2,$3,$4,$5,$9,$10} }' *sprime_results.txt | sort | uniq > sprime_Nmatch_all_results.txt
awk '{ if ($7 == "match") { print }}' sprime_Nmatch_all_results.txt > sprime_NmatchDmatch_all_results.txt
awk '{ if ($7 == "mismatch") { print }}' sprime_Nmatch_all_results.txt > sprime_NmatchDmis_all_results.txt
# clean intermediate files
rm *sprime_results.txt

# generate bed file where each SNP is a single region
awk '{print $1,$2,$2,$3}' sprime_Nmatch_all_results.txt | sed 's/^/chr/g' | sed 's/ /\t/g' > ${working_dir}/bedfiles/NeanderthalMatchSNPs_sprime_Browning2018.bed

#-------------------------------------------
# Neanderthal introgressed SNPs
# ref: Simonti et al. 2016
# ???
#-------------------------------------------