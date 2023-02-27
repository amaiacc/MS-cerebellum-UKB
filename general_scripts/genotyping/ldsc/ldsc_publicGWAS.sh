#!/bin/bash
#--------------------------------------
# donwnload and organize the summary statistics from publicly available GWASes
#--------------------------------------
gwas_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/downloaded_data/
working_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/ldsc/
cd ${gwas_dir}


#--------------------------------------
# convert to sumstats for LDSC
#--------------------------------------
# load modules and define tools/directories with resources
source activate ldsc
PATH=/export/home/acarrion/acarrion/projects/resources/github/ldsc/:${PATH}
resource_dir=/export/home/acarrion/acarrion/projects/resources/
## info on LD score formats/ estimation
# https://github.com/bulik/ldsc/wiki/LD-File-Formats
# https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
ldscores_dir=${resource_dir}/LDscores/

mkdir -p ${working_dir}
## sumstats format: https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format
# required columns:
##    SNP -- SNP identifier (e.g., rs number)
##  N -- sample size (which may vary from SNP to SNP).
##  Z -- z-score. Sign with respect to A1 (warning, possible gotcha)
##  A1 -- first allele (effect allele)
##  A2 -- second allele (other allele)
# Note that ldsc filters out all variants that are not SNPs and strand-ambiguous SNPs.

# do it one by one, formats may be different
#------------------------------------------------------
#------------------------------------------------------
cd ${gwas_dir}Lee2018_EA3
#1. GWAS_EA_excl23andMe.txt - Educational attainment (EA) GWAS meta-analysis of all discovery cohorts except 23andMe.
#2. GWAS_CP_all.txt - Cognitive performance (CP) GWAS meta-analysis of all discovery cohorts. 
##Sample Size = 257,828.
cd ${gwas_dir}Lee2018_EA3
p=Lee2018_CP
N1=257828 # taken from readme file
p_file=GWAS_CP_all.txt
head GWAS_CP_all.txt
#MarkerName      CHR     POS     A1      A2      EAF     Beta    SE      Pval
#rs2352974       3       49890613        T       C       0.5     -0.0319 0.00285 5.19e-29
if [ ! -f ${working_dir}/${p}_h2.log ]
 then
 munge_sumstats.py \
   --sumstats ${p_file} \
   --N ${N1} \
   --p Pval \
   --out ${working_dir}/${p} \
   --a1 A1 --a2 A2 \
   --merge-alleles ${ldscores_dir}/eur_w_ld_chr/w_hm3.snplist
 #---------------
 echo 'Estimate heritability for' ${p}
 ldsc.py \
   --h2 ${working_dir}/${p}.sumstats.gz \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${working_dir}/${p}_h2
fi

#------------------------------------------------------
cd ${gwas_dir}PGC/SCZ/SCZ2022
## primary dataset: all datasets
p=PGC3_SCZ_primary
p_file=PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz
#zless ${p_file}
#CHROM   ID      POS     A1      A2      FCAS    FCON    IMPINFO BETA    SE      PVAL    NGT     DIRE    NCAS    NCON    NEFFDIV2
#8       rs117278216     100516008       T       C       0.983   0.983   0.904   -0.0104039337153929     0.0366  0.7762  0       -??     67390   94015   78506.5
#8       rs62513865      101592213       C       T       0.937   0.932   0.949   0.0141001243787816      0.0169  0.4053  2       ---     74776   101023  85057.37

# remove header and fitler based on imputation quality
#zcat ${p_file} | grep -v "#" > tmp.tsv
if [ ! -f ${working_dir}/${p}_h2.log ]
then
 munge_sumstats.py \
   --sumstats tmp.tsv \
   --snp ID\
   --N-cas-col NCAS --N-con-col NCON \
   --info IMPINFO \
   --signed-sumstats BETA,0 \
   --p PVAL \
   --out ${working_dir}/${p} \
   --a1 A1 --a2 A2 \
   --merge-alleles ${ldscores_dir}/eur_w_ld_chr/w_hm3.snplist
  #---------------
 echo 'Estimate heritability for' ${p}
 ldsc.py \
   --h2 ${working_dir}/${p}.sumstats.gz \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${working_dir}/${p}_h2
fi
#------------------------------------------------------
cd ${gwas_dir}PGC/ASD  
p=Grove_ASD
p_file=iPSYCH-PGC_ASD_Nov2017.gz 
#iPSYCH-PGC_ASD_Nov2017.gz: Full ASD GWAS meta-analysis of samples of European ancestry
#(18,382 cases, 27,969 controls, see manuscript for methods)
N1=18382
N2=27969

#zless ${p_file}
#CHR     SNP     BP      A1      A2      INFO    OR      SE      P
#8       rs62513865      101592213       T       C       0.949   1.00652 0.027   0.8086
if [ ! -f ${working_dir}/${p}_h2.log ]
then
 munge_sumstats.py \
   --sumstats ${p_file} \
   --N-cas ${N1} --N-con ${N2} \
   --p p \
   --out ${working_dir}/${p} \
   --a1 A1 --a2 A2 \
   --merge-alleles ${ldscores_dir}/eur_w_ld_chr/w_hm3.snplist
  #---------------
 echo 'Estimate heritability for' ${p}
 ldsc.py \
   --h2 ${working_dir}/${p}.sumstats.gz \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${working_dir}/${p}_h2
fi

#------------------------------------------------------
## use GWAS files that have been QC'd for PGS analysis, following script: ../QC_sumstats/prepare_publicGWAS_QC.sh
#------------------------------------------------------
gwas_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/
cd ${gwas_dir}

## Handedness
p=LeftHandedness_CuellarPartida2020
p_file=${p}.QC.gz

#zless ${p_file}.
#rsid    a_1     a_0     Freq1   FreqSE  MinFreq MaxFreq N       Zscore  P       Direction       HetISq  HetChiSq        HetDf   HetPVal chr     pos
#rs2326918       a       g       0.8478  0.0045  0.8454  0.8560  206399  -0.763  0.4456  -+      0.0     0.731   1       0.3926  6       130840091
  
if [ ! -f ${working_dir}/${p}_h2.log ]
then
 munge_sumstats.py \
   --sumstats ${p_file} \
   --p P \
   --N-col Weight \
   --out ${working_dir}/${p} \
   --a1 a_1 --a2 a_0 \
   --merge-alleles ${ldscores_dir}/eur_w_ld_chr/w_hm3.snplist
 #---------------
 echo 'Estimate heritability for' ${p}
 ldsc.py \
   --h2 ${working_dir}/${p}.sumstats.gz \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${working_dir}/${p}_h2
fi

#------------------------------------------------------
gwas_dir=/export/home/acarrion/acarrion/projects/resources/datasets/GWAS_sumstats/QC/
cd ${gwas_dir}

p=Ambidextrous_CuellarPartida2020
p_file=${p}.QC.gz

#zless ${p_file}.

if [ ! -f ${working_dir}/${p}_h2.log ]
then
 munge_sumstats.py \
   --sumstats ${p_file} \
   --N-col N \
   --p P \
   --out ${working_dir}/${p} \
   --a1 a_1 --a2 a_0 \
   --merge-alleles ${ldscores_dir}/eur_w_ld_chr/w_hm3.snplist
 #---------------
 echo 'Estimate heritability for' ${p}
 ldsc.py \
   --h2 ${working_dir}/${p}.sumstats.gz \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${working_dir}/${p}_h2
fi


#------------------------------------------------------
# rest if any left

for p_file in $(ls *QC.gz | grep -v PGC)
do

#p_file=${p}.QC.gz
p=$(echo $p_file | sed 's/.QC.gz//g')

#chr     rsid    pos     a_1     a_0     a_1_Frequency   beta    se      P       N
#10      rs35418599      69083   T       C       0.7317  -0.0037 0.0018  0.0382120203183735      270059
#10      rs185642176     90127   T       C       0.0821  -0.0025 0.0017  0.145894862822189       270059

if [ ! -f ${working_dir}/${p}_h2.log ]
then
 munge_sumstats.py \
   --sumstats ${p_file} \
   --N-col N \
   --p P \
   --out ${working_dir}/${p} \
   --a1 a_1 --a2 a_0 \
   --merge-alleles ${ldscores_dir}/eur_w_ld_chr/w_hm3.snplist
 #---------------
 echo 'Estimate heritability for' ${p}
 ldsc.py \
   --h2 ${working_dir}/${p}.sumstats.gz \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${working_dir}/${p}_h2
fi
done
#------------------------------------------------------


#------------------------------------------------------
# summary table  h2
h2_file=$(ls -1 *h2.log )
paste <(echo 'File') <(echo 'h2 (se)') <(echo 'Lambda GC') <(echo 'Mean Chi^2') <(echo 'Intercept (se)') <(echo 'Ratio (se)')> summary_h2_ldsc.table
paste <(ls -1 ${h2_file}) \
      <(grep "h2:" ${h2_file} | awk '{print $5,$6}' ) \
      <(grep "Lambda GC:" ${h2_file} | awk '{print $3}') \
      <(grep "Mean Chi^2" ${h2_file} | awk '{print $3}') \
      <(grep "Intercept" ${h2_file} | awk '{print $2,$3}') \
      <(grep "Ratio" ${h2_file} | awk '{print $2,$3}') >> summary_h2_ldsc.table

#------------------------------------------------------
