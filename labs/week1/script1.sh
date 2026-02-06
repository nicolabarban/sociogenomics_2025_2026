#!/bin/bash


mkdir Sociogenomics
mkdir Sociogenomics/Data
mkdir Sociogenomics/Results
mkdir Sociogenomics/Software
mkdir Sociogenomics/Scripts

wget http://nicolabarban.com/sociogenomics2022/week1/data/week1.zip
mv week1.zip Sociogenomics/Data/

cd Sociogenomics/Data
unzip week1.zip 

head hapmap1.map 

head -10 hapmap1.map 

wc -l hapmap1.map  
wc -l hapmap1.ped


echo 'Hello' \
 'world'

echo 'Hello world' | wc

grep rs7540009 hapmap1

grep rs75 hapmap1

grep rs75 hapmap1 | wc -l



     
wget -O sumstatsUKB_height.tsv.gz https://www.dropbox.com/s/ou12jm89v74k55e/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz?dl= 

gunzip -d sumstatsUKB_height.tsv.gz

head sumstatsUKB_height.tsv
wc sumstatsUKB_height.tsv

more sumstatsUKB_height.tsv

cat BMI_pheno.txt
awk '{print}' BMI_pheno.txt


awk '{print $1, $2  }'  BMI_pheno.txt

awk ' /3:49860854/ {print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv
awk ' {print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv | grep 3:49860854


awk '{ print  NR, $1, $2  }'  BMI_pheno.txt

awk '{print $1,$NF}' BMI_pheno.txt

awk 'NR==3, NR==6 {print NR,$0}'  BMI_pheno.txt

awk '{print NR "- " $1 }' BMI_pheno.txt

awk '{ if(NR==1) print $1, $2  }'  BMI_pheno.txt

awk '{ if(NR<10) print $1, $2  }'  BMI_pheno.txt


awk '{ if(NR<10) print   }'  sumstatsUKB_height.tsv

awk '{ if($11<5e-08) print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv

awk '{ if($11<5e-08) print  $1, $2, $5, log($11) }'  sumstatsUKB_height.tsv

 awk 'END { print NR }' sumstatsUKB_height.tsv 

awk 'FS=":" {print $1, $2, $3}' sumstatsUKB_height.tsv | head
awk 'FS=":", OFS="-" {print $1, $2, $3}' sumstatsUKB_height.tsv | head

awk -F '[:"\t"]' '{ if(NR>1) print $1, $2, $3 , $4, $13}' sumstatsUKB_height.tsv | head

awk '{ if($11<5e-08) print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv > ../Results/sign_variants_UKB.txt

awk '{ if($3>0.1) print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv > ../Results/common_variants_UKB.txt

head  ../Results/common_variants_UKB.txt

awk -F '[:"\t"]' '{ if($6 >0.1 && $1==1) print $1, $2, $3, $4,  $6}' sumstatsUKB_height.tsv | wc -l

awk -F '[:"\t"]' '{ if($6 >0.1 && $1==21) print $1, $2, $3, $4,  $6}' sumstatsUKB_height.tsv | wc -l

cd $HOME
cd Sociogenomics/Software

wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip
unzip plink_linux_x86_64_20210606.zip

chmod +x plink

./plink --help 


cd $HOME/Sociogenomics
 ln -s Software/plink

./plink --help 

./plink --file   Data/hapmap1 --freq --out Results/test

cat Results/test.frq | more
