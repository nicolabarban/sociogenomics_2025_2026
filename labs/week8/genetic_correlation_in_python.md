
## This  part of the lab needs python 2.7 installed in your computer. We recommend the anaconda distribution

### The next commands show how to install it in google shell
First, some data cleaning to free some space on the google shell account

```
cd /home/nicola_barban/Sociogenomics/Data/
rm *.*

cd /home/nicola_barban/Sociogenomics/Results/
rm *.*
```

This requires time and installs miniconda (which is a python distribution) in the system
```
cd /home/nicola_barban/Sociogenomics/Software
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
chmod +x miniconda.sh
./miniconda.sh  -b
rm ~/miniconda.sh
source $HOME/miniconda3/bin/activate

```


This installs LDSC in your system
```
cd /home/nicola_barban/Sociogenomics/Software
git clone https://github.com/bulik/ldsc.git


cd ldsc
```

Congratulations you have installed LDSC
now you have to install some packeges that LDSC uses for its calculations

```
conda env create --file environment.yml
source activate ldsc
```

Now, let's check if everything works
```
./ldsc.py -h
./munge_sumstats.py -h
```

This is the LDSC command to estimate heritability.
We need the following files:
 1. Giant_Height2018.txt.gz
 2. w_hm3.snplist



```
cd /home/nicola_barban/Sociogenomics

```
First we need to "munge" summary statistics, that is to make ti comparable with reference data that have info on LD scores

##**This takes a lot of time. you can skip this and find the prepared data in the folder**
```
python Software/ldsc/munge_sumstats.py \
  --sumstats Data/Giant_Height2018.txt.gz \
  --snp SNP \
  --N-col N \
  --a1 Tested_Allele \
  --a2 Other_Allele \
  --p P \
  --frq Freq_Tested_Allele_in_HRS \
  --out Results/height \
  --merge-alleles Data/w_hm3.snplist

```

The folder *eur_w_ld_chr/* contains the LD scores that are merged with our summary statistics


```
python Software/ldsc/ldsc.py \
 --h2 Data/height.sumstats.gz \
 --ref-ld-chr Data/eur_w_ld_chr/ \
 --w-ld-chr Data/eur_w_ld_chr/ \
 --out Results/height_h2

```

Now we can try to estimate genetic correlation between height and educational attainment.

First, we need to prepare data for Educational Attainment

##**This takes a lot of time. you can skip this and find the prepared data in the folder**

```
python Software/ldsc/munge_sumstats.py \
   --sumstats Data/EA2_results.txt.gz \
   --snp SNP\
   --N 500000 \
   --a1 A1 \
   --a2 A2 \
   --p Pval \
   --frq EAF \
   --out Reults/Educ \
   --merge-alleles Data/w_hm3.snplist

```

And this is **finally** the command to estimate genetic correlation

```
python Software/ldsc/ldsc.py \
 --rg Data/height.sumstats.gz,Data/Educ.sumstats.gz  \
 --ref-ld-chr Data/eur_w_ld_chr/ \
 --w-ld-chr Data/eur_w_ld_chr/ \
 --out Results/height_educ

```
