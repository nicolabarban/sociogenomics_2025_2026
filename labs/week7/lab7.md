# Lab week 6. Polygenic scores



## Description
In this lab we will learn:
* Calculate a "monogenic" score
* Calculate PGS using the software PRSice
* Use the calculated PGS in R

## Data


```
cd $HOME
 wget -O lab_week7.zip https://www.dropbox.com/s/kwciw2cb19gkrzy/week6.zip?dl=0
unzip lab_week7.zip 
mv week6/*.*  ./
rm -r week6/
```

## Calculate a "monogenic" score

We start from a "monogenic score" based only on 1 SNP rs9930506
the file score_rs9930506.txt  includes the weights 

```
./plink 	--bfile 1kg_hm3_qc \
 		--score score_rs9930506.txt 1 2 3 \
 		--pheno BMI_pheno.txt \
 		--out 1kg_FTOscore
```

This is our "monogenic score"

```
head 1kg_FTOscore.profile
```

## Install R in google cloud

To calculate PGS we need additional software:

1. PRSice (https://choishingwan.github.io/PRSice/)
2. R 

```

wget https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip
unzip PRSice_linux.zip
```
Install R on the virtual machine
```
sudo apt-get install r-base r-base-dev
```

let's see if it works
```
R
install.packages(c("ggplot2", "data.table"))
```

```
./PRSice_linux -h
```




Calculating a Polygenic score on BMI based on BMI.txt
```
Rscript PRSice.R --dir . \
    --prsice ./PRSice_linux \
    --base BMI.txt \
    --target 1kg_hm3_qc \
    --snp MarkerName \
    --A1 A1 \
    --A2 A2 \
    --stat Beta \
    --pvalue Pval \
    --pheno-file BMI_pheno.txt \
    --bar-levels 1 \
    --fastscore \
    --binary-target F \
    --out BMI_score_all
```
See the ouput. there are 9 duplicated SNPs that can be removed with the command extract
```
Rscript PRSice.R --dir . \
    --prsice ./PRSice_linux \
    --base BMI.txt \
    --target 1kg_hm3_qc \
    --snp MarkerName \
    --A1 A1 \
    --A2 A2 \
    --stat Beta \
    --pvalue Pval \
    --pheno-file BMI_pheno.txt \
    --bar-levels 1 \
    --fastscore \
    --binary-target F \
    --extract BMI_score_all.valid \
    --out BMI_score_all
```


Calculating different scores at various pvalues
```
Rscript PRSice.R --dir . \
    --prsice ./PRSice_linux \
    --base BMI.txt \
    --target 1kg_hm3_qc \
    --thread 1 \
    --snp MarkerName \
    --A1 A1 \
    --A2 A2 \
    --stat Beta \
    --pvalue Pval \
    --bar-levels 5e-08,5e-07,5e-06,5e-05,5e-04,5e-03,5e-02,5e-01 \
    --fastscore \
    --all-score \
    --no-regress \
    --binary-target F \
    --extract BMI_score_all.valid \
    --out BMIscore_thresholds  
```
See results
```
head BMIscore_thresholds.all_score
```


PGS on binary trait
```
Rscript PRSice.R --dir . \
    --prsice ./PRSice_linux \
    --base BMI.txt \
    --target 1kg_hm3_qc \
    --thread 1 \
    --snp MarkerName \
    --A1 A1 \
    --A2 A2 \
    --stat Beta \
    --pvalue Pval \
    --no-clump F \
    --pheno-file Obesity_pheno.txt \
    --interval 0.00005 \
    --lower 0.0001 \
    --quantile 5 \
    --all-score \
    --binary-target T \
    --extract BMI_score_all.valid \
    --out Obesity_score_graphics

```

#THIS NEED TO BE EXECUTED IN R 

```
#import external data
setwd("/home/nicola_barban/Sociogenomics/")

data<-read.table("BMIscore_thresholds.all_score", header=T)

# show first rows of data
head(data)

# calculate new standardised variable
data$PGS=(data$Pt_1-mean(data$Pt_1))/sd(data$Pt_1, na.rm=T)

#plot histogram of the polygenic score 
 hist(data$PGS)

# import external data with phenotype
pheno_BMI<-read.table("BMI_pheno.txt", header=T)

# merge the two datasets
data.with.pheno<-merge(data,pheno_BMI, by="IID")
```


```

# run a linear regression model
mod1<-(lm(BMI~PGS, data=data.with.pheno))
summary(mod1)


summary(mod1)$r.square
```


```
# create a vector with column names
columns=c("FID", "IID", "pca1", "pca2", "pca3", "pca4", "pca5", "pca6", "pca7", "pca8", "pca9", "pca10")

# read external data with PCAs
pca<- read.table(file="../Data/pca.eigenvec", sep = "", header=F, col.names=columns)[,c(2:12)]

# merge file with covariates with the rest of the file
data.with.covars<-merge(data.with.pheno,pca, by="IID")

# Estimate new linear regression model
mod2<-lm(BMI~PGS+pca1+pca2+pca3+pca3+pca5+pca6+pca7+pca8+pca9+pca10, data=data.with.covars)
```

```
# Calculate R-Square
summary(mod2)$r.square

# Estimate linear regression model
mod2.no.pgs<-update(mod2,  ~ . - PGS)

# Estimate difference in R-squared
print(summary(mod2)$r.square-summary(mod2.no.pgs)$r.square)
# install new package in R to perform bootstrap

# install new package in R to perform bootstrap

install.packages("boot")
library(boot)
set.seed(12345)

```

```

# define new function to calculate differential r-squared.
# The following commands define the new function
rsq <- function(formula, PGS=PGS, data, indices) {
  d <- data[indices,] 
  fit1 <- lm(formula,  data=d)
  fit2 <- update(fit1,  ~ . + PGS)
  
  return(summary(fit2)$r.square-summary(fit1)$r.square)
} 
```


```
# bootstrapping with 1000 replications
results <- boot(data=data.with.covars, statistic=rsq, 
  	R=1000, formula=BMI~pca1+pca2+pca3+pca3+pca5+pca6+pca7+pca8+pca9+pca10, PGS=PGS)
boot.ci(results, type="norm")      
```
