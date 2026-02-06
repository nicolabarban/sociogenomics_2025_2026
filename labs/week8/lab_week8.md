# Lab week 8. Working with PGS in R





```
#import external data
setwd("/home/nicola_barban/Sociogenomics/Results")

data<-read.table("BMIscore_thresholds.all_score", header=T)

# show first rows of data
head(data)

# calculate new standardised variable
data$PGS=(data$Pt_1-mean(data$Pt_1))/sd(data$Pt_1, na.rm=T)

#plot histogram of the polygenic score 
 hist(data$PGS)

# import external data with phenotype
pheno_BMI<-read.table("../Data/BMI_pheno.txt", header=T)

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
