# Where to find summary statistics

1. [The PGC (Psychiatric Genomics Consortium)](https://www.med.unc.edu/pgc/results-and-downloads), has analyzed common psychiatric disorders (MDD, Schizophrenia, ADHD, OCD, Bipolar Disorder and more)

2. [The SSGAC (Social Sciences Genetic Association Consortium)](https://www.thessgac.org/data) performs genome wide association studies of a variety of social and psychological traits like education, personality, and reproductive behavior.

3. [The Nealelab](http://www.nealelab.is/uk-biobank) quickly ran and published online GWAS of >4000 traits that were measured as part of the UK Biobank. These traits include many disease (ICD-10 diagnostic codes, both self reported and based on hospital data), social traits (e.g. social deprivation), personality traits (e.g. neuroticism), cognition (e.g. memory) and many more (from snoring to the propensity to drive to fast). The Nealelab ran these GWAS very quickly and as a service to the field. Their GWAS of case/control traits use linear regression (linear probability model). Please read their extensive read me which describes their GWAS analysis in detail.

4. [The CCACE (Centre for Cognitive Ageing and Cognitive Epidemiology)](http://www.ccace.ed.ac.uk/node/335) has published GWAS on assorted personality traits, cognitive traits, and tiredness.


5. [The GPC (Genetics of Personality Consortium)](http://www.tweelingenregister.org/GPC/) published several, slightly dated, GWAS on the "Big 5" personality scales.

6. [The EGG (Early Growth Genetics)](https://egg-consortium.org/) Consortium performs GWAS of traits related to early growth.

7. [The GIANT consortium](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files) publishes GWAS, mainly about antropomorpic traits.

8. [GWAS Catalog](https://www.ebi.ac.uk/gwas/)
9. [Biobank Japan](https://pheweb.jp/)
10. [FInnGEn](https://www.finngen.fi/en/access_results)
11. [Pan UKKB](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=268241601) Pan Ancestry GWAS analysis of UK Biobank
12. 




# calculate heritability and genetic correlation using ldsc


This code provide examples on how to calculate h2 and rg using R


First, we need to install the library GenomicSEM


```
install.packages("devtools")
install.packages("vctrs")

library(devtools)
install_github("GenomicSEM/GenomicSEM", force=TRUE, dependencies = TRUE)


require(GenomicSEM)
```

We will base the analysis on European data.
The first example will use data from summary statistics of BMI and Height (only on chromosome 1)


```
files<-c("panUKB_BMI_EUR_chr1.txt", "Yengo_Height_EUR_chr1.txt")
```



Let's define the reference file being used to allign alleles across summary stats
 using hapmap3

```
hm3<-"eur_w_ld_chr/w_hm3.snplist"

```

Let's give a name to the traits
```
trait.names<-c("BMI","Height")
```


Let's list the sample sizes. Since Yengo has sample size in the file but pan UKB does not,  we provide that directly

```
N<-c(419163, NA)
```


Let's define the imputation quality filter and MAF filter 

```
info.filter=0.9
maf.filter=0.01
```

Let's run munge
```
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)
```


###Step 2: run LDSC


* traits = the name of the .sumstats.gz traits produced by munge
* ld = folder of LD scores 
* wld = folder of LD scores
* sample.prev = the proportion of cases to total sample size. For quantitative traits list NA
* population.prev = the population lifetime prevalence of the traits. For quantitative traits list NA
* trait.names = optional fifth argument to list trait names so they can be named in your model



```
traits<-c("BMI.sumstats.gz","Height.sumstats.gz")

ld <- "eur_w_ld_chr/"

wld <- "eur_w_ld_chr/"

sample.prev<-c(NA,NA)

population.prev<-c(NA,NA)

trait.names<-c("BMI", "Height")

```

Let's run LDSC

```
LDSC_EUR <- ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev, ld=ld, wld=wld,trait.names=trait.names)
save(LDSC_EUR, file = "LDSC_EUR.RData")
```

Some output from LDSC:

* genetic covariance matrix
* genetic correlation matrix

```
LDSC_EUR$S
cov2cor(LDSC_EUR$S)
```


Let's grab the standard errors from the V matrix

```
k<-nrow(LDSC_EUR$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSC_EUR$V))
```

matrix of Z-stats

```
Z<-LDSC_EUR$S/SE
```

matrix of p-values
```
P<-2*pnorm(Z,lower.tail=FALSE)
```

create a genetic heatmap 


```
require(corrplot)
rownames(LDSC_EUR$S)<-colnames(LDSC_EUR$S)

corrplot(corr = cov2cor(LDSC_EUR$S),
         method = "color",
         addCoef.col = "dark grey",
         add = F,
         bg = "white",
         diag = T,
         outline = T,
         mar = c(0,0,2,0),
         number.cex=2,
         cl.pos = "b",
         cl.ratio = 0.125,
         cl.align.text = "l",
         cl.offset = 0.2,
         tl.srt=45,
         tl.pos = "lt",
         tl.offset=0.2,
         tl.col = "black",
         pch.col = "white",
         addgrid.col = "black",
         xpd = T,
         tl.cex=1.3,
         is.corr=TRUE,
         title = "European")
```


