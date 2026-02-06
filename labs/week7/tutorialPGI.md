# Tutorial on PGI calculation
based on https://choishingwan.github.io/PRS-Tutorial/base/










 Download file 
 Attention very big!
 
```
wget "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90275001-GCST90276000/GCST90275127/harmonised/GCST90275127.h.tsv.gz"

gunzip -c GCST90275127.h.tsv.gz | head
````


 renaming the file into Parkinson_gwas
```
mv GCST90275127.h.tsv.gz Parkinson_gwas.tsv.gz
```
### Filtering variants with MAF<1% 


The bash code above does the following:

Decompresses and reads the Parkinson_gwas.tsv.gz file
Prints the header line (NR==1)
Prints any line with MAF above 0.05 ($7 because the eleventh column of the file contains the MAF information)
Compresses and writes the results to parkinson.gz

```
gunzip -c Parkinson_gwas.tsv.gz |\
awk 'NR==1 || ($7 > 0.05)  {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' |\
gzip  > parkinson.gz
```

# Start from here to avoid large download!!



### Removing duplicate SNPs


```
gunzip -c parkinson.gz |\
awk '{seen[$9]++; if(seen[$9]==1){ print}}' |\
gzip - > parkinson.nodup.gz
```

The above command does the following:

Decompresses and reads the parkinson.gz file
Count number of time SNP ID was observed, assuming the third column contian the SNP ID (seen[$9]++). If this is the first time seeing this SNP ID, print it.
Compresses and writes the results to parkinson.nodup.gz


### Removing ambiguous SNPS

```
gunzip -c parkinson.nodup.gz |\
awk '!( ($3=="A" && $4=="T") || \
        ($3=="T" && $4=="A") || \
        ($3=="G" && $4=="C") || \
        ($3=="C" && $4=="G")) {print}' |\
    gzip > parkinson.QC.gz
```

## Working on the target Data

selecting EUR samples
```
awk 'NR>1 && ($6=="EUR")  {print 0, $1}' 1kg_samples.txt >EUR_sample.txt

```
Quality control for Base file

```
./plink \
    --bfile 1kg_hm3 \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
	--keep EUR_sample.txt \
    --write-snplist \
    --make-just-fam \
    --out EUR.QC
```

## Calculating PGI with plink


### Clumping
Linkage disequilibrium, which corresponds to the correlation between the genotypes of genetic variants across the genome, makes identifying the contribution from causal independent genetic variants extremely challenging. One way of approximately capturing the right level of causal signal is to perform clumping, which removes SNPs in ways that only weakly correlated SNPs are retained but preferentially retaining the SNPs most associated with the phenotype under study. Clumping can be performed using the following command in plink:

```
./plink \
    --bfile 1kg_hm3 \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump parkinson.QC.gz \
    --clump-snp-field rsid \
    --clump-field p_value \
    --out EUR
	
```


Extracting from the file EUR.clumped the list of SNPs


This will generate EUR.clumped, containing the index SNPs after clumping is performed. We can extract the index SNP ID by performing the following command:
$3 because the third column contains the SNP ID


```
awk 'NR!=1{print $3}' EUR.clumped >  EUR.valid.snp
```


We read from the parkinson.QC.gz file, assuming that the 9th column is the SNP ID; 3rd column is the effective allele information; the 5th column is the effect size estimate; and that the file contains a header




```
gunzip -c parkinson.QC.gz > parkinson.QC.tsv
head parkinson.QC.tsv


./plink \
    --bfile 1kg_hm3 \
    --score parkinson.QC.tsv 9 3 5 header \
    --extract EUR.valid.snp \
    --out EUR
	
head EUR.profile 	
	
```

# PGI with PRSice

Calculate scores using all SNPs

```

Rscript PRSice.R \
    --prsice PRSice_mac \
    --base parkinson.QC.gz \
    --target 1kg_hm3 \
    --base-maf effect_allele_frequency:0.05 \
    --thread 1 \
    --snp rsid \
    --A1 effect_allele \
    --A2 other_allele \
    --stat beta \
    --pvalue p_value \
    --fastscore \
    --bar-levels 1 \
    --no-regress \
    --binary-target F \
    --extract EUR.valid.snp \
    --out Parkinson_score_all
```


here using different threesholds

```

Rscript PRSice.R \
    --prsice PRSice_mac \
    --base parkinson.QC.gz \
    --target 1kg_hm3 \
    --base-maf effect_allele_frequency:0.05 \
    --thread 1 \
    --snp rsid \
    --A1 effect_allele \
    --A2 other_allele \
    --stat beta \
    --pvalue p_value \
	--bar-levels 5e-08,5e-07,5e-06,5e-05,5e-04,5e-03,5e-02,5e-01 \
    --fastscore \
    --all-score \
    --no-regress \
    --binary-target F \
    --extract EUR.valid.snp \
    --out Parkinson_score_thresholds
```



in R

```

data_PRSice<-read.table("Parkinson_score_thresholds.all_score", header=T)
head(data_PRSice)


data_Plink<-read.table("EUR.profile", header=T)


data_merge<-merge(data_PRSice,data_Plink, by=c("FID", "IID"))
 cor(data_merge$SCORE, data_merge$Pt_1)
 
 cor(data_merge$SCORE, data_merge$Pt_5e.08)
 
 
 
```