# Lab Week 5. Sociogenomics

[Google Shell](https://cloud.google.com/shell/docs/launching-cloud-shell?hl=en)


## Description
In this lab we will learn:
* Calculate IBS and relatedness
* Calculate and visualize PCA
* Quality control 
* GWAS with covariates


YOu can download data from 1000Genome from week3 tutorial

```
cd $HOME

wget -O week3.zip https://www.dropbox.com/scl/fi/kvsdtvsl3m4gl19omle1y/week3.zip?rlkey=3fyj402e77jsvo97iwz8ke7sc&e=1&st=k1x60x1z&dl=0
unzip  week3.zip 
mv week3/*.*  ./
rm -r week3/
rm -r __MACOSX
```


## **Step 1: Compute IBD Using PLINK**
PLINK estimates the proportion of IBD alleles between pairs of individuals. To do this, use the `--genome` flag:

```
./plink --bfile 1kg_hm3 --genome --out ibd_results

```


This will generate an output file (ibd_results.genome) with pairwise relatedness estimates.

## **Step 2: Understanding the IBD Output**
The output file (`ibd_results.genome`) contains several columns. The most relevant for identifying relatives are:

- **PI_HAT**: Estimated proportion of alleles shared identical by descent (IBD).
- **Z0**: Probability that the pair shares **zero** alleles IBD.
- **Z1**: Probability that the pair shares **one** allele IBD.
- **Z2**: Probability that the pair shares **two** alleles IBD.
- **IBS0**: Number of SNPs where both individuals have different alleles (used to detect sample swaps).

### **Interpreting PI_HAT Values**
The **PI_HAT** value helps classify relationships between individuals. Use the following thresholds as a guide:

| **PI_HAT Value** | **Relationship** |
|------------------|------------------------------|
| > 0.9           | Duplicate samples / Identical twins |
| 0.4 - 0.6       | Full siblings |
| 0.2 - 0.4       | Half-siblings / Grandparent-grandchild |
| 0.1 - 0.2       | First cousins |
| < 0.1           | Unrelated individuals |

If PI_HAT is **greater than 0.9**, the individuals are likely **duplicates** or **identical twins**. If it falls between **0.4 and 0.6**, they are likely **full siblings**.

You can also check **Z0, Z1, and Z2** values:
- **Z2 ≈ 1** → The individuals are likely **monozygotic twins**.
- **Z1 ≈ 1** → Indicates **first-degree relatives** (e.g., siblings or parent-child).
- **Z0 ≈ 1** → The individuals are likely **unrelated**.


```
head ibd_results.genome 

```
# Identifying Relatives in a PLINK File Using `awk`

After computing **identity by descent (IBD)** with PLINK, you can use `awk` to filter relationships based on **PI_HAT** values from the `ibd_results.genome` file.

## **Filtering Relatives Using `awk`**

### **1. Find Possible Duplicates (PI_HAT > 0.9)**


```bash
awk '$10 > 0.9' ibd_results.genome
```


Identifies individuals who share more than 90% of their genome, likely duplicate samples or identical twins.

### **2.Find Full Siblings (0.4 ≤ PI_HAT ≤ 0.6)**


```bash
awk '$10 >= 0.4 && $10 <= 0.6' ibd_results.genome
```



### **3. Find Half-Siblings or Grandparent-Grandchild (0.2 ≤ PI_HAT ≤ 0.4)**


```bash
awk '$10 >= 0.2 && $10 <= 0.4' ibd_results.genome
```

### **4. Find First Cousins (0.1 ≤ PI_HAT ≤ 0.2)**


```bash
awk '$10 >= 0.1 && $10 <= 0.2' ibd_results.genome

```

### **5. Find Unrelated Individuals (PI_HAT < 0.1)**


```bash
awk '$10 < 0.1' ibd_results.genome

```
### **How many relatives**

we can use the command `wc -l`

```bash
awk '$10 > 0.9' ibd_results.genome | wc -l
awk '$10 >= 0.4 && $10 <= 0.6' ibd_results.genome | wc -l
awk '$10 >= 0.2 && $10 <= 0.4' ibd_results.genome | wc -l
awk '$10 >= 0.1 && $10 <= 0.2' ibd_results.genome | wc -l
awk '$10 < 0.1' ibd_results.genome | wc -l
```

## **Remove relatives Using PLINK’s `--rel-cutoff` Option**

PLINK provides a built-in option to **prune related individuals** based on their **PI_HAT** values. This ensures that only **unrelated individuals** remain in the dataset.

### **Command to Remove Related Individuals**
```bash
./plink --bfile 1kg_hm3 --rel-cutoff 0.1 --make-bed --out unrelated_samples

```
## Compare Allele Frequencies with a Reference Dataset

First, compute allele frequencies for both datasets using PLINK:

```bash
./plink --bfile 1kg_hm3  --freq --out 1kg_freq
./plink --bfile hapmap_CEU_QC  --freq --out reference_freq
```
Merge and visualize in R:

```R
library(ggplot2)
library(data.table)


freq_sample <- fread("1kg_freq.frq")
freq_ref <- fread("reference_freq.frq")

colnames(freq_sample) <- colnames(freq_ref) <- c("CHR", "SNP", "A1", "A2", "MAF", "N")

data_merged <- merge(freq_sample, freq_ref, by = "SNP", suffixes = c("_sample", "_ref"))

ggplot(data_merged, aes(x = MAF_sample, y = MAF_ref)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Allele Frequency Comparison", x = "Sample MAF", y = "Reference MAF") +
  theme_minimal()


```



## Compute Principal Components Analysis (PCA)

To detect population structure, compute PCA:


```bash
./plink --bfile  1kg_hm3 \
		--pca 10 \
		--out 1kg_pca
```	  

```R
library(ggplot2)
columns=c("fid", "Sample.name", "pca1", "pca2", "pca3", "pca4",   "pca5", "pca6", "pca7", "pca8", "pca9", "pca10")

pca<- read.table(file="1kg_pca.eigenvec", sep = "", header=F, col.names=columns)[,c(2:12)]

ggplot(pca, aes(x=pca1, y=pca2))+ geom_point()+ theme_bw()+  xlab("PC1") + ylab("PC2")
```	  

```R
geo <- read.table(file="1kg_samples.txt", sep = "\t",header=T)[,c(1,4,5,6,7)]


data <- merge(geo, pca, by="Sample.name")
ggplot(data, aes(x=pca1, y=pca2, col=Superpopulation.name))+
     	 geom_point()+ theme_bw()+
     	 xlab("PC1") + ylab("PC2")+
 	 labs(col = "")
```	 	 
	 
```R

ggplot(data, aes(x=pca1, y=pca2, col=Population.name))+
     	 geom_point()+ theme_bw()+
     	 xlab("PC1") + ylab("PC2")+
 	 labs(col = "")
```	 	 	 

## Run a GWAS after QC and adding covariates
**PLINK QC**. combine different commands in one go
```


./plink     --bfile 1kg_hm3 \
			--autosome \
			--snps-only \
       	  	--mind 0.03 \
       		--geno 0.05 \
       		--maf 0.05 \
    		--hwe 1e-06 \
			--rel-cutoff 0.1 \
		 	--keep 1kg_samples_EUR.txt \
        	--make-bed  --out 1kg_hm3_QC_CEU    
			
```
### Run a  GWAS
```
./plink    	--bfile 1kg_hm3_QC_CEU \
			--pheno height_sim2.phen \
       	 	--assoc \
      	 	--linear \
			--covar 1kg_pca.eigenvec \
      		--out height_gwas_covar
```		 