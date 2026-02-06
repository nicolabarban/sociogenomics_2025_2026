# Lab Week 4. Sociogenomics

[Google Shell](https://cloud.google.com/shell/docs/launching-cloud-shell?hl=en)


## Description
In this lab we will learn:
* Calculate IBS and relatedness
* Association analys

##  Managing files and directories. 
Let's have a look a the file. 

For this tutorial we will use data from hapmap and 1000 Genome. 
Download the data using this command.

```
wget https://www.nicolabarban.com/sociogenomics_lab/data/hapmap_CEU.bed --no-check-certificate
wget https://www.nicolabarban.com/sociogenomics_lab/data/hapmap_CEU.bim --no-check-certificate
wget https://www.nicolabarban.com/sociogenomics_lab/data/hapmap_CEU.fam --no-check-certificate

```
YOu can download data from 1000Genome from last's week tutorial

```
cd $HOME

wget -O week3.zip https://www.dropbox.com/scl/fi/kvsdtvsl3m4gl19omle1y/week3.zip?rlkey=3fyj402e77jsvo97iwz8ke7sc&e=1&st=k1x60x1z&dl=0
unzip  week3.zip 
mv week3/*.*  ./
rm -r week3/
rm -r __MACOSX
```


**PLINK QC**. combine different commands in one go
```


./plink     --bfile hapmap_CEU \
			--filter-founders \
			--autosome \
			--snps-only \
       	  	--mind 0.03 \
       		--geno 0.05 \
       		--maf 0.05 \
    		--hwe 0.00001 \
        	--make-bed  --out hapmap_CEU_QC     
			
```
## Select independent SNPS

Calculate linkage disequilibrium
```
./plink --bfile 1kg_hm3 \
	 	--ld rs1048488 rs3115850 \
		--out ld_example
```
Check also this: https://ldlink.nih.gov/?tab=home


## Calculate independent SNPs (Pruning)
```

./plink 	 --bfile 1kg_hm3  \
        	--indep-pairwise 50 5 0.2 \
        	--out  1kg_hm3
			
			
```


```
./plink --bfile 1kg_hm3 \
	 	--ld rs1048488 rs2519031 \
		--out ld_example2
```

Select from original sample independent SNPs
```
./plink		--bfile  1kg_hm3 \
			--extract 1kg_hm3.prune.in \
			--make-bed \
 			--out  1kg_hm3_indep

```


## Calculate IBS and relatedness

Calculate Identity By State matrix
```
./plink --bfile  1kg_hm3_indep \
		--keep 1kg_samples_EUR.txt \
		--distance --out ibs_matrix
```
Calculate relatedness matrix

```
./plink --bfile 1kg_hm3_indep --keep 1kg_samples_EUR.txt --make-rel --out rel_matrix

```

		
		
## GCTA

```
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.3-linux-kernel-3-x86_64.zip
unzip gcta-1.94.3-linux-kernel-3-x86_64.zip
mv gcta-1.94.3-linux-kernel-3-x86_64/gcta64 ./
rm -r gcta-1.94.3-linux-kernel-3-x86_64/
./gcta64

```

```
./gcta64 --bfile 1kg_hm3 		--keep 1kg_samples_EUR.txt 		 --make-grm --out 1kg_hm3_allSNPs

```

### Getting Phenotype

```
wget https://www.nicolabarban.com/sociogenomics_lab/data/height_sim2.phen --no-check-certificate
head height_sim2.phen


```


```
./gcta64  --reml  --grm 1kg_hm3_allSNPs	--pheno height_sim2.phen --grm-adj 0 --grm-cutoff 0.05 --out height_h2
head height_h2.hsq
```


## Association analys

### Linear additive model
```
./plink    	--bfile 1kg_hm3 \
			--pheno height_sim2.phen \
        	--snps rs9674439 \
       	 	--assoc \
      	 	--linear \
      		--out height_rs9674439
```


```
head height_rs9674439.assoc.linear
```


### Linear dominant analysis
```
./plink    	 --bfile 1kg_hm3 \
			--pheno height_sim2.phen \
        	 --snps rs9674439 \
       	 	--assoc \
      	 	--linear dominant \
      	 	--out height_rs9674439_dom
```		 
	
### All vaariants, a.k.a GWAS
```
./plink    	--bfile 1kg_hm3 \
			--pheno height_sim2.phen \
       	 	--assoc \
      	 	--linear \
      		--out height_gwas
```			