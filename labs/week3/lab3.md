# Lab Week 3. Sociogenomics

[Google Shell](https://cloud.google.com/shell/docs/launching-cloud-shell?hl=en)

## Description
In this lab we will learn:

* Compute descripitve statistics with PLINK on genetic files 
* Quality control procedure with plink




### Here are some basic Linux commands:

* `pwd`: Displays current directory
* `mkdir`: Creates directory echo: copies the command line to the screen 
* `ls` : lists the names of files in the current directory.
* `cat`: displays a text file
* `rm`: deletes a file
* `cp`: copies a file
* `mv`: changes the name of a file
* `grep`: searches for a string in a file
* `wc`:count the number of words/lines in a file
* `head`: displays the beginning of a file tail: displays the end of a file
* `sort`: displays the lines of a file in order
* `uniq`: Removes adjacent duplicate lines (file has to be sorted!)
* `gzip`, `gunzip`: Compress, decompress files with .gz extension.


### Remove all files from your directory


```
rm *.*
``` 

### Download the files we need for today's lecture
```
cd $HOME

wget -O week3.zip https://www.dropbox.com/scl/fi/kvsdtvsl3m4gl19omle1y/week3.zip?rlkey=3fyj402e77jsvo97iwz8ke7sc&e=1&st=k1x60x1z&dl=0
unzip -o week3.zip 
mv week3/*.*  ./
rm -r week3/
```

##  Managing files and directories.
Let's have a look a the file. 

```
head   1kg_hm3.fam
head   1kg_hm3.bim

```


```
 wc -l 1kg_hm3.fam
```

how many variants?

```
 wc -l 1kg_hm3.bim
```


##  Descriptive Statistics



### Allele frequency
We can calculate allele frequency

```
 ./plink --bfile 1kg_hm3  --freq --out 1kgAllele_Frequency
```

```
head 1kgAllele_Frequency.frq 
```


### Linkage disequilibrium
Calculate linkage disequilibrium
```
./plink --bfile 1kg_hm3 \
	 	--ld rs1048488 rs3115850 \
		--out ld_example
```

### Select individuals


```
./plink --bfile 1kg_hm3 \
	 	--keep 1kg_samples_EUR.txt \
		--make-bed \
		--out 1kg_hm3_EUR
```

```
./plink --bfile 1kg_hm3 \
	 	--remove 1kg_samples_EUR.txt \
		--make-bed \
		--out 1kg_hm3_NOT_EUR
```

Calculate linkage disequilibrium
```
./plink --bfile 1kg_hm3_EUR \
	 	--ld rs1048488 rs3115850 \
		--out ld_example_EUR
```

Calculate linkage disequilibrium
```
./plink --bfile 1kg_hm3_NOT_EUR \
	 	--ld rs1048488 rs3115850 \
		--out ld_example_NOT_EUR
```


Filter females
```

./plink     --bfile 1kg_hm3 \
            --filter-females \
            --make-bed \
       	 	--out 1kg_hm3_filter_females

```

## Missing values

### individuals
```


./plink --bfile 1kg_hm3 --missing --out 1kg_missing_data
```
### Sample missing data report

```
head 1kg_missing_data.imiss
```
### variant-based missing data report

```
1kg_missing_data.lmiss
```


### Select individuals with genotype at least 95% complete
We can select individuals based on the completness of their genotype
```
./plink --bfile 1kg_hm3 --make-bed --mind 0.05 --out 1kg_highgeno
```



### Add a phenotype into PLINK files

PLINK file can also store info on a phenotype


```

head 1kg_hm3.fam
```


This file contains info on BMI of different individuals
```

head BMI_pheno.txt
```

This is how we add a phenotipic information to a plink file
```


./plink      --bfile 1kg_hm3 \
             --pheno BMI_pheno.txt \
             --make-bed --out 1kg_BMI

```

```
 head 1kg_BMI.fam
```



## Quality control
```


./plink --bfile 1kg_BMI \
		--mind 0.05 \
		--make-bed \
		--out 1kg_hm3_mind005
```

Calculate heterozygocity
```

./plink --bfile 1kg_BMI \
		--het --out 1kg_BMI_het
```


Low call-rate SNPS
```

./plink --bfile 1kg_BMI \
		--geno 0.05 \
		--make-bed \
		--out 1kg_hm3_geno
```

Allele frequency
```

./plink --bfile 1kg_BMI \
	 	--maf 0.01 \
		--make-bed  --out 1kg_hm3_maf
```
deviation from HWE
```

./plink --bfile 1kg_BMI \
	 	--hwe 0.00001 \
		--make-bed  --out 1kg_hm3_hwe

```

## Combine different commands in one go

**PLINK QC**. combine different commands in one go
```


./plink     --bfile 1kg_BMI \
       	--mind 0.03 \
       	--geno 0.05 \
       	--maf 0.01 \
    	--hwe 0.00001 \
        --make-bed  --out 1kg_hm3_QC      
			
```
