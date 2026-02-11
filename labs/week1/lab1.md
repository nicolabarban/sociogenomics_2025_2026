# Lab 1.  Introducing the Unix shell

## Description

## Part I . Managing files and directories.


Create directories
```
mkdir data
mkdir labs
mkdir labs/week1
ls
ls -F
```

Upload files (upload them  if using Cloud or using GUI)

Navigate  directories
```
cd data
ls
cd ..
cd $HOME
```

Copying files and directories
```
cd .. 
cp data/hapmap1.map ../
ls 


cp -r labs data/
cd data
ls
```


Moving files and directories
```
mv data/hapmap1.map ../
ls 



```

Remove files and directories
```
rm data/hapmap1.map 
ls 


cd ..
rm -r data
rm -r labs

```

## Exercise 1

* Create one directory called Sociogenomics in your home
* Create the following subdirectories:		
	* Data
	* Scripts
	* Results
	* Software
* Move the files from Virtuale into the folders
* Remove the files used until now

---

## Looking at files

look at the first 6 rows
```
cd Sociogenomics/Data
head hapmap1.map 

```


look at the first x rows
```
head -10 hapmap1.map 

```

Count number of lines
```
wc -l hapmap1.map  
wc -l hapmap1.ped
```


Echo  allows you to send text into the terminal  the symbol (\) "backslash" allows you to write your command in several 

```
echo 'Hello' \
 'world'
```

It is possible to combine several commands by using pipilines
```
echo 'Hello world' | wc

```

Grep allows you to search for a pattern in a file

```
grep rs7540009 hapmap1.map

grep rs75 hapmap1.map

grep rs75 hapmap1.map | wc -l


```


Download GWAS summary statistics from UK Biobank
http://www.nealelab.is/uk-biobank

```
curl -L -o sumstatsUKB_height.tsv.gz "https://www.dropbox.com/s/ou12jm89v74k55e/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz?dl=1"

gunzip -d sumstatsUKB_height.tsv.gz
```

Have a look at the data
```
head sumstatsUKB_height.tsv
wc sumstatsUKB_height.tsv

more sumstatsUKB_height.tsv

```


## AWK
Awk is a scripting language used for manipulating data and generating reports. The awk command programming language requires no compiling and allows the user to use variables, numeric functions, string functions, and logical operators. 

Awk is a utility that enables a programmer to write tiny but effective programs in the form of statements that define text patterns that are to be searched for in each line of a document and the action that is to be taken when a match is found within a line. Awk is mostly used for pattern scanning and processing. It searches one or more files to see if they contain lines that matches with the specified patterns and then perform the associated actions. 

Awk is abbreviated from the names of the developers – Aho, Weinberger, and Kernighan. 

```
cat BMI_pheno.txt
awk '{print}' BMI_pheno.txt


```

```
awk '{print $1, $2  }'  BMI_pheno.txt

```



find a pattern on the data
```
awk ' /3:49860854/ {print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv
awk ' {print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv | grep 3:49860854

```

Built-In Variables In Awk

Awk’s built-in variables include the field variables—$1, $2, $3, and so on ($0 is the entire line) — that break a line of text into individual words or pieces called fields. 

NR: NR command keeps a current count of the number of input records. Remember that records are usually lines. Awk command performs the pattern/action statements once for each record in a file. 
NF: NF command keeps a count of the number of fields within the current input record. 
FS: FS command contains the field separator character which is used to divide fields on the input line. The default is “white space”, meaning space and tab characters. FS can be reassigned to another character (typically in BEGIN) to change the field separator. 

OFS: OFS command stores the output field separator, which separates the fields when Awk prints them. The default is a blank space. Whenever print has several parameters separated with commas, it will print the value of OFS in between each parameter. 



Print Number of row, and column number 1 and 2
```
awk '{ print  NR, $1, $2  }'  BMI_pheno.txt

```

Print number of fields
```
awk '{print $1,$NF}' BMI_pheno.txt
```

Print case 3-6
```
awk 'NR==3, NR==6 {print NR,$0}'  BMI_pheno.txt
```

Print symbols
```
awk '{print NR "- " $1 }' BMI_pheno.txt

```

Print specific rows
```
awk '{ if(NR==1) print $1, $2  }'  BMI_pheno.txt

awk '{ if(NR<10) print $1, $2  }'  BMI_pheno.txt

```


Print first 10 rows
```

awk '{ if(NR<10) print   }'  sumstatsUKB_height.tsv

```

print only rows with p-value>5e-08
```
awk '{ if($11<5e-08) print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv

```

transforming columns
```
awk '{ if($11<5e-08) print  $1, $2, $5, log($11) }'  sumstatsUKB_height.tsv

```


 To count the lines in a file:  
```
$ awk 'END { print NR }' sumstatsUKB_height.tsv 
```

Change delimiter character
```
awk 'FS=":" {print $1, $2, $3}' sumstatsUKB_height.tsv | head
awk 'FS=":", OFS="-" {print $1, $2, $3}' sumstatsUKB_height.tsv | head
```

multiple field separators

```
awk -F '[:"\t"]' '{ if(NR>1) print $1, $2, $3 , $4, $13}' sumstatsUKB_height.tsv | head
```

Redirecting output to new file
```
awk '{ if($11<5e-08) print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv > ../Results/sign_variants_UKB.txt
```

Selecting SNPs with MAF  greater than 10%
```
awk '{ if($3>0.1) print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv > ../Results/common_variants_UKB.txt

head  ../Results/common_variants_UKB.txt
```

Counting SNPS with MAF>10% pvalue>5e-08 in Chromosome 1
```
awk -F '[:"\t"]' '{ if($6 >0.1 && $1==1) print $1, $2, $3, $4,  $6}' sumstatsUKB_height.tsv | wc -l

```
Counting SNPS with MAF>10% pvalue>5e-08 in Chromosome 21

```
awk -F '[:"\t"]' '{ if($6 >0.1 && $1==21) print $1, $2, $3, $4,  $6}' sumstatsUKB_height.tsv | wc -l

```
