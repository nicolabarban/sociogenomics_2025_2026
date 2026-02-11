# Lab 1. Introduction to the Unix Shell

In this lab you will learn the basics of working in a Linux terminal using **Google Cloud Shell**. By the end you will be able to navigate the file system, organise files into directories, inspect genetic data files, and use `awk` for data manipulation.

---

## 0. Getting started — Google Cloud Shell

Open [Google Cloud Shell](https://shell.cloud.google.com/) in your browser. You get a free Linux virtual machine with a persistent home directory.

### Clone the course repository

```
git clone https://github.com/nicolabarban/sociogenomics_2025_2026.git
```

This downloads all lab materials and data into `~/sociogenomics_2025_2026/`.

### Uploading and downloading files

Cloud Shell has a built-in file transfer feature. Click the **three-dot menu (⋮)** in the top-right corner of the terminal:

- **Upload file** — transfer a file from your computer to the current directory.
- **Download file** — enter the path to a file (e.g. `Sociogenomics/Results/output.txt`) to save it to your computer.

You can also click **Open Editor** to browse and edit files visually.

---

## Part I. Navigating the file system

### Where am I?

```
pwd
```

This prints your **working directory** (usually `/home/your_username`).

### Listing files

```
ls
ls -l
ls -lh
```

`-l` shows a detailed list (permissions, size, date). `-h` makes sizes human-readable.

### Creating directories

Create a project directory with subdirectories for organising your work:

```
mkdir -p Sociogenomics/Data Sociogenomics/Results Sociogenomics/Scripts
```

Verify the structure:

```
ls Sociogenomics
```

### Navigating directories

```
cd Sociogenomics
pwd
cd Data
pwd
cd ..
cd ~
```

`cd ..` moves up one level. `cd ~` (or just `cd`) takes you back to your home directory.

### Copying files

Copy the HapMap data from the course repository into your project:

```
cp ~/sociogenomics_2025_2026/data/hapmap1.map ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/hapmap1.ped ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/BMI_pheno.txt ~/Sociogenomics/Data/
```

### Moving and renaming files

```
cp ~/Sociogenomics/Data/hapmap1.map ~/Sociogenomics/Data/hapmap1_backup.map
mv ~/Sociogenomics/Data/hapmap1_backup.map ~/Sociogenomics/Results/
ls ~/Sociogenomics/Results/
```

`mv` moves or renames a file.

### Removing files

```
rm ~/Sociogenomics/Results/hapmap1_backup.map
```

Use `rm -r` to remove a directory and its contents (be careful!).

### Exercise 1

1. Create a directory called `Sociogenomics/Temp`.
2. Copy `hapmap1.ped` into `Temp/`.
3. Rename it to `test.ped` using `mv`.
4. Remove the `Temp/` directory and everything inside it.

---

## Part II. Exploring the HapMap data

The HapMap project catalogued common genetic variants across human populations. We have two files:

- **`hapmap1.map`** — one row per SNP with 4 columns: chromosome, SNP ID, genetic distance, base-pair position.
- **`hapmap1.ped`** — one row per individual with 6 fixed columns (family ID, individual ID, father ID, mother ID, sex, phenotype) followed by two allele columns per SNP.

Navigate to the Data directory:

```
cd ~/Sociogenomics/Data
```

### Looking at files

View the first 10 lines of the map file:

```
head hapmap1.map
```

View the first 5 lines:

```
head -5 hapmap1.map
```

View the last 5 lines:

```
tail -5 hapmap1.map
```

The ped file is very wide, so `head` will only show the beginning of each row:

```
head -3 hapmap1.ped
```

### Counting lines

```
wc -l hapmap1.map
wc -l hapmap1.ped
```

**Question:** How many SNPs are in the dataset? How many individuals?

### Searching with grep

Search for a specific SNP by its rs ID:

```
grep rs7540009 hapmap1.map
```

Search for all SNPs whose ID starts with `rs75`:

```
grep rs75 hapmap1.map
```

Count the matches:

```
grep -c rs75 hapmap1.map
```

### Searching on a specific chromosome

Find all SNPs on chromosome 22:

```
grep "^22" hapmap1.map | head
```

Count them:

```
grep -c "^22" hapmap1.map
```

`^` means "start of line", so `^22` matches rows where the chromosome column is 22.

### Combining commands with pipes

Pipes (`|`) send the output of one command into another:

```
grep "^1	" hapmap1.map | wc -l
```

This counts SNPs on chromosome 1. The tab character after `1` avoids matching chromosomes 10–19.

### Exercise 2

1. How many SNPs are on chromosome 6?
2. Find the SNP `rs4558854` — which chromosome is it on?
3. How many SNPs have IDs starting with `rs10`?

---

## Part III. Data manipulation with AWK

`awk` is a powerful tool for processing column-based text files. It reads a file line by line and lets you select, filter, and transform columns.

### Basic syntax

```
awk '{ action }' filename
```

### Printing columns

The `.map` file has 4 columns: `$1` (chromosome), `$2` (SNP ID), `$3` (genetic distance), `$4` (position).

```
awk '{ print $1, $2 }' hapmap1.map | head
```

Print with custom labels:

```
awk '{ print "chr"$1, $2, "pos:"$4 }' hapmap1.map | head
```

### Built-in variables

- **NR** — current row number
- **NF** — number of fields in the current row
- **$0** — the entire line

```
awk '{ print NR, $1, $2 }' hapmap1.map | head
```

How many columns does the ped file have?

```
awk '{ print NF; exit }' hapmap1.ped
```

**Question:** The `.ped` file has 6 fixed columns followed by two allele columns per SNP. Given the number of fields, how many SNPs does this confirm?

### Filtering rows

Print only SNPs on chromosome 2:

```
awk '$1 == 2 { print $1, $2, $4 }' hapmap1.map | head
```

Print rows 100 to 105:

```
awk 'NR >= 100 && NR <= 105' hapmap1.map
```

### Pattern matching

Find a SNP by exact ID:

```
awk '$2 == "rs7540009" { print $0 }' hapmap1.map
```

Find all SNPs matching a pattern:

```
awk '$2 ~ /rs75/ { print $1, $2, $4 }' hapmap1.map
```

### Counting SNPs per chromosome

```
awk '{ count[$1]++ } END { for (chr in count) print "chr"chr, count[chr] }' hapmap1.map | sort -n -k1.4
```

### Working with the BMI phenotype file

```
head BMI_pheno.txt
```

This file has 3 tab-separated columns: FID, IID, BMI.

Print the first 5 individuals:

```
awk 'NR <= 6 { print $1, $2, $3 }' BMI_pheno.txt
```

Find individuals with BMI above 30:

```
awk 'NR > 1 && $3 > 30 { print $1, $2, $3 }' BMI_pheno.txt | head
```

Count them:

```
awk 'NR > 1 && $3 > 30' BMI_pheno.txt | wc -l
```

### Changing delimiters

Use tab as output separator:

```
awk 'BEGIN { OFS="\t" } { print $1, $2, $4 }' hapmap1.map | head
```

### Redirecting output to a file

Save all chromosome 1 SNPs to a new file:

```
awk '$1 == 1 { print $2, $4 }' hapmap1.map > ../Results/chr1_snps.txt
head ../Results/chr1_snps.txt
wc -l ../Results/chr1_snps.txt
```

You can download this file to your computer using the Cloud Shell menu: **⋮ → Download file**, then type `Sociogenomics/Results/chr1_snps.txt`.

### Counting lines with AWK

```
awk 'END { print NR }' hapmap1.map
```

### Exercise 3

1. Use `awk` to count how many SNPs are on each of chromosomes 1, 2, and 22.
2. Extract all SNPs on chromosome 6 with their positions and save to `Results/chr6_snps.txt`.
3. How many individuals in `BMI_pheno.txt` have a BMI below 20?
4. Compute the average BMI across all individuals. Hint: accumulate a sum and count in the main block, then print the result in the `END` block.
