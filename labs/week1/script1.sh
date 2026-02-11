#!/bin/bash
# =============================================================================
# Lab 1 - Introduction to the Unix Shell
# Sociogenomics 2025/2026
# Compatible with Google Cloud Shell
# =============================================================================

# --- 0. Getting started ---

git clone https://github.com/nicolabarban/sociogenomics_2025_2026.git

# --- Part I. Navigating the file system ---

pwd

ls
ls -l
ls -lh

mkdir -p Sociogenomics/Data Sociogenomics/Results Sociogenomics/Scripts

ls Sociogenomics

cd Sociogenomics
pwd
cd Data
pwd
cd ..
cd ~

# Copy data from the course repo
cp ~/sociogenomics_2025_2026/data/hapmap1.map ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/hapmap1.ped ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/BMI_pheno.txt ~/Sociogenomics/Data/

# Moving and renaming
cp ~/Sociogenomics/Data/hapmap1.map ~/Sociogenomics/Data/hapmap1_backup.map
mv ~/Sociogenomics/Data/hapmap1_backup.map ~/Sociogenomics/Results/
ls ~/Sociogenomics/Results/

# Removing files
rm ~/Sociogenomics/Results/hapmap1_backup.map

# --- Part II. Exploring the HapMap data ---

cd ~/Sociogenomics/Data

# Looking at files
head hapmap1.map
head -5 hapmap1.map
tail -5 hapmap1.map
head -3 hapmap1.ped

# Counting lines
wc -l hapmap1.map
wc -l hapmap1.ped

# Searching with grep
grep rs7540009 hapmap1.map
grep rs75 hapmap1.map
grep -c rs75 hapmap1.map

# Searching on a specific chromosome
grep "^22" hapmap1.map | head
grep -c "^22" hapmap1.map

# Combining commands with pipes
grep "^1" hapmap1.map | wc -l

# --- Part III. Data manipulation with AWK ---

# Printing columns
awk '{ print $1, $2 }' hapmap1.map | head
awk '{ print "chr"$1, $2, "pos:"$4 }' hapmap1.map | head

# Built-in variables
awk '{ print NR, $1, $2 }' hapmap1.map | head
awk '{ print NF; exit }' hapmap1.ped

# Filtering rows
awk '$1 == 2 { print $1, $2, $4 }' hapmap1.map | head
awk 'NR >= 100 && NR <= 105' hapmap1.map

# Pattern matching
awk '$2 == "rs7540009" { print $0 }' hapmap1.map
awk '$2 ~ /rs75/ { print $1, $2, $4 }' hapmap1.map

# Counting SNPs per chromosome
awk '{ count[$1]++ } END { for (chr in count) print "chr"chr, count[chr] }' hapmap1.map | sort -n -k1.4

# Working with the BMI phenotype file
head BMI_pheno.txt
awk 'NR <= 6 { print $1, $2, $3 }' BMI_pheno.txt
awk 'NR > 1 && $3 > 30 { print $1, $2, $3 }' BMI_pheno.txt | head
awk 'NR > 1 && $3 > 30' BMI_pheno.txt | wc -l

# Changing delimiters
awk 'BEGIN { OFS="\t" } { print $1, $2, $4 }' hapmap1.map | head

# Redirecting output to a file
awk '$1 == 1 { print $2, $4 }' hapmap1.map > ../Results/chr1_snps.txt
head ../Results/chr1_snps.txt
wc -l ../Results/chr1_snps.txt

# Counting lines with AWK
awk 'END { print NR }' hapmap1.map
