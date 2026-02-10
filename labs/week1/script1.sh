#!/bin/bash
set -e  # Stop the script if any command fails

# =============================================================================
# Lab 1 - Introducing the Unix Shell
# Sociogenomics 2025/2026
# Compatible with Google Cloud Shell
# =============================================================================

# --- PART I: Setting up the project directory structure ---

# Create the main project directory and subdirectories for organizing files
mkdir Sociogenomics
mkdir Sociogenomics/Data       # Raw data files (genotype data, phenotype files, GWAS summary stats)
mkdir Sociogenomics/Results    # Output files from analyses
mkdir Sociogenomics/Software   # Bioinformatics tools (e.g., PLINK)
mkdir Sociogenomics/Scripts    # Shell scripts and analysis code

# --- Download and extract lab data ---

# Download the week 1 data archive and move it to the Data folder
# Using curl (available by default on Google Cloud Shell)
curl -L -o week1.zip http://nicolabarban.com/sociogenomics2022/week1/data/week1.zip
mv week1.zip Sociogenomics/Data/

# Navigate to the Data directory and extract the zip file
cd Sociogenomics/Data
unzip week1.zip

# --- Exploring files with head and wc ---

# Display the first 10 lines of the hapmap genetic map file (default for head)
head hapmap1.map

# Display the first 10 lines explicitly
head -10 hapmap1.map

# Count the number of lines in the .map file (SNP information)
wc -l hapmap1.map
# Count the number of lines in the .ped file (genotype data per individual)
wc -l hapmap1.ped

# --- Using echo and pipelines ---

# echo prints text to the terminal; backslash (\) continues the command on the next line
echo 'Hello' \
 'world'

# Pipe (|) sends the output of echo into wc, which counts lines, words, and characters
echo 'Hello world' | wc

# --- Searching with grep ---

# Search for a specific SNP (rs7540009) in the hapmap file
grep rs7540009 hapmap1

# Search for all SNPs whose ID starts with "rs75"
grep rs75 hapmap1

# Count how many SNPs match the pattern "rs75"
grep rs75 hapmap1 | wc -l

# --- Downloading GWAS summary statistics from UK Biobank ---

# Download height GWAS summary statistics (Neale Lab, UK Biobank)
# -L follows redirects, -o specifies the output filename
curl -L -o sumstatsUKB_height.tsv.gz "https://www.dropbox.com/s/ou12jm89v74k55e/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz?dl=1"

# Decompress the gzipped file
gunzip -d sumstatsUKB_height.tsv.gz

# --- Inspecting the GWAS summary statistics ---

# View the first 10 lines (header + first rows) of the summary statistics file
head sumstatsUKB_height.tsv
# Count lines, words, and characters in the file
wc sumstatsUKB_height.tsv

# Preview the first 20 lines of the file (use 'less' or 'more' for interactive browsing)
head -20 sumstatsUKB_height.tsv

# --- Introduction to AWK ---

# Print the entire contents of the BMI phenotype file (equivalent to cat)
cat BMI_pheno.txt
awk '{print}' BMI_pheno.txt

# Print only the first two columns (e.g., FID and IID)
awk '{print $1, $2  }'  BMI_pheno.txt

# --- Searching for a specific SNP in GWAS results ---

# Find a specific variant (chr3, position 49860854) using awk pattern matching
# Print columns: variant ID ($1), position ($2), allele ($5), p-value ($11)
awk ' /3:49860854/ {print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv
# Same result using awk + grep pipeline
awk ' {print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv | grep 3:49860854

# --- AWK built-in variables: NR (row number) and NF (number of fields) ---

# Print the row number (NR) along with columns 1 and 2
awk '{ print  NR, $1, $2  }'  BMI_pheno.txt

# Print the first column and the last column ($NF = last field)
awk '{print $1,$NF}' BMI_pheno.txt

# Print rows 3 through 6, with their row numbers
awk 'NR==3, NR==6 {print NR,$0}'  BMI_pheno.txt

# Print row number followed by a dash and the first column
awk '{print NR "- " $1 }' BMI_pheno.txt

# --- Conditional printing with if statements ---

# Print columns 1 and 2 only for the first row (header)
awk '{ if(NR==1) print $1, $2  }'  BMI_pheno.txt

# Print columns 1 and 2 for the first 9 rows
awk '{ if(NR<10) print $1, $2  }'  BMI_pheno.txt

# Print the first 9 rows of the GWAS summary statistics (all columns)
awk '{ if(NR<10) print   }'  sumstatsUKB_height.tsv

# --- Filtering by p-value ---

# Print SNPs with genome-wide significant p-value (< 5e-08)
awk '{ if($11<5e-08) print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv

# Same filter but transform p-values to natural log scale
awk '{ if($11<5e-08) print  $1, $2, $5, log($11) }'  sumstatsUKB_height.tsv

# Count total number of rows (lines) in the file
 awk 'END { print NR }' sumstatsUKB_height.tsv

# --- Changing field delimiters ---

# Use colon (:) as field separator (useful for chr:pos format in variant column)
awk 'FS=":" {print $1, $2, $3}' sumstatsUKB_height.tsv | head
# Use colon as input separator and dash as output separator
awk 'FS=":", OFS="-" {print $1, $2, $3}' sumstatsUKB_height.tsv | head

# Use multiple field separators (colon, double-quote, tab) to parse complex fields
awk -F '[:"\t"]' '{ if(NR>1) print $1, $2, $3 , $4, $13}' sumstatsUKB_height.tsv | head

# --- Redirecting output to files ---

# Save genome-wide significant SNPs to a results file
awk '{ if($11<5e-08) print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv > ../Results/sign_variants_UKB.txt

# Save common variants (MAF > 10%) to a results file
awk '{ if($3>0.1) print  $1, $2, $5, $11 }'  sumstatsUKB_height.tsv > ../Results/common_variants_UKB.txt

# Verify the output file
head  ../Results/common_variants_UKB.txt

# --- Combining multiple conditions ---

# Count common SNPs (MAF > 10%) on Chromosome 1
awk -F '[:"\t"]' '{ if($6 >0.1 && $1==1) print $1, $2, $3, $4,  $6}' sumstatsUKB_height.tsv | wc -l

# Count common SNPs (MAF > 10%) on Chromosome 21
awk -F '[:"\t"]' '{ if($6 >0.1 && $1==21) print $1, $2, $3, $4,  $6}' sumstatsUKB_height.tsv | wc -l

# --- PART II: Installing PLINK (genetic analysis tool) ---

# Navigate to the Software directory
cd $HOME
cd Sociogenomics/Software

# Download and extract PLINK (v1.9, Linux x86_64 - compatible with Google Cloud Shell)
curl -L -o plink_linux_x86_64_20210606.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip
unzip plink_linux_x86_64_20210606.zip

# Make the PLINK binary executable
chmod +x plink

# Test that PLINK works
./plink --help

# --- Creating a symbolic link and running PLINK ---

# Create a symbolic link to plink in the main project directory
# This allows running plink from the project root without specifying the full path
cd $HOME/Sociogenomics
 ln -s Software/plink

# Verify plink is accessible via the symlink
./plink --help

# Run PLINK to compute allele frequencies from the hapmap data
# --file: input file prefix (reads hapmap1.ped and hapmap1.map)
# --freq: calculate minor allele frequencies
# --out: output file prefix
./plink --file   Data/hapmap1 --freq --out Results/test

# View the allele frequency results
head -20 Results/test.frq
