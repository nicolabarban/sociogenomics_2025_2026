# Lab 2. Getting Started with PLINK

In this lab you will install PLINK on Google Cloud Shell and learn how to work with genetic data in PLINK format. By the end you will be able to convert between file formats, compute summary statistics, and filter SNPs and individuals.

We use the **HapMap Phase III** dataset (`hapmap3`), which contains 1,184 individuals and ~1.4 million SNPs.

---

## 0. Getting started

Open [Google Cloud Shell](https://shell.cloud.google.com/) in your browser.

First, install Git LFS (needed to download the data files):

```
sudo apt-get install -y git-lfs
git lfs install
```

If you have not already cloned the course repository (from Lab 1), run:

```
cd $HOME
git clone https://github.com/nicolabarban/sociogenomics_2025_2026.git
```

If you already cloned it, pull the latest updates:

```
cd ~/sociogenomics_2025_2026
git pull
git lfs pull
cd $HOME
```

> **Note:** The data files are stored with Git LFS. The `git lfs pull` command downloads the actual file contents. If you skip this step, the data files will contain only small text pointers instead of real data.

Make sure your project directories from Lab 1 still exist:

```
mkdir -p ~/Sociogenomics/Data ~/Sociogenomics/Results ~/Sociogenomics/Scripts
```

Copy the HapMap data:

```
cp ~/sociogenomics_2025_2026/data/hapmap3.bed ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/hapmap3.bim ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/hapmap3.fam ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/BMI_pheno.txt ~/Sociogenomics/Data/
```

---

## 1. Installing PLINK

PLINK is a free, open-source tool for analysing genetic data. It can handle data management, quality control, and association analysis. We use **PLINK 1.9**, available from [cog-genomics.org/plink](https://www.cog-genomics.org/plink/).

Run the setup script from the course repository:

```
bash ~/sociogenomics_2025_2026/scripts/setup_plink19.sh
```

This downloads PLINK, makes it executable, and adds it to your PATH. Reload your shell so the PATH update takes effect:

```
source ~/.bashrc
```

Verify the installation:

```
plink --help
```

You should see a help message listing PLINK commands and options. If you get `command not found`, close and reopen your Cloud Shell terminal, then try again.

---

## Part I. PLINK file formats

PLINK works with two main file-format families:

| Format | Files | Description |
|--------|-------|-------------|
| **Text (PED/MAP)** | `.ped` + `.map` | Human-readable but large and slow |
| **Binary (BED/BIM/FAM)** | `.bed` + `.bim` + `.fam` | Compact and fast; `.bed` is not human-readable |

Our `hapmap3` dataset is already in **binary format**. Let's explore each file.

### The BIM file

The `.bim` file is an extended map file with one row per SNP and 6 columns: chromosome, SNP ID, genetic distance (cM), base-pair position, allele 1, allele 2.

```
cd ~/Sociogenomics/Data
head hapmap3.bim
```

### The FAM file

The `.fam` file has one row per individual with 6 columns: Family ID, Individual ID, Father ID, Mother ID, Sex (1=male, 2=female), Phenotype (-9=missing).

```
head hapmap3.fam
```

### The BED file

The `.bed` file stores genotype data in compressed binary. It is **not** human-readable:

```
head -c 20 hapmap3.bed
```

You will see garbled characters. This is expected.

### Counting individuals and SNPs

```
wc -l hapmap3.fam
wc -l hapmap3.bim
```

**Question:** How many individuals and how many SNPs are in the dataset?

### Converting binary to text format

You can convert binary files to PED/MAP with `--recode`:

```
plink --bfile hapmap3 --chr 22 --recode --out hapmap3_chr22_text
```

We add `--chr 22` to convert only chromosome 22 (the full dataset would produce a very large PED file).

```
head hapmap3_chr22_text.map
head -2 hapmap3_chr22_text.ped | cut -c1-80
```

The `.map` file has 4 columns: chromosome, SNP ID, genetic distance, base-pair position.

The `.ped` file has one row per individual. The first 6 columns are: Family ID, Individual ID, Father ID, Mother ID, Sex, Phenotype. The remaining columns contain two alleles per SNP. We use `cut` to show only the first 80 characters because each row is very wide.

### Converting text back to binary

Convert back to binary with `--make-bed`:

```
plink --file hapmap3_chr22_text --make-bed --out hapmap3_chr22_binary
ls -lh hapmap3_chr22_binary.bed hapmap3_chr22_binary.bim hapmap3_chr22_binary.fam
```

Notice how much smaller the binary files are compared to the text files.

### Exercise 1

1. How many individuals are in the dataset? How many are male and how many are female? Hint: use `awk` on the `.fam` file to count by the sex column.
2. How many SNPs are in `hapmap3.bim`?
3. Use `awk` to count how many SNPs in `hapmap3.bim` are on chromosome 1.
4. Use `grep` to find the SNP `rs9930506` in `hapmap3.bim`. What chromosome is it on and what are its two alleles?

---

## Part II. Summary statistics with PLINK

PLINK can compute a range of useful summary statistics directly from binary files.

### Allele frequencies

Calculate the minor allele frequency (MAF) for every SNP:

```
plink --bfile hapmap3 --freq --out allele_freq
```

Inspect the output:

```
head allele_freq.frq
```

The `.frq` file has columns: CHR, SNP, A1 (minor allele), A2 (major allele), MAF, NCHROBS (number of allele observations).

Find the frequency of a specific SNP:

```
grep rs9930506 allele_freq.frq
```

### Missing data rates

Compute per-SNP and per-individual missing rates:

```
plink --bfile hapmap3 --missing --out missing_report
```

This creates two files:

- `missing_report.imiss` — per-**i**ndividual missing rates
- `missing_report.lmiss` — per-**l**ocus (SNP) missing rates

```
head missing_report.imiss
head missing_report.lmiss
```

The key column is `F_MISS` — the fraction of missing genotypes.

Find individuals with the highest missingness:

```
sort -k6 -n -r missing_report.imiss | head
```

### Hardy-Weinberg equilibrium

Test each SNP for deviation from Hardy-Weinberg equilibrium (HWE):

```
plink --bfile hapmap3 --hardy --out hwe_report
```

```
head hwe_report.hwe
```

The output shows observed and expected genotype counts and a p-value. SNPs with very low p-values may indicate genotyping errors.

Find SNPs with HWE p-value below 1e-6:

```
awk '$9 < 1e-6' hwe_report.hwe | head
```

### Exercise 2

1. What is the minor allele frequency of SNP `rs9930506`? Use the `.frq` file to find it.
2. How many SNPs have a MAF below 0.05? Hint: use `awk` on the `.frq` file to filter by the MAF column and count with `wc -l`.
3. Which individual has the highest rate of missing genotypes? What is their missing rate?

---

## Part III. Filtering and extracting data

PLINK makes it easy to create subsets of your data by filtering on SNPs, individuals, or chromosomes.

### Extract a single SNP

Extract only the SNP `rs9930506`:

```
plink --bfile hapmap3 --snp rs9930506 --make-bed --out rs9930506_only
```

Check the result:

```
wc -l rs9930506_only.bim
wc -l rs9930506_only.fam
```

You should see 1 SNP and the same number of individuals as the original file.

### Extract a range of SNPs

Extract SNPs between `rs3751813` and `rs8044769` (inclusive, based on genomic order within a chromosome):

```
plink --bfile hapmap3 --from rs3751813 --to rs8044769 --make-bed --out snp_range
wc -l snp_range.bim
```

### Extract SNPs from a list

Create a text file with one SNP ID per line:

```
echo -e "rs9930506\nrs1048488\nrs7520934" > snp_list.txt
cat snp_list.txt
```

Use `--extract` to keep only these SNPs:

```
plink --bfile hapmap3 --extract snp_list.txt --make-bed --out selected_snps
wc -l selected_snps.bim
```

### Filter by chromosome

Keep only SNPs on chromosome 22:

```
plink --bfile hapmap3 --chr 22 --make-bed --out chr22_only
wc -l chr22_only.bim
```

### Filter by minor allele frequency

Keep only common SNPs (MAF >= 0.05):

```
plink --bfile hapmap3 --maf 0.05 --make-bed --out common_snps
wc -l common_snps.bim
```

Compare with the total number of SNPs:

```
wc -l hapmap3.bim
```

**Question:** How many SNPs were removed by the MAF filter?

### Filter by missingness

Remove SNPs with more than 2% missing data and individuals with more than 5% missing data:

```
plink --bfile hapmap3 --geno 0.02 --mind 0.05 --make-bed --out qc_filtered
```

Check the log to see how many SNPs and individuals were removed.

### Combine filters

You can combine multiple filters in a single command:

```
plink --bfile hapmap3 --chr 1 --maf 0.01 --geno 0.02 --make-bed --out chr1_clean
wc -l chr1_clean.bim
wc -l chr1_clean.fam
```

### Exercise 3

1. Extract all SNPs on chromosome 6 with a MAF above 0.10 into a new binary file called `chr6_common`. How many SNPs remain?
2. Create a text file with three SNP IDs of your choice (look at `hapmap3.bim` for options). Use `--extract` to create a new dataset with only those SNPs. Verify the number of SNPs in the output `.bim` file.
3. Apply the following quality control filters to the full dataset and create a cleaned file called `hapmap3_qc`:
   - Remove SNPs with more than 5% missing data (`--geno 0.05`)
   - Remove individuals with more than 10% missing data (`--mind 0.10`)
   - Remove SNPs with MAF below 0.01 (`--maf 0.01`)
   - Remove SNPs that deviate from Hardy-Weinberg equilibrium with p < 1e-6 (`--hwe 1e-6`)

   How many SNPs and individuals remain after QC? Check the PLINK log output.
