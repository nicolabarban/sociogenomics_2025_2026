# Lab 4. Genetic Tests and Principal Component Analysis

In this lab you will learn how to explore genetic data in R, run statistical tests, compute Principal Component Analysis (PCA) to detect population structure, and use PCA scores to predict ancestry. By the end you will be able to visualise genotype quality distributions, test Hardy-Weinberg equilibrium, run association tests, produce annotated PCA plots that reveal ancestry differences, and classify individuals into continental ancestry groups using machine learning.

We continue using the **HapMap Phase III** dataset (`hapmap3`) from the course repository.

---

## 0. Getting started

Open [Google Cloud Shell](https://shell.cloud.google.com/) in your browser.

Update the course repository and pull the data files:

```bash
cd ~/sociogenomics_2025_2026
git pull
git lfs pull
cd $HOME
```

Make sure your project directories exist:

```bash
mkdir -p ~/Sociogenomics/Data ~/Sociogenomics/Results ~/Sociogenomics/Scripts
```

Copy the data files (if not already present from Lab 3):

```bash
cp ~/sociogenomics_2025_2026/data/hapmap3.bed ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/hapmap3.bim ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/hapmap3.fam ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/BMI_pheno.txt ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/1kg_samples.txt ~/Sociogenomics/Data/
```

Move to your working data directory:

```bash
cd ~/Sociogenomics/Data
```

Make sure PLINK is available:

```bash
plink --version
```

If not, reinstall it:

```bash
bash ~/sociogenomics_2025_2026/scripts/setup_plink19.sh
source ~/.bashrc
```

Apply the standard QC filters from Lab 3 to get a clean dataset (skip this if `hapmap3_qc` already exists):

```bash
plink --bfile hapmap3 \
      --mind 0.05 \
      --geno 0.02 \
      --maf  0.01 \
      --hwe  1e-6 \
      --make-bed \
      --out hapmap3_qc
```

---

## Part I. Exploring Genetic Data in R

Before running any statistical tests, it is good practice to load the data into R and inspect its basic properties. PLINK binary files cannot be read directly in R, but PLINK can export summary tables that we can then explore.

### 1.0 Generate summary statistics from PLINK

```bash
cd ~/Sociogenomics/Data

# Compute per-SNP missingness, allele frequencies, and HWE
plink --bfile hapmap3_qc \
      --missing \
      --freq \
      --hardy \
      --out hapmap3_summary
```

Now start R:

```bash
R
```

### 1.0a Load and inspect summary statistics

```R
library(data.table)
library(ggplot2)

# --- SNP missingness ---
lmiss <- fread("hapmap3_summary.lmiss")
cat("SNPs in dataset:", nrow(lmiss), "\n")
cat("Range of SNP missingness:", range(lmiss$F_MISS), "\n")

ggplot(lmiss, aes(x = F_MISS)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = "white") +
  theme_bw() +
  xlab("Per-SNP missing rate") +
  ylab("Number of SNPs") +
  ggtitle("Distribution of SNP missingness")
ggsave("hist_snp_miss.png", width = 7, height = 5)

# --- Individual missingness ---
imiss <- fread("hapmap3_summary.imiss")
cat("Individuals in dataset:", nrow(imiss), "\n")
cat("Range of individual missingness:", range(imiss$F_MISS), "\n")

ggplot(imiss, aes(x = F_MISS)) +
  geom_histogram(bins = 40, fill = "coral", colour = "white") +
  theme_bw() +
  xlab("Per-individual missing rate") +
  ylab("Number of individuals") +
  ggtitle("Distribution of individual missingness")
ggsave("hist_ind_miss.png", width = 7, height = 5)
```

### 1.0b Minor allele frequency distribution

```R
# --- MAF distribution ---
frq <- fread("hapmap3_summary.frq")
cat("Mean MAF:", round(mean(frq$MAF), 4), "\n")
cat("Median MAF:", round(median(frq$MAF), 4), "\n")

ggplot(frq, aes(x = MAF)) +
  geom_histogram(bins = 50, fill = "darkgreen", colour = "white") +
  theme_bw() +
  xlab("Minor allele frequency") +
  ylab("Number of SNPs") +
  ggtitle("MAF distribution across all SNPs") +
  geom_vline(xintercept = 0.05, colour = "red",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 0.06, y = Inf, vjust = 1.5,
           label = "MAF = 0.05", colour = "red", size = 3.5)
ggsave("hist_maf.png", width = 7, height = 5)
```

The red dashed line shows the common MAF = 0.05 filter. SNPs to the left of it are rare variants.

### 1.0c Hardy-Weinberg equilibrium distribution

```R
# --- HWE p-value distribution ---
hwe <- fread("hapmap3_summary.hwe")
# Keep only "ALL" rows (not case/control split)
hwe_all <- hwe[TEST == "ALL"]

# -log10(p) plot
hwe_all[, log10p := -log10(P)]

ggplot(hwe_all, aes(x = log10p)) +
  geom_histogram(bins = 60, fill = "purple", colour = "white") +
  theme_bw() +
  xlab(expression(-log[10](p))) +
  ylab("Number of SNPs") +
  ggtitle("Distribution of HWE test statistics") +
  geom_vline(xintercept = 6, colour = "red",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 6.2, y = Inf, vjust = 1.5,
           label = "p = 1e-6 threshold", colour = "red", size = 3.5)
ggsave("hist_hwe.png", width = 7, height = 5)

cat("SNPs failing HWE (p < 1e-6):", sum(hwe_all$P < 1e-6, na.rm = TRUE), "\n")
```

### 1.0d Observed vs expected heterozygosity per individual

Individuals with extreme heterozygosity deviate from the population average, which may indicate contamination or inbreeding.

```R
# Compute heterozygosity per individual from the HWE file
# O(HET) and E(HET) are in columns 7 and 8
het_snp <- hwe_all[, .(SNP, `O(HET)`, `E(HET)`)]

# For per-individual heterozygosity we use the .het file (requires --het flag)
# Run from bash: plink --bfile hapmap3_qc --het --out hapmap3_het
# Then load:
het <- fread("hapmap3_het.het")
# F = inbreeding coefficient; negative = excess heterozygosity
het[, obs_het := (`O(HOM)` - `E(HOM)`) / `N(NM)`]

ggplot(het, aes(x = F)) +
  geom_histogram(bins = 50, fill = "orange", colour = "white") +
  geom_vline(xintercept = c(-0.15, 0.15), colour = "red",
             linetype = "dashed", linewidth = 0.8) +
  theme_bw() +
  xlab("Inbreeding coefficient F") +
  ylab("Number of individuals") +
  ggtitle("Per-individual inbreeding coefficient") +
  annotate("text", x = -0.16, y = Inf, vjust = 1.5, hjust = 1,
           label = "-0.15", colour = "red", size = 3.5) +
  annotate("text", x = 0.16, y = Inf, vjust = 1.5, hjust = 0,
           label = "+0.15", colour = "red", size = 3.5)
ggsave("hist_inbreeding.png", width = 7, height = 5)

outliers <- het[F < -0.15 | F > 0.15]
cat("Heterozygosity outliers:", nrow(outliers), "\n")
if (nrow(outliers) > 0) print(outliers[, .(FID, IID, F)])

q()
```

To generate the `.het` file used above, run this PLINK command before starting R:

```bash
plink --bfile hapmap3_qc --het --out hapmap3_het
```

---

## Part II. Statistical Tests on Genetic Data



### 2.1 Hardy-Weinberg Equilibrium (HWE) test

**Hardy-Weinberg Equilibrium** (HWE) states that in a large random-mating population with no selection, mutation, or migration, allele and genotype frequencies remain constant across generations. For a biallelic SNP with allele frequencies $p$ (allele A) and $q = 1 - p$ (allele B), the expected genotype frequencies are:

$$P(AA) = p^2, \quad P(AB) = 2pq, \quad P(BB) = q^2$$

**Why test for HWE?** Deviations from HWE in a control sample usually indicate genotyping errors — differential allelic dropout, probe failure, or sample contamination — rather than a true biological signal. SNPs with $p < 10^{-6}$ in controls are typically removed in QC.

HWE is tested using a chi-squared test with 1 degree of freedom (or an exact test for small samples):

$$\chi^2 = \frac{(O_{AA} - E_{AA})^2}{E_{AA}} + \frac{(O_{AB} - E_{AB})^2}{E_{AB}} + \frac{(O_{BB} - E_{BB})^2}{E_{BB}}$$

Compute HWE statistics for all SNPs:

```bash
plink --bfile hapmap3_qc \
      --hardy \
      --out hwe_results
```

Inspect the output:

```bash
head hwe_results.hwe
```

The columns are: CHR, SNP, TEST (ALL / AFF / UNAFF), A1, A2, GENO (observed counts: homozygous A1 / heterozygous / homozygous A2), O(HET) observed heterozygosity, E(HET) expected heterozygosity, P (p-value).

Count SNPs at different significance thresholds:

```bash
# How many SNPs deviate at p < 0.05?
awk 'NR>1 && $9 < 0.05' hwe_results.hwe | wc -l

# How many at p < 1e-6 (typical QC threshold)?
awk 'NR>1 && $9 < 1e-6' hwe_results.hwe | wc -l

# How many total SNPs?
awk 'NR>1' hwe_results.hwe | wc -l
```

Find the SNPs with the strongest HWE deviation:

```bash
awk 'NR>1 {print $2, $9}' hwe_results.hwe | sort -k2 -n | head -10
```

**Interpreting observed vs expected heterozygosity:**

An excess of observed heterozygosity (`O(HET)` >> `E(HET)`) often indicates sample contamination (two individuals' DNA mixed together). A deficit of heterozygosity suggests inbreeding or genotyping errors causing allele dropout.

Examine the heterozygosity for the most deviant SNP:

```bash
# Get the most deviant SNP name
TOP_SNP=$(awk 'NR>1 {print $2, $9}' hwe_results.hwe | sort -k2 -n | head -1 | awk '{print $1}')
echo "Most deviant SNP: $TOP_SNP"
grep "$TOP_SNP" hwe_results.hwe
```

### 2.2 Allele frequency statistics

Before running association tests, it is useful to understand the allele frequency distribution of your dataset.

Compute allele frequencies:

```bash
plink --bfile hapmap3_qc \
      --freq \
      --out allele_freq
```

Inspect the output:

```bash
head allele_freq.frq
```

The columns are: CHR, SNP, A1 (minor allele), A2 (major allele), MAF (minor allele frequency), NCHROBS (number of allele observations).

Summarise the frequency distribution:

```bash
awk 'NR>1 {
    if ($5 < 0.01) bin1++
    else if ($5 < 0.05) bin2++
    else if ($5 < 0.10) bin3++
    else if ($5 < 0.20) bin4++
    else bin5++
} END {
    print "MAF < 0.01  (rare):", bin1
    print "MAF 0.01-0.05:", bin2
    print "MAF 0.05-0.10:", bin3
    print "MAF 0.10-0.20:", bin4
    print "MAF >= 0.20 (common):", bin5
}' allele_freq.frq
```

### 2.3 Chi-squared association test

The simplest genetic association test asks: is a SNP's allele frequency different between cases and controls?

For a case-control study, PLINK computes a 2×2 contingency table (two alleles × two groups) and applies a chi-squared test with 1 degree of freedom:

|         | Allele A | Allele B |
|---------|----------|----------|
| Cases   | a        | b        |
| Controls| c        | d        |

$$\chi^2 = \frac{N(ad - bc)^2}{(a+b)(c+d)(a+c)(b+d)}$$

For a **quantitative trait** (like height or BMI), a linear regression is used instead:

$$y = \beta_0 + \beta_1 x + \epsilon$$

where $y$ is the phenotype, $x$ is the SNP genotype (coded 0, 1, 2 copies of the effect allele), and $\beta_1$ is the additive effect estimate.

Run a basic association test with BMI as a quantitative phenotype:

```bash
plink --bfile hapmap3_qc \
      --pheno BMI_pheno.txt \
      --assoc \
      --linear \
      --out bmi_assoc
```

Inspect the output:

```bash
head bmi_assoc.assoc.linear
```

The key columns are: CHR, SNP, BP (position), A1 (effect allele), TEST, NMISS, BETA (effect size), STAT (test statistic), P (p-value).

Find the most significant associations:

```bash
# Skip header, sort by p-value ascending
awk 'NR>1' bmi_assoc.assoc.linear | sort -k9 -n | head -20
```

Count genome-wide significant hits (p < 5×10⁻⁸, the standard GWAS threshold):

```bash
awk 'NR>1 && $9 < 5e-8' bmi_assoc.assoc.linear | wc -l
```

The GWAS significance threshold of $p < 5 \times 10^{-8}$ is a Bonferroni correction for approximately 1 million independent tests across the genome. It controls the family-wise error rate at 5%.

### 2.4 Test a single SNP

To examine a specific SNP in detail:

```bash
# Test a specific SNP (replace with a SNP of interest from your results)
TOP_SNP=$(awk 'NR>1' bmi_assoc.assoc.linear | sort -k9 -n | head -1 | awk '{print $2}')
echo "Top SNP: $TOP_SNP"

plink --bfile hapmap3_qc \
      --pheno BMI_pheno.txt \
      --snp $TOP_SNP \
      --assoc \
      --linear \
      --out top_snp_test

cat top_snp_test.assoc.linear
```

Test the same SNP under a dominant model (one copy of the effect allele is enough):

```bash
plink --bfile hapmap3_qc \
      --pheno BMI_pheno.txt \
      --snp $TOP_SNP \
      --assoc \
      --linear dominant \
      --out top_snp_dominant

cat top_snp_dominant.assoc.linear
```

Test under a recessive model (two copies required):

```bash
plink --bfile hapmap3_qc \
      --pheno BMI_pheno.txt \
      --snp $TOP_SNP \
      --assoc \
      --linear recessive \
      --out top_snp_recessive

cat top_snp_recessive.assoc.linear
```

**Question:** Does the p-value change across the three models? Which model fits best?

### Exercise 1

1. How many SNPs in `hapmap3_qc` deviate significantly from HWE at $p < 10^{-6}$? What could cause this?
2. Look at the 10 SNPs with the strongest HWE deviation. Do their observed and expected heterozygosities suggest excess or deficit heterozygosity?
3. Run the BMI association and report the top 5 SNPs by p-value. What are their chromosomal positions?

---

## Part III. Principal Component Analysis (PCA)

### 3.1 Why PCA in genetics?

PCA is used in genomics to:
1. **Detect population stratification** — systematic ancestry differences between cases and controls that can confound association tests
2. **Visualise ancestry** — reveal continental and sub-continental genetic clusters
3. **Create covariates** — PC scores are added to GWAS regression models to control for stratification

PCA works by decomposing the genotype matrix $\mathbf{G}$ (individuals × SNPs) into orthogonal components that capture the directions of maximum variance. Each PC summarises a linear combination of allele frequencies. Individuals with similar ancestry will cluster together in PC space because they share similar allele frequencies genome-wide.

**Genomic PCA requires LD-pruned SNPs.** Regions in high LD would over-represent certain loci, distorting the PC directions. We first prune to an approximately independent SNP set.

### 3.2 Prepare the pruned SNP set

```bash
# LD pruning: remove one SNP from each pair with r² > 0.2
plink --bfile hapmap3_qc \
      --indep-pairwise 50 5 0.2 \
      --out hapmap3_pruned

# Extract the pruned SNP set
plink --bfile hapmap3_qc \
      --extract hapmap3_pruned.prune.in \
      --make-bed \
      --out hapmap3_pruned_set
```

How many SNPs remain after pruning?

```bash
wc -l hapmap3_pruned.prune.in
```

### 3.3 Compute PCA with PLINK

```bash
plink --bfile hapmap3_pruned_set \
      --pca 20 \
      --out hapmap3_pca
```

The `--pca 20` flag computes the top 20 principal components. PLINK produces two output files:
- `hapmap3_pca.eigenvec` — PC scores for each individual (eigenvectors)
- `hapmap3_pca.eigenval` — eigenvalues (variance explained by each PC)

Inspect the outputs:

```bash
head hapmap3_pca.eigenvec
head hapmap3_pca.eigenval
```

The `eigenvec` file has columns: FID, IID, PC1, PC2, ..., PC20.

### 3.4 Variance explained by each PC

The eigenvalues tell you how much variance each PC explains. Compute the proportion of variance explained:

```bash
# Total variance = sum of all eigenvalues
# Proportion for each PC = eigenvalue / total
awk 'BEGIN{sum=0} {val[NR]=$1; sum+=$1}
     END{
       for(i=1; i<=NR; i++)
         printf "PC%d: %.2f%%\n", i, val[i]/sum*100
     }' hapmap3_pca.eigenval
```

Typically, the first 2-5 PCs capture the major axes of continental ancestry, while later PCs capture finer-scale structure or noise.

### 3.5 Visualise PCA in R

Start R:

```bash
R
```

**Plot PC1 vs PC2 (all individuals, no colour):**

```R
library(ggplot2)
library(data.table)

# Read PC scores
pca <- fread("hapmap3_pca.eigenvec", header = FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:20))

# Basic scatter plot
ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.6, size = 1.5) +
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  ggtitle("PCA of HapMap3 — all individuals")

ggsave("pca_basic.png", width = 7, height = 6)
```

**Add population labels from the 1000 Genomes sample file:**

```R
# Read sample information
geo <- fread("1kg_samples.txt", sep = "\t", header = TRUE)

# Rename the sample column to match the PCA file
# The sample file uses "Sample name" as column 1
colnames(geo)[1] <- "IID"

# Merge on individual ID
data <- merge(pca, geo[, c("IID", "Population code",
                            "Population name",
                            "Superpopulation code",
                            "Superpopulation name")],
              by = "IID")

# Plot coloured by superpopulation
ggplot(data, aes(x = PC1, y = PC2,
                 colour = `Superpopulation name`)) +
  geom_point(alpha = 0.7, size = 1.5) +
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "Superpopulation",
       title = "PCA coloured by continental ancestry") +
  theme(legend.position = "right")

ggsave("pca_superpop.png", width = 8, height = 6)
```

You should see clearly separated clusters corresponding to the five 1000 Genomes superpopulations: **AFR** (African), **EUR** (European), **EAS** (East Asian), **SAS** (South Asian), **AMR** (Admixed American).

**Plot coloured by sub-population:**

```R
ggplot(data, aes(x = PC1, y = PC2,
                 colour = `Population name`)) +
  geom_point(alpha = 0.7, size = 1.5) +
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "Population",
       title = "PCA coloured by sub-population") +
  theme(legend.position = "right",
        legend.text = element_text(size = 7))

ggsave("pca_population.png", width = 10, height = 6)
```

**Plot PC1 vs PC3 to reveal a third axis of variation:**

```R
ggplot(data, aes(x = PC1, y = PC3,
                 colour = `Superpopulation name`)) +
  geom_point(alpha = 0.7, size = 1.5) +
  theme_bw() +
  xlab("PC1") + ylab("PC3") +
  labs(colour = "Superpopulation",
       title = "PC1 vs PC3")

ggsave("pca_pc1_pc3.png", width = 8, height = 6)
```

**Scree plot — variance explained by each PC:**

```R
eigenval <- fread("hapmap3_pca.eigenval", header = FALSE)
colnames(eigenval) <- "eigenvalue"
eigenval$PC <- 1:nrow(eigenval)
eigenval$pct <- eigenval$eigenvalue / sum(eigenval$eigenvalue) * 100

ggplot(eigenval, aes(x = PC, y = pct)) +
  geom_col(fill = "steelblue") +
  geom_line(aes(group = 1)) +
  geom_point() +
  theme_bw() +
  xlab("Principal Component") +
  ylab("Variance explained (%)") +
  ggtitle("Scree plot") +
  scale_x_continuous(breaks = 1:20)

ggsave("pca_scree.png", width = 8, height = 5)
```

Exit R:

```R
q()
```

### 3.6 Population-specific PCA

You can run PCA within a single population to reveal finer-scale structure. Create population-specific sample lists (reuse from Lab 3 if already present):

```bash
awk -F'\t' 'NR>1 && $6 == "EUR" {print $1, $1}' 1kg_samples.txt > samples_EUR.txt
awk -F'\t' 'NR>1 && $6 == "AFR" {print $1, $1}' 1kg_samples.txt > samples_AFR.txt
awk -F'\t' 'NR>1 && $6 == "EAS" {print $1, $1}' 1kg_samples.txt > samples_EAS.txt
```

Run PCA within Europeans:

```bash
plink --bfile hapmap3_pruned_set \
      --keep samples_EUR.txt \
      --pca 10 \
      --out pca_EUR
```

Visualise in R:

```bash
R
```

```R
library(ggplot2)
library(data.table)

pca_eur <- fread("pca_EUR.eigenvec", header = FALSE)
colnames(pca_eur) <- c("FID", "IID", paste0("PC", 1:10))

geo <- fread("1kg_samples.txt", sep = "\t", header = TRUE)
colnames(geo)[1] <- "IID"

data_eur <- merge(pca_eur, geo[, c("IID", "Population name")], by = "IID")

ggplot(data_eur, aes(x = PC1, y = PC2,
                     colour = `Population name`)) +
  geom_point(alpha = 0.8, size = 2) +
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "European population",
       title = "PCA within European populations") +
  theme(legend.position = "right")

ggsave("pca_EUR_subpop.png", width = 8, height = 6)

q()
```

Within European populations, PCs often separate Northern Europeans (Finnish, British, Scandinavian) from Southern Europeans (Iberian, Tuscan) and from Admixed Americans with European ancestry.

### 3.7 PCA as covariates in GWAS

Population stratification can cause spurious associations if cases and controls differ in ancestry. The standard solution is to include PC scores as covariates in the regression model.

Run a stratification-corrected association test using the first 10 PCs:

```bash
plink --bfile hapmap3_qc \
      --pheno BMI_pheno.txt \
      --assoc \
      --linear \
      --covar hapmap3_pca.eigenvec \
      --covar-number 1-10 \
      --out bmi_assoc_pca_corrected
```

Compare the number of significant hits before and after PC correction:

```bash
# Before correction
awk 'NR>1 && $9 < 5e-8' bmi_assoc.assoc.linear | wc -l

# After correction
awk 'NR>1 && $9 < 5e-8' bmi_assoc_pca_corrected.assoc.linear | wc -l
```

**Genomic inflation factor (lambda):** A quick check for residual stratification after PC correction is to compute the genomic inflation factor $\lambda_{GC}$. It is the ratio of the median observed chi-squared statistic to the expected median under the null (0.4549):

$$\lambda_{GC} = \frac{\text{median}(\chi^2_{\text{observed}})}{0.4549}$$

A value of $\lambda_{GC} \approx 1.00$ indicates no inflation. Values much above 1 suggest residual stratification or polygenicity.

Compute $\lambda_{GC}$ in R:

```bash
R
```

```R
library(data.table)

# Uncorrected
res_raw <- fread("bmi_assoc.assoc.linear")
colnames(res_raw)[9] <- "P"
res_raw <- res_raw[!is.na(P)]
chisq_raw <- qchisq(res_raw$P, df = 1, lower.tail = FALSE)
lambda_raw <- median(chisq_raw, na.rm = TRUE) / 0.4549
cat("Lambda (uncorrected):", round(lambda_raw, 3), "\n")

# PC-corrected
res_corr <- fread("bmi_assoc_pca_corrected.assoc.linear")
res_corr <- res_corr[TEST == "ADD"]
colnames(res_corr)[9] <- "P"
res_corr <- res_corr[!is.na(P)]
chisq_corr <- qchisq(res_corr$P, df = 1, lower.tail = FALSE)
lambda_corr <- median(chisq_corr, na.rm = TRUE) / 0.4549
cat("Lambda (PC-corrected):", round(lambda_corr, 3), "\n")

q()
```

### Exercise 2

1. Examine the scree plot. How many PCs are needed to capture the main axes of genetic variation? Where does the curve flatten out?
2. In the PC1 vs PC2 plot coloured by superpopulation, which two superpopulations are most separated along PC1? Along PC2?
3. Within Europeans, which populations are most separated from each other in PC1 vs PC2? What historical events might explain this?
4. How does the genomic inflation factor $\lambda_{GC}$ change before and after adding 10 PCs as covariates? What does this tell you about population stratification in this dataset?

---

## Part IV. Detecting and Removing Population Outliers

In a homogeneous cohort study (e.g., a study of British individuals), individuals who cluster far from the main group in PCA space are likely to have a different ancestry background. They should be identified and removed before association analysis to avoid stratification.

### 4.1 Identify European individuals within a multi-ancestry dataset

Suppose we want to keep only individuals who cluster near the EUR superpopulation in PC space.

```bash
R
```

```R
library(data.table)
library(ggplot2)

pca <- fread("hapmap3_pca.eigenvec", header = FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:20))

geo <- fread("1kg_samples.txt", sep = "\t", header = TRUE)
colnames(geo)[1] <- "IID"
data <- merge(pca, geo[, c("IID", "Superpopulation code")], by = "IID")

# Compute EUR centroid
eur_mean_pc1 <- mean(data[`Superpopulation code` == "EUR", PC1])
eur_mean_pc2 <- mean(data[`Superpopulation code` == "EUR", PC2])
eur_sd_pc1   <- sd(data[`Superpopulation code` == "EUR", PC1])
eur_sd_pc2   <- sd(data[`Superpopulation code` == "EUR", PC2])

cat("EUR centroid: PC1 =", round(eur_mean_pc1, 4),
    ", PC2 =", round(eur_mean_pc2, 4), "\n")

# Flag individuals within 3 SD of EUR centroid on PC1 and PC2
data[, eur_like := abs(PC1 - eur_mean_pc1) < 3 * eur_sd_pc1 &
                   abs(PC2 - eur_mean_pc2) < 3 * eur_sd_pc2]

cat("Individuals within 3 SD of EUR on PC1+PC2:", sum(data$eur_like), "\n")

# Save the EUR-like sample list
eur_keep <- data[eur_like == TRUE, .(FID, IID)]
fwrite(eur_keep, "samples_EUR_like.txt", sep = " ", col.names = FALSE)

# Visualise the selection
ggplot(data, aes(x = PC1, y = PC2,
                 colour = `Superpopulation code`,
                 shape = eur_like)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_shape_manual(values = c(1, 16),
                     labels = c("Excluded", "EUR-like")) +
  theme_bw() +
  labs(colour = "Superpopulation", shape = "Selection",
       title = "EUR-like individuals (within 3 SD of EUR centroid)")

ggsave("pca_eur_selection.png", width = 9, height = 6)
q()
```

Apply the filter in PLINK:

```bash
plink --bfile hapmap3_qc \
      --keep samples_EUR_like.txt \
      --make-bed \
      --out hapmap3_EUR

wc -l hapmap3_EUR.fam
```

---

## Part V. Ancestry Prediction from PCA

PCA scores can be used to **predict the ancestry** of individuals of unknown origin by comparing their PC coordinates to reference populations with known ancestry. This is the principle behind many commercial direct-to-consumer genetic ancestry tests.

The approach here uses the 1000 Genomes individuals as a **labelled reference panel** and treats any sample we want to classify as a **query**. We train a classifier on the reference PC scores and predict superpopulation labels.

### 5.1 Concept: projection onto reference PCs

The gold-standard approach in research is:
1. Compute PCs on the **reference panel** (known ancestry)
2. **Project** query samples onto the reference PC axes (so the axes are not distorted by the query)
3. Assign ancestry based on proximity in PC space

For this exercise, because all our samples are already in the same PCA run, we can directly use the PC scores from `hapmap3_pca.eigenvec`.

### 5.2 k-Nearest Neighbours (k-NN) classifier

The simplest classifier for ancestry assignment is **k-Nearest Neighbours (k-NN)**. For each individual, it finds the $k$ closest individuals (by Euclidean distance in PC space) among the labelled reference and assigns the majority label.

Open R:

```bash
R
```

```R
library(data.table)
library(ggplot2)
library(class)   # for knn()

# Load PC scores (already computed)
pca <- fread("hapmap3_pca.eigenvec", header = FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:20))

# Load population labels
geo <- fread("1kg_samples.txt", sep = "\t", header = TRUE)
colnames(geo)[1] <- "IID"

# Merge labels onto PCA
data <- merge(pca, geo[, c("IID", "Superpopulation code")], by = "IID")
setnames(data, "Superpopulation code", "superpop")

cat("Labelled individuals:", nrow(data), "\n")
cat("Superpopulation counts:\n")
print(table(data$superpop))
```

### 5.3 Train / test split

Split the dataset: 80% training (reference), 20% testing (held out to evaluate performance):

```R
set.seed(42)
n <- nrow(data)
train_idx <- sample(1:n, size = floor(0.8 * n), replace = FALSE)
test_idx  <- setdiff(1:n, train_idx)

# Use the first 10 PCs as features
pc_cols <- paste0("PC", 1:10)

train_X <- as.matrix(data[train_idx, ..pc_cols])
test_X  <- as.matrix(data[test_idx,  ..pc_cols])
train_y <- data$superpop[train_idx]
test_y  <- data$superpop[test_idx]

cat("Training set size:", nrow(train_X), "\n")
cat("Test set size:", nrow(test_X), "\n")
```

### 5.4 Run k-NN and evaluate accuracy

```R
# k-NN classifier with k = 5
predicted <- knn(train = train_X,
                 test  = test_X,
                 cl    = train_y,
                 k     = 5)

# Confusion matrix
conf_mat <- table(Predicted = predicted, True = test_y)
print(conf_mat)

# Overall accuracy
accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
cat("Overall accuracy:", round(accuracy * 100, 1), "%\n")

# Per-population accuracy
per_pop <- diag(conf_mat) / colSums(conf_mat)
cat("\nPer-population accuracy:\n")
print(round(per_pop * 100, 1))
```

You should see very high overall accuracy (typically > 95%) reflecting the clear separation of continental populations in PC space.

### 5.5 Try different values of k

```R
# Test accuracy for k = 1, 3, 5, 10, 20
results <- data.frame(k = c(1, 3, 5, 10, 20), accuracy = NA)

for (i in seq_len(nrow(results))) {
  k_val <- results$k[i]
  pred  <- knn(train_X, test_X, train_y, k = k_val)
  cm    <- table(pred, test_y)
  results$accuracy[i] <- sum(diag(cm)) / sum(cm) * 100
}

print(results)

ggplot(results, aes(x = k, y = accuracy)) +
  geom_line() +
  geom_point(size = 3, colour = "steelblue") +
  theme_bw() +
  xlab("k (number of neighbours)") +
  ylab("Classification accuracy (%)") +
  ggtitle("k-NN ancestry classification accuracy") +
  ylim(90, 100)
ggsave("knn_accuracy.png", width = 7, height = 5)
```

### 5.6 Visualise predictions in PC space

```R
# Predict ancestry for ALL individuals using k=5
predicted_all <- knn(train = train_X,
                     test  = as.matrix(data[, ..pc_cols]),
                     cl    = train_y,
                     k     = 5)

data$predicted_superpop <- as.character(predicted_all)
data$correct <- data$predicted_superpop == data$superpop

# Plot: colour = true superpop, shape = correct/wrong
ggplot(data, aes(x = PC1, y = PC2,
                 colour = superpop,
                 shape  = correct)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_shape_manual(values = c(4, 16),
                     labels = c("Misclassified", "Correct")) +
  theme_bw() +
  labs(colour = "True superpopulation",
       shape  = "Classification",
       title  = "k-NN ancestry predictions (k = 5)") +
  theme(legend.position = "right")
ggsave("knn_predictions.png", width = 9, height = 6)

# Where do misclassifications occur?
misclass <- data[correct == FALSE, .(IID, superpop, predicted_superpop, PC1, PC2)]
cat("Misclassified individuals:", nrow(misclass), "\n")
print(head(misclass, 20))
```

Misclassifications are most common among **AMR** (Admixed American) individuals, who have mixed European, Native American, and sometimes African ancestry, making them genetically intermediate between clusters.

### 5.7 Random Forest classifier (bonus)

Random Forests typically outperform k-NN for ancestry classification because they capture non-linear boundaries.

```R
# install.packages("randomForest")  # run once if needed
library(randomForest)

# Prepare data frames
train_df <- data.frame(train_X, superpop = factor(train_y))
test_df  <- data.frame(test_X)

# Train Random Forest with 500 trees
rf_model <- randomForest(superpop ~ ., data = train_df,
                         ntree = 500, importance = TRUE)

# Predict on test set
rf_pred <- predict(rf_model, newdata = test_df)
rf_cm   <- table(Predicted = rf_pred, True = test_y)
rf_acc  <- sum(diag(rf_cm)) / sum(rf_cm) * 100

cat("Random Forest accuracy:", round(rf_acc, 1), "%\n")
print(rf_cm)

# Variable importance: which PCs matter most?
importance_df <- data.frame(
  PC       = rownames(importance(rf_model)),
  MeanDecGini = importance(rf_model)[, "MeanDecreaseGini"]
)

ggplot(importance_df, aes(x = reorder(PC, MeanDecGini),
                           y = MeanDecGini)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw() +
  xlab("Principal Component") +
  ylab("Mean Decrease in Gini") +
  ggtitle("PC importance for ancestry classification (Random Forest)")
ggsave("rf_importance.png", width = 7, height = 5)

q()
```

The variable importance plot shows which PCs contribute most to separating populations. PC1 and PC2, which capture major continental ancestry axes, typically dominate.

### Exercise 3

1. What is the k-NN classification accuracy for $k = 5$ using 10 PCs? Which superpopulation is hardest to classify and why?
2. Repeat the k-NN analysis using only PC1 and PC2 (two features instead of ten). How much does accuracy drop?
3. Look at the misclassified individuals. Which true superpopulation do they belong to, and which superpopulation are they assigned to? What does this pattern tell you about genetic admixture?
4. Compare k-NN vs Random Forest accuracy. Which classifier performs better, and by how much?

---

## Summary

Here is the full pipeline in one block for reference:

```bash
# 1. QC
plink --bfile hapmap3 \
      --mind 0.05 --geno 0.02 --maf 0.01 --hwe 1e-6 \
      --make-bed --out hapmap3_qc

# 2. HWE test
plink --bfile hapmap3_qc --hardy --out hwe_results

# 3. Association test (no covariates)
plink --bfile hapmap3_qc \
      --pheno BMI_pheno.txt --assoc --linear \
      --out bmi_assoc

# 4. LD pruning for PCA
plink --bfile hapmap3_qc \
      --indep-pairwise 50 5 0.2 \
      --out hapmap3_pruned
plink --bfile hapmap3_qc \
      --extract hapmap3_pruned.prune.in \
      --make-bed --out hapmap3_pruned_set

# 5. PCA
plink --bfile hapmap3_pruned_set --pca 20 --out hapmap3_pca

# 6. Association test corrected for 10 PCs
plink --bfile hapmap3_qc \
      --pheno BMI_pheno.txt --assoc --linear \
      --covar hapmap3_pca.eigenvec --covar-number 1-10 \
      --out bmi_assoc_pca_corrected
```

### Key concepts from this lab

| Concept | What it measures | Typical threshold / use |
|---------|-----------------|-------------------------|
| SNP missingness | Data completeness per SNP | Remove SNPs with F\_MISS > 0.02 |
| Individual missingness | Data completeness per sample | Remove individuals with F\_MISS > 0.05 |
| Inbreeding coefficient F | Heterozygosity deviation | Flag \|F\| > 0.15 |
| HWE p-value | Genotyping quality | Remove SNPs with $p < 10^{-6}$ |
| MAF | Minor allele frequency | Remove SNPs with MAF $< 0.01$ |
| GWAS p-value | Association significance | $p < 5 \times 10^{-8}$ |
| $\lambda_{GC}$ | Genomic inflation / stratification | Should be close to 1.0 |
| PC scores | Ancestry axes | Use top 10 as GWAS covariates |
| k-NN / RF classifier | Ancestry prediction | >95% accuracy for continental groups |
