# Lab 4. Genetic Tests and Principal Component Analysis

In this lab you will learn how to explore genetic data in R, run statistical tests, compute Principal Component Analysis (PCA) to detect population structure, and use PCA scores to predict ancestry. By the end you will be able to visualise genotype quality distributions, test Hardy-Weinberg equilibrium, run association tests, produce annotated PCA plots that reveal ancestry differences, and classify individuals into continental ancestry groups using machine learning.

We continue using the **HapMap Phase III** dataset (`hapmap3`) from the course repository.

**Tools used in this lab:**
- **[Google Cloud Shell](https://shell.cloud.google.com/)** — for all PLINK commands (bash)
- **[Google Colab (R)](https://colab.to/r)** — for data exploration and visualisation (R)

---

## 0. Getting started in Cloud Shell

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

Apply the standard QC filters from Lab 3 (skip this if `hapmap3_qc` already exists):

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

## How to use R in Google Colab

Google Colab supports R natively. To open an R notebook:

1. Go to **[https://colab.to/r](https://colab.to/r)** — this opens a new Colab notebook with an **R runtime** already selected.
2. Alternatively: open [colab.research.google.com](https://colab.research.google.com), create a new notebook, then go to **Runtime → Change runtime type → R**.

All R code in this lab runs in Colab cells exactly as written. The standard packages (`ggplot2`, `data.table`, `class`, `randomForest`) are pre-installed.

### Transferring files from Cloud Shell to Colab

PLINK output files are plain text. After running the PLINK commands, download them from Cloud Shell and upload them to Colab.

**Step 1 — In Cloud Shell:** package all output files into a zip:

```bash
cd ~/Sociogenomics/Data
zip lab4_results.zip \
    hapmap3_summary.lmiss \
    hapmap3_summary.imiss \
    hapmap3_summary.frq \
    hapmap3_summary.hwe \
    hapmap3_het.het \
    hapmap3_pca.eigenvec \
    hapmap3_pca.eigenval \
    pca_EUR.eigenvec \
    bmi_assoc.assoc.linear \
    bmi_assoc_pca_corrected.assoc.linear \
    1kg_samples.txt
```

**Step 2 — Download from Cloud Shell:** click the three-dot menu (⋮) at the top right → **Download** → type `~/Sociogenomics/Data/lab4_results.zip` → Download.

**Step 3 — Upload to Colab:** in your Colab R notebook, click the **folder icon** in the left sidebar to open the Files panel. Click the **upload icon** (↑) and select `lab4_results.zip`. Then unzip it with this R cell:

```r
unzip("lab4_results.zip")
list.files()   # confirm files are present
```

All subsequent R code assumes files are in the Colab working directory (`/content/`).

---

## Part I. Exploring Genetic Data in R (Colab)

### 1.0 Generate summary statistics in Cloud Shell

Run these PLINK commands in **Cloud Shell** first:

```bash
cd ~/Sociogenomics/Data

# Per-SNP and per-individual missingness
plink --bfile hapmap3_qc --missing --out hapmap3_summary

# Allele frequencies
plink --bfile hapmap3_qc --freq --out hapmap3_summary

# Hardy-Weinberg equilibrium
plink --bfile hapmap3_qc --hardy --out hapmap3_summary

# Per-individual inbreeding coefficient
plink --bfile hapmap3_qc --het --out hapmap3_het
```

Then transfer the files to Colab as described above.

### 1.0a Load and inspect summary statistics

Open your **Colab R notebook** and run the following cells.

**Install / load libraries** (first cell):

```r
library(data.table)
library(ggplot2)
```

**SNP and individual missingness:**

```r
# --- SNP missingness ---
lmiss <- fread("hapmap3_summary.lmiss")
cat("SNPs in dataset:", nrow(lmiss), "\n")
cat("SNP missingness range:", range(lmiss$F_MISS), "\n")

p1 <- ggplot(lmiss, aes(x = F_MISS)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = "white") +
  theme_bw() +
  xlab("Per-SNP missing rate") + ylab("Number of SNPs") +
  ggtitle("SNP missingness")

# --- Individual missingness ---
imiss <- fread("hapmap3_summary.imiss")
cat("Individuals in dataset:", nrow(imiss), "\n")
cat("Individual missingness range:", range(imiss$F_MISS), "\n")

p2 <- ggplot(imiss, aes(x = F_MISS)) +
  geom_histogram(bins = 40, fill = "coral", colour = "white") +
  theme_bw() +
  xlab("Per-individual missing rate") + ylab("Number of individuals") +
  ggtitle("Individual missingness")

library(patchwork)
p1 + p2
```

### 1.0b Minor allele frequency distribution

```r
frq <- fread("hapmap3_summary.frq")
cat("Mean MAF:  ", round(mean(frq$MAF), 4), "\n")
cat("Median MAF:", round(median(frq$MAF), 4), "\n")

ggplot(frq, aes(x = MAF)) +
  geom_histogram(bins = 50, fill = "darkgreen", colour = "white") +
  geom_vline(xintercept = 0.05, colour = "red", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 0.07, y = Inf, vjust = 2,
           label = "MAF = 0.05", colour = "red", size = 4) +
  theme_bw() +
  xlab("Minor allele frequency") + ylab("Number of SNPs") +
  ggtitle("MAF distribution across all SNPs")
```

The red dashed line shows the common MAF = 0.05 filter. SNPs to the left of it are rare variants.

### 1.0c Hardy-Weinberg equilibrium distribution

```r
hwe <- fread("hapmap3_summary.hwe")
hwe_all <- hwe[TEST == "ALL"]
hwe_all[, log10p := -log10(P)]

ggplot(hwe_all, aes(x = log10p)) +
  geom_histogram(bins = 60, fill = "purple", colour = "white") +
  geom_vline(xintercept = 6, colour = "red", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 6.3, y = Inf, vjust = 2,
           label = "p = 1e-6", colour = "red", size = 4) +
  theme_bw() +
  xlab(expression(-log[10](p))) + ylab("Number of SNPs") +
  ggtitle("Distribution of HWE test statistics")

cat("SNPs failing HWE (p < 1e-6):", sum(hwe_all$P < 1e-6, na.rm = TRUE), "\n")
cat("\nTop 10 most deviant SNPs:\n")
print(hwe_all[order(P)][1:10, .(SNP, CHR, `O(HET)`, `E(HET)`, P)])
```

### 1.0d Per-individual inbreeding coefficient

Individuals with extreme inbreeding coefficients (F) may indicate contamination (very negative F) or inbreeding (very positive F).

```r
het <- fread("hapmap3_het.het")

ggplot(het, aes(x = F)) +
  geom_histogram(bins = 50, fill = "orange", colour = "white") +
  geom_vline(xintercept = c(-0.15, 0.15), colour = "red",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = -0.17, y = Inf, vjust = 2, hjust = 1,
           label = "-0.15", colour = "red", size = 4) +
  annotate("text", x =  0.17, y = Inf, vjust = 2, hjust = 0,
           label = "+0.15", colour = "red", size = 4) +
  theme_bw() +
  xlab("Inbreeding coefficient F") + ylab("Number of individuals") +
  ggtitle("Per-individual inbreeding coefficient")

outliers <- het[F < -0.15 | F > 0.15]
cat("Heterozygosity outliers:", nrow(outliers), "\n")
if (nrow(outliers) > 0) print(outliers[, .(FID, IID, F)])
```

---

## Part II. Statistical Tests on Genetic Data (Cloud Shell)

All commands in this part run in **Cloud Shell**.

### 2.1 Hardy-Weinberg Equilibrium (HWE) test

**Hardy-Weinberg Equilibrium** (HWE) states that in a large random-mating population with no selection, mutation, or migration, allele and genotype frequencies remain constant across generations. For a biallelic SNP with allele frequencies $p$ (allele A) and $q = 1 - p$ (allele B), the expected genotype frequencies are:

$$P(AA) = p^2, \quad P(AB) = 2pq, \quad P(BB) = q^2$$

**Why test for HWE?** Deviations from HWE usually indicate genotyping errors — differential allelic dropout, probe failure, or sample contamination. SNPs with $p < 10^{-6}$ in controls are typically removed in QC.

HWE is tested using a chi-squared test with 1 degree of freedom:

$$\chi^2 = \frac{(O_{AA} - E_{AA})^2}{E_{AA}} + \frac{(O_{AB} - E_{AB})^2}{E_{AB}} + \frac{(O_{BB} - E_{BB})^2}{E_{BB}}$$

```bash
plink --bfile hapmap3_qc --hardy --out hwe_results
head hwe_results.hwe
```

The columns are: CHR, SNP, TEST (ALL / AFF / UNAFF), A1, A2, GENO, O(HET), E(HET), P.

Count SNPs at different thresholds:

```bash
awk 'NR>1 && $9 < 0.05' hwe_results.hwe | wc -l    # p < 0.05
awk 'NR>1 && $9 < 1e-6' hwe_results.hwe | wc -l    # p < 1e-6 (QC threshold)
```

Find the most deviant SNPs:

```bash
awk 'NR>1 {print $2, $9}' hwe_results.hwe | sort -k2 -n | head -10
```

### 2.2 Allele frequency statistics

```bash
plink --bfile hapmap3_qc --freq --out allele_freq
head allele_freq.frq
```

Summarise by frequency bin:

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

### 2.3 Association test

For a **quantitative trait** (like BMI), PLINK fits a linear regression:

$$y = \beta_0 + \beta_1 x + \epsilon$$

where $y$ is the phenotype, $x$ is the genotype (0, 1, 2 copies of the effect allele), and $\beta_1$ is the additive effect estimate.

```bash
plink --bfile hapmap3_qc \
      --pheno BMI_pheno.txt \
      --assoc --linear \
      --out bmi_assoc

head bmi_assoc.assoc.linear
```

The key columns are: CHR, SNP, BP, A1, TEST, NMISS, BETA, STAT, P.

```bash
# Top associations
awk 'NR>1' bmi_assoc.assoc.linear | sort -k9 -n | head -20

# Genome-wide significant hits (p < 5e-8)
awk 'NR>1 && $9 < 5e-8' bmi_assoc.assoc.linear | wc -l
```

### 2.4 Test a single SNP under different genetic models

```bash
TOP_SNP=$(awk 'NR>1' bmi_assoc.assoc.linear | sort -k9 -n | head -1 | awk '{print $2}')
echo "Top SNP: $TOP_SNP"

# Additive model
plink --bfile hapmap3_qc --pheno BMI_pheno.txt \
      --snp $TOP_SNP --assoc --linear --out top_snp_additive

# Dominant model (one copy sufficient)
plink --bfile hapmap3_qc --pheno BMI_pheno.txt \
      --snp $TOP_SNP --assoc --linear dominant --out top_snp_dominant

# Recessive model (two copies required)
plink --bfile hapmap3_qc --pheno BMI_pheno.txt \
      --snp $TOP_SNP --assoc --linear recessive --out top_snp_recessive

cat top_snp_additive.assoc.linear
cat top_snp_dominant.assoc.linear
cat top_snp_recessive.assoc.linear
```

**Question:** Does the p-value change across models? Which model fits best?

### Exercise 1

1. How many SNPs deviate significantly from HWE at $p < 10^{-6}$? What could cause this?
2. Do the most deviant SNPs show excess or deficit heterozygosity?
3. Report the top 5 BMI-associated SNPs. What are their chromosomal positions?

---

## Part III. Principal Component Analysis (PCA)

### 3.1 Why PCA in genetics?

PCA is used in genomics to:
1. **Detect population stratification** — systematic ancestry differences between cases and controls that confound association tests
2. **Visualise ancestry** — reveal continental and sub-continental genetic clusters
3. **Create covariates** — PC scores are added to GWAS regression models to control for stratification

PCA decomposes the genotype matrix $\mathbf{G}$ (individuals × SNPs) into orthogonal components that capture the directions of maximum variance. Individuals with similar ancestry cluster together in PC space because they share similar allele frequencies genome-wide.

**Genomic PCA requires LD-pruned SNPs.** Regions in high LD would over-represent certain loci, distorting the PC directions.

### 3.2 Compute PCA in Cloud Shell

```bash
# LD pruning
plink --bfile hapmap3_qc \
      --indep-pairwise 50 5 0.2 \
      --out hapmap3_pruned

plink --bfile hapmap3_qc \
      --extract hapmap3_pruned.prune.in \
      --make-bed \
      --out hapmap3_pruned_set

wc -l hapmap3_pruned.prune.in   # how many SNPs remain?

# Compute PCA (top 20 components)
plink --bfile hapmap3_pruned_set \
      --pca 20 \
      --out hapmap3_pca

head hapmap3_pca.eigenvec
head hapmap3_pca.eigenval
```

Proportion of variance explained by each PC:

```bash
awk 'BEGIN{sum=0} {val[NR]=$1; sum+=$1}
     END{
       for(i=1; i<=NR; i++)
         printf "PC%d: %.2f%%\n", i, val[i]/sum*100
     }' hapmap3_pca.eigenval
```

Within-European PCA (for section 3.5):

```bash
awk -F'\t' 'NR>1 && $6 == "EUR" {print $1, $1}' 1kg_samples.txt > samples_EUR.txt

plink --bfile hapmap3_pruned_set \
      --keep samples_EUR.txt \
      --pca 10 \
      --out pca_EUR
```

### 3.3 Scree plot in Colab (R)

Switch to your **Colab R notebook**. Upload `lab4_results.zip` if not already done.

```r
library(data.table)
library(ggplot2)

eigenval <- fread("hapmap3_pca.eigenval", header = FALSE, col.names = "eigenvalue")
eigenval[, PC  := seq_len(.N)]
eigenval[, pct := eigenvalue / sum(eigenvalue) * 100]

ggplot(eigenval, aes(x = PC, y = pct)) +
  geom_col(fill = "steelblue", colour = "white") +
  geom_line(aes(group = 1)) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:20) +
  theme_bw() +
  xlab("Principal Component") + ylab("Variance explained (%)") +
  ggtitle("Scree plot")

print(eigenval[, .(PC, pct = round(pct, 2))])
```

### 3.4 PCA plots coloured by ancestry in Colab (R)

**Load PCA scores and population labels:**

```r
pc_cols <- c("FID", "IID", paste0("PC", 1:20))
pca <- fread("hapmap3_pca.eigenvec", header = FALSE, col.names = pc_cols)

geo <- fread("1kg_samples.txt", sep = "\t", header = TRUE)
setnames(geo, "Sample name", "IID")

data <- merge(pca, geo[, .(`IID`, `Population code`, `Population name`,
                            `Superpopulation code`, `Superpopulation name`)],
              by = "IID")

cat("Individuals with labels:", nrow(data), "\n")
print(table(data$`Superpopulation name`))
```

**PC1 vs PC2 coloured by superpopulation:**

```r
ggplot(data, aes(x = PC1, y = PC2, colour = `Superpopulation name`)) +
  geom_point(alpha = 0.7, size = 1.5) +
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "Superpopulation",
       title = "PCA coloured by continental ancestry")
```

You should see clearly separated clusters: **AFR** (African), **EUR** (European), **EAS** (East Asian), **SAS** (South Asian), **AMR** (Admixed American).

**PC1 vs PC2 coloured by sub-population:**

```r
ggplot(data, aes(x = PC1, y = PC2, colour = `Population name`)) +
  geom_point(alpha = 0.7, size = 1.5) +
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "Population",
       title = "PCA coloured by sub-population") +
  theme(legend.text = element_text(size = 7))
```

**PC1 vs PC3:**

```r
ggplot(data, aes(x = PC1, y = PC3, colour = `Superpopulation name`)) +
  geom_point(alpha = 0.7, size = 1.5) +
  theme_bw() +
  xlab("PC1") + ylab("PC3") +
  labs(colour = "Superpopulation", title = "PC1 vs PC3")
```

### 3.5 Within-population PCA (Europeans)

```r
pc_cols_eur <- c("FID", "IID", paste0("PC", 1:10))
pca_eur <- fread("pca_EUR.eigenvec", header = FALSE, col.names = pc_cols_eur)

data_eur <- merge(pca_eur, geo[, .(`IID`, `Population name`)], by = "IID")

ggplot(data_eur, aes(x = PC1, y = PC2, colour = `Population name`)) +
  geom_point(alpha = 0.8, size = 2) +
  theme_bw() +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "European population",
       title = "PCA within European populations")
```

Within Europeans, PCs often separate Northern Europeans (Finnish, British) from Southern Europeans (Iberian, Tuscan).

### 3.6 PCA as covariates in GWAS (Cloud Shell)

```bash
plink --bfile hapmap3_qc \
      --pheno BMI_pheno.txt \
      --assoc --linear \
      --covar hapmap3_pca.eigenvec \
      --covar-number 1-10 \
      --out bmi_assoc_pca_corrected

# Compare hits before and after
awk 'NR>1 && $9 < 5e-8' bmi_assoc.assoc.linear | wc -l
awk 'NR>1 && $9 < 5e-8' bmi_assoc_pca_corrected.assoc.linear | wc -l
```

### 3.7 Genomic inflation factor λ_GC in Colab (R)

The genomic inflation factor $\lambda_{GC}$ measures residual stratification:

$$\lambda_{GC} = \frac{\text{median}(\chi^2_{\text{observed}})}{0.4549}$$

A value close to 1.0 indicates no inflation.

```r
compute_lambda <- function(pvals) {
  pvals <- pvals[!is.na(pvals)]
  chisq <- qchisq(pvals, df = 1, lower.tail = FALSE)
  median(chisq) / 0.4549
}

# Uncorrected
res_raw <- fread("bmi_assoc.assoc.linear")
lam_raw <- compute_lambda(res_raw$P)
cat("Lambda (uncorrected): ", round(lam_raw, 3), "\n")

# PC-corrected (keep only the ADD test rows)
res_corr <- fread("bmi_assoc_pca_corrected.assoc.linear")
res_corr <- res_corr[TEST == "ADD"]
lam_corr <- compute_lambda(res_corr$P)
cat("Lambda (PC-corrected):", round(lam_corr, 3), "\n")
```

### Exercise 2

1. Examine the scree plot. How many PCs are needed to capture the main axes of variation?
2. In the PC1 vs PC2 plot, which superpopulations are most separated along PC1? Along PC2?
3. Within Europeans, which populations are most separated? What historical events might explain this?
4. How does $\lambda_{GC}$ change before and after adding 10 PCs as covariates?

---

## Part IV. Detecting and Removing Population Outliers

In a homogeneous cohort study, individuals who cluster far from the main group in PCA space likely have different ancestry and should be removed before association analysis.

### 4.1 Identify EUR-like individuals in Colab (R)

```r
# Compute the EUR centroid and standard deviations
eur_mean_pc1 <- mean(data[`Superpopulation code` == "EUR", PC1])
eur_mean_pc2 <- mean(data[`Superpopulation code` == "EUR", PC2])
eur_sd_pc1   <- sd(data[`Superpopulation code` == "EUR", PC1])
eur_sd_pc2   <- sd(data[`Superpopulation code` == "EUR", PC2])

cat("EUR centroid: PC1 =", round(eur_mean_pc1, 4),
    ", PC2 =", round(eur_mean_pc2, 4), "\n")

# Flag individuals within 3 SD of the EUR centroid
data[, eur_like := abs(PC1 - eur_mean_pc1) < 3 * eur_sd_pc1 &
                   abs(PC2 - eur_mean_pc2) < 3 * eur_sd_pc2]

cat("Individuals within 3 SD of EUR centroid:", sum(data$eur_like), "\n")

# Visualise the selection
ggplot(data, aes(x = PC1, y = PC2,
                 colour = `Superpopulation code`,
                 shape  = eur_like)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_shape_manual(values = c(4, 16),
                     labels = c("Excluded", "EUR-like (kept)")) +
  theme_bw() +
  labs(colour = "Superpopulation", shape = "Selection",
       title = "EUR-like individuals (within 3 SD of EUR centroid)")

# Save the sample list for PLINK
eur_keep <- data[eur_like == TRUE, .(FID, IID)]
fwrite(eur_keep, "samples_EUR_like.txt", sep = " ", col.names = FALSE)
cat("Saved", nrow(eur_keep), "EUR-like individuals to samples_EUR_like.txt\n")
```

**Download `samples_EUR_like.txt` from Colab:** in the Files panel (left sidebar), right-click the file → **Download**. Then upload it to Cloud Shell using the ⋮ menu → **Upload**, and apply the filter:

```bash
plink --bfile hapmap3_qc \
      --keep samples_EUR_like.txt \
      --make-bed \
      --out hapmap3_EUR

wc -l hapmap3_EUR.fam
```

---

## Part V. Ancestry Prediction from PCA (Colab — R)

PCA scores can be used to **predict the ancestry** of individuals of unknown origin. This is the principle behind commercial genetic ancestry tests (e.g. 23andMe, AncestryDNA).

We use the 1000 Genomes individuals as a **labelled reference panel** and classify them using their PC coordinates.

### 5.1 Concept: projection onto reference PCs

The gold-standard research approach is:
1. Compute PCs on the **reference panel** (known ancestry)
2. **Project** query samples onto the reference PC axes
3. Assign ancestry based on proximity in PC space

For this exercise, all samples are already in the same PCA run, so we can use the PC scores from `hapmap3_pca.eigenvec` directly.

### 5.2 k-Nearest Neighbours (k-NN) classifier

The simplest classifier for ancestry assignment is **k-Nearest Neighbours (k-NN)**. For each individual it finds the $k$ closest individuals (by Euclidean distance in PC space) among the labelled reference and assigns the majority label.

```r
library(class)   # for knn()

# Features: first 10 PCs
pc_features <- paste0("PC", 1:10)

X <- as.matrix(data[, ..pc_features])
y <- data$`Superpopulation code`

cat("Total individuals:", nrow(X), "\n")
print(table(y))
```

### 5.3 Train / test split

```r
set.seed(42)
n         <- nrow(data)
train_idx <- sample(n, size = floor(0.8 * n), replace = FALSE)
test_idx  <- setdiff(seq_len(n), train_idx)

train_X <- X[train_idx, ]
test_X  <- X[test_idx,  ]
train_y <- y[train_idx]
test_y  <- y[test_idx]

cat("Training set:", nrow(train_X), "individuals\n")
cat("Test set:    ", nrow(test_X),  "individuals\n")
```

### 5.4 Run k-NN and evaluate accuracy

```r
predicted <- knn(train = train_X, test = test_X, cl = train_y, k = 5)

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

You should see overall accuracy > 95%, reflecting the clear separation of continental populations in PC space.

### 5.5 Effect of k on accuracy

```r
k_values   <- c(1, 3, 5, 10, 20)
accuracies <- numeric(length(k_values))

for (i in seq_along(k_values)) {
  pred           <- knn(train_X, test_X, train_y, k = k_values[i])
  cm             <- table(pred, test_y)
  accuracies[i]  <- sum(diag(cm)) / sum(cm) * 100
}

results <- data.frame(k = k_values, accuracy = accuracies)
print(results)

ggplot(results, aes(x = k, y = accuracy)) +
  geom_line() +
  geom_point(size = 3, colour = "steelblue") +
  theme_bw() +
  xlab("k (number of neighbours)") +
  ylab("Classification accuracy (%)") +
  ggtitle("k-NN ancestry classification accuracy") +
  ylim(90, 100)
```

### 5.6 Visualise predictions in PC space

```r
# Predict ancestry for ALL individuals (k = 5)
predicted_all <- knn(train = train_X,
                     test  = X,
                     cl    = train_y,
                     k     = 5)

data[, predicted   := as.character(predicted_all)]
data[, correct     := predicted == `Superpopulation code`]

ggplot(data, aes(x = PC1, y = PC2,
                 colour = `Superpopulation code`,
                 shape  = correct)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_shape_manual(values = c(4, 16),
                     labels = c("Misclassified", "Correct")) +
  theme_bw() +
  labs(colour = "True superpopulation",
       shape  = "Classification",
       title  = "k-NN ancestry predictions (k = 5)")

misclass <- data[correct == FALSE, .(`IID`, `Superpopulation code`, predicted, PC1, PC2)]
cat("Misclassified individuals:", nrow(misclass), "\n")
print(head(misclass, 20))
```

Misclassifications are most common among **AMR** (Admixed American) individuals, who have mixed European, Native American, and African ancestry.

### 5.7 Random Forest classifier (bonus)

Random Forests typically outperform k-NN because they capture non-linear decision boundaries.

```r
library(randomForest)

train_df <- data.frame(train_X, superpop = factor(train_y))
test_df  <- data.frame(test_X)

rf_model <- randomForest(superpop ~ ., data = train_df,
                         ntree = 500, importance = TRUE)

rf_pred <- predict(rf_model, newdata = test_df)
rf_cm   <- table(Predicted = rf_pred, True = test_y)
rf_acc  <- sum(diag(rf_cm)) / sum(rf_cm) * 100

cat("Random Forest accuracy:", round(rf_acc, 1), "%\n")
cat("k-NN accuracy:         ", round(accuracy * 100, 1), "%\n")
print(rf_cm)

# Variable importance: which PCs matter most?
imp <- data.frame(
  PC          = rownames(importance(rf_model)),
  MeanDecGini = importance(rf_model)[, "MeanDecreaseGini"]
)

ggplot(imp, aes(x = reorder(PC, MeanDecGini), y = MeanDecGini)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw() +
  xlab("Principal Component") +
  ylab("Mean Decrease in Gini") +
  ggtitle("PC importance for ancestry classification\n(Random Forest)")
```

PC1 and PC2, which capture major continental ancestry axes, typically dominate the importance ranking.

### Exercise 3

1. What is the k-NN accuracy for $k = 5$ using 10 PCs? Which superpopulation is hardest to classify and why?
2. Repeat using only PC1 and PC2 (`pc_features <- c("PC1", "PC2")`). How much does accuracy drop?
3. Look at the misclassified individuals. Which true superpopulation do they belong to, and which are they assigned to?
4. Compare k-NN vs Random Forest accuracy. Which performs better?

---

## Summary

### Cloud Shell — PLINK pipeline

```bash
# 1. QC
plink --bfile hapmap3 \
      --mind 0.05 --geno 0.02 --maf 0.01 --hwe 1e-6 \
      --make-bed --out hapmap3_qc

# 2. Summary statistics
plink --bfile hapmap3_qc --missing --freq --hardy --out hapmap3_summary
plink --bfile hapmap3_qc --het --out hapmap3_het

# 3. Association test
plink --bfile hapmap3_qc --pheno BMI_pheno.txt --assoc --linear --out bmi_assoc

# 4. LD pruning + PCA
plink --bfile hapmap3_qc --indep-pairwise 50 5 0.2 --out hapmap3_pruned
plink --bfile hapmap3_qc --extract hapmap3_pruned.prune.in \
      --make-bed --out hapmap3_pruned_set
plink --bfile hapmap3_pruned_set --pca 20 --out hapmap3_pca

# 5. PCA-corrected association
plink --bfile hapmap3_qc --pheno BMI_pheno.txt --assoc --linear \
      --covar hapmap3_pca.eigenvec --covar-number 1-10 \
      --out bmi_assoc_pca_corrected

# 6. Package files for Colab
zip lab4_results.zip hapmap3_summary.* hapmap3_het.het \
    hapmap3_pca.eigenvec hapmap3_pca.eigenval pca_EUR.eigenvec \
    bmi_assoc.assoc.linear bmi_assoc_pca_corrected.assoc.linear \
    1kg_samples.txt
```

### Colab (R) — key libraries

| Library | Purpose |
|---------|---------|
| `data.table` | Fast reading of large PLINK output files |
| `ggplot2` | Plotting |
| `patchwork` | Combining multiple plots |
| `class` | k-NN classifier (`knn()`) |
| `randomForest` | Random Forest classifier |

### Key concepts from this lab

| Concept | What it measures | Typical threshold / use |
|---------|-----------------|-------------------------|
| SNP missingness | Data completeness per SNP | Remove F\_MISS > 0.02 |
| Individual missingness | Data completeness per sample | Remove F\_MISS > 0.05 |
| Inbreeding coefficient F | Heterozygosity deviation | Flag \|F\| > 0.15 |
| HWE p-value | Genotyping quality | Remove $p < 10^{-6}$ |
| MAF | Minor allele frequency | Remove MAF $< 0.01$ |
| GWAS p-value | Association significance | $p < 5 \times 10^{-8}$ |
| $\lambda_{GC}$ | Genomic inflation / stratification | Should be close to 1.0 |
| PC scores | Ancestry axes | Use top 10 as GWAS covariates |
| k-NN / RF classifier | Ancestry prediction | >95% accuracy for continental groups |
