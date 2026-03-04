# Lab 4. Genetic Tests and Principal Component Analysis

In this lab you will learn how to explore genetic data in R, run statistical tests, compute Principal Component Analysis (PCA) to detect population structure, and use PCA scores to predict ancestry. By the end you will be able to visualise genotype quality distributions, test Hardy-Weinberg equilibrium, run association tests, produce annotated PCA plots that reveal ancestry differences, and classify individuals into continental ancestry groups using machine learning.

We continue using the **HapMap Phase III** dataset (`hapmap3`) from the course repository.

**Tools used in this lab:**
- **[Google Cloud Shell](https://shell.cloud.google.com/)** — for all PLINK commands (bash)
- **Google Colab (R)** — for data exploration and visualisation

### Open the R notebook in Google Colab

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/nicolabarban/sociogenomics_2025_2026/blob/gh_pages/labs/week4/lab4_colab.ipynb)

> The notebook uses an **R runtime**. When it opens, Colab may prompt you to switch to an R kernel — click **Yes**. Data files are downloaded automatically from the course repository.

Alternatively, download the **[R script](lab4.R)** and run it locally in RStudio. It downloads the data automatically.

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

All R code in this lab runs in Colab cells exactly as written. The standard packages (`ggplot2`, `class`) are pre-installed.

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
    pca_AFR.eigenvec \
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

**SNP and individual missingness:**

```r
lmiss <- read.table("hapmap3_summary.lmiss", header = TRUE)
imiss <- read.table("hapmap3_summary.imiss", header = TRUE)

ggplot(lmiss, aes(x = F_MISS)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = "white") +
  xlab("Per-SNP missing rate") + ylab("Number of SNPs") +
  ggtitle("SNP missingness")

ggplot(imiss, aes(x = F_MISS)) +
  geom_histogram(bins = 40, fill = "coral", colour = "white") +
  xlab("Per-individual missing rate") + ylab("Number of individuals") +
  ggtitle("Individual missingness")
```

### 1.0b Minor allele frequency distribution

```r
frq <- read.table("hapmap3_summary.frq", header = TRUE)
cat("Mean MAF:  ", round(mean(frq$MAF), 4), "\n")
cat("Median MAF:", round(median(frq$MAF), 4), "\n")

ggplot(frq, aes(x = MAF)) +
  geom_histogram(bins = 50, fill = "darkgreen", colour = "white") +
  geom_vline(xintercept = 0.05, colour = "red", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 0.07, y = Inf, vjust = 2,
           label = "MAF = 0.05", colour = "red", size = 4) +
  xlab("Minor allele frequency") + ylab("Number of SNPs") +
  ggtitle("MAF distribution across all SNPs")
```

The red dashed line shows the common MAF = 0.05 filter. SNPs to the left of it are rare variants.

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

Within-African PCA (for section 3.5):

```bash
awk -F'\t' 'NR>1 && $6 == "AFR" {print $1, $1}' 1kg_samples.txt > samples_AFR.txt

plink --bfile hapmap3_pruned_set \
      --keep samples_AFR.txt \
      --pca 10 \
      --out pca_AFR
```

### 3.3 Scree plot in Colab (R)

Switch to your **Colab R notebook**. Data files are loaded automatically from the course repository.

```r
eigenval      <- read.table("hapmap3_pca.eigenval", header = FALSE)
colnames(eigenval) <- "eigenvalue"
eigenval$PC   <- seq_len(nrow(eigenval))
eigenval$pct  <- eigenval$eigenvalue / sum(eigenval$eigenvalue) * 100

ggplot(eigenval, aes(x = PC, y = pct)) +
  geom_col(fill = "steelblue", colour = "white") +
  geom_line(aes(group = 1)) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:20) +
  xlab("Principal Component") + ylab("Variance explained (%)") +
  ggtitle("Scree plot")

print(eigenval[, c("PC", "pct")])
```

### 3.4 PCA plots coloured by ancestry in Colab (R)

**Load PCA scores and population labels:**

```r
# Note: read.table converts spaces in column names to dots
# e.g. "Sample name" → "Sample.name", "Superpopulation name" → "Superpopulation.name"
pc_cols <- c("FID", "IID", paste0("PC", 1:20))
pca <- read.table("hapmap3_pca.eigenvec", header = FALSE, col.names = pc_cols)
geo <- read.table("1kg_samples.txt", sep = "\t", header = TRUE)

data <- merge(pca, geo[, c("Sample.name", "Population.code", "Population.name",
                             "Superpopulation.code", "Superpopulation.name")],
              by.x = "IID", by.y = "Sample.name")

cat("Individuals with labels:", nrow(data), "\n")
print(table(data$Superpopulation.name))
```

**PC1 vs PC2 coloured by superpopulation:**

```r
ggplot(data, aes(x = PC1, y = PC2, colour = Superpopulation.name)) +
  geom_point(alpha = 0.7, size = 1.5) +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "Superpopulation",
       title = "PCA coloured by continental ancestry")
```

You should see clearly separated clusters: **AFR** (African), **EUR** (European), **EAS** (East Asian), **SAS** (South Asian), **AMR** (Admixed American).

**PC1 vs PC2 coloured by sub-population:**

```r
ggplot(data, aes(x = PC1, y = PC2, colour = Population.name)) +
  geom_point(alpha = 0.7, size = 1.5) +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "Population",
       title = "PCA coloured by sub-population") +
  theme(legend.text = element_text(size = 7))
```

**PC1 vs PC3:**

```r
ggplot(data, aes(x = PC1, y = PC3, colour = Superpopulation.name)) +
  geom_point(alpha = 0.7, size = 1.5) +
  xlab("PC1") + ylab("PC3") +
  labs(colour = "Superpopulation", title = "PC1 vs PC3")
```

### 3.5 Within-population PCA (Africans)

```r
pc_cols_afr <- c("FID", "IID", paste0("PC", 1:10))
pca_afr <- read.table("pca_AFR.eigenvec", header = FALSE, col.names = pc_cols_afr)

data_afr <- merge(pca_afr, geo[, c("Sample.name", "Population.name")],
                  by.x = "IID", by.y = "Sample.name")

ggplot(data_afr, aes(x = PC1, y = PC2, colour = Population.name)) +
  geom_point(alpha = 0.8, size = 2) +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "African population",
       title = "PCA within African populations")
```

The three HapMap3 African groups — **YRI** (Yoruba, Nigeria), **LWK** (Luhya, Kenya), and **ASW** (African Americans, SW USA) — show distinct clustering. ASW individuals often appear intermediate between YRI/LWK and other continents due to admixture.

### Exercise 2

1. Examine the scree plot. How many PCs are needed to capture the main axes of variation?
2. In the PC1 vs PC2 plot, which superpopulations are most separated along PC1? Along PC2?
3. Within Africans, which populations are most separated? How does the ASW cluster position reflect their history of admixture?

---

## Part IV. Detecting and Removing Population Outliers

In a homogeneous cohort study, individuals who cluster far from the main group in PCA space likely have different ancestry and should be removed before association analysis.

### 4.1 Identify EUR-like individuals in Colab (R)

```r
eur_mean_pc1 <- mean(data$PC1[data$Superpopulation.code == "EUR"])
eur_mean_pc2 <- mean(data$PC2[data$Superpopulation.code == "EUR"])
eur_sd_pc1   <- sd(data$PC1[data$Superpopulation.code == "EUR"])
eur_sd_pc2   <- sd(data$PC2[data$Superpopulation.code == "EUR"])

data$eur_like <- abs(data$PC1 - eur_mean_pc1) < 3 * eur_sd_pc1 &
                 abs(data$PC2 - eur_mean_pc2) < 3 * eur_sd_pc2

cat("Individuals within 3 SD of EUR centroid:", sum(data$eur_like), "\n")

ggplot(data, aes(x = PC1, y = PC2,
                 colour = Superpopulation.code,
                 shape  = eur_like)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_shape_manual(values = c(4, 16),
                     labels = c("Excluded", "EUR-like (kept)")) +
  labs(colour = "Superpopulation", shape = "Selection",
       title = "EUR-like individuals (within 3 SD of EUR centroid)")

eur_keep <- data[data$eur_like, c("FID", "IID")]
write.table(eur_keep, "samples_EUR_like.txt",
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
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

pc_features <- paste0("PC", 1:10)
X <- as.matrix(data[, pc_features])
y <- data$Superpopulation.code

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
predicted_all  <- knn(train = train_X, test = X, cl = train_y, k = 5)
data$predicted <- as.character(predicted_all)
data$correct   <- data$predicted == data$Superpopulation.code

ggplot(data, aes(x = PC1, y = PC2,
                 colour = Superpopulation.code,
                 shape  = correct)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_shape_manual(values = c(4, 16),
                     labels = c("Misclassified", "Correct")) +
  labs(colour = "True superpopulation",
       shape  = "Classification",
       title  = "k-NN ancestry predictions (k = 5)")
```

Misclassifications are most common among **AMR** (Admixed American) individuals, who have mixed European, Native American, and African ancestry.

### Exercise 3

1. What is the k-NN accuracy for $k = 5$ using 10 PCs? Which superpopulation is hardest to classify and why?
2. Repeat using only PC1 and PC2 (`pc_features <- c("PC1", "PC2")`). How much does accuracy drop?
3. Look at the misclassified individuals. Which true superpopulation do they belong to, and which are they assigned to?

---

## Summary

### Cloud Shell — PLINK pipeline

```bash
# 1. QC
plink --bfile hapmap3 \
      --mind 0.05 --geno 0.02 --maf 0.01 --hwe 1e-6 \
      --make-bed --out hapmap3_qc

# 2. Summary statistics
plink --bfile hapmap3_qc --missing --freq --out hapmap3_summary

# 3. LD pruning + PCA
plink --bfile hapmap3_qc --indep-pairwise 50 5 0.2 --out hapmap3_pruned
plink --bfile hapmap3_qc --extract hapmap3_pruned.prune.in \
      --make-bed --out hapmap3_pruned_set
plink --bfile hapmap3_pruned_set --pca 20 --out hapmap3_pca

# 5. Package files for Colab
zip lab4_results.zip hapmap3_summary.lmiss hapmap3_summary.imiss \
    hapmap3_summary.frq hapmap3_pca.eigenvec hapmap3_pca.eigenval \
    pca_AFR.eigenvec 1kg_samples.txt
```

### Colab (R) — key libraries

| Library | Purpose |
|---------|---------|
| `ggplot2` | Plotting |
| `class` | k-NN classifier (`knn()`) |

### Key concepts from this lab

| Concept | What it measures | Typical threshold / use |
|---------|-----------------|-------------------------|
| SNP missingness | Data completeness per SNP | Remove F\_MISS > 0.02 |
| Individual missingness | Data completeness per sample | Remove F\_MISS > 0.05 |
| Inbreeding coefficient F | Heterozygosity deviation | Flag \|F\| > 0.15 |
| HWE p-value | Genotyping quality | Remove $p < 10^{-6}$ |
| MAF | Minor allele frequency | Remove MAF $< 0.01$ |
| GWAS p-value | Association significance | $p < 5 \times 10^{-8}$ |
| PC scores | Ancestry axes | Use top 10 as GWAS covariates |
| k-NN / RF classifier | Ancestry prediction | >95% accuracy for continental groups |
