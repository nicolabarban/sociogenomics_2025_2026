# Lab 4. Genetic Tests and Principal Component Analysis

In this lab you will learn how to explore genetic data in Python, run statistical tests, compute Principal Component Analysis (PCA) to detect population structure, and use PCA scores to predict ancestry. By the end you will be able to visualise genotype quality distributions, test Hardy-Weinberg equilibrium, run association tests, produce annotated PCA plots that reveal ancestry differences, and classify individuals into continental ancestry groups using machine learning.

We continue using the **HapMap Phase III** dataset (`hapmap3`) from the course repository.

**Tools used in this lab:**
- **[Google Cloud Shell](https://shell.cloud.google.com/)** — for all PLINK commands (bash)
- **[Google Colab](https://colab.research.google.com/)** — for data exploration and visualisation (Python)

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

## How to transfer files from Cloud Shell to Google Colab

All PLINK output files are plain text. After running the PLINK commands below, you will download them to your computer and upload them to Colab.

**Step 1 — In Cloud Shell:** package all results into a zip file:

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

**Step 2 — Download from Cloud Shell:** click the three-dot menu (⋮) at the top right of the Cloud Shell panel → **Download** → type `~/Sociogenomics/Data/lab4_results.zip` → click Download.

**Step 3 — In Google Colab:** upload the zip and extract it:

```python
# Run this cell first in every Colab session
from google.colab import files
import zipfile, os

uploaded = files.upload()          # select lab4_results.zip from your computer
with zipfile.ZipFile("lab4_results.zip", "r") as z:
    z.extractall(".")

print("Files available:", os.listdir("."))
```

All subsequent Python code in this lab assumes files are in the current Colab working directory.

---

## Part I. Exploring Genetic Data in Python (Colab)

### 1.0 Generate summary statistics in Cloud Shell

Run these PLINK commands in **Cloud Shell** first, then transfer the output files to Colab.

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

### 1.0a Load and inspect summary statistics in Colab

Open [Google Colab](https://colab.research.google.com/), create a new notebook, and run the following cells.

**Install and import libraries** (first cell):

```python
# All libraries below come pre-installed in Colab
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

sns.set_theme(style="whitegrid", font_scale=1.2)
```

**SNP and individual missingness:**

```python
# --- SNP missingness ---
lmiss = pd.read_csv("hapmap3_summary.lmiss", sep=r"\s+")
print(f"SNPs in dataset: {len(lmiss)}")
print(f"SNP missingness range: {lmiss['F_MISS'].min():.4f} – {lmiss['F_MISS'].max():.4f}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].hist(lmiss["F_MISS"], bins=50, color="steelblue", edgecolor="white")
axes[0].set_xlabel("Per-SNP missing rate")
axes[0].set_ylabel("Number of SNPs")
axes[0].set_title("SNP missingness")

# --- Individual missingness ---
imiss = pd.read_csv("hapmap3_summary.imiss", sep=r"\s+")
print(f"Individuals in dataset: {len(imiss)}")
print(f"Individual missingness range: {imiss['F_MISS'].min():.4f} – {imiss['F_MISS'].max():.4f}")

axes[1].hist(imiss["F_MISS"], bins=40, color="coral", edgecolor="white")
axes[1].set_xlabel("Per-individual missing rate")
axes[1].set_ylabel("Number of individuals")
axes[1].set_title("Individual missingness")

plt.tight_layout()
plt.savefig("hist_missingness.png", dpi=150)
plt.show()
```

### 1.0b Minor allele frequency distribution

```python
frq = pd.read_csv("hapmap3_summary.frq", sep=r"\s+")
print(f"Mean MAF:   {frq['MAF'].mean():.4f}")
print(f"Median MAF: {frq['MAF'].median():.4f}")

plt.figure(figsize=(8, 5))
plt.hist(frq["MAF"], bins=50, color="darkgreen", edgecolor="white")
plt.axvline(0.05, color="red", linestyle="--", linewidth=1.5, label="MAF = 0.05")
plt.xlabel("Minor allele frequency")
plt.ylabel("Number of SNPs")
plt.title("MAF distribution across all SNPs")
plt.legend()
plt.tight_layout()
plt.savefig("hist_maf.png", dpi=150)
plt.show()

# Frequency bins
bins = [0, 0.01, 0.05, 0.10, 0.20, 0.51]
labels = ["<0.01 (rare)", "0.01–0.05", "0.05–0.10", "0.10–0.20", "≥0.20 (common)"]
frq["bin"] = pd.cut(frq["MAF"], bins=bins, labels=labels, right=False)
print("\nSNPs per MAF bin:")
print(frq["bin"].value_counts().sort_index())
```

### 1.0c Hardy-Weinberg equilibrium distribution

```python
hwe = pd.read_csv("hapmap3_summary.hwe", sep=r"\s+")
hwe_all = hwe[hwe["TEST"] == "ALL"].copy()

hwe_all["log10p"] = -np.log10(hwe_all["P"])

plt.figure(figsize=(8, 5))
plt.hist(hwe_all["log10p"], bins=60, color="purple", edgecolor="white")
plt.axvline(6, color="red", linestyle="--", linewidth=1.5, label="p = 1e-6 threshold")
plt.xlabel(r"$-\log_{10}(p)$")
plt.ylabel("Number of SNPs")
plt.title("Distribution of HWE test statistics")
plt.legend()
plt.tight_layout()
plt.savefig("hist_hwe.png", dpi=150)
plt.show()

n_fail = (hwe_all["P"] < 1e-6).sum()
print(f"SNPs failing HWE (p < 1e-6): {n_fail}")
print(f"Top 10 most deviant SNPs:")
print(hwe_all.nsmallest(10, "P")[["SNP", "CHR", "O(HET)", "E(HET)", "P"]])
```

### 1.0d Per-individual inbreeding coefficient

Individuals with extreme inbreeding coefficients (F) deviate from the population average, which may indicate contamination (F very negative) or inbreeding (F very positive).

```python
het = pd.read_csv("hapmap3_het.het", sep=r"\s+")
# F = inbreeding coefficient (already in the file)

plt.figure(figsize=(8, 5))
plt.hist(het["F"], bins=50, color="orange", edgecolor="white")
plt.axvline(-0.15, color="red", linestyle="--", linewidth=1.5, label="±0.15 threshold")
plt.axvline( 0.15, color="red", linestyle="--", linewidth=1.5)
plt.xlabel("Inbreeding coefficient F")
plt.ylabel("Number of individuals")
plt.title("Per-individual inbreeding coefficient")
plt.legend()
plt.tight_layout()
plt.savefig("hist_inbreeding.png", dpi=150)
plt.show()

outliers = het[(het["F"] < -0.15) | (het["F"] > 0.15)]
print(f"Heterozygosity outliers: {len(outliers)}")
if len(outliers) > 0:
    print(outliers[["FID", "IID", "F"]])
```

---

## Part II. Statistical Tests on Genetic Data (Cloud Shell)

All commands in this part run in **Cloud Shell**.

### 2.1 Hardy-Weinberg Equilibrium (HWE) test

**Hardy-Weinberg Equilibrium** (HWE) states that in a large random-mating population with no selection, mutation, or migration, allele and genotype frequencies remain constant across generations. For a biallelic SNP with allele frequencies $p$ (allele A) and $q = 1 - p$ (allele B), the expected genotype frequencies are:

$$P(AA) = p^2, \quad P(AB) = 2pq, \quad P(BB) = q^2$$

**Why test for HWE?** Deviations from HWE in a control sample usually indicate genotyping errors — differential allelic dropout, probe failure, or sample contamination — rather than a true biological signal. SNPs with $p < 10^{-6}$ in controls are typically removed in QC.

HWE is tested using a chi-squared test with 1 degree of freedom:

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

The columns are: CHR, SNP, TEST (ALL / AFF / UNAFF), A1, A2, GENO (observed counts: homozygous A1 / heterozygous / homozygous A2), O(HET), E(HET), P.

Count SNPs at different significance thresholds:

```bash
# How many SNPs deviate at p < 0.05?
awk 'NR>1 && $9 < 0.05' hwe_results.hwe | wc -l

# How many at p < 1e-6 (typical QC threshold)?
awk 'NR>1 && $9 < 1e-6' hwe_results.hwe | wc -l
```

Find the SNPs with the strongest HWE deviation:

```bash
awk 'NR>1 {print $2, $9}' hwe_results.hwe | sort -k2 -n | head -10
```

### 2.2 Allele frequency statistics

Compute allele frequencies:

```bash
plink --bfile hapmap3_qc \
      --freq \
      --out allele_freq
head allele_freq.frq
```

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

For a **quantitative trait** (like BMI), a linear regression is used:

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
awk 'NR>1' bmi_assoc.assoc.linear | sort -k9 -n | head -20
```

Count genome-wide significant hits ($p < 5 \times 10^{-8}$, the standard GWAS threshold):

```bash
awk 'NR>1 && $9 < 5e-8' bmi_assoc.assoc.linear | wc -l
```

The GWAS significance threshold of $p < 5 \times 10^{-8}$ is a Bonferroni correction for approximately 1 million independent tests across the genome.

### 2.4 Test a single SNP under different genetic models

```bash
# Identify the top SNP
TOP_SNP=$(awk 'NR>1' bmi_assoc.assoc.linear | sort -k9 -n | head -1 | awk '{print $2}')
echo "Top SNP: $TOP_SNP"

# Additive model
plink --bfile hapmap3_qc --pheno BMI_pheno.txt \
      --snp $TOP_SNP --assoc --linear \
      --out top_snp_additive
cat top_snp_additive.assoc.linear

# Dominant model (one copy of the effect allele is sufficient)
plink --bfile hapmap3_qc --pheno BMI_pheno.txt \
      --snp $TOP_SNP --assoc --linear dominant \
      --out top_snp_dominant
cat top_snp_dominant.assoc.linear

# Recessive model (two copies required)
plink --bfile hapmap3_qc --pheno BMI_pheno.txt \
      --snp $TOP_SNP --assoc --linear recessive \
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

PCA works by decomposing the genotype matrix $\mathbf{G}$ (individuals × SNPs) into orthogonal components that capture the directions of maximum variance. Individuals with similar ancestry cluster together in PC space because they share similar allele frequencies genome-wide.

**Genomic PCA requires LD-pruned SNPs.** Regions in high LD would over-represent certain loci, distorting the PC directions.

### 3.2 Compute PCA in Cloud Shell

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

# How many SNPs remain after pruning?
wc -l hapmap3_pruned.prune.in

# Compute PCA (top 20 components)
plink --bfile hapmap3_pruned_set \
      --pca 20 \
      --out hapmap3_pca
```

PLINK produces two output files:
- `hapmap3_pca.eigenvec` — PC scores for each individual
- `hapmap3_pca.eigenval` — eigenvalues (variance explained by each PC)

Inspect the output:

```bash
head hapmap3_pca.eigenvec
head hapmap3_pca.eigenval
```

Compute the proportion of variance explained by each PC:

```bash
awk 'BEGIN{sum=0} {val[NR]=$1; sum+=$1}
     END{
       for(i=1; i<=NR; i++)
         printf "PC%d: %.2f%%\n", i, val[i]/sum*100
     }' hapmap3_pca.eigenval
```

Also run PCA within Europeans (for section 3.6):

```bash
awk -F'\t' 'NR>1 && $6 == "EUR" {print $1, $1}' 1kg_samples.txt > samples_EUR.txt

plink --bfile hapmap3_pruned_set \
      --keep samples_EUR.txt \
      --pca 10 \
      --out pca_EUR
```

### 3.3 Visualise PCA in Colab

Switch to **Google Colab**. Make sure you have uploaded `lab4_results.zip` (see the file transfer section at the top).

**Scree plot — variance explained:**

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid", font_scale=1.2)

# Read eigenvalues
eigenval = pd.read_csv("hapmap3_pca.eigenval", header=None, names=["eigenvalue"])
eigenval["PC"] = range(1, len(eigenval) + 1)
eigenval["pct"] = eigenval["eigenvalue"] / eigenval["eigenvalue"].sum() * 100

plt.figure(figsize=(9, 5))
plt.bar(eigenval["PC"], eigenval["pct"], color="steelblue", edgecolor="white")
plt.plot(eigenval["PC"], eigenval["pct"], "o-", color="navy", markersize=5)
plt.xlabel("Principal Component")
plt.ylabel("Variance explained (%)")
plt.title("Scree plot")
plt.xticks(range(1, 21))
plt.tight_layout()
plt.savefig("pca_scree.png", dpi=150)
plt.show()

print(eigenval[["PC", "pct"]].to_string(index=False))
```

**Load PCA scores and population labels:**

```python
# Read PC scores
pc_cols = ["FID", "IID"] + [f"PC{i}" for i in range(1, 21)]
pca = pd.read_csv("hapmap3_pca.eigenvec", sep=r"\s+", header=None, names=pc_cols)

# Read population metadata
geo = pd.read_csv("1kg_samples.txt", sep="\t")
geo = geo.rename(columns={"Sample name": "IID"})

# Merge
data = pca.merge(geo[["IID", "Population code", "Population name",
                        "Superpopulation code", "Superpopulation name"]],
                 on="IID", how="inner")

print(f"Individuals with labels: {len(data)}")
print(data["Superpopulation name"].value_counts())
```

**PC1 vs PC2 coloured by superpopulation:**

```python
fig, ax = plt.subplots(figsize=(9, 6))

for superpop, grp in data.groupby("Superpopulation name"):
    ax.scatter(grp["PC1"], grp["PC2"], label=superpop,
               alpha=0.7, s=15, linewidths=0)

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("PCA — continental ancestry")
ax.legend(title="Superpopulation", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.savefig("pca_superpop.png", dpi=150, bbox_inches="tight")
plt.show()
```

You should see clearly separated clusters for: **AFR** (African), **EUR** (European), **EAS** (East Asian), **SAS** (South Asian), **AMR** (Admixed American).

**PC1 vs PC2 coloured by sub-population:**

```python
fig, ax = plt.subplots(figsize=(11, 7))

populations = data["Population name"].unique()
palette = sns.color_palette("tab20", len(populations))

for i, pop in enumerate(sorted(populations)):
    grp = data[data["Population name"] == pop]
    ax.scatter(grp["PC1"], grp["PC2"], label=pop,
               color=palette[i], alpha=0.7, s=15, linewidths=0)

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("PCA — sub-population")
ax.legend(title="Population", bbox_to_anchor=(1.02, 1), loc="upper left",
          fontsize=7, ncol=2)
plt.tight_layout()
plt.savefig("pca_subpop.png", dpi=150, bbox_inches="tight")
plt.show()
```

**PC1 vs PC3:**

```python
fig, ax = plt.subplots(figsize=(9, 6))

for superpop, grp in data.groupby("Superpopulation name"):
    ax.scatter(grp["PC1"], grp["PC3"], label=superpop,
               alpha=0.7, s=15, linewidths=0)

ax.set_xlabel("PC1")
ax.set_ylabel("PC3")
ax.set_title("PC1 vs PC3")
ax.legend(title="Superpopulation", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.savefig("pca_pc1_pc3.png", dpi=150, bbox_inches="tight")
plt.show()
```

### 3.4 Within-population PCA (Europeans)

```python
pc_cols_eur = ["FID", "IID"] + [f"PC{i}" for i in range(1, 11)]
pca_eur = pd.read_csv("pca_EUR.eigenvec", sep=r"\s+", header=None, names=pc_cols_eur)

data_eur = pca_eur.merge(geo[["IID", "Population name"]], on="IID", how="inner")

fig, ax = plt.subplots(figsize=(9, 6))
for pop, grp in data_eur.groupby("Population name"):
    ax.scatter(grp["PC1"], grp["PC2"], label=pop, alpha=0.8, s=25, linewidths=0)

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("PCA within European populations")
ax.legend(title="Population", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
plt.tight_layout()
plt.savefig("pca_EUR_subpop.png", dpi=150, bbox_inches="tight")
plt.show()
```

Within European populations, PCs often separate Northern Europeans (Finnish, British) from Southern Europeans (Iberian, Tuscan).

### 3.5 PCA as covariates in GWAS (Cloud Shell)

Population stratification can cause spurious associations. The standard solution is to include PC scores as covariates.

```bash
plink --bfile hapmap3_qc \
      --pheno BMI_pheno.txt \
      --assoc \
      --linear \
      --covar hapmap3_pca.eigenvec \
      --covar-number 1-10 \
      --out bmi_assoc_pca_corrected
```

Compare hits before and after correction:

```bash
# Before
awk 'NR>1 && $9 < 5e-8' bmi_assoc.assoc.linear | wc -l

# After
awk 'NR>1 && $9 < 5e-8' bmi_assoc_pca_corrected.assoc.linear | wc -l
```

### 3.6 Genomic inflation factor λ_GC in Colab

The genomic inflation factor $\lambda_{GC}$ measures residual stratification:

$$\lambda_{GC} = \frac{\text{median}(\chi^2_{\text{observed}})}{0.4549}$$

A value of $\lambda_{GC} \approx 1.00$ indicates no inflation.

```python
from scipy import stats

def compute_lambda(pvals):
    """Compute genomic inflation factor from a series of p-values."""
    pvals = pvals.dropna()
    chisq = stats.chi2.ppf(1 - pvals, df=1)
    return np.median(chisq) / 0.4549

# Uncorrected
res_raw = pd.read_csv("bmi_assoc.assoc.linear", sep=r"\s+")
lam_raw = compute_lambda(res_raw["P"])
print(f"Lambda (uncorrected):   {lam_raw:.3f}")

# PC-corrected (keep only the ADD test rows)
res_corr = pd.read_csv("bmi_assoc_pca_corrected.assoc.linear", sep=r"\s+")
res_corr = res_corr[res_corr["TEST"] == "ADD"]
lam_corr = compute_lambda(res_corr["P"])
print(f"Lambda (PC-corrected):  {lam_corr:.3f}")
```

### Exercise 2

1. Examine the scree plot. How many PCs are needed to capture the main axes of genetic variation?
2. In the PC1 vs PC2 plot, which two superpopulations are most separated along PC1? Along PC2?
3. Within Europeans, which populations are most separated in PC space? What historical events might explain this?
4. How does $\lambda_{GC}$ change before and after adding 10 PCs as covariates?

---

## Part IV. Detecting and Removing Population Outliers

In a homogeneous cohort study, individuals who cluster far from the main group in PCA space likely have different ancestry. We identify and remove them.

### 4.1 Identify EUR-like individuals in Colab

```python
# Compute the EUR centroid and standard deviations
eur = data[data["Superpopulation code"] == "EUR"]
eur_mean = eur[["PC1", "PC2"]].mean()
eur_std  = eur[["PC1", "PC2"]].std()

print(f"EUR centroid: PC1 = {eur_mean['PC1']:.4f}, PC2 = {eur_mean['PC2']:.4f}")

# Flag individuals within 3 SD of the EUR centroid on PC1 and PC2
data["eur_like"] = (
    (np.abs(data["PC1"] - eur_mean["PC1"]) < 3 * eur_std["PC1"]) &
    (np.abs(data["PC2"] - eur_mean["PC2"]) < 3 * eur_std["PC2"])
)

print(f"Individuals within 3 SD of EUR centroid: {data['eur_like'].sum()}")

# Visualise
fig, ax = plt.subplots(figsize=(9, 6))
colours = {True: "steelblue", False: "lightgrey"}
shapes  = {True: "o", False: "x"}
labels  = {True: "EUR-like (kept)", False: "Excluded"}

for superpop, grp in data.groupby("Superpopulation code"):
    for keep, sub in grp.groupby("eur_like"):
        ax.scatter(sub["PC1"], sub["PC2"],
                   c=colours[keep], marker=shapes[keep],
                   alpha=0.6, s=15, linewidths=0.5,
                   label=f"{superpop} — {labels[keep]}" if keep else None)

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("EUR-like selection (within 3 SD of EUR centroid)")
plt.tight_layout()
plt.savefig("pca_eur_selection.png", dpi=150)
plt.show()

# Save the EUR-like sample list for PLINK
eur_keep = data[data["eur_like"]][["FID", "IID"]]
eur_keep.to_csv("samples_EUR_like.txt", sep=" ", index=False, header=False)
print(f"Saved {len(eur_keep)} EUR-like individuals to samples_EUR_like.txt")
```

**Download `samples_EUR_like.txt` from Colab** (Files panel → right-click → Download), then upload it to Cloud Shell and apply the filter:

```bash
# In Cloud Shell — upload samples_EUR_like.txt first using Cloud Shell upload button
plink --bfile hapmap3_qc \
      --keep samples_EUR_like.txt \
      --make-bed \
      --out hapmap3_EUR

wc -l hapmap3_EUR.fam
```

---

## Part V. Ancestry Prediction from PCA (Colab)

PCA scores can be used to **predict the ancestry** of individuals of unknown origin. This is the principle behind commercial direct-to-consumer genetic ancestry tests (e.g. 23andMe, AncestryDNA).

We use the 1000 Genomes individuals as a **labelled reference panel** and classify them using their PC coordinates.

### 5.1 Prepare features and labels

```python
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, classification_report
import matplotlib.pyplot as plt
import seaborn as sns

# Features: top 10 PCs
pc_features = [f"PC{i}" for i in range(1, 11)]

X = data[pc_features].values
y = data["Superpopulation code"].values

print(f"Total individuals: {len(X)}")
print("Superpopulation counts:")
for pop, n in zip(*np.unique(y, return_counts=True)):
    print(f"  {pop}: {n}")
```

### 5.2 Train / test split

```python
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

print(f"Training set: {len(X_train)} individuals")
print(f"Test set:     {len(X_test)} individuals")
```

### 5.3 k-Nearest Neighbours classifier

```python
knn = KNeighborsClassifier(n_neighbors=5)
knn.fit(X_train, y_train)
y_pred_knn = knn.predict(X_test)

# Overall accuracy
acc_knn = (y_pred_knn == y_test).mean()
print(f"k-NN accuracy (k=5): {acc_knn * 100:.1f}%")

# Detailed report
print("\nClassification report:")
print(classification_report(y_test, y_pred_knn))

# Confusion matrix
cm = confusion_matrix(y_test, y_pred_knn, labels=np.unique(y))
fig, ax = plt.subplots(figsize=(7, 5))
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
            xticklabels=np.unique(y), yticklabels=np.unique(y), ax=ax)
ax.set_xlabel("Predicted")
ax.set_ylabel("True")
ax.set_title("k-NN confusion matrix (k=5)")
plt.tight_layout()
plt.savefig("knn_confusion.png", dpi=150)
plt.show()
```

### 5.4 Effect of k on accuracy

```python
k_values = [1, 3, 5, 10, 20, 50]
accuracies = []

for k in k_values:
    knn_k = KNeighborsClassifier(n_neighbors=k)
    knn_k.fit(X_train, y_train)
    acc = knn_k.score(X_test, y_test) * 100
    accuracies.append(acc)
    print(f"k={k:3d}: {acc:.1f}%")

plt.figure(figsize=(7, 5))
plt.plot(k_values, accuracies, "o-", color="steelblue", markersize=8)
plt.xlabel("k (number of neighbours)")
plt.ylabel("Classification accuracy (%)")
plt.title("k-NN ancestry classification accuracy")
plt.ylim(90, 100)
plt.tight_layout()
plt.savefig("knn_accuracy.png", dpi=150)
plt.show()
```

### 5.5 Visualise predictions in PC space

```python
# Predict ancestry for all individuals
y_pred_all = knn.predict(X)
data["predicted"] = y_pred_all
data["correct"]   = data["predicted"] == data["Superpopulation code"]

fig, ax = plt.subplots(figsize=(9, 6))
for superpop, grp in data.groupby("Superpopulation code"):
    correct   = grp[grp["correct"]]
    incorrect = grp[~grp["correct"]]
    ax.scatter(correct["PC1"],   correct["PC2"],   alpha=0.6, s=12, linewidths=0, label=superpop)
    ax.scatter(incorrect["PC1"], incorrect["PC2"], marker="x", s=60, linewidths=1.5, color="black")

ax.scatter([], [], marker="x", color="black", s=60, label="Misclassified")
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("k-NN ancestry predictions (k=5)\nBlack × = misclassified")
ax.legend(title="True superpopulation", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.savefig("knn_predictions.png", dpi=150, bbox_inches="tight")
plt.show()

misclass = data[~data["correct"]][["IID", "Superpopulation code", "predicted", "PC1", "PC2"]]
print(f"\nMisclassified individuals: {len(misclass)}")
print(misclass.head(20).to_string(index=False))
```

Misclassifications are most common among **AMR** (Admixed American) individuals, who have mixed European, Native American, and sometimes African ancestry.

### 5.6 Random Forest classifier

Random Forests typically outperform k-NN for ancestry classification because they capture non-linear decision boundaries.

```python
rf = RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1)
rf.fit(X_train, y_train)
y_pred_rf = rf.predict(X_test)

acc_rf = (y_pred_rf == y_test).mean()
print(f"Random Forest accuracy: {acc_rf * 100:.1f}%")
print(f"k-NN accuracy:          {acc_knn * 100:.1f}%")

print("\nClassification report (Random Forest):")
print(classification_report(y_test, y_pred_rf))

# Confusion matrix
cm_rf = confusion_matrix(y_test, y_pred_rf, labels=np.unique(y))
fig, ax = plt.subplots(figsize=(7, 5))
sns.heatmap(cm_rf, annot=True, fmt="d", cmap="Greens",
            xticklabels=np.unique(y), yticklabels=np.unique(y), ax=ax)
ax.set_xlabel("Predicted")
ax.set_ylabel("True")
ax.set_title("Random Forest confusion matrix")
plt.tight_layout()
plt.savefig("rf_confusion.png", dpi=150)
plt.show()
```

**Variable importance — which PCs matter most?**

```python
importance = pd.DataFrame({
    "PC": pc_features,
    "importance": rf.feature_importances_
}).sort_values("importance", ascending=True)

plt.figure(figsize=(7, 5))
plt.barh(importance["PC"], importance["importance"], color="steelblue")
plt.xlabel("Mean Decrease in Impurity")
plt.title("PC importance for ancestry classification\n(Random Forest)")
plt.tight_layout()
plt.savefig("rf_importance.png", dpi=150)
plt.show()
```

PC1 and PC2 typically dominate because they capture the major continental ancestry axes.

### Exercise 3

1. What is the k-NN classification accuracy for $k = 5$ using 10 PCs? Which superpopulation is hardest to classify and why?
2. Repeat the k-NN analysis using only PC1 and PC2 (set `pc_features = ["PC1", "PC2"]`). How much does accuracy drop?
3. Look at the misclassified individuals. Which true superpopulation do they belong to, and which are they assigned to? What does this tell you about genetic admixture?
4. Compare k-NN vs Random Forest accuracy. Which performs better?

---

## Summary

### Cloud Shell pipeline (PLINK commands)

```bash
# 1. QC
plink --bfile hapmap3 \
      --mind 0.05 --geno 0.02 --maf 0.01 --hwe 1e-6 \
      --make-bed --out hapmap3_qc

# 2. Summary statistics for Colab exploration
plink --bfile hapmap3_qc --missing --freq --hardy --out hapmap3_summary
plink --bfile hapmap3_qc --het --out hapmap3_het

# 3. Association test (no covariates)
plink --bfile hapmap3_qc \
      --pheno BMI_pheno.txt --assoc --linear \
      --out bmi_assoc

# 4. LD pruning + PCA
plink --bfile hapmap3_qc --indep-pairwise 50 5 0.2 --out hapmap3_pruned
plink --bfile hapmap3_qc --extract hapmap3_pruned.prune.in \
      --make-bed --out hapmap3_pruned_set
plink --bfile hapmap3_pruned_set --pca 20 --out hapmap3_pca

# 5. PCA-corrected association
plink --bfile hapmap3_qc \
      --pheno BMI_pheno.txt --assoc --linear \
      --covar hapmap3_pca.eigenvec --covar-number 1-10 \
      --out bmi_assoc_pca_corrected

# 6. Package results for Colab
zip lab4_results.zip hapmap3_summary.* hapmap3_het.het \
    hapmap3_pca.eigenvec hapmap3_pca.eigenval pca_EUR.eigenvec \
    bmi_assoc.assoc.linear bmi_assoc_pca_corrected.assoc.linear \
    1kg_samples.txt
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
