# Lab 7. Polygenic Scores with PRSice-2

In this lab you will build and evaluate **polygenic scores (PGS)** using **PRSice-2**.

All the tools (PLINK, PRSice) and data live in a single directory `~/Sociogenomics`. A zip file with all required data is provided.

---

## 0. Setup

### 0.1 Download and unpack data

```bash
mkdir -p ~/Sociogenomics/Data ~/Sociogenomics/Results
cd ~/Sociogenomics/Data
wget -O lab7_data.zip "https://github.com/nicolabarban/sociogenomics_2025_2026/raw/gh_pages/labs/week7/lab7_data.zip"
unzip -o lab7_data.zip
ls
```

You should see: `1kg_hm3.*`, `1kg_hm3_QC_CEU.*`, `Trait2.ma`, `1kg.Trait2.phen`, `1kg_pca.eigenvec`, `EUR.id`, and others.

### 0.2 Install PLINK

```bash
cd ~/Sociogenomics
wget -q https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip -O plink.zip
unzip -o plink.zip -d .
chmod +x plink
./plink --version
```

### 0.3 Install PRSice-2

```bash
cd ~/Sociogenomics
wget -q https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip
unzip -o PRSice_linux.zip -d .
chmod +x PRSice_linux
./PRSice_linux --help | head -3
```

### 0.4 Install R and required packages

```bash
sudo apt-get update && sudo apt-get install -y r-base
R -e 'install.packages(c("ggplot2","data.table","optparse"), repos="https://cloud.r-project.org", Ncpus=4)'
Rscript --version
```

> **Note:** This step takes a few minutes and must be repeated each Cloud Shell session. Alternatively, use **Posit Cloud** ([posit.cloud](https://posit.cloud)) where R is pre-installed.

### 0.5 Check everything works

```bash
cd ~/Sociogenomics
./plink --version
./PRSice_linux --help | head -1
Rscript --version
ls Data/1kg_hm3_QC_CEU.bed Data/Trait2.ma Data/1kg.Trait2.phen
```

All four commands should succeed. If any fails, fix it before proceeding.

---

## Part I. Explore the data

```bash
cd ~/Sociogenomics/Data
```

### Summary statistics (discovery GWAS)

```bash
head Trait2.ma
wc -l Trait2.ma
```

**Question:** How many SNPs? How many columns? What do `BETA` and `P` represent?

Count genome-wide significant SNPs:

```bash
awk 'NR>1 && $8 < 5e-8' Trait2.ma | wc -l
```

### Target genotype data (European subset, QC'd)

```bash
wc -l 1kg_hm3_QC_CEU.bim   # SNPs
wc -l 1kg_hm3_QC_CEU.fam   # individuals
```

### Phenotype

```bash
head 1kg.Trait2.phen
```

Three columns: `FID IID Trait2` (continuous, $h^2 \approx 0.2$).

---

## Part II. Run PRSice-2

PRSice-2 implements the **Clumping + Thresholding (C+T)** method in a single
automated pipeline. Conceptually it does three things:

1. **Align** the base (GWAS sumstats) and target (your PLINK files) on SNP
   identifier, effect allele and strand;
2. **Clump** the base to keep one SNP per LD block (so we don't add the same
   signal twice);
3. **Score** each target individual at a grid of $p$-value thresholds and
   regress the phenotype on each resulting PGS to pick the best threshold.

Each command-line flag controls one of those three steps.

### 2.1 A first, well-annotated run

```bash
cd ~/Sociogenomics

Rscript PRSice.R --dir . \
    --prsice ./PRSice_linux \
    --base Data/Trait2.ma \
    --target Data/1kg_hm3_QC_CEU \
    --snp SNP \
    --A1 A1 \
    --A2 A2 \
    --stat BETA \
    --pvalue P \
    --beta \
    --pheno Data/1kg.Trait2.phen \
    --binary-target F \
    --bar-levels 5e-08,5e-06,5e-04,0.05,0.5,1 \
    --fastscore \
    --all-score \
    --out Results/Trait2_PRSice
```

### 2.2 What each flag does

**Wrapper / executables**

| Flag | Meaning |
|------|---------|
| `--dir .` | Working directory PRSice uses for intermediate files |
| `--prsice ./PRSice_linux` | Path to the compiled PRSice binary |
| `--out Results/Trait2_PRSice` | Prefix for every output file |

**Base (discovery GWAS)**

| Flag | Meaning |
|------|---------|
| `--base Data/Trait2.ma` | Summary-statistics file |
| `--snp SNP` | Column name with the SNP ID (rsID) |
| `--A1 A1` | Effect allele column |
| `--A2 A2` | Reference (non-effect) allele column |
| `--stat BETA` | Column with the effect size |
| `--beta` | Tells PRSice `BETA` is a linear regression coefficient. Use `--or` (odds ratio) for binary traits |
| `--pvalue P` | Column with the $p$-value |

*Optional but often useful base flags:* `--chr CHR`, `--bp BP`, `--se SE`,
`--info INFO,0.8` (filter on imputation quality), `--maf MAF,0.01` (filter on
allele frequency in the base).

**Target (your sample)**

| Flag | Meaning |
|------|---------|
| `--target Data/1kg_hm3_QC_CEU` | PLINK `.bed/.bim/.fam` prefix |
| `--pheno Data/1kg.Trait2.phen` | Phenotype file (FID IID pheno) |
| `--binary-target F` | `F` = continuous trait; `T` = case/control |

**Thresholding**

| Flag | Meaning |
|------|---------|
| `--bar-levels 5e-08,5e-06,5e-04,0.05,0.5,1` | Thresholds for the $R^2$ barplot |
| `--fastscore` | Only evaluate the bar levels above (fast) |
| *(alternative)* `--lower 5e-08 --upper 0.5 --interval 5e-05` | High-resolution scan over many thresholds |
| `--no-full` | Skip the $P_T = 1$ score |

**Output**

| Flag | Meaning |
|------|---------|
| `--all-score` | Write a PGS for every individual at every threshold (needed for Part III) |
| `--print-snp` | Write the list of SNPs kept at each threshold |
| `--quantile 10` | Also produce decile-based quantile plots |

**Clumping (defaults shown)**

| Flag | Default | Meaning |
|------|---------|---------|
| `--clump-kb` | `250` | Physical window (kb) for clumping |
| `--clump-r2` | `0.1` | LD threshold; SNPs with $r^2$ above this are dropped |
| `--clump-p` | `1` | Only SNPs with base $p$ below this enter clumping |
| `--no-clump` | off | Turn clumping off altogether (don't do this on unclumped GWAS) |

### 2.3 Inspect the output

```bash
cat Results/Trait2_PRSice.summary
cat Results/Trait2_PRSice.prsice
ls Results/Trait2_PRSice*
```

| File | Content |
|------|---------|
| `*.summary` | Best threshold and $R^2$ |
| `*.prsice` | $R^2$ at each threshold |
| `*.all_score` | PGS for each individual at each threshold |
| `*.snp` | Selected SNPs with weights |
| `*_BARPLOT_*.png` | Barplot of $R^2$ by threshold |

### Exercise 1

1. Which $p$-value threshold gives the best $R^2$?
2. How many SNPs are included at that threshold?
3. Does $R^2$ increase monotonically with threshold? Why or why not?

---

## Part II-bis. Try different PRSice options

Now that the baseline run is in place, repeat PRSice while changing **one set
of options at a time**. Save each run under a different `--out` prefix so you
can compare the results side-by-side.

### A. High-resolution $p$-value scan

`--fastscore` only tests the thresholds listed in `--bar-levels`. To find the
true optimum, use a dense grid via `--lower / --upper / --interval`:

```bash
Rscript PRSice.R --dir . \
    --prsice ./PRSice_linux \
    --base Data/Trait2.ma \
    --target Data/1kg_hm3_QC_CEU \
    --snp SNP --A1 A1 --A2 A2 --stat BETA --pvalue P --beta \
    --pheno Data/1kg.Trait2.phen --binary-target F \
    --lower 5e-08 --upper 0.5 --interval 5e-05 \
    --all-score \
    --out Results/Trait2_PRSice_hires
```

Look at `Trait2_PRSice_hires.prsice`: there is now one row per tested
threshold. The `*_HIGH-RES_*.png` plot shows $R^2$ vs $P_T$ as a smooth curve
with the optimum flagged.

**Compare:** does the best threshold found by the fine scan agree with the
coarse `--fastscore` run?

### B. Loose vs strict clumping

Clumping parameters control **how much LD we allow between the SNPs that enter
the score**. A stricter window + lower $r^2$ keeps fewer but more independent
SNPs.

```bash
# Stricter: 500 kb window, r^2 < 0.01 (very independent SNPs)
Rscript PRSice.R --dir . --prsice ./PRSice_linux \
    --base Data/Trait2.ma --target Data/1kg_hm3_QC_CEU \
    --snp SNP --A1 A1 --A2 A2 --stat BETA --pvalue P --beta \
    --pheno Data/1kg.Trait2.phen --binary-target F \
    --clump-kb 500 --clump-r2 0.01 \
    --bar-levels 5e-08,5e-06,5e-04,0.05,0.5,1 --fastscore --all-score \
    --out Results/Trait2_PRSice_strict

# Looser: 100 kb, r^2 < 0.5 (more correlated SNPs allowed)
Rscript PRSice.R --dir . --prsice ./PRSice_linux \
    --base Data/Trait2.ma --target Data/1kg_hm3_QC_CEU \
    --snp SNP --A1 A1 --A2 A2 --stat BETA --pvalue P --beta \
    --pheno Data/1kg.Trait2.phen --binary-target F \
    --clump-kb 100 --clump-r2 0.5 \
    --bar-levels 5e-08,5e-06,5e-04,0.05,0.5,1 --fastscore --all-score \
    --out Results/Trait2_PRSice_loose
```

**Compare** the number of SNPs kept (`wc -l Results/*_strict.snp` vs
`*_loose.snp`) and the peak $R^2$. Stricter clumping usually keeps fewer SNPs
but reduces the risk of double-counting LD-correlated signals.

### C. Adjusting for covariates inside PRSice

Instead of standardising the PGS in R and regressing in a second step, you can
let PRSice fit the full model internally. Covariates must be a tab/space
separated file with a header, first two columns `FID IID`.

Make a small covariate file from the PCs:

```bash
awk 'BEGIN{OFS="\t"; print "FID","IID","PC1","PC2","PC3","PC4","PC5"}
     {print $1,$2,$3,$4,$5,$6,$7}' Data/1kg_pca.eigenvec \
   > Data/pcs.cov
head Data/pcs.cov
```

Run PRSice with the covariates:

```bash
Rscript PRSice.R --dir . --prsice ./PRSice_linux \
    --base Data/Trait2.ma --target Data/1kg_hm3_QC_CEU \
    --snp SNP --A1 A1 --A2 A2 --stat BETA --pvalue P --beta \
    --pheno Data/1kg.Trait2.phen --binary-target F \
    --cov Data/pcs.cov --cov-col PC1,PC2,PC3,PC4,PC5 \
    --bar-levels 5e-08,5e-06,5e-04,0.05,0.5,1 --fastscore --all-score \
    --out Results/Trait2_PRSice_cov
```

Now `*.summary` reports $R^2$ as the **incremental** $R^2$ over the
covariates-only model. Compare it with the raw $R^2$ from the baseline run.

### D. Scoring method: average vs sum

`--score` controls how per-SNP contributions are aggregated:

| Value | Formula |
|-------|---------|
| `avg` (default) | $\mathrm{PGS}_i = \tfrac{1}{M}\sum_j \beta_j\,x_{ij}$ |
| `sum` | $\mathrm{PGS}_i = \sum_j \beta_j\,x_{ij}$ |
| `std` | Standardises genotypes before summing |

```bash
Rscript PRSice.R --dir . --prsice ./PRSice_linux \
    --base Data/Trait2.ma --target Data/1kg_hm3_QC_CEU \
    --snp SNP --A1 A1 --A2 A2 --stat BETA --pvalue P --beta \
    --pheno Data/1kg.Trait2.phen --binary-target F \
    --score sum \
    --bar-levels 5e-08,5e-06,5e-04,0.05,0.5,1 --fastscore --all-score \
    --out Results/Trait2_PRSice_sum
```

$R^2$ should be **identical** to the default run — only the scale of the PGS
changes. Plot both PGS distributions in R to confirm.

### E. Decile / quantile plot

`--quantile 10` produces a "decile plot": the sample is split into 10 equal
bins of PGS, and the mean phenotype is plotted per bin. It is the classic way
to visualise risk gradients.

```bash
Rscript PRSice.R --dir . --prsice ./PRSice_linux \
    --base Data/Trait2.ma --target Data/1kg_hm3_QC_CEU \
    --snp SNP --A1 A1 --A2 A2 --stat BETA --pvalue P --beta \
    --pheno Data/1kg.Trait2.phen --binary-target F \
    --bar-levels 5e-08,5e-06,5e-04,0.05,0.5,1 --fastscore --all-score \
    --quantile 10 \
    --out Results/Trait2_PRSice_q10
```

Look for `*_QUANTILES_*.png`.

### F. Restrict the base by MAF or INFO

Very rare or poorly imputed SNPs are the first suspects when a PGS looks
noisy. Filter them out at the source:

```bash
Rscript PRSice.R --dir . --prsice ./PRSice_linux \
    --base Data/Trait2.ma --target Data/1kg_hm3_QC_CEU \
    --snp SNP --A1 A1 --A2 A2 --stat BETA --pvalue P --beta \
    --maf MAF,0.05 \
    --pheno Data/1kg.Trait2.phen --binary-target F \
    --bar-levels 5e-08,5e-06,5e-04,0.05,0.5,1 --fastscore --all-score \
    --out Results/Trait2_PRSice_maf05
```

(Only works if `Trait2.ma` has a `MAF` column — otherwise PRSice will tell
you.) Replace with `--info INFO,0.8` if the base has imputation quality.

### Exercise 1-bis

Run at least **two** of the variants above (suggested: high-res scan + strict
vs loose clumping) and fill in the comparison table:

| Run | Best $P_T$ | # SNPs | Best $R^2$ |
|-----|------------|--------|------------|
| Baseline (`Trait2_PRSice`) | | | |
| High-res (`_hires`) | | | |
| Strict clumping (`_strict`) | | | |
| Loose clumping (`_loose`) | | | |
| With covariates (`_cov`) | | | |

Write 2-3 sentences on what moves the $R^2$ the most: threshold, clumping, or
covariates?

---

## Part III. Analyse the PGS in R

Start R:

```bash
cd ~/Sociogenomics
R
```

### 3.1 Load data

```r
scores <- read.table("Results/Trait2_PRSice.all_score", header = TRUE)
pheno  <- read.table("Data/1kg.Trait2.phen", header = FALSE,
                     col.names = c("FID", "IID", "Trait2"))
pca    <- read.table("Data/1kg_pca.eigenvec", header = FALSE,
                     col.names = c("FID", "IID", paste0("PC", 1:10)))

d <- merge(merge(scores, pheno, by = c("FID", "IID")),
           pca, by = c("FID", "IID"))

# Standardise the best PGS (adjust column name to your best threshold)
d$PGS <- scale(d[["Pt_0.05"]])
```

### 3.2 PGS distribution

```r
hist(d$PGS, breaks = 30, col = "steelblue", border = "white",
     main = "PGS distribution (EUR)", xlab = "Standardised PGS")
```

**Question:** Is it approximately normal? Why?

### 3.3 Incremental $R^2$

```r
mod0 <- lm(Trait2 ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = d)
mod1 <- lm(Trait2 ~ PGS+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = d)

summary(mod1)$coefficients["PGS", ]

delta_r2 <- summary(mod1)$r.squared - summary(mod0)$r.squared
cat("R2 PCs only: ", round(summary(mod0)$r.squared, 4), "\n")
cat("R2 PCs+PGS:  ", round(summary(mod1)$r.squared, 4), "\n")
cat("Incremental: ", round(delta_r2, 4), "\n")
```

### 3.4 Compare all thresholds

```r
thresholds <- c("Pt_5e.08", "Pt_5e.06", "Pt_5e.04", "Pt_0.05", "Pt_0.5", "Pt_1")
results <- data.frame(threshold = thresholds, delta_r2 = NA)

for (i in seq_along(thresholds)) {
  d$tmp <- scale(d[[thresholds[i]]])
  mod <- lm(Trait2 ~ tmp+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = d)
  results$delta_r2[i] <- summary(mod)$r.squared - summary(mod0)$r.squared
}

print(results)
barplot(results$delta_r2, names.arg = results$threshold,
        col = "steelblue", border = NA, las = 2,
        main = "Incremental R2 by threshold",
        ylab = expression(Delta ~ R^2))
```

### 3.5 Bootstrap 95% CI

```r
library(boot)
set.seed(2026)

boot_fn <- function(data, idx) {
  ds <- data[idx, ]
  m0 <- lm(Trait2 ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = ds)
  m1 <- lm(Trait2 ~ PGS+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = ds)
  summary(m1)$r.squared - summary(m0)$r.squared
}

b <- boot(d, boot_fn, R = 1000)
boot.ci(b, type = "norm")
```

### Exercise 2

1. What is the incremental $R^2$ of the best PGS?
2. What is the 95% bootstrap CI?
3. Which threshold is optimal and why?

---

## Part IV. Cross-ancestry PGS portability

The PGS was trained on EUR GWAS. Does it predict in other populations?

### 4.1 Score all 1,092 individuals

Back in the terminal (quit R with `q()`):

```bash
cd ~/Sociogenomics

./plink --bfile Data/1kg_hm3 \
        --score Results/Trait2_PRSice.snp 1 2 4 header \
        --out Results/Trait2_pgs_all_pops
```

### 4.2 Analyse in R

```r
pgs_all <- read.table("Results/Trait2_pgs_all_pops.profile", header = TRUE)
pheno   <- read.table("Data/1kg.Trait2.phen", header = FALSE,
                      col.names = c("FID", "IID", "Trait2"))
pop     <- read.table("Data/1kg-sample-2504-phased.txt", header = TRUE)

da <- merge(pgs_all[, c("FID","IID","SCORE")], pheno, by = c("FID","IID"))
da <- merge(da, pop[, c("sample","super_pop")], by.x = "IID", by.y = "sample")

cat("R2 by super-population:\n")
for (p in c("EUR", "EAS", "AFR", "AMR")) {
  sub <- da[da$super_pop == p, ]
  r2 <- cor(sub$SCORE, sub$Trait2)^2
  cat(sprintf("  %s: R2 = %.4f (N = %d)\n", p, r2, nrow(sub)))
}

da$PGS <- scale(da$SCORE)
boxplot(PGS ~ super_pop, data = da,
        col = c("tomato","gold","skyblue","plum"),
        main = "PGS by super-population",
        xlab = "Super-population", ylab = "Standardised PGS")
```

### Exercise 3

1. In which population is the PGS most predictive? Least?
2. Why does the PGS transfer poorly to non-EUR populations?
3. What are possible solutions? (multi-ancestry GWAS, PRS-CSx, ...)

---

