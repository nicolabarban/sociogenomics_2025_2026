# Lab 7. Polygenic Scores with PRSice-2

In this lab you will build, evaluate, and interpret **polygenic scores (PGS)** using **PRSice-2**, the most widely used tool for the clumping + thresholding (C+T) method.

The lab covers:

* Setting up R and PRSice-2 on Google Cloud Shell (persistent install via conda)
* Running a full C+T pipeline in one command
* Comparing PGS $R^2$ at different $p$-value thresholds
* Incremental $R^2$ and bootstrap confidence intervals in R
* Cross-ancestry PGS portability
* Gene $\times$ Environment interaction with PGS

---

## 0. Setup

### 0.1 Install Git LFS and clone the repo

Open **Google Cloud Shell**:

```bash
sudo apt-get update && sudo apt-get install -y git-lfs
git lfs install
```

Clone or refresh the course repo:

```bash
cd ~
if [ -d ~/sociogenomics_2025_2026/.git ]; then
  cd ~/sociogenomics_2025_2026 && git pull
else
  rm -rf ~/sociogenomics_2025_2026
  git clone https://github.com/nicolabarban/sociogenomics_2025_2026.git
fi
cd ~/sociogenomics_2025_2026 && git lfs pull
```

### 0.2 Directories and data

```bash
mkdir -p ~/Sociogenomics/Data ~/Sociogenomics/Results ~/Sociogenomics/Software
rm -rf ~/Sociogenomics/Data/*

cp ~/sociogenomics_2025_2026/data/1kg_hm3.* \
   ~/sociogenomics_2025_2026/data/1kg_hm3_QC_CEU.* \
   ~/sociogenomics_2025_2026/data/1kg_pca.eigenvec \
   ~/sociogenomics_2025_2026/data/1kg_pca.eigenval \
   ~/sociogenomics_2025_2026/data/Trait2.ma \
   ~/sociogenomics_2025_2026/data/1kg.Trait2.phen \
   ~/sociogenomics_2025_2026/data/1kg-sample-2504-phased.txt \
   ~/sociogenomics_2025_2026/data/EUR.id \
   ~/Sociogenomics/Data/
```

Verify:

```bash
ls ~/Sociogenomics/Data/1kg_hm3_QC_CEU.* ~/Sociogenomics/Data/Trait2.ma
```

### 0.3 Install PLINK

```bash
plink --version 2>/dev/null || { bash ~/sociogenomics_2025_2026/scripts/setup_plink19.sh && source ~/.bashrc; }
```

### 0.4 Install R via conda (persistent --- no reinstall next session)

Google Cloud Shell resets `apt` packages between sessions, so `sudo apt-get install r-base` must be repeated every time. **Conda installs R into `$HOME`**, which persists.

**First time only** (~3 minutes):

```bash
if [ ! -d ~/miniforge3 ]; then
  curl -L -o ~/miniforge.sh \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
  bash ~/miniforge.sh -b -p ~/miniforge3
  ~/miniforge3/bin/conda init bash
  source ~/.bashrc
fi
```

Create an environment with R and all packages PRSice needs:

```bash
if ! conda env list | grep -q sociogen; then
  mamba create -y -n sociogen -c conda-forge \
    r-base r-ggplot2 r-data.table r-optparse r-boot
fi
```

**Every subsequent session** (2 seconds):

```bash
conda activate sociogen
R --version | head -1
```

> **Alternative:** you can also use **Posit Cloud** ([posit.cloud](https://posit.cloud)), where R is pre-installed. Clone the repo via *New Project → Git Repository*, then run the PLINK/PRSice setup from the RStudio Terminal.

### 0.5 Install PRSice-2

```bash
cd ~/Sociogenomics/Software
wget -q https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip
unzip -o PRSice_linux.zip
chmod +x PRSice_linux
```

Verify:

```bash
./PRSice_linux --help | head -3
Rscript --version
```

Both commands should print version info without errors. If `Rscript` is not found, activate conda: `conda activate sociogen`.

Create symlinks for convenience:

```bash
cd ~/Sociogenomics
ln -sf ~/Sociogenomics/Software/PRSice.R .
ln -sf ~/Sociogenomics/Software/PRSice_linux .
```

---

## Part I. Explore the data

```bash
cd ~/Sociogenomics/Data
```

### Summary statistics (base sample)

```bash
head Trait2.ma
wc -l Trait2.ma
```

The `.ma` format has columns including `SNP A1 A2 BETA P N`. This is the **discovery** GWAS.

**Question:** How many SNPs are in the summary statistics? How many reach genome-wide significance ($p < 5 \times 10^{-8}$)?

```bash
awk 'NR>1 && $8 < 5e-8' Trait2.ma | wc -l
```

### Target genotype data

```bash
wc -l 1kg_hm3_QC_CEU.bim   # SNPs after QC
wc -l 1kg_hm3_QC_CEU.fam   # European individuals
```

### Phenotype

```bash
head 1kg.Trait2.phen
```

Three columns: `FID IID Trait2` (continuous, $h^2 \approx 0.2$).

---

## Part II. Run PRSice-2

### 2.1 Run at fixed $p$-value thresholds

```bash
cd ~/Sociogenomics

Rscript PRSice.R --dir . \
    --prsice ./PRSice_linux \
    --base ~/Sociogenomics/Data/Trait2.ma \
    --target ~/Sociogenomics/Data/1kg_hm3_QC_CEU \
    --snp SNP \
    --A1 A1 \
    --A2 A2 \
    --stat BETA \
    --pvalue P \
    --beta \
    --pheno ~/Sociogenomics/Data/1kg.Trait2.phen \
    --binary-target F \
    --bar-levels 5e-08,5e-06,5e-04,0.05,0.5,1 \
    --fastscore \
    --all-score \
    --out ~/Sociogenomics/Results/Trait2_PRSice
```

### 2.2 What PRSice does under the hood

1. Reads the GWAS summary statistics (base)
2. Aligns SNPs between base and target (matching alleles)
3. **Clumps** (LD $r^2 < 0.1$ in 250 kb windows) to select independent SNPs
4. **Computes PGS** at each $p$-value threshold
5. **Regresses** phenotype on PGS and reports $R^2$
6. Produces plots

### 2.3 Inspect the output

```bash
ls ~/Sociogenomics/Results/Trait2_PRSice*
```

| File | Content |
|------|---------|
| `*.summary` | Best threshold and $R^2$ |
| `*.prsice` | $R^2$ at each threshold |
| `*.all_score` | PGS for each individual at each threshold |
| `*.snp` | Selected SNPs with weights |
| `*_BARPLOT_*.png` | Barplot of $R^2$ by threshold |

```bash
cat ~/Sociogenomics/Results/Trait2_PRSice.summary
cat ~/Sociogenomics/Results/Trait2_PRSice.prsice
```

**Question:** Which $p$-value threshold gives the best $R^2$? How many SNPs are included?

### 2.4 The barplot

Download the barplot image to your computer:

```bash
ls ~/Sociogenomics/Results/Trait2_PRSice_BARPLOT_*.png
```

> In Cloud Shell, click the **three-dot menu (⋮)** → **Download file** and enter the path.

**Question:** Does $R^2$ increase monotonically with the threshold? Why or why not?

### Exercise 1

1. What is the best $p$-value threshold for Trait2?
2. How many SNPs are included at that threshold?
3. What is the $R^2$ of the best PGS?

---

## Part III. Analyse the PGS in R

Open R (in Cloud Shell: `conda activate sociogen && R`; or in Posit Cloud).

### 3.1 Load the PGS scores

```r
all_scores <- read.table("~/Sociogenomics/Results/Trait2_PRSice.all_score",
                         header = TRUE)
head(all_scores)
```

Columns: `FID IID Pt_5e.08 Pt_5e.06 Pt_0.0005 Pt_0.05 Pt_0.5 Pt_1` --- one PGS per threshold.

> **Note:** R converts `-` to `.` in column names, so `Pt_5e-08` becomes `Pt_5e.08`.

### 3.2 Load phenotype and PCs

```r
pheno <- read.table("~/Sociogenomics/Data/1kg.Trait2.phen",
                    header = FALSE,
                    col.names = c("FID", "IID", "Trait2"))

pca <- read.table("~/Sociogenomics/Data/1kg_pca.eigenvec",
                  header = FALSE,
                  col.names = c("FID", "IID", paste0("PC", 1:10)))

d <- merge(all_scores, pheno, by = c("FID", "IID"))
d <- merge(d, pca, by = c("FID", "IID"))

# Standardise the best PGS (replace Pt_0.05 with your best threshold)
d$PGS <- scale(d[["Pt_0.05"]])
```

### 3.3 Plot the PGS distribution

```r
hist(d$PGS, breaks = 30, col = "steelblue", border = "white",
     main = "PGS distribution (EUR, best threshold)",
     xlab = "Standardised PGS")
```

**Question:** Is it approximately normal? Why? (Hint: Central Limit Theorem.)

### 3.4 Incremental $R^2$

```r
# Baseline model (PCs only)
mod0 <- lm(Trait2 ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = d)

# Full model (PCs + PGS)
mod1 <- lm(Trait2 ~ PGS+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = d)

# Compare
summary(mod1)$coefficients["PGS", ]

delta_r2 <- summary(mod1)$r.squared - summary(mod0)$r.squared
cat("R2 baseline (PCs):  ", round(summary(mod0)$r.squared, 4), "\n")
cat("R2 full (PCs + PGS):", round(summary(mod1)$r.squared, 4), "\n")
cat("Incremental R2:     ", round(delta_r2, 4), "\n")
```

### 3.5 Compare all thresholds

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
        main = "Incremental R2 by p-value threshold",
        ylab = expression(Delta ~ R^2))
```

### 3.6 Bootstrap 95% confidence interval

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
3. Does including more SNPs always improve $R^2$? What is the optimal threshold?

---

## Part IV. Cross-ancestry PGS portability

The PGS was trained on EUR GWAS. How well does it predict in other populations?

### 4.1 Compute PGS for all 1,092 individuals

Back in the terminal:

```bash
cd ~/Sociogenomics

plink --bfile ~/Sociogenomics/Data/1kg_hm3 \
      --score ~/Sociogenomics/Results/Trait2_PRSice.snp 1 2 4 header \
      --out ~/Sociogenomics/Results/Trait2_pgs_all_pops
```

### 4.2 Analyse in R

```r
pgs_all <- read.table("~/Sociogenomics/Results/Trait2_pgs_all_pops.profile",
                      header = TRUE)
pheno_all <- read.table("~/Sociogenomics/Data/1kg.Trait2.phen",
                        header = FALSE,
                        col.names = c("FID", "IID", "Trait2"))
pop <- read.table("~/Sociogenomics/Data/1kg-sample-2504-phased.txt",
                  header = TRUE)

da <- merge(pgs_all[, c("FID", "IID", "SCORE")], pheno_all,
            by = c("FID", "IID"))
da <- merge(da, pop[, c("sample", "super_pop")],
            by.x = "IID", by.y = "sample")

# R2 by super-population
cat("R2 by super-population:\n")
for (p in c("EUR", "EAS", "AFR", "AMR")) {
  sub <- da[da$super_pop == p, ]
  r2 <- cor(sub$SCORE, sub$Trait2)^2
  cat(sprintf("  %s: R2 = %.4f (N = %d)\n", p, r2, nrow(sub)))
}
```

### 4.3 Boxplot by population

```r
da$PGS <- scale(da$SCORE)

boxplot(PGS ~ super_pop, data = da,
        col = c("tomato", "gold", "skyblue", "plum"),
        main = "PGS distribution by super-population",
        xlab = "Super-population", ylab = "Standardised PGS")
```

### Exercise 3

1. In which population is the PGS most predictive? Least predictive?
2. Is there a systematic shift in PGS means across populations? Why?
3. Why is the PGS less accurate in non-European populations?
4. What are possible solutions? (multi-ancestry GWAS, PRS-CSx, ...)

---

## Part V. Gene $\times$ Sex interaction with PGS

Does the PGS effect differ between males and females?

### 5.1 Prepare the data

```r
# Use the EUR data from Part III (object d)
fam <- read.table("~/Sociogenomics/Data/1kg_hm3_QC_CEU.fam",
                  col.names = c("FID", "IID", "fa", "mo", "sex", "ph"))
d <- merge(d, fam[, c("FID", "IID", "sex")], by = c("FID", "IID"))
d$female <- ifelse(d$sex == 2, 1, 0)
```

### 5.2 Fit the interaction model

$$Y_i = \beta_0 + \beta_1\,\text{PGS}_i + \beta_2\,\text{Female}_i + \beta_3\,(\text{PGS}_i \times \text{Female}_i) + \text{PCs} + \varepsilon_i$$

```r
mod_gxe <- lm(Trait2 ~ PGS * female +
              PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = d)
summary(mod_gxe)
```

### 5.3 Visualise

```r
library(ggplot2)
d$Sex <- ifelse(d$female == 1, "Female", "Male")

ggplot(d, aes(x = PGS, y = Trait2, colour = Sex)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "PGS x Sex interaction",
       x = "PGS (standardised)", y = "Trait2") +
  theme_minimal()
```

### 5.4 Interpretation

* $\beta_3 > 0$: **reinforcing** --- PGS effect is larger in females
* $\beta_3 < 0$: **compensating** --- PGS effect is attenuated in females
* $\beta_3 \approx 0$: no interaction --- genetic effect is the same in both sexes

### Exercise 4

1. Report $\hat{\beta}_3$, its SE and $p$-value. Is the interaction significant?
2. Do the regression lines in the plot have different slopes?
3. Why is testing GxE important from a public health perspective?
4. **Bonus:** repeat with a continuous moderator (e.g.\ PC1 as a proxy for ancestry gradient within EUR).

---

## Summary

| Step | Tool | What you learned |
|------|------|------------------|
| Inspect data | bash | Structure of GWAS summary stats |
| Run PRSice-2 | `Rscript PRSice.R` | Automated C+T pipeline |
| Analyse PGS | R | Distribution, incremental $R^2$, threshold comparison |
| Bootstrap CI | R (`boot`) | Uncertainty quantification |
| Cross-ancestry | PLINK + R | PGS portability limitations |
| GxE interaction | R | PGS $\times$ environment |

**Key messages:**

* PRSice-2 automates C+T in one command: clumping, scoring, regression, and plots
* Always control for ancestry PCs and report **incremental** $R^2$
* The optimal $p$-value threshold is trait-specific
* PGS trained in EUR transfer poorly to non-EUR populations
* GxE interactions reveal whether genetic effects depend on the environment
