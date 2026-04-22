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

### 2.1 Run at fixed $p$-value thresholds

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

### 2.2 Inspect the output

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

## Part V. Gene $\times$ Sex interaction

Does the PGS effect differ between males and females?

```r
# Use the EUR data from Part III (object d)
fam <- read.table("Data/1kg_hm3_QC_CEU.fam",
                  col.names = c("FID","IID","fa","mo","sex","ph"))
d <- merge(d, fam[, c("FID","IID","sex")], by = c("FID","IID"))
d$female <- ifelse(d$sex == 2, 1, 0)

# Interaction model
mod_gxe <- lm(Trait2 ~ PGS * female +
              PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = d)
summary(mod_gxe)

# Plot
library(ggplot2)
d$Sex <- ifelse(d$female == 1, "Female", "Male")
ggplot(d, aes(x = PGS, y = Trait2, colour = Sex)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "PGS x Sex interaction",
       x = "PGS (standardised)", y = "Trait2") +
  theme_minimal()
```

**Interpretation:**
* $\beta_3 > 0$: PGS effect **larger** in females (reinforcing)
* $\beta_3 < 0$: PGS effect **attenuated** in females (compensating)
* $\beta_3 \approx 0$: no interaction

### Exercise 4

1. Report $\hat{\beta}_3$, SE, $p$-value. Is the interaction significant?
2. Do the regression lines have different slopes?
3. **Bonus:** repeat replacing sex with PC1 (ancestry gradient within EUR).

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `Rscript: command not found` | Run `sudo apt-get install -y r-base` |
| `./plink: Permission denied` | Run `chmod +x plink` |
| `./PRSice_linux: Permission denied` | Run `chmod +x PRSice_linux` |
| R package install fails | Try `install.packages("pkg", repos="https://cloud.r-project.org")` |
| `Error: cannot open connection` (file not found) | Check you are in `~/Sociogenomics` with `pwd` |
| PRSice produces empty output | Check `Trait2.ma` column names match `--snp SNP --A1 A1 --pvalue P` |
| Barplot PNG not visible | Download via Cloud Shell menu (⋮ → Download file) |
| Cloud Shell session expired, R gone | Run `sudo apt-get install -y r-base` again |

---

## Summary

| Step | Tool | What you learned |
|------|------|------------------|
| Run PRSice-2 | `Rscript PRSice.R` | Automated C+T pipeline |
| Analyse PGS | R | Distribution, incremental $R^2$, threshold comparison |
| Bootstrap CI | R (`boot`) | Uncertainty quantification |
| Cross-ancestry | PLINK + R | PGS portability limitations |
| GxE interaction | R + ggplot2 | PGS $\times$ environment |

**Key messages:**

* PRSice-2 automates C+T in one command
* Always control for ancestry PCs and report **incremental** $R^2$
* The optimal $p$-value threshold is trait-specific
* PGS trained in EUR transfer poorly to non-EUR populations
* GxE interactions reveal whether genetic effects depend on the environment
