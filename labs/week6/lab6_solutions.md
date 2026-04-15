# Lab 6. Polygenic Scores — Solutions

Solutions and expected outputs for the exercises in [Lab 6](lab6.md). Numbers may vary slightly with different random seeds or PLINK versions.

---

## Part II — Reusing QC files from Week 5

**Question:** How many SNPs and how many European individuals are in the QC'd file?

```bash
wc -l ~/Sociogenomics/Data/1kg_hm3_QC_CEU.bim
# 761583
wc -l ~/Sociogenomics/Data/1kg_hm3_QC_CEU.fam
# 374
```

$\approx$ **761 K SNPs**, **374 European individuals**.

---

## Part III — PGS with PRSice-2 (Exercise 2)

After running PRSice, the file `Trait2_PRSice.summary` reports:

| Field | Value |
|---|---|
| Best threshold | **$p < 0.5$** |
| $R^2$ (PRSice) | **0.144** (14.4%) |
| $N$ SNPs | **2,089** |
| $p$-value of coefficient | $\sim 9 \times 10^{-14}$ |

And `Trait2_PRSice.prsice` shows the full grid:

| Threshold | $R^2$ | $N$ SNPs |
|---|---|---|
| $5 \times 10^{-8}$ | 0.118 | 738 |
| $5 \times 10^{-6}$ | 0.131 | 943 |
| $5 \times 10^{-4}$ | 0.139 | 1,235 |
| $0.05$ | 0.142 | 1,695 |
| **$0.5$** | **0.144** | **2,089** |
| $1$ | 0.143 | 2,250 |

### 1. What is the best $p$-value threshold for Trait2?

**$p < 0.5$**.

### 2. How many SNPs are included at that threshold?

**2,089**.

### 3. What is the $R^2$ of the best PGS?

**$\approx 14.4\%$** (unadjusted --- the $R^2$ reported by PRSice already controls for the intercept but not for PCs).

### 4. Compare the mean and standard deviation of the PGS from PRSice and PLINK

After computing both, the two PGS columns are almost perfectly correlated:

```r
cor(prsice_pgs$Pt_0.5, plink_pgs$SCORE)
# 0.99989
```

They differ only by a scaling factor (PLINK divides by the number of non-missing alleles per individual, PRSice averages differently), so the means differ by about 10%, but the regression coefficient and the $R^2$ are identical.

---

## Part IV — Analyse the PGS in R (Exercise 3)

### 1. What is the incremental $R^2$ of the best PGS?

```r
mod0 <- lm(Trait2 ~ PC1 + PC2 + ... + PC10, data = d)
mod1 <- lm(Trait2 ~ PGS_best + PC1 + PC2 + ... + PC10, data = d)
summary(mod1)$r.squared - summary(mod0)$r.squared
# 0.147
```

$\Delta R^2 \approx$ **14.7%** (slightly higher than PRSice's 14.4% because we also adjust for ancestry PCs).

### 2. What is the 95% bootstrap CI?

```r
set.seed(12345)
boot.ci(boot_results, type = "norm")
# 95% CI: (0.091, 0.211)
```

$\Delta R^2 \in$ **(9.1%, 21.1%)** with 95% confidence.

### 3. Does including more SNPs always improve $R^2$? What is the optimal threshold?

The incremental $R^2$ at each threshold (controlling for PCs):

| Threshold | $\Delta R^2$ |
|---|---|
| $5 \times 10^{-8}$ | 0.121 |
| $5 \times 10^{-6}$ | 0.136 |
| $5 \times 10^{-4}$ | 0.145 |
| $0.05$ | 0.145 |
| **$0.5$** | **0.147** |
| $1$ | 0.147 |

The $R^2$ **plateaus** at $p < 0.5$. Adding all remaining SNPs ($p \leq 1$) gives a tiny additional gain ($<$0.1%), consistent with the idea that `Trait2` is polygenic: most of the heritability is captured by SNPs with $p < 0.5$, not just genome-wide significant ones.

### 4. Repeat the analysis for `Trait1`. Why might the result be different?

`Trait1` shows almost no association: summary stats have no SNPs at $p < 5 \times 10^{-8}$, and PGS $R^2 \approx 0$ at any threshold.

Why? Because `Trait1` is simulated as a **null trait** (no genetic signal) --- it is included in the SISG 2024 dataset precisely to show what happens when the true causal effect is zero. The PGS picks up only noise.

---

## Part V — Cross-ancestry portability (Exercise 4, optional)

After computing the PGS on all 2,504 individuals and comparing $R^2$ by super-population:

| Population | $R^2$ | $N$ | Relative accuracy (vs.\ EUR) |
|---|---|---|---|
| **EUR** | 0.145 | 364 | 100% |
| EAS | 0.080 | 276 | 55% |
| AMR | 0.079 | 172 | 54% |
| AFR | 0.046 | 221 | 32% |
| SAS | NA | 0 | (not in this subset) |

### 1. Boxplot of PGS values by super-population --- any systematic shift?

```r
library(ggplot2)
ggplot(da, aes(x = super_pop, y = SCORE, fill = super_pop)) +
  geom_boxplot() + theme_minimal()
```

Yes: the **mean PGS differs across populations**, because allele frequencies differ. This is a classic manifestation of population stratification --- the PGS values for African-ancestry individuals are systematically shifted relative to Europeans, even before looking at the phenotype.

### 2. Why is the PGS less accurate in African-ancestry populations?

Several compounding reasons:

* **LD differences:** SNPs tagging the causal variant in Europeans may not tag it (or may tag a different one) in Africans. LD blocks are shorter in African populations.
* **Allele frequencies differ**, weakening the fit.
* **Possible causal effect heterogeneity** across populations (GxE).
* **Imperfect portability of the discovery GWAS**: effects estimated in Europeans may not transfer.

### 3. What are possible solutions?

* **Multi-ancestry GWAS:** larger GWAS in non-European populations (All of Us, PAGE, BioBank Japan, GenomeAsia).
* **Trans-ethnic meta-analysis** combining ancestries.
* **Multi-ancestry PRS methods:** PRS-CSx, DendroPRS, Bayesian-MTAG.
* **LD-aware methods** with population-specific reference panels.
* **Long-term goal:** diverse and representative genomic datasets so that PGS does not exacerbate health inequalities.

---

## Summary of expected results

| Metric | Value |
|---|---|
| SNPs after QC (EUR) | $\sim$761,000 |
| Individuals (EUR, QC) | 374 |
| Best $p$-value threshold | 0.5 |
| SNPs at best threshold | 2,089 |
| PGS $R^2$ (PRSice) | 14.4% |
| Incremental $\Delta R^2$ (with PCs) | 14.7% |
| 95% Bootstrap CI | (9.1%, 21.1%) |
| Cross-ancestry $R^2$ (AFR / EUR) | 32% |
