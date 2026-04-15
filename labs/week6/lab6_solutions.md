# Lab 6. Polygenic Scores — Solutions

Solutions and expected outputs for the exercises in [Lab 6](lab6.md). Numbers may vary slightly with different random seeds or PLINK versions.

---

## Part II — Reusing QC files from Week 5

**Question:** How many SNPs and how many European individuals are in the QC'd file?

```bash
wc -l ~/Sociogenomics/Data/1kg_hm3_QC_CEU.bim   # 761 583
wc -l ~/Sociogenomics/Data/1kg_hm3_QC_CEU.fam   #    374
```

$\approx$ **761 K SNPs**, **374 European individuals**.

---

## Part III — Monogenic FTO score

Genotype distribution of the FTO score:

| SCORE | A-allele count | $N$ |
|---|---|---|
| 0.0 | 0 | 87 |
| 0.2 | 1 | 160 |
| 0.4 | 2 | 127 |

Each individual's `SCORE` is simply $0.4 \times$ their A-allele count at rs9930506. There is no regression exercise here --- this section is just a warm-up to see `--score` in action on a single SNP.

---

## Part IV — C+T step by step with PLINK (Exercise 1)

### 1. How many ambiguous SNPs were removed?

**2,355** out of 32,260 SNPs in `Trait2.ma` have A/T or C/G alleles and were filtered out.

### 2. How many SNPs required strand flipping?

**Zero** in this dataset. Both the base summary statistics (UK Biobank) and the target genotypes (1000 Genomes) are reported on the forward strand from the same reference build. In real analyses --- especially when combining data from different consortia and platforms --- a non-trivial fraction usually requires flipping.

### 3. How many clumps did PLINK produce?

**$\sim$2,249 clumps** from $\sim$17,100 variants that match between target and cleaned base.

```bash
wc -l ~/Sociogenomics/Results/Trait2_clumped.clumped
# 2252 (includes header + blank line)
```

### 4. How many SNPs pass $p < 5 \times 10^{-8}$ among the clumped SNPs?

```bash
wc -l ~/Sociogenomics/Results/score_5e8.txt
# 737
```

**737 independent genome-wide significant SNPs**.

### 5. Does a looser threshold change the PGS distribution?

Yes. At $p < 0.05$ there are $\sim$2,000 SNPs, the PGS variance is larger, and the correlation with `Trait2` is typically higher. This is the classical polygenic signal --- many small effects together explain more variance.

---

## Part V — PGS with PRSice-2 (Exercise 2)

After running PRSice, the file `Trait2_PRSice.summary` reports:

| Field | Value |
|---|---|
| Best threshold | **$p < 0.5$** |
| $R^2$ (PRSice) | **0.144** (14.4%) |
| $N$ SNPs | **2,089** |
| $p$-value of coefficient | $\sim 9 \times 10^{-14}$ |

The full grid (`Trait2_PRSice.prsice`):

| Threshold | $R^2$ | $N$ SNPs |
|---|---|---|
| $5 \times 10^{-8}$ | 0.118 | 738 |
| $5 \times 10^{-6}$ | 0.131 | 943 |
| $5 \times 10^{-4}$ | 0.139 | 1,235 |
| $0.05$ | 0.142 | 1,695 |
| **$0.5$** | **0.144** | **2,089** |
| $1$ | 0.143 | 2,250 |

### 1. Best $p$-value threshold

**$p < 0.5$**.

### 2. Number of SNPs at that threshold

**2,089**.

### 3. $R^2$ of the best PGS

**$\approx 14.4\%$**.

### 4. PRSice vs.\ manual PLINK `--score` at the same threshold

The PGS from PRSice at $p < 5 \times 10^{-8}$ (`Pt_5e.08` column in `*.all_score`) and the PLINK `--score` output using `score_5e8.txt` are **nearly perfectly correlated**:

```r
cor(prsice_pgs$Pt_5e.08, plink_pgs$SCORE)
# ~ 1.00
```

Mean values differ by about 10% because PLINK normalises by allele count per individual; PRSice sums without the same normalisation. The regression coefficient and $R^2$ are identical.

---

## Part VI — Analyse the PGS in R (Exercise 3)

### 1. Incremental $R^2$

```r
mod0 <- lm(Trait2 ~ PC1 + ... + PC10, data = d)
mod1 <- lm(Trait2 ~ PGS_best + PC1 + ... + PC10, data = d)
summary(mod1)$r.squared - summary(mod0)$r.squared
# 0.147
```

$\Delta R^2 \approx$ **14.7%** (slightly higher than PRSice's 14.4% because we also adjust for PCs).

### 2. Bootstrap 95% CI

```r
set.seed(12345)
boot.ci(boot_results, type = "norm")
# (0.091, 0.211)
```

$\Delta R^2 \in$ **(9.1%, 21.1%)**.

### 3. Does including more SNPs always improve $R^2$?

No. The incremental $R^2$ **plateaus** around $p < 0.5$:

| Threshold | $\Delta R^2$ |
|---|---|
| $5 \times 10^{-8}$ | 0.121 |
| $5 \times 10^{-6}$ | 0.136 |
| $5 \times 10^{-4}$ | 0.145 |
| $0.05$ | 0.145 |
| **$0.5$** | **0.147** |
| $1$ | 0.147 |

The optimal threshold here is $p < 0.5$. Adding *all* remaining SNPs ($p \leq 1$) gives a negligible gain --- consistent with polygenicity.

### 4. `Trait1`: why so different?

`Trait1` is simulated as a **null trait** --- no genetic signal. Summary stats show no SNPs at $p < 5 \times 10^{-8}$ and PGS $R^2 \approx 0$ at any threshold. It is included as a negative control: with a true causal effect of zero, the PGS can only pick up noise.

---

## Part VII — Cross-ancestry portability (Exercise 4, optional)

Applying the score file to all ancestries and computing $R^2$ by super-population:

| Population | $R^2$ | $N$ | Relative accuracy (vs.\ EUR) |
|---|---|---|---|
| **EUR** | 0.145 | 364 | 100% |
| EAS | 0.080 | 276 | 55% |
| AMR | 0.079 | 172 | 54% |
| AFR | 0.046 | 221 | 32% |
| SAS | NA | 0 | (not in this subset) |

### 1. Systematic shift in PGS by super-population?

Yes. Allele frequencies differ across populations $\to$ the mean PGS differs too. A boxplot by `super_pop` shows clear shifts, even *before* looking at phenotype.

### 2. Why is the PGS less accurate in African-ancestry populations?

* **LD differences:** shorter LD blocks in African populations; SNPs tagging the causal variant in Europeans may not tag it in Africans.
* **Allele frequencies differ**, weakening the match.
* **Possible effect-size heterogeneity** (GxE).
* **Discovery GWAS is European-ancestry** $\to$ biased weights.

### 3. Possible solutions

* **Multi-ancestry GWAS** (All of Us, PAGE, BBJ, GenomeAsia)
* **Trans-ethnic meta-analysis**
* **Multi-ancestry PRS methods** (PRS-CSx, DendroPRS)
* **Population-specific LD reference panels**
* **Long-term goal:** diverse datasets $\to$ PGS that does not exacerbate health disparities.

---

## Summary of expected results

| Metric | Value |
|---|---|
| SNPs after QC (EUR) | $\sim$761,000 |
| Individuals (EUR, QC) | 374 |
| Clumps (PLINK `--clump`) | $\sim$2,249 |
| SNPs at $p < 5 \times 10^{-8}$ (after clumping) | 737 |
| Best $p$-value threshold (PRSice) | 0.5 |
| SNPs at best threshold | 2,089 |
| PGS $R^2$ (PRSice) | 14.4% |
| Incremental $\Delta R^2$ (with PCs) | 14.7% |
| 95% bootstrap CI | (9.1%, 21.1%) |
| Cross-ancestry $R^2$ (AFR / EUR) | 32% |
