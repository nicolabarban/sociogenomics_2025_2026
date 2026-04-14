# Lab 7. Gene × Environment Interaction

In this lab we will learn how to detect and analyse **Gene × Environment (GxE) interactions** --- situations in which the effect of a genetic variant on a trait depends on an environmental exposure.

We will use:

* `bpdata.csv` (1,000 individuals, 11 candidate SNPs, blood pressure, sex, BMI)
* The PGS for `Trait2` computed in Lab 6
* R (Google Colab or RStudio)

The lab covers:

* Testing **single-SNP × environment** interactions
* Testing **PGS × environment** interactions
* Interpreting and visualising interactions
* Recognising the **diathesis-stress** and **differential susceptibility** patterns

---

## 0. Getting started

Update the course repository:

```bash
cd ~/sociogenomics_2025_2026
git pull
git lfs pull
```

Make sure the working directories exist:

```bash
cd $HOME
mkdir -p ~/Sociogenomics/Data ~/Sociogenomics/Results ~/Sociogenomics/Scripts
```

Copy the data we need from the course repo:

```bash
cp ~/sociogenomics_2025_2026/data/week7/bpdata.csv ~/Sociogenomics/Data/

# If not already done in Lab 6, copy also the phenotype and population files
cp ~/sociogenomics_2025_2026/data/week6/1kg.Trait2.phen ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/week6/1kg-sample-2504-phased.txt ~/Sociogenomics/Data/
```

> **Note:** Most of this lab is in R. You can use Google Colab (R runtime) or RStudio.

---

## Part I. The GxE model

### The statistical model

Gene $\times$ Environment interaction is tested with a regression model of the form:

$$y_i = \beta_0 + \beta_G\,G_i + \beta_E\,E_i + \beta_{GE}\,G_i \times E_i + \varepsilon_i$$

Where:

* $y_i$: phenotype
* $G_i$: genotype (0, 1, 2 for SNP count; or continuous for PGS)
* $E_i$: environmental exposure (sex, BMI, SES, \ldots)
* $\beta_{GE}$: the **interaction coefficient** --- the quantity of interest

A significant $\beta_{GE}$ means the genetic effect **varies** across environments.

### Two theoretical frames

* **Diathesis-stress:** genes amplify vulnerability to adverse environments
* **Differential susceptibility** (Belsky): genes increase sensitivity to *both* adverse *and* positive environments --- the ``for better and for worse'' model

We will see both patterns in today's lab.

---

## Part II. Load and inspect the data

### Load `bpdata.csv`

In R:

```r
# Load the blood pressure data
bp <- read.csv("~/Sociogenomics/Data/bpdata.csv",
               na.strings = c("NA", ""))
head(bp)
dim(bp)
summary(bp)
```

The dataset contains:

| Variable | Description |
|----------|-------------|
| `sex` | MALE / FEMALE |
| `sbp` | systolic blood pressure (mmHg) |
| `dbp` | diastolic blood pressure (mmHg) |
| `bmi` | body mass index (kg/m$^2$) |
| `snp1`--`snp11` | candidate SNPs (character genotypes) |

### Recode genotypes to numeric (allele counts)

Most analyses need genotypes as 0/1/2:

```r
# Function to recode a single SNP to 0/1/2 based on the minor allele
recode_snp <- function(geno) {
  g <- na.omit(geno)
  all_alleles <- unlist(strsplit(g, ""))
  tab <- sort(table(all_alleles))
  minor <- names(tab)[1]   # least frequent allele
  sapply(geno, function(x) {
    if (is.na(x)) return(NA)
    sum(strsplit(x, "")[[1]] == minor)
  })
}

snp_cols <- grep("^snp", names(bp), value = TRUE)
for (s in snp_cols) {
  bp[[paste0(s, "_n")]] <- recode_snp(bp[[s]])
}

# Check the first recoded SNP
head(bp[, c("snp1", "snp1_n")])
```

**Question:** What are the genotype counts (0/1/2) for `snp1`? What is its minor allele frequency?

```r
table(bp$snp1_n)
mean(bp$snp1_n, na.rm = TRUE) / 2   # MAF
```

---

## Part III. Single-SNP × Environment interactions

### Marginal SNP effects

Let's first check which SNPs are associated with systolic blood pressure, controlling for sex and BMI:

```r
# Loop over 11 SNPs, fit linear model sbp ~ snp + bmi + sex
results <- data.frame(snp = character(), beta = numeric(),
                      se = numeric(), p = numeric())

for (s in snp_cols) {
  snp_n <- paste0(s, "_n")
  mod <- lm(sbp ~ bp[[snp_n]] + bmi + sex, data = bp)
  coef_snp <- summary(mod)$coefficients[2, ]
  results <- rbind(results,
                   data.frame(snp = s,
                              beta = coef_snp[1],
                              se = coef_snp[2],
                              p = coef_snp[4]))
}

print(results[order(results$p), ])
```

**Question:** Which SNP(s) show the strongest marginal association with SBP?

### Interaction with sex

Test whether the effect of the top SNP differs between males and females:

```r
top_snp <- results$snp[which.min(results$p)]
top_snp_n <- paste0(top_snp, "_n")

mod_gxe_sex <- lm(sbp ~ bp[[top_snp_n]] * sex + bmi, data = bp)
summary(mod_gxe_sex)
```

**Question:** Is the interaction coefficient (`bp[[top_snp_n]]:sexMALE`) significant?

### Interaction with BMI (continuous environment)

```r
# Center BMI for easier interpretation
bp$bmi_c <- bp$bmi - mean(bp$bmi, na.rm = TRUE)

mod_gxe_bmi <- lm(sbp ~ bp[[top_snp_n]] * bmi_c + sex, data = bp)
summary(mod_gxe_bmi)
```

The interaction coefficient tells us: *by how much does the SNP effect change per 1 kg/m$^2$ increase in BMI?*

### Visualise the interaction

```r
library(ggplot2)

# Rename the SNP column to something predictable for the plot
bp$snp <- bp[[top_snp_n]]

# Prediction grid at BMI = -5, 0, +5 (centred)
grid <- expand.grid(snp = 0:2,
                    bmi_c = c(-5, 0, 5),
                    sex = "FEMALE")
grid$pred <- predict(lm(sbp ~ snp * bmi_c + sex, data = bp),
                     newdata = grid)

ggplot(grid, aes(x = snp, y = pred, colour = factor(bmi_c))) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  labs(x = "Allele count", y = "Predicted SBP (mmHg)",
       colour = "BMI (centred)",
       title = paste("GxE:", top_snp, "× BMI interaction")) +
  theme_minimal()
```

### Exercise 1

1. Loop over all 11 SNPs and test each for interaction with BMI. Report the top 3 SNPs by interaction $p$-value.
2. What is the Bonferroni-corrected significance threshold for 11 SNPs $\times$ 2 environments (sex, BMI)?
3. Repeat the analysis for `dbp` (diastolic blood pressure). Are the results consistent with `sbp`?

---

## Part IV. PGS × Environment interactions

Testing individual SNPs $\times$ environment is underpowered. A more powerful approach uses the **polygenic score** as a summary of genetic predisposition.

### Load the PGS from Lab 6

If you have the `Trait2_PRSice.all_score` file from Lab 6, load it:

```r
# Load all scores at different thresholds
all_scores <- read.table("~/Sociogenomics/Results/Trait2_PRSice.all_score",
                         header = TRUE)
head(all_scores)
```

Load the phenotype and population info:

```r
pheno <- read.table("~/Sociogenomics/Data/1kg.Trait2.phen",
                    header = FALSE,
                    col.names = c("FID", "IID", "Trait2"))

pop <- read.table("~/Sociogenomics/Data/1kg-sample-2504-phased.txt",
                  header = TRUE)

# Merge
d <- merge(all_scores, pheno, by = c("FID", "IID"))
d <- merge(d, pop[, c("sample", "super_pop")],
           by.x = "IID", by.y = "sample")

# Standardise PGS at the "best" threshold from Lab 6
d$PGS <- scale(d$Pt_0.5)   # adjust threshold if different
```

### Main effect of PGS

```r
mod_main <- lm(Trait2 ~ PGS, data = d)
summary(mod_main)$coefficients
summary(mod_main)$r.squared
```

### PGS × super-population interaction

The environmental variable here is super-population (ancestry). Does the PGS effect vary across ancestries?

```r
mod_pgs_pop <- lm(Trait2 ~ PGS * super_pop, data = d)
summary(mod_pgs_pop)

# Compare nested models (F-test on the interaction)
anova(mod_main, mod_pgs_pop)
```

**Question:** Is there significant variation in the PGS effect across super-populations? What does this tell us about **PGS portability**?

### Visualise the interaction

```r
library(ggplot2)

ggplot(d, aes(x = PGS, y = Trait2, colour = super_pop)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Standardised PGS", y = "Trait2",
       colour = "Super-population",
       title = "PGS × Ancestry interaction") +
  theme_minimal()
```

### Slope estimates by population

```r
library(dplyr)

slopes <- d |>
  group_by(super_pop) |>
  summarise(beta = coef(lm(Trait2 ~ PGS))[2],
            se = summary(lm(Trait2 ~ PGS))$coefficients[2, 2],
            n = n())

print(slopes)
```

**Question:** In which population is the PGS effect strongest? Weakest? How does this compare to the cross-ancestry portability discussed in the lecture?

### Exercise 2

1. Fit the same PGS × super-population model on `Trait1`. Are the interactions as strong?
2. Plot the PGS distribution by super-population. Is it centred at 0 in all groups?
3. How does failing to detect interaction across ancestries affect clinical use of PGS?

---

## Part V. A common pitfall: scale dependence

A significant GxE can appear (or disappear!) simply by transforming the outcome.

### Compare linear and log-transformed outcome

```r
# Original scale
mod1 <- lm(sbp ~ bp[[top_snp_n]] * bmi_c + sex, data = bp)
coef_orig <- summary(mod1)$coefficients

# Log-transformed outcome
mod2 <- lm(log(sbp) ~ bp[[top_snp_n]] * bmi_c + sex, data = bp)
coef_log <- summary(mod2)$coefficients

print(coef_orig)
print(coef_log)
```

**Question:** Does the interaction coefficient have the same sign on both scales? The same significance?

### Takeaway

* GxE results **depend on the scale** of the outcome
* Always justify the choice of scale (linear, log, rank-normalised, \ldots)
* Report interactions on the scale that has substantive meaning
* **Replication** in an independent dataset is the only way to trust a GxE finding

---

## Part VI. Diathesis-stress vs differential susceptibility (optional)

Suppose your top GxE result shows that the slope of the SNP on SBP is:

| BMI group | Slope |
|-----------|-------|
| Low | 0 or small |
| Mid | moderate |
| High | large |

This is a **diathesis-stress** pattern: the SNP only matters in an adverse environment (high BMI).

In contrast, if low-BMI individuals show a *negative* slope, it's **differential susceptibility**: the genotype makes people more sensitive to both positive and negative environments.

### Stratified analysis

```r
# Stratify by BMI tertile
bp$bmi_tertile <- cut(bp$bmi,
                      breaks = quantile(bp$bmi,
                                        c(0, 1/3, 2/3, 1),
                                        na.rm = TRUE),
                      include.lowest = TRUE,
                      labels = c("Low", "Mid", "High"))

# Slope of top SNP by BMI tertile
for (t in c("Low", "Mid", "High")) {
  sub <- bp[bp$bmi_tertile == t, ]
  mod <- lm(sbp ~ sub[[top_snp_n]] + sex, data = sub)
  beta <- coef(mod)[2]
  se   <- summary(mod)$coefficients[2, 2]
  cat(sprintf("%-4s BMI: beta = %.3f (SE = %.3f)\n", t, beta, se))
}
```

**Question:** Which theoretical model does your result match? Diathesis-stress or differential susceptibility?

### Exercise 3 (optional)

1. Apply the same tertile-based analysis to the PGS $\times$ super-population interaction. What pattern do you see?
2. Propose a sociological example where GxE might follow differential susceptibility (e.g., education PGS × school quality).
3. What sample size would you need to reliably detect a GxE interaction of $R^2 = 0.2\%$? (Hint: power calculations assume $\alpha = 0.05$, power = 0.8).

---

## Summary

In this lab we learned how to:

| Step | Tool | Command / concept |
|------|------|-------------------|
| Recode character genotypes | R | `strsplit` + allele counts |
| Marginal SNP test | R | `lm(y ~ snp + cov)` |
| GxE regression | R | `lm(y ~ snp * E + cov)` |
| Interaction plot | ggplot2 | `geom_smooth(method = "lm")` |
| PGS $\times$ E | R | `lm(y ~ PGS * E + cov)` |
| Stratified slopes | R | loop over groups |

**Key messages:**

* GxE requires **very large samples** --- single-study results are often not replicable
* The **scale** of the outcome matters
* **PGS × ancestry** interactions reveal the portability problem
* Always Bonferroni-correct for the number of SNPs $\times$ environments
* Replication is the gold standard

---

## Further reading

* Thompson et al.\ (2022) *The Genetic Lottery* --- chapters on GxE
* Belsky \& Pluess (2009) --- differential susceptibility
* Manuck \& McCaffery (2014) *Ann Rev Psychol* --- GxE methodology
* Duncan \& Keller (2011) *AJP* --- why candidate-gene GxE often fails
