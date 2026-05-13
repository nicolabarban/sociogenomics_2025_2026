---
title: "Lab Recap & Final Report Pipeline"
---

# Lab Recap and Final Report Pipeline

We use the `1kg_height` PLINK fileset (1.092 individuals from 1000G, 851k SNPs, with simulated height in the phenotype column of the `.fam`) and `1kg_samples.txt` for the super-population labels.

## Part A — Ancestry and population structure (Lab 4)

### Step 1. LD-prune SNPs (for a clean PCA)

```bash
plink --bfile 1kg_height \
      --maf 0.01 --geno 0.05 \
      --indep-pairwise 200 50 0.2 \
      --out work_prune
```

### Step 2. Run PCA on the pruned set

```bash
plink --bfile 1kg_height \
      --extract work_prune.prune.in \
      --pca 10 \
      --out work_pca

head -3 work_pca.eigenvec       # FID IID PC1 PC2 ... PC10
cat work_pca.eigenval           # eigenvalues
```

### Step 3. Classify ancestry in R

```r
library(ggplot2)
library(dplyr)

ev <- read.table("work_pca.eigenvec",
                 header = FALSE,
                 col.names = c("FID", "IID", paste0("PC", 1:10)))

# 1000G super-population labels
geo <- read.table("1kg_samples.txt", sep = "\t", header = TRUE)

# Merge: 1kg_samples.txt has 'Sample.name' — match on IID
ev <- merge(ev, geo[, c("Sample.name", "Superpopulation.code")],
            by.x = "IID", by.y = "Sample.name", all.x = TRUE)
names(ev)[names(ev) == "Superpopulation.code"] <- "superpop"
ev$is_ref <- !is.na(ev$superpop)

# PCA cloud coloured by known super-population
ggplot(ev, aes(PC1, PC2, colour = superpop)) +
  geom_point(alpha = 0.6) +
  scale_colour_manual(values = c("EUR" = "#1F3A5F", "AFR" = "#C0392B",
                                 "EAS" = "#2A9D8F", "AMR" = "#E9A23B",
                                 "SAS" = "#8E44AD")) +
  labs(title = "PCA with 1000G super-population labels",
       x = "PC1", y = "PC2") +
  theme_minimal()

# k-means clustering on PC1..PC5 with 4 clusters
# (1kg_height contains AFR, AMR, EAS, EUR — no SAS)
pc_cols <- paste0("PC", 1:5)

set.seed(42)
km <- kmeans(ev[, pc_cols], centers = 4, nstart = 25)
ev$cluster <- km$cluster

# Label each cluster with the most-frequent known super-population.
cluster_to_pop <- ev %>%
  filter(is_ref) %>%
  count(cluster, superpop) %>%
  group_by(cluster) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  select(cluster, superpop) %>%
  rename(pred_pop = superpop)

ev <- merge(ev, cluster_to_pop, by = "cluster", all.x = TRUE)
table(ev$pred_pop)

# Save EUR-like individuals in PLINK keep format
eur_keep <- subset(ev, pred_pop == "EUR", select = c("FID", "IID"))
write.table(eur_keep, "samples_EUR_like.txt",
            sep = " ", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
cat("EUR retained:", nrow(eur_keep), "\n")

# Classified scatter for the report
ggplot(ev, aes(PC1, PC2, colour = pred_pop)) +
  geom_point(alpha = 0.6) +
  scale_colour_manual(values = c("EUR" = "#1F3A5F", "AFR" = "#C0392B",
                                 "EAS" = "#2A9D8F", "AMR" = "#E9A23B",
                                 "SAS" = "#8E44AD")) +
  labs(title = "Predicted super-population (k-means)", colour = NULL) +
  theme_minimal()
```

**Required deliverable:** a PCA plot coloured by 1000G super-population, the k-means-predicted classification scatter, a count of how many individuals you classified as EUR, the file `samples_EUR_like.txt`, and 3–5 sentences on *why* a single-ancestry GWAS matters.

---

## Part B — Quality control (Lab 3)

Apply the canonical filters on **the EUR subset only**:

| Filter | PLINK flag | Purpose |
|--------|------------|---------|
| Minor allele frequency > 0.05 | `--maf 0.05` | Drop rare SNPs (low power, unstable allele frequencies) |
| Individual call rate > 95% | `--mind 0.05` | Drop badly-genotyped individuals |
| SNP call rate > 95% | `--geno 0.05` | Drop badly-genotyped SNPs |
| HWE p > 10⁻⁵ | `--hwe 1e-5` | Drop SNPs with implausible allele frequencies (genotyping errors) |
| Relatedness ϕ̂ < 0.1 | `--rel-cutoff 0.1` | Drop one of each related pair (independence assumption) |

`--rel-cutoff` cannot be combined with `--make-bed` in plink 1.9, so we do two passes:

```bash
# Pass 1: write the list of unrelated individuals
plink --bfile 1kg_height \
      --keep samples_EUR_like.txt \
      --maf 0.05 --mind 0.05 --geno 0.05 --hwe 1e-5 \
      --rel-cutoff 0.1 \
      --out eur_unrel

# Pass 2: produce the QC'd binary fileset
plink --bfile 1kg_height \
      --keep eur_unrel.rel.id \
      --maf 0.05 --mind 0.05 --geno 0.05 --hwe 1e-5 \
      --make-bed --out eur_qc

wc -l eur_qc.{bim,fam}
```

> **Common mistake:** running QC on the full multi-ancestry cohort and *then* subsetting to EUR. The HWE filter is the giveaway: mixed-ancestry samples look out-of-equilibrium even when each ancestry is fine. Always **first ancestry, then QC**.

### PCA on the QC'd EUR subset (for GWAS covariates)

The PCs we computed in Part A separate the *super-populations*. As covariates inside a single ancestry they are nearly collinear. Compute a fresh PCA on `eur_qc` and use **those** PCs in the GWAS:

```bash
plink --bfile eur_qc \
      --indep-pairwise 200 50 0.2 \
      --out eur_prune

plink --bfile eur_qc \
      --extract eur_prune.prune.in \
      --pca 10 \
      --out eur_pca
```

---

## Build a covariate file

`1kg_height` has the height phenotype in the `.fam` file (column 6) and sex in column 5. We need to add **age** as a covariate — `1kg_height` does not ship with age, so we synthesise a uniform age.

```r
fam <- read.table("eur_qc.fam", header = FALSE,
                  col.names = c("FID","IID","PAT","MAT","SEX","height"))

# Use sex from the .fam (1 = male, 2 = female; drop unknowns)
fam <- subset(fam, SEX %in% c(1, 2))
fam$sex <- fam$SEX

# Synthesise a plausible adult age (uniform 25-75)
set.seed(2026)
fam$age <- runif(nrow(fam), 25, 75)

# Merge in the within-EUR PCs
pcs <- read.table("eur_pca.eigenvec",
                  header = FALSE,
                  col.names = c("FID", "IID", paste0("PC", 1:10)))
cov_pc <- merge(fam[, c("FID","IID","age","sex","height")],
                pcs[, c("IID", paste0("PC", 1:10))], by = "IID")

write.table(cov_pc[, c("FID","IID","age","sex", paste0("PC", 1:10))],
            "covariates_with_pcs.txt",
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

cat("Covariate file rows:", nrow(cov_pc), "\n")
```

---

## Part C — GWAS (Lab 4)

Height is already in `eur_qc.fam` (column 6), so plink picks it up automatically as the phenotype.

### Step 1. GWAS without covariates

```bash
plink --bfile eur_qc \
      --chr 1-22 \
      --linear \
      --out gwas_nocov
```

### Step 2. GWAS with covariates

```bash
plink --bfile eur_qc \
      --chr 1-22 \
      --linear hide-covar \
      --covar covariates_with_pcs.txt \
      --covar-name age,sex,PC1-PC10 \
      --out gwas_cov
```

### Step 3. Top-10 SNPs in each model

```bash
# In a plink 1.9 .assoc.linear the P-value is column 9
echo "TOP 10 -- NO COVARIATES"
( head -1 gwas_nocov.assoc.linear ; \
  tail -n +2 gwas_nocov.assoc.linear | sort -g -k9,9 | head -10 ) | column -t

echo "TOP 10 -- WITH COVARIATES"
( head -1 gwas_cov.assoc.linear ; \
  tail -n +2 gwas_cov.assoc.linear | sort -g -k9,9 | head -10 ) | column -t
```

### Step 4. Manhattan plot

```r
library(qqman)

g <- read.table("gwas_cov.assoc.linear", header = TRUE)
g <- subset(g, !is.na(P))

cat("SNPs tested:", nrow(g), "\n")
cat("Genome-wide significant (p<5e-8):", sum(g$P < 5e-8), "\n")
cat("Suggestive (p<1e-5):", sum(g$P < 1e-5), "\n")

png("fig_manhattan.png", width = 1400, height = 600, res = 130)
manhattan(g, chr = "CHR", bp = "BP", snp = "SNP", p = "P",
          col = c("#1F3A5F", "#2A9D8F"),
          suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-8),
          main = "GWAS Manhattan -- height (with covariates)")
dev.off()
```

---

## Part D — Polygenic score with PRSice-2 (Lab 6)

```bash
PRSICE_BIN="$HOME/PRSice/PRSice_linux"
PRSICE_R="$HOME/PRSice/PRSice.R"

Rscript "$PRSICE_R" \
    --prsice "$PRSICE_BIN" \
    --base Height_GWAS_sumstats.txt \
    --target eur_qc \
    --binary-target F \
    --cov covariates_with_pcs.txt \
    --cov-col age,sex,PC1-PC10 \
    --clump-kb 250 --clump-r2 0.1 \
    --bar-levels 0.005,0.05,0.5,1 \
    --fastscore \
    --out prsice_height
```

PRSice picks up the phenotype from `eur_qc.fam` automatically. It writes:

- `prsice_height.summary` — winning p-threshold and R² across thresholds
- `prsice_height.best` — per-individual PGS at the best threshold
- `prsice_height.prsice` — R² at every tested threshold

Incremental R² of the PGS over a covariates-only model:

```r
prs <- read.table("prsice_height.best",         header = TRUE)  # FID IID In_Regression PRS
cov <- read.table("covariates_with_pcs.txt",    header = TRUE)
fam <- read.table("eur_qc.fam", header = FALSE,
                  col.names = c("FID","IID","PAT","MAT","SEX","height"))

d <- Reduce(function(a, b) merge(a, b, by = c("FID","IID")),
            list(prs, cov, fam[, c("FID","IID","height")]))
d$PRS_z <- as.numeric(scale(d$PRS))

m0 <- lm(height ~ age + sex + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
         data = d)
m1 <- lm(height ~ PRS_z + age + sex + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
         data = d)
inc_r2 <- summary(m1)$r.squared - summary(m0)$r.squared

cat(sprintf("Incremental R^2 of PGS-height on height = %.4f\n", inc_r2))
```

---

## Part E — G×E interaction (Labs 7–8)

The required regression is:

$$Y_i = \beta_0 + \beta_1\,\mathrm{PGS}_i + \beta_2\,\mathrm{Gender}_i + \beta_3\,(\mathrm{PGS}_i \cdot \mathrm{Gender}_i) + \mathrm{age}_i + \mathrm{PC1}_i + \dots + \mathrm{PC10}_i + \varepsilon_i$$

**Always standardise the PGS** before interacting it with another variable, otherwise the interaction coefficient is not interpretable.

```r
fit <- lm(height ~ PRS_z * factor(sex) + age +
                   PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
          data = d)
summary(fit)
# Report beta_3, its SE and p-value: the term `PRS_z:factor(sex)2`.
```

### Interaction plot

Predicted height vs PGS by gender, holding age at the sample mean and PCs at 0:

```r
grid <- expand.grid(
  PRS_z = seq(-3, 3, by = 0.1),
  sex   = c(1, 2),
  age   = mean(d$age),
  PC1 = 0, PC2 = 0, PC3 = 0, PC4 = 0, PC5 = 0,
  PC6 = 0, PC7 = 0, PC8 = 0, PC9 = 0, PC10 = 0
)
pred <- predict(fit, newdata = grid, interval = "confidence")
grid <- cbind(grid, pred)

grid$sex_label <- factor(grid$sex,
                         levels = c(1, 2),
                         labels = c("Male", "Female"))

p_gxe <- ggplot(grid, aes(PRS_z, fit, colour = sex_label, fill = sex_label)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("Male" = "#1F3A5F", "Female" = "#C0392B")) +
  scale_fill_manual  (values = c("Male" = "#1F3A5F", "Female" = "#C0392B")) +
  labs(x = "PGS (standardised)", y = "Predicted height (cm)",
       colour = NULL, fill = NULL,
       title = "G x E: PGS x gender on height",
       subtitle = "Predicted height at mean age, PCs set to 0") +
  theme_minimal(base_size = 11)

ggsave("fig_gxe.png", p_gxe, width = 6.5, height = 4.5, dpi = 150)
```

### Reading the sign of β₃

| Sign of β₃ | Interpretation | Direction |
|------------|----------------|-----------|
| Positive | PGS effect is **larger** in the reference-coded-2 group | reinforcing |
| Negative | PGS effect is **smaller** in the reference-coded-2 group | compensating |
| ~0 | PGS effect does not depend on the moderator | null |
| Slope sign flips across groups | **cross-over** interaction | |

Make sure you know which level of `sex` is coded 1 vs 2 before reading the sign.
