---
title: "Lab Recap & Final Report Pipeline"
---

# Lab Recap and Final Report Pipeline

We use:

- the `1kg_hm3` PLINK fileset (1.092 individuals from 1000G, 851k SNPs, rs-IDs),
- `height_cm.phen` — the simulated height phenotype (FID IID HEIGHT, no header),
- `1kg_samples.txt` — 1000G super-population labels.

## Part A — Ancestry and population structure (Lab 4)

### Step 1. LD-prune SNPs (for a clean PCA)

```bash
plink --bfile 1kg_hm3 \
      --maf 0.01 --geno 0.05 \
      --indep-pairwise 200 50 0.2 \
      --out work_prune
```

### Step 2. Run PCA on the pruned set

```bash
plink --bfile 1kg_hm3 \
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
# (1kg_hm3 contains AFR, AMR, EAS, EUR — no SAS)
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
plink --bfile 1kg_hm3 \
      --keep samples_EUR_like.txt \
      --maf 0.05 --mind 0.05 --geno 0.05 --hwe 1e-5 \
      --rel-cutoff 0.1 \
      --out eur_unrel

# Pass 2: produce the QC'd binary fileset
plink --bfile 1kg_hm3 \
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

`1kg_hm3` has sex in the `.fam` (column 5) but no age and no phenotype. We synthesise an age, read height from `height_cm.phen`, and merge everything together with the within-EUR PCs.

```r
fam <- read.table("eur_qc.fam", header = FALSE,
                  col.names = c("FID","IID","PAT","MAT","SEX","PHENO"))

# Use sex from the .fam (1 = male, 2 = female; drop unknowns)
fam <- subset(fam, SEX %in% c(1, 2))
fam$sex <- fam$SEX

# Synthesise a plausible adult age (uniform 25-75)
set.seed(2026)
fam$age <- runif(nrow(fam), 25, 75)

# Read height from the phenotype file
phe <- read.table("height_cm.phen", header = FALSE,
                  col.names = c("FID","IID","height"))

# Merge in the within-EUR PCs
pcs <- read.table("eur_pca.eigenvec",
                  header = FALSE,
                  col.names = c("FID", "IID", paste0("PC", 1:10)))

cov_pc <- Reduce(function(a, b) merge(a, b, by = c("FID","IID")),
                 list(fam[, c("FID","IID","age","sex")],
                      phe,
                      pcs[, c("FID","IID", paste0("PC", 1:10))]))

write.table(cov_pc[, c("FID","IID","age","sex", paste0("PC", 1:10))],
            "covariates_with_pcs.txt",
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
```

---

## Part C — GWAS (Lab 4)

The phenotype lives in `height_cm.phen`; we pass it to plink with `--pheno`.

### Step 1. GWAS without covariates

```bash
plink --bfile eur_qc \
      --chr 1-22 \
      --pheno height_cm.phen \
      --linear \
      --out gwas_nocov
```

### Step 2. GWAS with covariates

```bash
plink --bfile eur_qc \
      --chr 1-22 \
      --pheno height_cm.phen \
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

PRSice is shipped in `data/PRSice/`:

- `PRSice_mac` — the Mac binary (x86_64, runs on Apple Silicon via Rosetta)
- `PRSice.R` — the R wrapper

```bash
PRSICE_BIN="../data/PRSice/PRSice_mac"      # adjust to where you launch from
PRSICE_R="../data/PRSice/PRSice.R"

Rscript "$PRSICE_R" \
    --prsice "$PRSICE_BIN" \
    --base Height_GWAS_sumstats.txt \
    --target eur_qc \
    --pheno height_cm.phen --pheno-col height \
    --binary-target F \
    --cov covariates_with_pcs.txt \
    --cov-col age,sex,@PC[1-10] \
    --clump-kb 250 --clump-r2 0.1 \
    --bar-levels 0.005,0.05,0.5,1 \
    --fastscore \
    --out prsice_height
```

Notes on the syntax:

- `@PC[1-10]` is PRSice's shorthand for `PC1,PC2,...,PC10` — PLINK's `PC1-PC10` range is **not** understood by PRSice.
- The phenotype file `height_cm.phen` has no header, so we name the trait `height` here.

PRSice writes:

- `prsice_height.summary` — winning p-threshold and R² across thresholds
- `prsice_height.best` — per-individual PGS at the best threshold
- `prsice_height.prsice` — R² at every tested threshold

Incremental R² of the PGS over a covariates-only model:

```r
prs <- read.table("prsice_height.best",         header = TRUE)  # FID IID In_Regression PRS
cov <- read.table("covariates_with_pcs.txt",    header = TRUE)
phe <- read.table("height_cm.phen", header = FALSE,
                  col.names = c("FID","IID","height"))

d <- Reduce(function(a, b) merge(a, b, by = c("FID","IID")),
            list(prs, cov, phe))
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

Same idiom as Lab 8 — use `interact_plot()` from the **interactions** package:

```r
# install.packages(c("interactions", "jtools"))
library(interactions)

interact_plot(fit,
              pred = "PRS_z", modx = "sex",
              modx.values = c(1, 2),
              modx.labels = c("Male", "Female"),
              interval = TRUE, int.width = 0.95,
              x.label = "PGS (SD)",
              y.label = "Predicted height (cm)",
              legend.main = "Gender")
```

### Reading the sign of β₃

| Sign of β₃ | Interpretation | Direction |
|------------|----------------|-----------|
| Positive | PGS effect is **larger** in the reference-coded-2 group | reinforcing |
| Negative | PGS effect is **smaller** in the reference-coded-2 group | compensating |
| ~0 | PGS effect does not depend on the moderator | null |
| Slope sign flips across groups | **cross-over** interaction | |

Make sure you know which level of `sex` is coded 1 vs 2 before reading the sign.
