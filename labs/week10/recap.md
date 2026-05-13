---
title: "Lab Recap & Final Report Pipeline"
---

# Lab Recap and Final Report Pipeline

## Part A — Ancestry and population structure (Lab 4)

### Step 1. LD-prune SNPs (for a clean PCA)

```bash
plink2 --bfile 1kg_hm3_anon \
       --maf 0.01 --geno 0.05 \
       --indep-pairwise 200 50 0.2 \
       --out work_prune
```

### Step 2. Run PCA on the pruned set

```bash
plink2 --bfile 1kg_hm3_anon \
       --extract work_prune.prune.in \
       --pca 10 \
       --out work_pca

head -3 work_pca.eigenvec       # FID IID PC1 PC2 ... PC10
cat work_pca.eigenval           # % variance per PC
```

### Step 3. Classify ancestry in R

```r
library(ggplot2)
library(dplyr)

ev  <- read.table("work_pca.eigenvec",
                  header = FALSE,
                  col.names = c("FID", "IID", paste0("PC", 1:10)))
ref <- read.table("reference_panel.txt",
                  header = TRUE)          # FID IID superpop (40 labels)

# Tag the 40 reference individuals by their super-population
ev <- merge(ev, ref[, c("IID", "superpop")], by = "IID", all.x = TRUE)
ev$is_ref <- !is.na(ev$superpop)

# Quick look at the cloud
ggplot(ev, aes(PC1, PC2)) +
  geom_point(aes(colour = is.na(superpop)), alpha = 0.4) +
  geom_point(data = subset(ev, !is.na(superpop)),
             aes(colour = superpop), size = 3) +
  scale_colour_manual(values = c("EUR" = "#1F3A5F", "AFR" = "#C0392B",
                                 "EAS" = "#2A9D8F", "AMR" = "#E9A23B",
                                 "TRUE" = "grey80", "FALSE" = "grey80")) +
  labs(title = "PCA: all individuals + 40 labelled references",
       x = "PC1", y = "PC2") +
  theme_minimal()

# k-means clustering on PC1..PC5 with 4 clusters (we expect 4 super-pops)
pc_cols <- paste0("PC", 1:5)

set.seed(42)
km <- kmeans(ev[, pc_cols], centers = 4, nstart = 25)
ev$cluster <- km$cluster

# Label each cluster with the most-frequent super-population among
# the reference individuals that fell into it.
cluster_to_pop <- ev %>%
  filter(is_ref) %>%
  count(cluster, superpop) %>%
  group_by(cluster) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  select(cluster, superpop) %>%
  rename(pred_pop = superpop)

ev <- merge(ev, cluster_to_pop, by = "cluster", all.x = TRUE)

table(ev$pred_pop)

# Save EUR IDs in PLINK keep format
my_eur <- subset(ev, pred_pop == "EUR", select = c("FID", "IID"))
write.table(my_eur, "my_eur_ids.txt",
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
cat("EUR retained:", nrow(my_eur), "\n")

# Classified scatter for the report
ggplot(ev, aes(PC1, PC2, colour = pred_pop)) +
  geom_point(alpha = 0.6) +
  scale_colour_manual(values = c("EUR" = "#1F3A5F", "AFR" = "#C0392B",
                                 "EAS" = "#2A9D8F", "AMR" = "#E9A23B")) +
  labs(title = "Predicted super-population", colour = NULL) +
  theme_minimal()
```

**Required deliverable:** a PCA plot with the 40 labelled references overlaid, the predicted super-population scatter, a count of how many individuals you classified as EUR, the file `my_eur_ids.txt`, and 3–5 sentences on *why* a single-ancestry GWAS matters.

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

```bash
plink2 --bfile 1kg_hm3_anon \
       --keep my_eur_ids.txt \
       --maf 0.05 \
       --mind 0.05 \
       --geno 0.05 \
       --hwe 1e-5 \
       --rel-cutoff 0.1 \
       --make-bed --out eur_qc

wc -l eur_qc.{bim,fam}
```

> **Common mistake:** running QC on the full multi-ancestry cohort and *then* subsetting to EUR. The HWE filter is the giveaway: mixed-ancestry samples look out-of-equilibrium even when each ancestry is fine. Always **first ancestry, then QC**.

---

## Simulate a phenotype

Before running the GWAS, we need a phenotype. For this tutorial we **simulate a random continuous trait** $Y_{sim}$ for every individual in the QC'd subset.

```r
library(dplyr)

fam <- read.table("eur_qc.fam", header = FALSE,
                  col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))
cov <- read.table("covariates.txt", header = TRUE)

# Merge to bring age and sex alongside FID/IID
fam <- merge(fam[, c("FID", "IID")], cov, by = c("FID", "IID"))

# Simulate Y_sim = 0.30 * age + 0.50 * sex + Gaussian noise
set.seed(2026)
fam$Ysim <- 0.30 * fam$age + 0.50 * fam$sex + rnorm(nrow(fam), mean = 0, sd = 5)

# Save in PLINK format
write.table(fam[, c("FID", "IID", "Ysim")], "pheno_sim.txt",
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

summary(fam$Ysim)
```

`pheno_sim.txt` now contains `FID  IID  Ysim` for every QC'd individual. We use it for the rest of the pipeline.

---

## Part C — GWAS (Lab 4)

### Step 1. GWAS without covariates

```bash
plink2 --bfile eur_qc \
       --pheno pheno_sim.txt --pheno-name Ysim \
       --glm \
       --out gwas_nocov
```

### Step 2. Merge PCs into the covariate file

```bash
# Pair eigenvec rows with covariates by IID
awk 'NR==FNR{a[$2]=$0; next} $2 in a{print $0,a[$2]}' \
    work_pca.eigenvec covariates.txt > _cov_with_pcs.txt

{
  echo -e "FID\tIID\tage\tsex\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10"
  awk '{print $1,$2,$3,$4,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' OFS='\t' \
      _cov_with_pcs.txt
} > covariates_with_pcs.txt
rm _cov_with_pcs.txt
```

### Step 3. GWAS with covariates

```bash
plink2 --bfile eur_qc \
       --pheno pheno_sim.txt --pheno-name Ysim \
       --covar covariates_with_pcs.txt \
       --covar-name age,sex,PC1-PC10 \
       --glm hide-covar \
       --out gwas_cov
```

### Step 4. Top-10 SNPs in each model

```bash
echo "TOP 10 -- NO COVARIATES"
( head -1 gwas_nocov.Ysim.glm.linear ; \
  tail -n +2 gwas_nocov.Ysim.glm.linear | sort -gk12 | head -10 ) | column -t

echo "TOP 10 -- WITH COVARIATES"
( head -1 gwas_cov.Ysim.glm.linear ; \
  tail -n +2 gwas_cov.Ysim.glm.linear | sort -gk12 | head -10 ) | column -t
```

### Step 5. Manhattan plot

```r
library(qqman)

g <- read.table("gwas_cov.Ysim.glm.linear",
                header = TRUE, comment.char = "",
                check.names = FALSE)
# Rename PLINK2 columns to the friendly names qqman expects
names(g)[names(g) == "#CHROM"] <- "CHR"
names(g)[names(g) == "POS"]    <- "BP"
names(g)[names(g) == "ID"]     <- "SNP"

g <- subset(g, !is.na(P))

cat("SNPs tested:", nrow(g), "\n")
cat("Genome-wide significant (p<5e-8):", sum(g$P < 5e-8), "\n")

png("fig_manhattan.png", width = 1400, height = 600, res = 130)
manhattan(g, chr = "CHR", bp = "BP", snp = "SNP", p = "P",
          col = c("#1F3A5F", "#2A9D8F"),
          suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-8),
          main = "GWAS Manhattan -- Ysim (with covariates)")
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
    --pheno pheno_sim.txt \
    --pheno-col Ysim \
    --binary-target F \
    --cov covariates_with_pcs.txt \
    --cov-col age,sex,PC1-PC10 \
    --clump-kb 250 --clump-r2 0.1 \
    --bar-levels 0.005,0.05,0.5,1 \
    --fastscore \
    --out prsice_sim
```

PRSice writes three files:

- `prsice_sim.summary` — winning p-threshold and R² across thresholds
- `prsice_sim.best` — per-individual PGS at the best threshold
- `prsice_sim.prsice` — R² at every tested threshold

Post-process in R to compute the **incremental R²** of the PGS over a covariates-only model:

```r
prs <- read.table("prsice_sim.best",           header = TRUE)  # FID IID In_Regression PRS
phe <- read.table("pheno_sim.txt",             header = TRUE)  # FID IID Ysim
cov <- read.table("covariates_with_pcs.txt",   header = TRUE)  # FID IID age sex PC1..PC10

d <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID")),
            list(prs, phe, cov))
d$PRS_z <- as.numeric(scale(d$PRS))               # standardise to mean 0, SD 1

m0 <- lm(Ysim ~ age + sex + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
         data = d)
m1 <- lm(Ysim ~ PRS_z + age + sex + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
         data = d)
inc_r2 <- summary(m1)$r.squared - summary(m0)$r.squared

cat(sprintf("Incremental R^2 of PGS-height on Ysim = %.4f\n", inc_r2))
```

> **Honest interpretation.** Because `Ysim` is simulated independently of the genome, expect a *small* incremental R². The PGS for height is used here to show you the full pipeline, **not** to predict the trait accurately.

---

## Part E — G×E interaction (Labs 7–8)

The required regression is:

$$Y_i = \beta_0 + \beta_1\,\mathrm{PGS}_i + \beta_2\,\mathrm{Gender}_i + \beta_3\,(\mathrm{PGS}_i \cdot \mathrm{Gender}_i) + \mathrm{age}_i + \mathrm{PC1}_i + \dots + \mathrm{PC10}_i + \varepsilon_i$$

**Always standardise the PGS** before interacting it with another variable, otherwise the interaction coefficient is not interpretable.

```r
fit <- lm(Ysim ~ PRS_z * factor(sex) + age +
                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
          data = d)
summary(fit)
# Report beta_3, its SE and p-value: the term `PRS_z:factor(sex)2`.
```

### Interaction plot

Predicted Y vs PGS by gender, holding age at the sample mean and PCs at 0:

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
  labs(x = "PGS (standardised)", y = "Predicted Ysim",
       colour = NULL, fill = NULL,
       title = "G x E: PGS x gender on Ysim",
       subtitle = "Predicted Ysim at mean age, PCs set to 0") +
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
