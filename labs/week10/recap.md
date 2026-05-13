---
title: "Lab Recap & Final Report Pipeline"
---

# Lab Recap and Final Report Pipeline

This page is a single-stop reference for the **Final Report**. It bundles all the commands you have learnt across labs 1–9, organised in the order you will need them to complete the report. Read it from top to bottom: it is also the script of the in-class recap lecture.

> **Deadline:** Monday 8 June 2026 — 10 pages max, code in appendix.

---

## 1. The big picture

| Lab | Topic | Tools | Report part |
|-----|-------|-------|-------------|
| 1 | Unix shell, files, `awk` | Cloud Shell, bash | all |
| 2 | PLINK basics | PLINK 1.9, file formats | all |
| 3 | QC, relatedness, LD | PLINK, KING | **B** |
| 4 | Genetic tests + PCA + ancestry | R + PLINK `--pca` | **A, C** |
| 5 | Phenotype simulation | R, bash | — |
| 6 | PGS with PRSice-2 | PRSice-2, R | **D** |
| 7 | G×E interaction | R `lm()`, ggplot | **E** |
| 8 | Replicating Walter (2016) | R, HRS subset | **E** |
| 9 | PGI tutorial (pre-computed PGIs) | R, scoring | **D** |

One pipeline: **raw genotypes** → **cleaned data** → **ancestry-restricted subset** → **GWAS** → **PGS** → **G×E**.

---

## 2. What the report contains

| Part | Output | Labs | Pts |
|------|--------|------|-----|
| A. Ancestry & population structure | PCA plots, EUR classification, `my_eur_ids.txt` | 1, 4 | 25 |
| B. Quality control | Filtered EUR subset, table of pass/fail counts | 1–3 | 20 |
| C. GWAS | Top-10 SNPs, Manhattan, QQ plot (with/without covs) | 4 | 20 |
| D. Polygenic score | PRSice bar-plot, best threshold, incremental R² | 6, 9 | 20 |
| E. G×E interaction | PGS×Gender regression, interaction plot | 7, 8 | 15 |

Group of 2–3 students per phenotype:

| Group | Phenotype | File |
|-------|-----------|------|
| 1, 2 | Height | `pheno_height.txt` |
| 3, 4 | BMI | `pheno_BMI.txt` |
| 5, 6 | Systolic Blood Pressure | `pheno_SBP.txt` |
| 7, 8 | Diastolic Blood Pressure | `pheno_DBP.txt` |
| 9, 10 | Heart Rate | `pheno_HR.txt` |

---

## 3. Setup

```bash
# Move into the report folder distributed by the instructor
cd report_2026
ls -la 1kg_hm3_anon.{bed,bim,fam}
wc -l 1kg_hm3_anon.fam          # number of individuals
wc -l 1kg_hm3_anon.bim          # number of variants
head -3 1kg_hm3_anon.fam        # FID IID PAT MAT SEX PHENO
head -3 1kg_hm3_anon.bim        # CHR SNP CM BP A1 A2
head -3 covariates.txt          # FID IID age gender
head -3 pheno_DBP.txt           # FID IID DBP  (replace with your group's pheno)
head -3 reference_panel.txt     # FID IID superpop  (40 labels)

# How many SNPs per chromosome?
awk '{print $1}' 1kg_hm3_anon.bim | sort | uniq -c | sort -n
```

You need:

- **PLINK 1.9** (`plink`) and **PLINK 2.0** (`plink2`) — the latter is much faster.
- **R** with packages `data.table`, `ggplot2`, `dplyr`, `qqman`, `class`.
- **PRSice-2** for Part D (the binary `PRSice_linux` plus `PRSice.R`).

---

## 4. Part A — Ancestry and population structure (Lab 4)

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
library(data.table)
library(ggplot2)
library(class)

ev  <- fread("work_pca.eigenvec",
             col.names = c("FID", "IID", paste0("PC", 1:10)))
ref <- fread("reference_panel.txt")          # FID IID superpop (40 labels)

# Merge: tag the 40 reference individuals by their super-population
ev <- merge(ev, ref[, .(IID, superpop)], by = "IID", all.x = TRUE)
ev[, is_ref := !is.na(superpop)]

# Quick look at the cloud
ggplot(ev, aes(PC1, PC2)) +
  geom_point(aes(colour = is.na(superpop)), alpha = 0.4) +
  geom_point(data = ev[!is.na(superpop)],
             aes(colour = superpop), size = 3) +
  scale_colour_manual(values = c("EUR" = "#1F3A5F", "AFR" = "#C0392B",
                                 "EAS" = "#2A9D8F", "AMR" = "#E9A23B",
                                 "TRUE" = "grey80", "FALSE" = "grey80")) +
  labs(title = "PCA: all individuals + 40 labelled references",
       x = "PC1", y = "PC2") +
  theme_minimal()

# k-NN classification using the 40 labelled samples on PC1-PC5
ref_mat <- as.matrix(ev[is_ref == TRUE,  paste0("PC", 1:5), with = FALSE])
unk_mat <- as.matrix(ev[is_ref == FALSE, paste0("PC", 1:5), with = FALSE])
ref_lab <- ev[is_ref == TRUE, superpop]

set.seed(42)
ev[is_ref == FALSE, pred_pop := as.character(knn(ref_mat, unk_mat, ref_lab, k = 3))]
ev[is_ref == TRUE,  pred_pop := superpop]

table(ev$pred_pop)

# Save EUR IDs in PLINK keep format
my_eur <- ev[pred_pop == "EUR", .(FID, IID)]
fwrite(my_eur, "my_eur_ids.txt", sep = "\t", col.names = FALSE)
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

## 5. Part B — Quality control (Lab 3)

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

### Build a filter-by-filter table (for the report)

```bash
for filt in "--maf 0.05" "--mind 0.05" "--geno 0.05" "--hwe 1e-5"; do
  echo ">>> $filt"
  plink2 --bfile 1kg_hm3_anon --keep my_eur_ids.txt $filt \
         --make-just-bim --out _tmp 2>/dev/null
  wc -l _tmp.bim
done
rm -f _tmp.*
```

> **Common mistake:** running QC on the full multi-ancestry cohort and *then* subsetting to EUR. The HWE filter is the giveaway: mixed-ancestry samples look out-of-equilibrium even when each ancestry is fine. Always **first ancestry, then QC**.

---

## 6. Part C — GWAS (Lab 4)

### Step 1. GWAS without covariates

```bash
plink2 --bfile eur_qc \
       --pheno pheno_DBP.txt --pheno-name DBP \
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
       --pheno pheno_DBP.txt --pheno-name DBP \
       --covar covariates_with_pcs.txt \
       --covar-name age,sex,PC1-PC10 \
       --glm hide-covar \
       --out gwas_cov
```

### Step 4. Top-10 SNPs in each model

```bash
echo "TOP 10 -- NO COVARIATES"
( head -1 gwas_nocov.DBP.glm.linear ; \
  tail -n +2 gwas_nocov.DBP.glm.linear | sort -gk12 | head -10 ) | column -t

echo "TOP 10 -- WITH COVARIATES"
( head -1 gwas_cov.DBP.glm.linear ; \
  tail -n +2 gwas_cov.DBP.glm.linear | sort -gk12 | head -10 ) | column -t
```

### Step 5. Manhattan and QQ plots

Two recipes — quick (`qqman`) and from scratch (`ggplot2`). Use whichever you prefer for the report.

```r
library(data.table)
library(ggplot2)

g <- fread("gwas_cov.DBP.glm.linear")
setnames(g, old = c("#CHROM", "POS", "ID"),
            new = c("CHR", "BP", "SNP"),
            skip_absent = TRUE)
g <- g[!is.na(P)]

cat("SNPs tested:", nrow(g), "\n")
cat("Genome-wide significant (p<5e-8):", sum(g$P < 5e-8), "\n")
```

#### Option A — Quick with `qqman`

```r
library(qqman)

png("fig_manhattan_qqman.png", width = 1400, height = 600, res = 130)
manhattan(g, chr = "CHR", bp = "BP", snp = "SNP", p = "P",
          col = c("#1F3A5F", "#2A9D8F"),
          suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-8),
          main = "GWAS Manhattan -- DBP (with covariates)")
dev.off()

png("fig_qq_qqman.png", width = 600, height = 600, res = 130)
qq(g$P, main = "QQ plot -- DBP (with covariates)")
dev.off()
```

#### Option B — From scratch with `ggplot2` (more control)

A clean Manhattan needs four ingredients:

1. `-log10(p)` for each SNP
2. A **cumulative position** so chromosomes do not pile up at x = 0
3. **Alternating colours** per chromosome for readability
4. Optional: labels on the top hits

```r
# (1) -log10(p)
g[, logp := -log10(P)]

# (2) cumulative BP across chromosomes
chr_len <- g[, .(chr_len = max(BP)), by = CHR][order(CHR)]
chr_len[, chr_start := cumsum(as.numeric(chr_len)) - chr_len]
g <- merge(g, chr_len[, .(CHR, chr_start)], by = "CHR")
g[, bp_cum := BP + chr_start]
axis_df <- g[, .(centre = mean(bp_cum)), by = CHR][order(CHR)]

# (3) thin the cloud for speed: keep all p<0.05, subsample 10% of the rest
set.seed(1)
plot_dat <- rbind(g[P < 0.05], g[P >= 0.05][sample(.N, .N * 0.10)])

# (4) label the top 5
top <- g[order(P)][1:5]

p_man <- ggplot(plot_dat,
                aes(x = bp_cum, y = logp, colour = factor(CHR %% 2))) +
  geom_point(alpha = 0.7, size = 0.6) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", colour = "#C0392B") +
  geom_hline(yintercept = -log10(1e-5), linetype = "dotted", colour = "#E9A23B") +
  geom_text(data = top, aes(label = SNP), nudge_y = 0.4, size = 3,
            colour = "black", check_overlap = TRUE) +
  scale_colour_manual(values = c("0" = "#1F3A5F", "1" = "#2A9D8F"),
                      guide = "none") +
  scale_x_continuous(breaks = axis_df$centre, labels = axis_df$CHR) +
  labs(x = "Chromosome", y = expression(-log[10](p)),
       title = "GWAS Manhattan -- DBP (with covariates)",
       subtitle = "Dashed red = 5e-8 (GW), dotted gold = 1e-5 (suggestive)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("fig_manhattan_ggplot.png", p_man,
       width = 11, height = 4.5, dpi = 150)
```

QQ plot with a 95% confidence band and the genomic-inflation factor λ_GC printed on top:

```r
n <- nrow(g)
qq_dat <- data.table(
  obs = -log10(sort(g$P)),
  exp = -log10(ppoints(n))
)
qq_dat[, lo := -log10(qbeta(0.975, 1:n, n - 1:n + 1))]
qq_dat[, hi := -log10(qbeta(0.025, 1:n, n - 1:n + 1))]

lambda <- median(qchisq(1 - g$P, df = 1)) / qchisq(0.5, 1)
cat(sprintf("Genomic inflation lambda_GC = %.3f\n", lambda))

p_qq <- ggplot(qq_dat, aes(exp, obs)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") +
  geom_abline(slope = 1, intercept = 0,
              colour = "#C0392B", linetype = "dashed") +
  geom_point(colour = "#1F3A5F", size = 0.7, alpha = 0.7) +
  annotate("label", x = 0.5, y = max(qq_dat$obs) * 0.9,
           label = sprintf("lambda[GC] == %.3f", lambda), parse = TRUE) +
  labs(x = expression("Expected " * -log[10](p)),
       y = expression("Observed " * -log[10](p)),
       title = "QQ plot -- DBP (with covariates)") +
  theme_minimal(base_size = 11)

ggsave("fig_qq_ggplot.png", p_qq, width = 5.5, height = 5, dpi = 150)
```

**How to read the QQ plot:**

- Points hug the dashed red line for most of the range, then **lift away in the tail** → genuine polygenic signal.
- Points lift early (small expected `-log10(p)`) and never come back → **inflation** (stratification, cryptic relatedness, batch effects).
- λ_GC close to **1.0** means well-controlled inflation; values much above 1.05 are a warning sign.

---

## 7. Part D — Polygenic score with PRSice-2 (Lab 6)

```bash
PRSICE_BIN="$HOME/PRSice/PRSice_linux"
PRSICE_R="$HOME/PRSice/PRSice.R"

Rscript "$PRSICE_R" \
    --prsice "$PRSICE_BIN" \
    --base Height_GWAS_sumstats.txt \
    --target eur_qc \
    --pheno pheno_DBP.txt \
    --pheno-col DBP \
    --binary-target F \
    --cov covariates_with_pcs.txt \
    --cov-col age,sex,PC1-PC10 \
    --clump-kb 250 --clump-r2 0.1 \
    --bar-levels 0.005,0.05,0.5,1 \
    --fastscore \
    --out prsice_dbp
```

PRSice writes three files:

- `prsice_dbp.summary` — winning p-threshold and R² across thresholds
- `prsice_dbp.best` — per-individual PGS at the best threshold
- `prsice_dbp.prsice` — R² at every tested threshold

Post-process in R to compute the **incremental R²** of the PGS over a covariates-only model:

```r
library(data.table)

prs <- fread("prsice_dbp.best")           # FID IID In_Regression PRS
phe <- fread("pheno_DBP.txt")             # FID IID DBP
cov <- fread("covariates_with_pcs.txt")   # FID IID age sex PC1..PC10

d <- Reduce(function(a, b) merge(a, b, by = c("FID", "IID")),
            list(prs, phe, cov))
d[, PRS_z := scale(PRS)[, 1]]             # standardise to mean 0, SD 1

m0 <- lm(DBP ~ age + sex + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
         data = d)
m1 <- lm(DBP ~ PRS_z + age + sex + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
         data = d)
inc_r2 <- summary(m1)$r.squared - summary(m0)$r.squared

cat(sprintf("Incremental R^2 of PGS-height on DBP = %.4f\n", inc_r2))
```

> **Honest interpretation.** Unless your phenotype is itself height, expect a *small* incremental R². The PGS for height is used here to show you the full pipeline, **not** to predict your trait accurately.

---

## 8. Part E — G×E interaction (Labs 7–8)

The required regression is:

$$Y_i = \beta_0 + \beta_1\,\mathrm{PGS}_i + \beta_2\,\mathrm{Gender}_i + \beta_3\,(\mathrm{PGS}_i \cdot \mathrm{Gender}_i) + \mathrm{age}_i + \mathrm{PC1}_i + \dots + \mathrm{PC10}_i + \varepsilon_i$$

**Always standardise the PGS** before interacting it with another variable, otherwise the interaction coefficient is not interpretable.

```r
fit <- lm(DBP ~ PRS_z * factor(sex) + age +
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
  labs(x = "PGS (standardised)", y = "Predicted DBP",
       colour = NULL, fill = NULL,
       title = "G x E: PGS x gender on DBP",
       subtitle = "Predicted DBP at mean age, PCs set to 0") +
  theme_minimal(base_size = 11)

ggsave("fig_gxe.png", p_gxe, width = 6.5, height = 4.5, dpi = 150)
```

### Bonus — age terciles

```r
d[, age_tercile := cut(age,
                       breaks = quantile(age, c(0, 1/3, 2/3, 1)),
                       include.lowest = TRUE,
                       labels = c("Young", "Middle", "Old"))]

fit_age <- lm(DBP ~ PRS_z * age_tercile +
                    PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
              data = d)
summary(fit_age)
```

### Reading the sign of β₃

| Sign of β₃ | Interpretation | Direction |
|------------|----------------|-----------|
| Positive | PGS effect is **larger** in the reference-coded-2 group | reinforcing |
| Negative | PGS effect is **smaller** in the reference-coded-2 group | compensating |
| ~0 | PGS effect does not depend on the moderator | null |
| Slope sign flips across groups | **cross-over** interaction | |

Make sure you know which level of `sex` is coded 1 vs 2 before reading the sign.

---

## 9. Common pitfalls

1. **Running QC on the multi-ancestry cohort.** First classify ancestry, then QC the EUR subset.
2. **Forgetting PCs in the GWAS covariate set.** Inflated test statistics and false positives. Always include PC1–PC10.
3. **Not standardising the PGS** before the interaction. β₃ becomes uninterpretable.
4. **Reporting p-values without effect sizes.** A statistically significant β with no effect-size interpretation is not science.
5. **Submitting code that does not reproduce the figures** in the report.

---

## 10. Practical roadmap

| Day | Task |
|-----|------|
| 1 | Setup, download data, sanity-check `.fam` and phenotype files |
| 2 | Part A — PCA + ancestry classification → `my_eur_ids.txt` |
| 3 | Part B — QC on the EUR subset |
| 4 | Part C — GWAS, Manhattan + QQ |
| 5 | Part D — PRSice-2 and incremental R² |
| 6 | Part E — G×E regression and plots |
| 7 | Write-up, cross-check figures and code |

> **Group division:** one *analysis lead* on PLINK, one *visualisation lead* on R/ggplot, one *writing lead* on the narrative — but everyone reads everything before submission.

---

## 11. What we look for when grading

**Earns points**

- Clear narrative connecting all five parts
- Numbers reported at every step, not just final plots
- Honest discussion of *negative* or *noisy* results
- Correct interpretation of β₃ in Part E
- Reproducible, well-commented code

**Loses points**

- Plots without axis labels or sample sizes
- PGS regressions without standardisation
- Over-claiming causality from a cross-sectional R²
- Code that does not run / paths hard-coded to one machine
- Missing the page limit (10 pages text + figures)

---

## 12. Outputs you should attach to the report

- `fig_pca_overview.png` (Part A)
- `my_eur_ids.txt` (Part A)
- Per-filter QC table (Part B)
- `fig_manhattan_ggplot.png`, `fig_qq_ggplot.png` (Part C)
- PRSice bar-plot from `prsice_dbp_BARPLOT_*` and `prsice_dbp.summary` (Part D)
- `fig_gxe.png` (Part E)
- Code appendix (single `.R` and/or `.sh` file that reproduces everything above)

---

## 13. Submission & non-attending students

- **Deadline:** Monday 8 June 2026
- **Format:** one PDF for the report (max 10 pages) + a zip with code
- **Channel:** Virtuale / e-mail (follow course conventions)
- **Non-attending students:** substitute the lab project with an oral discussion of **two papers of your choice** from the [Week 8 G×E reading list](../../labs/week8/).

Good luck.
