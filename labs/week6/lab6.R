# =============================================================
# Lab 6 — Polygenic Scores (R analysis only)
# Sociogenomics 2025/2026 · University of Bologna · Nicola Barban
# =============================================================
#
# This script runs the R analysis component of Lab 6 on pre-computed
# PGS outputs. It downloads:
#   - lab6_results.zip   : PRSice and PLINK --score outputs
#   - 1kg.Trait2.phen    : phenotype
#   - 1kg_pca.eigenvec   : principal components
#   - 1kg-sample-2504-phased.txt : population info
#
# and reproduces the analysis from Part VI and Part VII of the lab:
#   - incremental R² controlling for ancestry PCs
#   - bootstrap 95% CI
#   - incremental R² across p-value thresholds
#   - cross-ancestry portability
#
# Run it in RStudio: Source, or line-by-line with Ctrl/Cmd+Enter.
# -------------------------------------------------------------

# ---- 0. Working directory and downloads ---------------------

# Create a clean working directory for the lab outputs
lab_dir <- file.path(tempdir(), "lab6")
dir.create(lab_dir, showWarnings = FALSE)
setwd(lab_dir)

# Download the all-in-one results zip (PGS outputs + phenotype + PCs + pop info)
zip_url <- "https://github.com/nicolabarban/sociogenomics_2025_2026/raw/gh_pages/labs/week6/lab6_results.zip"
download.file(zip_url, destfile = "lab6_results.zip", mode = "wb")
unzip("lab6_results.zip", overwrite = TRUE)

list.files()


# ---- Part III — Monogenic FTO score -------------------------

fto <- read.table("FTOscore.profile", header = TRUE)
head(fto)
table(fto$SCORE)


# ---- Part VI — Analyse the PGS in R -------------------------

# Load PGS at each p-value threshold (from PRSice)
all_scores <- read.table("Trait2_PRSice.all_score", header = TRUE)
head(all_scores)

# Columns: R replaces '-' with '.', so Pt_5e-08 becomes Pt_5e.08
colnames(all_scores)

# Load phenotype and principal components
pheno <- read.table("1kg.Trait2.phen", header = FALSE,
                    col.names = c("FID", "IID", "Trait2"))

pca_cols <- c("FID", "IID", paste0("PC", 1:10))
pca <- read.table("1kg_pca.eigenvec", header = FALSE, col.names = pca_cols)

# Merge
d <- merge(all_scores, pheno, by = c("FID", "IID"))
d <- merge(d, pca, by = c("FID", "IID"))
nrow(d)

# Standardise the PGS
d$PGS_best <- scale(d$Pt_0.5)
d$PGS_gw   <- scale(d$Pt_5e.08)
d$PGS_all  <- scale(d$Pt_1)

# PGS distribution (should be approximately normal)
hist(d$PGS_best, breaks = 30, col = "steelblue", border = "white",
     main = "Distribution of the standardised PGS",
     xlab = "PGS (z-score)")


# Incremental R²: baseline (PCs) vs. full (PCs + PGS)
mod0 <- lm(Trait2 ~ PC1 + PC2 + PC3 + PC4 + PC5 +
                     PC6 + PC7 + PC8 + PC9 + PC10, data = d)
mod1 <- lm(Trait2 ~ PGS_best + PC1 + PC2 + PC3 + PC4 + PC5 +
                    PC6 + PC7 + PC8 + PC9 + PC10, data = d)

delta_r2 <- summary(mod1)$r.squared - summary(mod0)$r.squared
cat("R2 baseline (PCs):  ", round(summary(mod0)$r.squared, 4), "\n")
cat("R2 full (PCs+PGS):  ", round(summary(mod1)$r.squared, 4), "\n")
cat("Incremental R2:     ", round(delta_r2, 4), "\n")
summary(mod1)$coefficients["PGS_best", ]


# Compare incremental R² across all thresholds
thresholds <- c("Pt_5e.08", "Pt_5e.06", "Pt_0.0005",
                "Pt_0.05",  "Pt_0.5",   "Pt_1")
labels     <- c("5e-8", "5e-6", "5e-4", "0.05", "0.5", "1")

r2_base <- summary(mod0)$r.squared
results <- data.frame(threshold = labels, delta_r2 = NA)

for (i in seq_along(thresholds)) {
  d$tmp <- scale(d[[thresholds[i]]])
  mod   <- lm(Trait2 ~ tmp + PC1 + PC2 + PC3 + PC4 + PC5 +
                        PC6 + PC7 + PC8 + PC9 + PC10, data = d)
  results$delta_r2[i] <- summary(mod)$r.squared - r2_base
}

print(results)

# Barplot
bp <- barplot(results$delta_r2, names.arg = results$threshold,
              col = "steelblue", border = NA,
              ylim = c(0, max(results$delta_r2) * 1.2),
              main = "Incremental R² by p-value threshold",
              xlab = "p-value threshold", ylab = "Incremental R²")
text(bp, results$delta_r2 + 0.005, round(results$delta_r2, 3), cex = 0.8)


# Bootstrap 95% CI for incremental R²
library(boot)
set.seed(12345)

rsq_fn <- function(data, indices) {
  ds <- data[indices, ]
  m0 <- lm(Trait2 ~ PC1 + PC2 + PC3 + PC4 + PC5 +
                    PC6 + PC7 + PC8 + PC9 + PC10, data = ds)
  m1 <- lm(Trait2 ~ PGS_best + PC1 + PC2 + PC3 + PC4 + PC5 +
                    PC6 + PC7 + PC8 + PC9 + PC10, data = ds)
  summary(m1)$r.squared - summary(m0)$r.squared
}

results_boot <- boot(data = d, statistic = rsq_fn, R = 1000)
boot.ci(results_boot, type = "norm")


# ---- Part VII — Cross-ancestry portability ------------------

pgs_all <- read.table("Trait2_pgs_all_pops.profile", header = TRUE)
pop     <- read.table("1kg-sample-2504-phased.txt", header = TRUE)

da <- merge(pgs_all[, c("FID", "IID", "SCORE")], pheno,
            by = c("FID", "IID"))
da <- merge(da, pop[, c("sample", "super_pop")],
            by.x = "IID", by.y = "sample")

table(da$super_pop)

# R² by super-population
for (p in c("EUR", "EAS", "SAS", "AFR", "AMR")) {
  sub <- da[da$super_pop == p, ]
  if (nrow(sub) == 0) {
    cat(sprintf("%s: no samples\n", p))
    next
  }
  r2 <- cor(sub$SCORE, sub$Trait2)^2
  cat(sprintf("%s: R² = %.4f  (N = %d)\n", p, r2, nrow(sub)))
}

# Boxplot of PGS distribution by super-population
boxplot(SCORE ~ super_pop, data = da,
        col = "lightblue",
        main = "PGS distribution by super-population",
        xlab = "Super-population", ylab = "Raw PGS (PLINK --score)")


# -------------------------------------------------------------
# End of Lab 6 R analysis.
# -------------------------------------------------------------
