# =============================================================
# Lab 6 (supplement) — LDpred2 with bigsnpr
# Sociogenomics 2025/2026 · University of Bologna · Nicola Barban
# =============================================================
#
# This script computes a PGS using LDpred2-auto (Bayesian) and compares
# the result with the C+T approach from the main Lab 6.
#
# Data: same Trait2 summary statistics (chromosome 20) and 1000 Genomes
# European target sample used in the main lab.
#
# Requirements: R >= 4.1, packages bigsnpr, bigstatsr, ggplot2
#
# Run in RStudio, Posit Cloud, or Google Colab (R kernel).
# Estimated run time: 5-10 minutes.
# -------------------------------------------------------------


# ---- 0. Install packages and download data ------------------

if (!requireNamespace("bigsnpr", quietly = TRUE)) {
  install.packages("bigsnpr", repos = "https://cloud.r-project.org")
}
library(bigsnpr)
library(bigstatsr)
library(ggplot2)

# Create a working directory
lab_dir <- file.path(tempdir(), "lab6_ldpred2")
dir.create(lab_dir, showWarnings = FALSE)
setwd(lab_dir)
cat("Working directory:", lab_dir, "\n")

# Download the all-in-one results zip (contains phenotype, PCs, etc.)
zip_url <- "https://github.com/nicolabarban/sociogenomics_2025_2026/raw/gh_pages/labs/week6/lab6_results.zip"
download.file(zip_url, destfile = "lab6_results.zip", mode = "wb")
unzip("lab6_results.zip", overwrite = TRUE)

# Download the PLINK binary files (bed/bim/fam) for the QC'd European sample
data_url <- "https://github.com/nicolabarban/sociogenomics_2025_2026/raw/gh_pages/data/"
for (ext in c("bed", "bim", "fam")) {
  f <- paste0("1kg_hm3_QC_CEU.", ext)
  download.file(paste0(data_url, f), destfile = f, mode = "wb", quiet = TRUE)
}

# Download the summary statistics
download.file(paste0(data_url, "Trait2.ma"), destfile = "Trait2.ma", quiet = TRUE)

list.files()


# ---- 1. Read the target genotype data -----------------------

# snp_readBed converts PLINK files to the bigsnpr backing-file format
# This only needs to be done once; afterwards use snp_attach()
snp_readBed("1kg_hm3_QC_CEU.bed")

obj.bigSNP <- snp_attach("1kg_hm3_QC_CEU.rds")

G   <- obj.bigSNP$genotypes
map <- obj.bigSNP$map
fam <- obj.bigSNP$fam

cat("Genotype matrix:", nrow(G), "individuals x", ncol(G), "SNPs\n")
head(map)


# ---- 2. Load and prepare the summary statistics -------------

sumstats <- read.table("Trait2.ma", header = TRUE)
head(sumstats)

# LDpred2 requires columns: chr, pos, a0, a1, beta, beta_se, n_eff
# Trait2.ma has:             CHR, SNP, A1, A2, AF1, BETA, SE, P, N

# Rename to the bigsnpr convention
sumstats_prep <- data.frame(
  chr     = sumstats$CHR,
  rsid    = sumstats$SNP,
  pos     = NA,  # will be filled by snp_match
  a1      = sumstats$A1,   # effect allele
  a0      = sumstats$A2,   # other allele
  beta    = sumstats$BETA,
  beta_se = sumstats$SE,
  n_eff   = sumstats$N
)

cat("Summary statistics:", nrow(sumstats_prep), "SNPs\n")


# ---- 3. Match summary statistics with target genotypes ------

# Prepare the target map in the bigsnpr format
map_bigsnpr <- setNames(map[, c(1, 2, 4, 6, 5)],
                        c("chr", "rsid", "pos", "a1", "a0"))

# snp_match aligns alleles and flips signs as needed
df_beta <- snp_match(sumstats_prep, map_bigsnpr, join_by_pos = FALSE)

cat("Matched SNPs:", nrow(df_beta), "\n")

# Restrict to chromosome 20 (where our summary stats come from)
df_beta <- df_beta[df_beta$chr == 20, ]
cat("Chr 20 SNPs after matching:", nrow(df_beta), "\n")


# ---- 4. Compute the LD correlation matrix -------------------

# Get the column indices in the genotype matrix for matched SNPs
ind.col <- df_beta[["_NUM_ID_"]]

# Compute LD using a 500-kb window (physical positions)
# (Ideally we would use genetic positions in cM, but these are
# not available in our .bim file, so we use a bp-based window)
cat("Computing LD matrix (this may take 1-2 minutes)...\n")

tmp <- tempfile(tmpdir = lab_dir)

corr0 <- snp_cor(
  G,
  ind.col  = ind.col,
  size     = 500 / 1000,    # 500 kb window (in units of 1000 bp)
  ncores   = 1
)

# Convert to Sparse Filled-in Banded Matrix (SFBM) for LDpred2
corr <- as_SFBM(corr0, tmp, compact = TRUE)

cat("LD matrix:", nrow(corr0), "x", ncol(corr0), "\n")


# ---- 5. Estimate h² with LD score regression ----------------

ld <- Matrix::colSums(corr0^2)

ldsc <- with(df_beta,
  snp_ldsc(
    ld,
    length(ld),
    chi2        = (beta / beta_se)^2,
    sample_size = n_eff,
    blocks      = NULL
  )
)

h2_est <- ldsc[["h2"]]
cat("LDSC h² estimate:", round(h2_est, 4), "\n")


# ---- 6. Run LDpred2-auto ------------------------------------

# LDpred2-auto jointly estimates h², polygenicity (p), and
# effect sizes — no tuning parameters needed.

cat("Running LDpred2-auto (this may take 2-5 minutes)...\n")

set.seed(42)

multi_auto <- snp_ldpred2_auto(
  corr,
  df_beta,
  h2_init         = pmin(pmax(h2_est, 0.01), 0.5),  # clamp to [0.01, 0.5]
  vec_p_init      = seq_log(1e-4, 0.5, length.out = 30),
  ncores          = 1,
  burn_in         = 200,
  num_iter        = 200,
  allow_jump_sign = FALSE,
  shrink_corr     = 0.5   # use 0.5 because our LD ref is small (N=374)
)

cat("LDpred2-auto completed:", length(multi_auto), "chains\n")


# ---- 7. QC: filter diverged chains --------------------------

range_corr <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep       <- which(range_corr > (0.95 * quantile(range_corr, 0.95, na.rm = TRUE)))
cat("Chains kept after QC:", length(keep), "out of", length(multi_auto), "\n")

# Posterior h² and p estimates from kept chains
if (length(keep) > 0) {
  all_h2 <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est, 200))
  all_p  <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est,  200))
  cat("Posterior h²: ", round(median(all_h2), 4),
      " (95% CI:", round(quantile(all_h2, 0.025), 4), "-",
      round(quantile(all_h2, 0.975), 4), ")\n")
  cat("Posterior p:  ", round(median(all_p), 4),
      " (95% CI:", round(quantile(all_p, 0.025), 4), "-",
      round(quantile(all_p, 0.975), 4), ")\n")
}


# ---- 8. Compute the PGS from posterior betas ----------------

# Average posterior effect sizes across kept chains
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))

# Compute PGS for all individuals
pgs_ldpred2 <- big_prodVec(G, beta_auto, ind.col = df_beta[["_NUM_ID_"]])

cat("PGS computed for", length(pgs_ldpred2), "individuals\n")
cat("Mean:", round(mean(pgs_ldpred2), 6), "  SD:", round(sd(pgs_ldpred2), 6), "\n")


# ---- 9. Evaluate: incremental R² ---------------------------

# Load phenotype and PCs
pheno <- read.table("1kg.Trait2.phen", header = FALSE,
                    col.names = c("FID", "IID", "Trait2"))
pca   <- read.table("1kg_pca.eigenvec", header = FALSE,
                    col.names = c("FID", "IID", paste0("PC", 1:10)))

# Build a data frame with individual IDs from the fam file
d <- data.frame(FID = fam$family.ID, IID = fam$sample.ID,
                PGS_ldpred2 = scale(pgs_ldpred2))
d <- merge(d, pheno, by = c("FID", "IID"))
d <- merge(d, pca,   by = c("FID", "IID"))
cat("Merged N:", nrow(d), "\n")

# Baseline model (PCs only)
mod0 <- lm(Trait2 ~ PC1 + PC2 + PC3 + PC4 + PC5 +
                     PC6 + PC7 + PC8 + PC9 + PC10, data = d)

# Full model (PCs + LDpred2 PGS)
mod1 <- lm(Trait2 ~ PGS_ldpred2 + PC1 + PC2 + PC3 + PC4 + PC5 +
                    PC6 + PC7 + PC8 + PC9 + PC10, data = d)

delta_r2_ldpred2 <- summary(mod1)$r.squared - summary(mod0)$r.squared
cat("\n=== LDpred2 results ===\n")
cat("R² baseline (PCs):      ", round(summary(mod0)$r.squared, 4), "\n")
cat("R² full (PCs + LDpred2):", round(summary(mod1)$r.squared, 4), "\n")
cat("Incremental R²:         ", round(delta_r2_ldpred2, 4), "\n")


# ---- 10. Compare with C+T from Lab 6 -----------------------

# Load the PRSice PGS from the results zip
prsice <- read.table("Trait2_PRSice.all_score", header = TRUE)
d2 <- merge(d, prsice[, c("FID", "IID", "Pt_0.5")], by = c("FID", "IID"))
d2$PGS_ct <- scale(d2$Pt_0.5)

mod_ct <- lm(Trait2 ~ PGS_ct + PC1 + PC2 + PC3 + PC4 + PC5 +
                       PC6 + PC7 + PC8 + PC9 + PC10, data = d2)
delta_r2_ct <- summary(mod_ct)$r.squared - summary(mod0)$r.squared

cat("\n=== Comparison ===\n")
cat("C+T (PRSice, best threshold):  incremental R² =", round(delta_r2_ct, 4), "\n")
cat("LDpred2-auto:                  incremental R² =", round(delta_r2_ldpred2, 4), "\n")
cat("Ratio LDpred2 / C+T:         ", round(delta_r2_ldpred2 / delta_r2_ct, 2), "\n")


# ---- 11. Visualise both PGS --------------------------------

# Scatter plot: C+T vs LDpred2
d2$PGS_ldpred2 <- d[match(d2$IID, d$IID), "PGS_ldpred2"]

cat("\nCorrelation C+T vs LDpred2:", round(cor(d2$PGS_ct, d2$PGS_ldpred2), 3), "\n")

ggplot(d2, aes(x = PGS_ct, y = PGS_ldpred2)) +
  geom_point(alpha = 0.5, colour = "steelblue") +
  geom_smooth(method = "lm", colour = "coral", se = FALSE) +
  labs(x = "PGS (C+T, PRSice best threshold)",
       y = "PGS (LDpred2-auto)",
       title = "Comparison of C+T and LDpred2 polygenic scores") +
  theme_minimal()

# Histogram comparison
par(mfrow = c(1, 2))
hist(d2$PGS_ct, breaks = 25, col = "steelblue", border = "white",
     main = "C+T (PRSice)", xlab = "Standardised PGS")
hist(d2$PGS_ldpred2, breaks = 25, col = "coral", border = "white",
     main = "LDpred2-auto", xlab = "Standardised PGS")
par(mfrow = c(1, 1))


# ---- Summary ------------------------------------------------

cat("\n")
cat("===============================================\n")
cat("  Method             | Delta R²  | N SNPs\n")
cat("-----------------------------------------------\n")
cat(sprintf("  C+T (PRSice)       | %.4f    | ~2,089\n", delta_r2_ct))
cat(sprintf("  LDpred2-auto       | %.4f    | %d\n", delta_r2_ldpred2, nrow(df_beta)))
cat("===============================================\n")
cat("\n")
cat("\nWHY DOES LDpred2 UNDERPERFORM HERE?\n")
cat("-----------------------------------\n")
cat("LDpred2 needs a large LD reference (N > 2,000) to estimate the\n")
cat("correlation matrix accurately. Our sample has only N = 374, which\n")
cat("makes the LD matrix very noisy. As a result:\n")
cat("  - The LDSC h² estimate is inflated (> 1, which is impossible)\n")
cat("  - The posterior effect sizes are poorly shrunk\n")
cat("  - The PGS prediction is worse than simple C+T\n")
cat("\n")
cat("IN PRACTICE (with real data):\n")
cat("  - Use a large LD reference panel (e.g., UK Biobank Europeans, N > 300K)\n")
cat("  - Use genome-wide SNPs (not just one chromosome)\n")
cat("  - LDpred2 then typically outperforms C+T by 10-30%\n")
cat("\n")
cat("This exercise illustrates both the LDpred2 WORKFLOW and the\n")
cat("importance of a proper LD reference for Bayesian PGS methods.\n")
