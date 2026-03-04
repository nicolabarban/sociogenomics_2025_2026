# =============================================================================
# Lab 4 — Genetic Tests and Principal Component Analysis
# Sociogenomics 2025/2026 · University of Bologna · Prof. Nicola Barban
# =============================================================================
# Run this script after downloading lab4_results.zip from the course repo:
#   https://github.com/nicolabarban/sociogenomics_2025_2026/raw/gh_pages/labs/week4/lab4_results.zip
#
# Then set your working directory to the folder where the files are extracted.
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

# Install missing packages (run once)
pkgs <- c("ggplot2", "data.table", "patchwork", "class")
missing_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(data.table)
library(ggplot2)
library(patchwork)

theme_set(theme_bw(base_size = 13))

# Download and unzip data (skip if files already present)
if (!file.exists("hapmap3_pca.eigenvec")) {
  url <- "https://github.com/nicolabarban/sociogenomics_2025_2026/raw/gh_pages/labs/week4/lab4_results.zip"
  download.file(url, destfile = "lab4_results.zip", mode = "wb")
  unzip("lab4_results.zip", overwrite = TRUE)
  cat("Files downloaded and extracted.\n")
}

# =============================================================================
# PART I — Exploring Genetic Data
# =============================================================================

# --- 1.1 SNP and individual missingness --------------------------------------

lmiss <- fread("hapmap3_summary.lmiss")
imiss <- fread("hapmap3_summary.imiss")

cat("SNPs in dataset:        ", nrow(lmiss), "\n")
cat("Individuals in dataset: ", nrow(imiss), "\n")
cat("SNP missingness range:  ", range(lmiss$F_MISS), "\n")
cat("Ind missingness range:  ", range(imiss$F_MISS), "\n")

p1 <- ggplot(lmiss, aes(x = F_MISS)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = "white") +
  xlab("Per-SNP missing rate") + ylab("Number of SNPs") +
  ggtitle("SNP missingness")

p2 <- ggplot(imiss, aes(x = F_MISS)) +
  geom_histogram(bins = 40, fill = "coral", colour = "white") +
  xlab("Per-individual missing rate") + ylab("Number of individuals") +
  ggtitle("Individual missingness")

print(p1 + p2)

# --- 1.2 Minor allele frequency distribution ---------------------------------

frq <- fread("hapmap3_summary.frq")
cat("\nMean MAF:  ", round(mean(frq$MAF), 4), "\n")
cat("Median MAF:", round(median(frq$MAF), 4), "\n")

frq[, bin := fcase(
  MAF < 0.01, "< 0.01 (rare)",
  MAF < 0.05, "0.01 - 0.05",
  MAF < 0.10, "0.05 - 0.10",
  MAF < 0.20, "0.10 - 0.20",
  default   = ">= 0.20 (common)"
)]
print(frq[, .N, by = bin][order(bin)])

p_maf <- ggplot(frq, aes(x = MAF)) +
  geom_histogram(bins = 50, fill = "darkgreen", colour = "white") +
  geom_vline(xintercept = 0.05, colour = "red",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 0.07, y = Inf, vjust = 2,
           label = "MAF = 0.05", colour = "red", size = 4) +
  xlab("Minor allele frequency") + ylab("Number of SNPs") +
  ggtitle("MAF distribution")

print(p_maf)

# --- 1.3 Hardy-Weinberg Equilibrium distribution -----------------------------

hwe     <- fread("hapmap3_summary.hwe")
hwe_all <- hwe[TEST == "ALL"]
hwe_all[, log10p := -log10(P)]

cat("\nSNPs failing HWE (p < 1e-6):", sum(hwe_all$P < 1e-6, na.rm = TRUE), "\n")
cat("Top 10 most deviant SNPs:\n")
print(hwe_all[order(P)][1:10, .(SNP, CHR, `O(HET)`, `E(HET)`, P)])

p_hwe <- ggplot(hwe_all, aes(x = log10p)) +
  geom_histogram(bins = 60, fill = "purple", colour = "white") +
  geom_vline(xintercept = 6, colour = "red",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 6.3, y = Inf, vjust = 2,
           label = "p = 1e-6", colour = "red", size = 4) +
  xlab(expression(-log[10](p))) + ylab("Number of SNPs") +
  ggtitle("HWE test statistics")

print(p_hwe)

# --- 1.4 Per-individual inbreeding coefficient -------------------------------

het <- fread("hapmap3_het.het")

outliers <- het[F < -0.15 | F > 0.15]
cat("\nHeterozygosity outliers:", nrow(outliers), "\n")
if (nrow(outliers) > 0) print(outliers[, .(FID, IID, F)])

p_het <- ggplot(het, aes(x = F)) +
  geom_histogram(bins = 50, fill = "orange", colour = "white") +
  geom_vline(xintercept = c(-0.15, 0.15), colour = "red",
             linetype = "dashed", linewidth = 0.8) +
  xlab("Inbreeding coefficient F") + ylab("Number of individuals") +
  ggtitle("Per-individual inbreeding coefficient")

print(p_het)

# =============================================================================
# PART II — Principal Component Analysis (PCA)
# =============================================================================

# --- 2.1 Scree plot ----------------------------------------------------------

eigenval       <- fread("hapmap3_pca.eigenval", header = FALSE, col.names = "eigenvalue")
eigenval[, PC  := seq_len(.N)]
eigenval[, pct := eigenvalue / sum(eigenvalue) * 100]

cat("\nVariance explained per PC:\n")
print(eigenval[, .(PC, pct = round(pct, 2))])

p_scree <- ggplot(eigenval, aes(x = PC, y = pct)) +
  geom_col(fill = "steelblue", colour = "white") +
  geom_line(aes(group = 1)) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:20) +
  xlab("Principal Component") + ylab("Variance explained (%)") +
  ggtitle("Scree plot")

print(p_scree)

# --- 2.2 Load PCA scores and population labels -------------------------------

pc_cols <- c("FID", "IID", paste0("PC", 1:20))
pca     <- fread("hapmap3_pca.eigenvec", header = FALSE, col.names = pc_cols)

geo <- fread("1kg_samples.txt", sep = "\t", header = TRUE)
setnames(geo, "Sample name", "IID")

data <- merge(
  pca,
  geo[, .(`IID`, `Population code`, `Population name`,
           `Superpopulation code`, `Superpopulation name`)],
  by = "IID"
)

cat("\nIndividuals with population labels:", nrow(data), "\n")
print(table(data$`Superpopulation name`))

# --- 2.3 PC1 vs PC2 — superpopulation ----------------------------------------

p_pca_super <- ggplot(data, aes(x = PC1, y = PC2,
                                colour = `Superpopulation name`)) +
  geom_point(alpha = 0.7, size = 1.5) +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "Superpopulation",
       title  = "PCA — continental ancestry (PC1 vs PC2)")

print(p_pca_super)

# --- 2.4 PC1 vs PC2 — sub-population -----------------------------------------

p_pca_sub <- ggplot(data, aes(x = PC1, y = PC2,
                               colour = `Population name`)) +
  geom_point(alpha = 0.7, size = 1.5) +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "Population",
       title  = "PCA — sub-population (PC1 vs PC2)") +
  theme(legend.text = element_text(size = 7))

print(p_pca_sub)

# --- 2.5 PC1 vs PC3 ----------------------------------------------------------

p_pca_pc3 <- ggplot(data, aes(x = PC1, y = PC3,
                               colour = `Superpopulation name`)) +
  geom_point(alpha = 0.7, size = 1.5) +
  xlab("PC1") + ylab("PC3") +
  labs(colour = "Superpopulation", title = "PC1 vs PC3")

print(p_pca_pc3)

# --- 2.6 Within-European PCA -------------------------------------------------

pc_cols_eur <- c("FID", "IID", paste0("PC", 1:10))
pca_eur     <- fread("pca_EUR.eigenvec", header = FALSE, col.names = pc_cols_eur)
data_eur    <- merge(pca_eur, geo[, .(`IID`, `Population name`)], by = "IID")

p_eur <- ggplot(data_eur, aes(x = PC1, y = PC2,
                               colour = `Population name`)) +
  geom_point(alpha = 0.8, size = 2) +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "European population",
       title  = "PCA within European populations")

print(p_eur)

# =============================================================================
# PART III — Detecting Population Outliers
# =============================================================================

eur_mean_pc1 <- mean(data[`Superpopulation code` == "EUR", PC1])
eur_mean_pc2 <- mean(data[`Superpopulation code` == "EUR", PC2])
eur_sd_pc1   <- sd(data[`Superpopulation code`   == "EUR", PC1])
eur_sd_pc2   <- sd(data[`Superpopulation code`   == "EUR", PC2])

data[, eur_like := abs(PC1 - eur_mean_pc1) < 3 * eur_sd_pc1 &
                   abs(PC2 - eur_mean_pc2) < 3 * eur_sd_pc2]

cat("\nIndividuals within 3 SD of EUR centroid:", sum(data$eur_like), "\n")

p_eur_sel <- ggplot(data, aes(x = PC1, y = PC2,
                               colour = `Superpopulation code`,
                               shape  = eur_like)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_shape_manual(values = c(4, 16),
                     labels = c("Excluded", "EUR-like (kept)")) +
  labs(colour = "Superpopulation", shape = "Selection",
       title  = "EUR-like individuals (within 3 SD of EUR centroid)")

print(p_eur_sel)

# Save sample list for PLINK
eur_keep <- data[eur_like == TRUE, .(FID, IID)]
fwrite(eur_keep, "samples_EUR_like.txt", sep = " ", col.names = FALSE)
cat("Saved", nrow(eur_keep), "EUR-like individuals to samples_EUR_like.txt\n")
cat("Upload this file to Cloud Shell and run:\n")
cat("  plink --bfile hapmap3_qc --keep samples_EUR_like.txt --make-bed --out hapmap3_EUR\n")

# =============================================================================
# BONUS — Ancestry Prediction with k-Nearest Neighbours
# =============================================================================

library(class)

pc_features <- paste0("PC", 1:10)
X <- as.matrix(data[, ..pc_features])
y <- data$`Superpopulation code`

set.seed(42)
n         <- nrow(data)
train_idx <- sample(n, size = floor(0.8 * n), replace = FALSE)
test_idx  <- setdiff(seq_len(n), train_idx)

train_X <- X[train_idx, ];  test_X <- X[test_idx, ]
train_y <- y[train_idx];    test_y <- y[test_idx]

predicted <- knn(train = train_X, test = test_X, cl = train_y, k = 5)

conf_mat <- table(Predicted = predicted, True = test_y)
accuracy <- sum(diag(conf_mat)) / sum(conf_mat)

cat("\nk-NN confusion matrix (k=5):\n")
print(conf_mat)
cat("\nOverall accuracy:", round(accuracy * 100, 1), "%\n")

per_pop <- diag(conf_mat) / colSums(conf_mat)
cat("Per-population accuracy:\n")
print(round(per_pop * 100, 1))

# Visualise predictions
predicted_all <- knn(train = train_X, test = X, cl = train_y, k = 5)
data[, predicted := as.character(predicted_all)]
data[, correct   := predicted == `Superpopulation code`]

p_knn <- ggplot(data, aes(x = PC1, y = PC2,
                           colour = `Superpopulation code`,
                           shape  = correct)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_shape_manual(values = c(4, 16),
                     labels = c("Misclassified", "Correct")) +
  labs(colour = "True superpopulation", shape = "Classification",
       title  = "k-NN ancestry predictions (k = 5)")

print(p_knn)

cat("\nDone! All plots displayed.\n")
