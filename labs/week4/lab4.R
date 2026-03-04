# =============================================================================
# Lab 4 — Genetic Tests and Principal Component Analysis
# Sociogenomics 2025/2026 · University of Bologna · Prof. Nicola Barban
# =============================================================================
# Requirements: R with ggplot2 installed.
# Run in RStudio or any R environment — data is downloaded automatically.
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)
theme_set(theme_bw(base_size = 13))

# Download and unzip data from the course repository (skip if already present)
if (!file.exists("hapmap3_pca.eigenvec")) {
  url <- "https://github.com/nicolabarban/sociogenomics_2025_2026/raw/gh_pages/labs/week4/lab4_results.zip"
  download.file(url, destfile = "lab4_results.zip", mode = "wb")
  unzip("lab4_results.zip", overwrite = TRUE)
  cat("Data downloaded and extracted.\n")
}

# =============================================================================
# PART I — Exploring Genetic Data
# =============================================================================

# --- 1.1 SNP and individual missingness --------------------------------------

lmiss <- read.table("hapmap3_summary.lmiss", header = TRUE)
imiss <- read.table("hapmap3_summary.imiss", header = TRUE)

ggplot(lmiss, aes(x = F_MISS)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = "white") +
  xlab("Per-SNP missing rate") + ylab("Number of SNPs") +
  ggtitle("SNP missingness")

ggplot(imiss, aes(x = F_MISS)) +
  geom_histogram(bins = 40, fill = "coral", colour = "white") +
  xlab("Per-individual missing rate") + ylab("Number of individuals") +
  ggtitle("Individual missingness")

# --- 1.2 Minor allele frequency distribution ---------------------------------

frq <- read.table("hapmap3_summary.frq", header = TRUE)
cat("\nMean MAF:  ", round(mean(frq$MAF), 4), "\n")
cat("Median MAF:", round(median(frq$MAF), 4), "\n")

ggplot(frq, aes(x = MAF)) +
  geom_histogram(bins = 50, fill = "darkgreen", colour = "white") +
  geom_vline(xintercept = 0.05, colour = "red",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 0.07, y = Inf, vjust = 2,
           label = "MAF = 0.05", colour = "red", size = 4) +
  xlab("Minor allele frequency") + ylab("Number of SNPs") +
  ggtitle("MAF distribution across all SNPs")

# =============================================================================
# PART II — Principal Component Analysis (PCA)
# =============================================================================

# --- 2.1 Scree plot ----------------------------------------------------------

eigenval      <- read.table("hapmap3_pca.eigenval", header = FALSE)
colnames(eigenval) <- "eigenvalue"
eigenval$PC   <- seq_len(nrow(eigenval))
eigenval$pct  <- eigenval$eigenvalue / sum(eigenval$eigenvalue) * 100

cat("\nVariance explained per PC:\n")
print(eigenval[, c("PC", "pct")])

ggplot(eigenval, aes(x = PC, y = pct)) +
  geom_col(fill = "steelblue", colour = "white") +
  geom_line(aes(group = 1)) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:20) +
  xlab("Principal Component") + ylab("Variance explained (%)") +
  ggtitle("Scree plot")

# --- 2.2 Load PCA scores and population labels -------------------------------

pc_cols <- c("FID", "IID", paste0("PC", 1:20))
pca     <- read.table("hapmap3_pca.eigenvec", header = FALSE, col.names = pc_cols)

# 1000 Genomes sample metadata
# Note: read.table converts spaces in column names to dots
# e.g. "Sample name" -> "Sample.name", "Superpopulation name" -> "Superpopulation.name"
geo <- read.table("1kg_samples.txt", sep = "\t", header = TRUE)

data <- merge(pca, geo[, c("Sample.name", "Population.code", "Population.name",
                             "Superpopulation.code", "Superpopulation.name")],
              by.x = "IID", by.y = "Sample.name")

cat("\nIndividuals with population labels:", nrow(data), "\n")
print(table(data$Superpopulation.name))

# --- 2.3 PC1 vs PC2 — superpopulation ----------------------------------------

ggplot(data, aes(x = PC1, y = PC2, colour = Superpopulation.name)) +
  geom_point(alpha = 0.7, size = 1.5) +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "Superpopulation",
       title  = "PCA — continental ancestry (PC1 vs PC2)")

# --- 2.4 PC1 vs PC2 — sub-population -----------------------------------------

ggplot(data, aes(x = PC1, y = PC2, colour = Population.name)) +
  geom_point(alpha = 0.7, size = 1.5) +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "Population",
       title  = "PCA — sub-population (PC1 vs PC2)") +
  theme(legend.text = element_text(size = 7))

# --- 2.5 PC1 vs PC3 ----------------------------------------------------------

ggplot(data, aes(x = PC1, y = PC3, colour = Superpopulation.name)) +
  geom_point(alpha = 0.7, size = 1.5) +
  xlab("PC1") + ylab("PC3") +
  labs(colour = "Superpopulation", title = "PC1 vs PC3")

# --- 2.6 Within-African PCA --------------------------------------------------

pc_cols_afr <- c("FID", "IID", paste0("PC", 1:10))
pca_afr     <- read.table("pca_AFR.eigenvec", header = FALSE, col.names = pc_cols_afr)
data_afr    <- merge(pca_afr, geo[, c("Sample.name", "Population.name")],
                     by.x = "IID", by.y = "Sample.name")

ggplot(data_afr, aes(x = PC1, y = PC2, colour = Population.name)) +
  geom_point(alpha = 0.8, size = 2) +
  xlab("PC1") + ylab("PC2") +
  labs(colour = "African population",
       title  = "PCA within African populations")

# =============================================================================
# PART III — Detecting Population Outliers
# =============================================================================

eur_mean_pc1 <- mean(data$PC1[data$Superpopulation.code == "EUR"])
eur_mean_pc2 <- mean(data$PC2[data$Superpopulation.code == "EUR"])
eur_sd_pc1   <- sd(data$PC1[data$Superpopulation.code == "EUR"])
eur_sd_pc2   <- sd(data$PC2[data$Superpopulation.code == "EUR"])

data$eur_like <- abs(data$PC1 - eur_mean_pc1) < 3 * eur_sd_pc1 &
                 abs(data$PC2 - eur_mean_pc2) < 3 * eur_sd_pc2

cat("\nIndividuals within 3 SD of EUR centroid:", sum(data$eur_like), "\n")

ggplot(data, aes(x = PC1, y = PC2,
                 colour = Superpopulation.code,
                 shape  = eur_like)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_shape_manual(values = c(4, 16),
                     labels = c("Excluded", "EUR-like (kept)")) +
  labs(colour = "Superpopulation", shape = "Selection",
       title  = "EUR-like individuals (within 3 SD of EUR centroid)")

# Save sample list for PLINK
eur_keep <- data[data$eur_like, c("FID", "IID")]
write.table(eur_keep, "samples_EUR_like.txt",
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Saved", nrow(eur_keep), "EUR-like individuals to samples_EUR_like.txt\n")

# =============================================================================
# BONUS — Ancestry Prediction with k-Nearest Neighbours
# =============================================================================

if (!requireNamespace("class", quietly = TRUE)) install.packages("class")
library(class)

pc_features <- paste0("PC", 1:10)
X <- as.matrix(data[, pc_features])
y <- data$Superpopulation.code

set.seed(42)
n         <- nrow(data)
train_idx <- sample(n, size = floor(0.8 * n), replace = FALSE)
test_idx  <- setdiff(seq_len(n), train_idx)

train_X <- X[train_idx, ];  test_X <- X[test_idx, ]
train_y <- y[train_idx];    test_y <- y[test_idx]

predicted <- knn(train = train_X, test = test_X, cl = train_y, k = 5)

conf_mat <- table(Predicted = predicted, True = test_y)
accuracy <- sum(diag(conf_mat)) / sum(conf_mat)

cat("\nk-NN confusion matrix (k = 5):\n")
print(conf_mat)
cat("\nOverall accuracy:", round(accuracy * 100, 1), "%\n")
cat("Per-population accuracy:\n")
print(round(diag(conf_mat) / colSums(conf_mat) * 100, 1))

# Visualise predictions
predicted_all    <- knn(train = train_X, test = X, cl = train_y, k = 5)
data$predicted   <- as.character(predicted_all)
data$correct     <- data$predicted == data$Superpopulation.code

ggplot(data, aes(x = PC1, y = PC2,
                 colour = Superpopulation.code,
                 shape  = correct)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_shape_manual(values = c(4, 16),
                     labels = c("Misclassified", "Correct")) +
  labs(colour = "True superpopulation", shape = "Classification",
       title  = "k-NN ancestry predictions (k = 5)")

cat("\nDone!\n")
