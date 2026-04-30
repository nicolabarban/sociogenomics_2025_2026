## ========================================================================
## Lab 8 (part B) -- G x E: replicating Walter et al. (2016, JAMA)
## Sociogenomics, Master Course 2025/2026, University of Bologna
## ------------------------------------------------------------------------
## Data: anonymised HRS subset (N ~ 8,450; PGS3 BMI, RAND HRS phenotypes)
## Source chapter: Felix Tropf, "Applying PGS: GxE applications"
## ========================================================================

## --- 0. Setup ----------------------------------------------------------

# install.packages(c("data.table", "ggplot2", "interactions", "jtools"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(interactions)
  library(jtools)
})

data_path <- "~/sociogenomics_2025_2026/labs/week8/data/hrs_lab8.csv"

## --- 1. Load and inspect -----------------------------------------------

d <- fread(data_path)

str(d)
summary(d[, .(BMI_AV, Age_AV, birth_year, pgs_bmi)])
table(d$sex)

ggplot(d, aes(birth_year)) +
  geom_histogram(binwidth = 2, fill = "steelblue") +
  labs(x = "Year of birth", y = "Count") +
  theme_minimal()

## --- 1.5 PC inspection: flag non-EUR individuals -----------------------

pcs <- paste0("pc", 1:10)

ggplot(d, aes(pc1, pc2)) +
  geom_point(alpha = 0.4, size = 0.6, colour = "steelblue") +
  labs(x = "PC1", y = "PC2",
       title = "Genetic PCs - HRS lab subset") +
  theme_minimal()

PC  <- as.matrix(d[, ..pcs])
mu  <- colMeans(PC)
S   <- cov(PC)
mhd <- mahalanobis(PC, center = mu, cov = S)

cutoff <- qchisq(0.999, df = length(pcs))
d[, pc_outlier := mhd > cutoff]
print(table(d$pc_outlier))

ggplot(d, aes(pc1, pc2, colour = pc_outlier)) +
  geom_point(alpha = 0.5, size = 0.7) +
  scale_colour_manual(values = c(`FALSE` = "steelblue", `TRUE` = "red"),
                      labels = c(`FALSE` = "EUR-like", `TRUE` = "outlier")) +
  labs(x = "PC1", y = "PC2", colour = NULL,
       title = "PC outliers via Mahalanobis distance",
       subtitle = sprintf("%d flagged of %d (chi-sq 99.9%% cutoff)",
                          sum(d$pc_outlier), nrow(d))) +
  theme_minimal()

# Sensitivity subset (use d_eur in place of d below if you want to check)
d_eur <- d[pc_outlier == FALSE]
cat("After dropping PC outliers: n =", nrow(d_eur), "\n")

## --- 2. Main effect ----------------------------------------------------

fmla_main <- as.formula(paste(
  "BMI_AV ~ pgs_bmi + birth_year + sex +",
  paste(pcs, collapse = " + ")
))

m_main <- lm(fmla_main, data = d)
summ(m_main, digits = 3)
coef(summary(m_main))["pgs_bmi", ]

## --- 3. G x birth-year interaction -------------------------------------

d[, by_c := birth_year - 1944]    # centre at the Walter pivot

fmla_gxe <- as.formula(paste(
  "BMI_AV ~ pgs_bmi * by_c + sex +",
  paste(pcs, collapse = " + ")
))

m_gxe <- lm(fmla_gxe, data = d)
summ(m_gxe, digits = 3)
coef(summary(m_gxe))[c("pgs_bmi", "by_c", "pgs_bmi:by_c"), ]

anova(m_main, m_gxe)
delta_r2 <- summary(m_gxe)$r.squared - summary(m_main)$r.squared
cat(sprintf("delta R^2 from G x E = %.5f\n", delta_r2))

## --- 4. Visualise the interaction --------------------------------------

# Predicted BMI by PGS, separate lines for four cohorts
interact_plot(m_gxe,
              pred = "pgs_bmi", modx = "by_c",
              modx.values = c(-30, -10, 10, 30),
              interval = TRUE, int.width = 0.95,
              x.label = "PGS-BMI (SD)",
              y.label = "Predicted BMI",
              legend.main = "Birth year - 1944")

# Conditional PGS slope at a grid of birth years
ss <- sim_slopes(m_gxe, pred = "pgs_bmi", modx = "by_c",
                 modx.values = seq(-40, 35, by = 5))
print(ss)

# Johnson-Neyman: where is the PGS slope significantly non-zero?
johnson_neyman(m_gxe, pred = "pgs_bmi", modx = "by_c", alpha = 0.05)

## --- 5. Walter design: pre-1944 vs post-1944 split ---------------------

d[, walter := factor(ifelse(birth_year < 1944, "pre1944", "post1944"),
                     levels = c("pre1944", "post1944"))]
table(d$walter)

run_walter <- function(group) {
  fit <- lm(fmla_main, data = d[walter == group])
  c(beta = unname(coef(fit)["pgs_bmi"]),
    se   = unname(coef(summary(fit))["pgs_bmi", "Std. Error"]),
    n    = nobs(fit))
}

walter_pre  <- run_walter("pre1944")
walter_post <- run_walter("post1944")

walter_df <- data.frame(
  group = factor(c("pre1944", "post1944"),
                 levels = c("pre1944", "post1944")),
  beta  = c(walter_pre["beta"],  walter_post["beta"]),
  se    = c(walter_pre["se"],    walter_post["se"]),
  n     = c(walter_pre["n"],     walter_post["n"])
)
print(walter_df)

ggplot(walter_df, aes(group, beta)) +
  geom_pointrange(aes(ymin = beta - 1.96 * se,
                      ymax = beta + 1.96 * se), size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  labs(x = NULL,
       y = "PGS-BMI slope (BMI per SD of PGS)",
       title = "Replication of Walter et al. (2016)",
       subtitle = paste("HRS, anonymised teaching subset, N =",
                        sum(walter_df$n))) +
  theme_minimal(base_size = 12)

## --- 6. Optional 3-way: PGS x cohort x sex (Herd-style) ---------------

fmla_3way <- as.formula(paste(
  "BMI_AV ~ pgs_bmi * by_c * sex +",
  paste(pcs, collapse = " + ")
))

m_3way <- lm(fmla_3way, data = d)
summ(m_3way, digits = 3)
anova(m_gxe, m_3way)

interact_plot(m_3way,
              pred = "pgs_bmi", modx = "by_c", mod2 = "sex",
              modx.values = c(-20, 20),
              interval = TRUE, int.width = 0.95,
              x.label = "PGS-BMI (SD)",
              y.label = "Predicted BMI",
              legend.main = "Birth year - 1944")
