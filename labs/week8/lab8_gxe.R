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

# Path to the dataset. Change to wherever you saved hrs_lab8.csv.
# If you downloaded it next to this script, just "hrs_lab8.csv" works.
data_path <- "hrs_lab8.csv"

## --- 1. Load and inspect -----------------------------------------------

d <- fread(data_path)

str(d)
summary(d[, .(BMI_AV, birth_year, pgs_bmi)])
table(d$sex)

ggplot(d, aes(birth_year)) +
  geom_histogram(binwidth = 2, fill = "steelblue") +
  labs(x = "Year of birth", y = "Count") +
  theme_minimal()

# Distribution of the polygenic index (PGI / PGS)
ggplot(d, aes(pgs_bmi)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(x = "PGS-BMI (z-scored)", y = "Count",
       title = "Polygenic index for BMI",
       subtitle = "After standardisation: mean 0, SD 1") +
  theme_minimal()

## --- 1.1 Scatter and binscatter: BMI vs PGS ---------------------------

# Raw scatter (all individuals): dense cloud, hard to read
ggplot(d, aes(pgs_bmi, BMI_AV)) +
  geom_point(alpha = 0.15, size = 0.5, colour = "steelblue") +
  geom_smooth(method = "lm", se = FALSE,
              colour = "firebrick", linewidth = 0.7) +
  labs(x = "PGS-BMI (z-scored)", y = "BMI",
       title = "Raw scatter: BMI by PGS-BMI") +
  theme_minimal()

# Binscatter: 25 equal-width PGS bins, mean BMI within each bin
ggplot(d, aes(pgs_bmi, BMI_AV)) +
  stat_summary_bin(fun = mean, bins = 25, geom = "point",
                   colour = "steelblue", size = 2) +
  geom_smooth(method = "lm", se = FALSE,
              colour = "firebrick", linewidth = 0.7) +
  labs(x = "PGS-BMI (z-scored)", y = "Mean BMI in bin",
       title = "Binscatter: BMI by PGS-BMI",
       subtitle = "25 equal-width bins, OLS fit on raw data") +
  theme_minimal()

## --- 1.5 PC inspection: visualise structure ----------------------------

pcs <- paste0("pc", 1:10)

ggplot(d, aes(pc1, pc2)) +
  geom_point(alpha = 0.4, size = 0.6, colour = "steelblue") +
  labs(x = "PC1", y = "PC2",
       title = "Genetic PCs - HRS lab subset") +
  theme_minimal()

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
              modx.labels = c("1914 cohort", "1934 cohort",
                              "1954 cohort", "1974 cohort"),
              interval = TRUE, int.width = 0.95,
              x.label = "PGS-BMI (SD)",
              y.label = "Predicted BMI",
              legend.main = "Birth cohort")

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
              modx.labels = c("1924 cohort", "1964 cohort"),
              interval = TRUE, int.width = 0.95,
              x.label = "PGS-BMI (SD)",
              y.label = "Predicted BMI",
              legend.main = "Birth cohort")

## --- 7. G x E with education as the moderator -------------------------
##  Compare three specifications against Walter's birth-year result.
##  None of these are exogenous like birth_year -- read descriptively.

## (a) Continuous, z-scored years of schooling
d[, edu_z := scale(raedyrs)]
fmla_edu <- as.formula(paste(
  "BMI_AV ~ pgs_bmi * edu_z + birth_year + sex +",
  paste(pcs, collapse = " + ")
))
m_edu <- lm(fmla_edu, data = d)
cat("\n--- (a) continuous edu (z) ---\n")
print(round(coef(summary(m_edu))[c("pgs_bmi","edu_z","pgs_bmi:edu_z"), ], 4))

interact_plot(m_edu, pred = "pgs_bmi", modx = "edu_z",
              modx.values = c(-1, 0, 1),
              interval = TRUE, int.width = 0.95,
              x.label = "PGS-BMI (SD)",
              y.label = "Predicted BMI",
              legend.main = "Education (SD)")

## (b) Three-level education (HS / some college / college+)
d[, edu3 := cut(raedyrs,
                breaks = c(-Inf, 12, 15, Inf),
                labels = c("HS_or_less", "some_coll", "coll_plus"))]
print(table(d$edu3))

fmla_e3 <- as.formula(paste(
  "BMI_AV ~ pgs_bmi * edu3 + birth_year + sex +",
  paste(pcs, collapse = " + ")
))
m_e3 <- lm(fmla_e3, data = d)
cat("\n--- (b) 3-level education ---\n")
co <- coef(summary(m_e3))
print(round(co[grep("pgs_bmi|edu3", rownames(co)), ], 4))

## (c) Parental education (childhood SES proxy, cleaner E)
d_par <- d[!is.na(rameduc) & !is.na(rafeduc)]
d_par[, par_edu := scale((rameduc + rafeduc) / 2)]

fmla_pe <- as.formula(paste(
  "BMI_AV ~ pgs_bmi * par_edu + birth_year + sex +",
  paste(pcs, collapse = " + ")
))
m_pe <- lm(fmla_pe, data = d_par)
cat("\n--- (c) parental education (z) ---\n")
print(round(coef(summary(m_pe))[c("pgs_bmi","par_edu","pgs_bmi:par_edu"), ], 4))
cat(sprintf("  n = %d (drop rows missing parental edu)\n", nrow(d_par)))

cat("\nContrast: birth_year G x E is significant (~0.025, p < 1e-7),\n",
    "but every education spec gives a near-null interaction.\n", sep = "")
