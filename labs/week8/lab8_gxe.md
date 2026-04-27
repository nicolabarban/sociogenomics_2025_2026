# Lab 8 (part B). G$\times$E: replicating Walter et al. (2016)

In this lab we put the G$\times$E theory from class to work and **replicate the central finding of Walter et al. (2016, *JAMA*)**: the polygenic score for body-mass index (PGS-BMI) predicts BMI much more strongly in cohorts born after the Second World War than in cohorts born before. The interpretation is that the *post-war obesogenic environment* "unlocked" the BMI PGS — a textbook example of how an environmental shift modulates the phenotypic expression of genetic predisposition.

We use anonymised data from the [Health and Retirement Study](https://hrs.isr.umich.edu/), distributed by the instructor under the HRS data-use agreement. The teaching subset includes the BMI PGS (PGS3 release, GIANT 2015 base GWAS), wave-averaged BMI and age, birth year, sex, and 10 genetic principal components for $N \approx 8{,}450$ European-ancestry participants born 1905–1980.

The lab covers:

* Loading and inspecting an HRS-style PGS dataset
* Estimating the **main effect** of a PGS on a phenotype with appropriate controls (age, sex, PCs)
* Adding a **G$\times$E interaction** with birth year and reading $\beta_{GE}$
* Visualising the interaction: marginal PGS slope by birth year, and predicted BMI by PGS$\times$cohort
* Replicating the Walter design with a **pre-1944 vs post-1944 split**
* Optional 3-way extension à la Herd et al.\ (2019): does the interaction differ by sex?

---

## 0. Getting started

This lab can be run on **Google Colab**, **Posit Cloud**, or any local **R** install (≥ 4.1). It does **not** need PLINK or PRSice — all the genetic work has already been done; we work directly with the pre-computed PGS.

### 0.1 Get the data

The dataset `hrs_lab8.csv` (≈ 1.3 MB) is **not** in the public GitHub repository because of the HRS data-use agreement. The instructor will share a private download link on the course Slack/Moodle. Save the file in a folder of your choice; the script below assumes:

```
~/sociogenomics_2025_2026/labs/week8/data/hrs_lab8.csv
```

If you are on Colab, upload the file with `files.upload()` and adjust the path accordingly.

### 0.2 Install R packages (one-off)

In R:

```r
install.packages(c("data.table", "ggplot2", "interactions", "jtools"))
```

`interactions::interact_plot` is what produces the marginal-slope plots. `jtools::summ` gives nicely formatted regression summaries. `data.table::fread` reads CSVs fast.

### 0.3 The companion script

A self-contained R script for this lab is at `lab8_gxe.R` in the same folder. The code blocks below are extracted from it — you can either copy/paste them into a fresh R session, or open the script and run it line-by-line.

---

## 1. Load and inspect the data

```r
library(data.table)
library(ggplot2)
library(interactions)

d <- fread("~/sociogenomics_2025_2026/labs/week8/data/hrs_lab8.csv")
str(d)
summary(d[, .(BMI_AV, Age_AV, birth_year, pgs_bmi)])
table(d$sex)
```

The variables are:

| name | meaning |
|---|---|
| `id` | anonymised subject identifier (`S00001`…) |
| `pgs_bmi` | BMI PGS, **standardised** to mean 0 / SD 1 (PGS3 / GIANT 2015) |
| `BMI_AV` | mean BMI across HRS waves 1–12 |
| `Age_AV` | mean age (years) across the same waves |
| `birth_year` | year of birth (1905–1980) |
| `sex` | factor: `male`, `female` |
| `pc1_5a` … `pc6_10e` | first 10 genetic PCs |

> **Why standardise the PGS?** With $z$-scored PGS, $\beta_G$ is the change in BMI for a one-SD increase in genetic predisposition. This makes effects comparable across PGS, traits, and papers.

> **Why include the PCs?** Population structure produces correlated allele frequencies and outcome differences. Without PCs the PGS coefficient is partly confounded by ancestry. Always include at least 5–10 PCs.

Quick sanity-check plot:

```r
ggplot(d, aes(birth_year)) + geom_histogram(binwidth = 2, fill = "steelblue") +
  labs(x = "Year of birth", y = "Count") + theme_minimal()
```

---

## 2. Main effect: does the PGS predict BMI?

A baseline OLS without interaction:

```r
pcs <- paste0(c("pc1_5", "pc6_10"), rep(letters[1:5], each = 2))
fmla_main <- as.formula(paste(
  "BMI_AV ~ pgs_bmi + birth_year + Age_AV + sex +",
  paste(pcs, collapse = " + ")
))

m_main <- lm(fmla_main, data = d)
coef(summary(m_main))["pgs_bmi", ]
```

> **Checkpoint.** You should find $\beta_{\text{PGS}} \approx 1.5$ BMI units per SD of PGS, with $p < 10^{-100}$. A 1-SD higher PGS is associated with $\approx 1.5$ kg/m² higher BMI on average.

What does this *not* tell us? Whether the PGS effect is the same for everyone — across cohorts, sexes, environments. That is a G$\times$E question.

---

## 3. The interaction model

Walter's hypothesis: the slope of BMI on PGS depends on birth cohort. Centring birth year at 1944 makes the main effect of PGS interpretable as the slope *at* the Walter pivot:

```r
d[, by_c := birth_year - 1944]

fmla_gxe <- as.formula(paste(
  "BMI_AV ~ pgs_bmi * by_c + Age_AV + sex +",
  paste(pcs, collapse = " + ")
))

m_gxe <- lm(fmla_gxe, data = d)
coef(summary(m_gxe))[c("pgs_bmi", "by_c", "pgs_bmi:by_c"), ]
```

Interpretation:

* $\beta_{\text{pgs\_bmi}}$ — slope of BMI on PGS for someone born **in 1944** (because we centred there).
* $\beta_{\text{by\_c}}$ — main effect of birth year on BMI, holding PGS fixed (BMI rises across cohorts on average — the secular obesity trend).
* $\beta_{\text{pgs\_bmi:by\_c}}$ — **the G$\times$E coefficient.** It tells you how much the PGS slope **changes per additional year of birth**.

> **Checkpoint.** Expect $\beta_{GE} \approx 0.025$, $p < 10^{-7}$. Each extra year of birth adds $\approx 0.025$ BMI units to the slope. Across the 75-year span of the data this is a substantial cumulative shift.

To test the interaction formally:

```r
anova(m_main, m_gxe)         # F-test on the interaction term
summary(m_gxe)$r.squared - summary(m_main)$r.squared   # delta R^2
```

The $\Delta R^2$ from the interaction is small in absolute terms but the F-test is highly significant. **Tiny variance, big mechanism**: G$\times$E rarely moves $R^2$ much, but it changes who the coefficient applies to.

---

## 4. Visualising the interaction

Two complementary plots.

### 4.1 Marginal slope of PGS by birth year

```r
interact_plot(m_gxe,
              pred = "pgs_bmi", modx = "by_c",
              modx.values = c(-30, -10, 10, 30),    # 1914, 1934, 1954, 1974
              interval = TRUE, int.width = 0.95,
              x.label = "PGS-BMI (SD)", y.label = "Predicted BMI",
              legend.main = "Birth year (centred at 1944)")
```

You should see a fan: lines for older cohorts are *flatter*, lines for younger cohorts are *steeper*. Same PGS, more BMI in younger cohorts.

### 4.2 Slope of PGS as a function of birth year (with CI band)

This is the iconic Walter-style plot:

```r
library(jtools)
sim_slopes(m_gxe, pred = "pgs_bmi", modx = "by_c",
           modx.values = seq(-40, 35, by = 5))
```

`sim_slopes` returns the conditional PGS slope at each birth-year value. To plot the trajectory, use `interactions::johnson_neyman`:

```r
johnson_neyman(m_gxe, pred = "pgs_bmi", modx = "by_c", alpha = 0.05)
```

The Johnson-Neyman plot shows where the PGS slope is statistically distinguishable from zero across cohorts. In our data, the slope is positive everywhere — the PGS always matters — but it is much *larger* in younger cohorts.

---

## 5. The Walter design: pre-1944 vs post-1944 split

The original paper splits the sample into two cohort groups and estimates the PGS effect separately in each. This is a more conservative test that does not impose linearity.

```r
d[, walter := factor(ifelse(birth_year < 1944, "pre1944", "post1944"),
                     levels = c("pre1944", "post1944"))]
table(d$walter)

run_walter <- function(group) {
  fit <- lm(fmla_main, data = d[walter == group])
  c(beta = coef(fit)["pgs_bmi"],
    se   = coef(summary(fit))["pgs_bmi", "Std. Error"],
    n    = nobs(fit))
}

walter_pre  <- run_walter("pre1944")
walter_post <- run_walter("post1944")
walter_pre; walter_post
```

> **Checkpoint.** Expect $\beta_{\text{pre}} \approx 1.20$, $\beta_{\text{post}} \approx 1.98$ — the PGS effect is roughly **65% larger in the post-1944 cohort**. This is the core Walter result.

Show it on one plot:

```r
walter_df <- rbind(
  data.frame(group = "pre1944",  beta = walter_pre["beta.pgs_bmi"],
             se = walter_pre["se"],  n = walter_pre["n"]),
  data.frame(group = "post1944", beta = walter_post["beta.pgs_bmi"],
             se = walter_post["se"], n = walter_post["n"])
)
walter_df$group <- factor(walter_df$group, levels = c("pre1944", "post1944"))

ggplot(walter_df, aes(group, beta)) +
  geom_pointrange(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se),
                  size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  labs(x = NULL, y = "PGS-BMI slope (BMI per SD)",
       title = "Replication of Walter et al. (2016)") +
  theme_minimal(base_size = 12)
```

---

## 6. Discussion

Answer the following before next class. Two-three sentences each.

1. **rGE.** Birth year is a "clean" environment in the sense that genotype cannot cause it. But could there be **selection** that makes our pre- vs post-1944 comparison non-causal? (Hint: HRS samples people who survive to old age — what does that imply for the pre-1944 group?)
2. **Population, not individual.** The post-1944 PGS slope is 1.98 BMI units per SD. Can we conclude that *for a specific person*, "65% of their BMI risk is genetic in the post-war environment"? Why or why not?
3. **Mechanism.** The interaction term $\beta_{GE}$ is positive but small per year. Suggest two *concrete* environmental changes between 1920 and 1980 that could plausibly drive this trend. How would you test which one is doing the work?
4. **Portability.** All participants in our subset are of European ancestry. If we re-ran this analysis in an East Asian sample, what would you predict, and why?

---

## 7. Optional extension: a 3-way interaction (Herd et al.\ 2019)

Herd and colleagues showed that the PGS-EA × cohort pattern was strongly moderated by **gender**: in older cohorts the PGS predicted education much less in women than in men, with the gap narrowing over time. The same logic — genetic potential gated by structural opportunity — could apply to BMI. Test it:

```r
fmla_3way <- as.formula(paste(
  "BMI_AV ~ pgs_bmi * by_c * sex + Age_AV +",
  paste(pcs, collapse = " + ")
))

m_3way <- lm(fmla_3way, data = d)
summary(m_3way)
anova(m_gxe, m_3way)
```

Interpret the three-way coefficient `pgs_bmi:by_c:sexfemale`. Does the obesogenic-environment effect on the BMI PGS look the same for men and women? Plot the predicted BMI by PGS, faceted by sex and split by birth-year tertile.

---

## References

* Walter, S., Mejía-Guevara, I., Estrada, K., Liu, S. Y., & Glymour, M. M. (2016). *Association of a Genetic Risk Score with Body Mass Index Across Different Birth Cohorts.* **JAMA**, 316(1), 63–69.
* Herd, P., Freese, J., Sicinski, K., Domingue, B. W., Mullan Harris, K., Wei, C., & Hauser, R. M. (2019). *Genes, Gender Inequality, and Educational Attainment.* **American Sociological Review**, 84(6), 1069–1098.
* Domingue, B. W., Conley, D., Fletcher, J., & Boardman, J. D. (2016). *Cohort effects in the genetic influence on smoking.* **Behavior Genetics**, 46(1), 31–42.
* Tropf, F. C. (2019). *Applying Polygenic Scores: G$\times$E applications.* In Conley & Fletcher (eds), *The Genome Factor / textbook companion materials* — chapter and code on which this lab is built.
* Health and Retirement Study (HRS). Public-use dataset. Polygenic Score release 3 (PGS3), Ware et al.\ 2018.
