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

This lab is designed to be run **locally in RStudio** (R ≥ 4.1). It does **not** need PLINK or PRSice — all the genetic work has already been done; we work directly with the pre-computed PGS.

### 0.1 Get the data

The dataset `hrs_lab8.csv` (≈ 1.4 MB) is **not** stored in the GitHub repository because of the HRS data-use agreement; the instructor distributes it via a private Dropbox link instead. Download it from:

<https://www.dropbox.com/scl/fi/6t4iqjec1l08eodssbk8d/hrs_lab8.csv?rlkey=7qm0daymlogvart7yl41ppaj2&dl=0>

Or, equivalently, fetch it from R (the `dl=1` suffix forces a direct download). Save it wherever you like — just point the snippet below to that path:

```r
download.file(
  "https://www.dropbox.com/scl/fi/6t4iqjec1l08eodssbk8d/hrs_lab8.csv?rlkey=7qm0daymlogvart7yl41ppaj2&dl=1",
  "hrs_lab8.csv",      # <-- change to wherever you want the file saved
  mode = "wb"
)
```

In the rest of the lab, replace the path in `fread(...)` with the location you chose. Working inside an RStudio Project keeps things tidy: put `hrs_lab8.csv` next to `lab8_gxe.R` (or in a `data/` subfolder) and reference it with a relative path.

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

d <- fread("hrs_lab8.csv")     # <-- adjust to your path
str(d)
summary(d[, .(BMI_AV, birth_year, pgs_bmi)])
table(d$sex)
```

The variables are:

| name | meaning |
|---|---|
| `id` | anonymised subject identifier (`S00001`…) |
| `pgs_bmi` | BMI PGS, **standardised** to mean 0 / SD 1 (PGS3 / GIANT 2015) |
| `BMI_AV` | mean BMI across HRS waves 1–12 |
| `birth_year` | year of birth (1905–1980) |
| `sex` | factor: `male`, `female` |
| `raedyrs` | respondent's years of completed schooling (0–17) |
| `rameduc` / `rafeduc` | mother's / father's years of completed schooling |
| `smoke_last` | 1 if currently smoking at the **last** wave the respondent was surveyed (else 0) |
| `drink_last` | 1 if currently drinks alcohol at the last wave |
| `shlt_last` | self-rated health at the last wave: 1 = excellent, 5 = poor |
| `pc1` … `pc10` | first 10 genetic PCs |

> **Why "last" and not the average?** Smoking, drinking, and self-rated health change over the life course. The last available observation gives you the most recent state before the respondent dropped out / was censored, which is closer to the BMI we observe at later waves than a long-run average. (For a more careful treatment you would model these as time-varying — out of scope here.)

> **Why standardise the PGS?** With $z$-scored PGS, $\beta_G$ is the change in BMI for a one-SD increase in genetic predisposition. This makes effects comparable across PGS, traits, and papers.

> **Why include the PCs?** Population structure produces correlated allele frequencies and outcome differences. Without PCs the PGS coefficient is partly confounded by ancestry. Always include at least 5–10 PCs.

Quick sanity-check plots:

```r
# Birth-year distribution
ggplot(d, aes(birth_year)) + geom_histogram(binwidth = 2, fill = "steelblue") +
  labs(x = "Year of birth", y = "Count") + theme_minimal()

# Distribution of the BMI polygenic index (PGI / PGS)
ggplot(d, aes(pgs_bmi)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  labs(x = "PGS-BMI (z-scored)", y = "Count",
       title = "Polygenic index for BMI",
       subtitle = "After standardisation: mean 0, SD 1") +
  theme_minimal()
```

The PGS histogram should look approximately Gaussian — that's the typical shape of a polygenic score (sum of many small allele effects, central limit theorem).

### 1.1 Scatter and binscatter: BMI vs PGS

Start with a **raw scatter** of every individual:

```r
ggplot(d, aes(pgs_bmi, BMI_AV)) +
  geom_point(alpha = 0.15, size = 0.5, colour = "steelblue") +
  geom_smooth(method = "lm", se = FALSE,
              colour = "firebrick", linewidth = 0.7) +
  labs(x = "PGS-BMI (z-scored)", y = "BMI",
       title = "Raw scatter: BMI by PGS-BMI") +
  theme_minimal()
```

With 8,000+ points the cloud is dense and the underlying signal is hard to read by eye — only the OLS line betrays the positive trend. A **binscatter** shows the same relationship cleanly: split the PGS into equal-width bins, average BMI within each bin, then plot the bin means.

```r
ggplot(d, aes(pgs_bmi, BMI_AV)) +
  stat_summary_bin(fun = mean, bins = 25, geom = "point",
                   colour = "steelblue", size = 2) +
  geom_smooth(method = "lm", se = FALSE,
              colour = "firebrick", linewidth = 0.7) +
  labs(x = "PGS-BMI (z-scored)", y = "Mean BMI in bin",
       title = "Binscatter: BMI by PGS-BMI",
       subtitle = "25 equal-width bins, OLS fit on raw data") +
  theme_minimal()
```

The 25 bin means trace an almost perfectly linear, positive relationship: each bin step up in PGS adds $\approx 1$–$1.5$ BMI units. This is the "main effect" we will quantify in § 2.

---

## 1.5 Inspect the genetic PCs and flag non-European individuals

The HRS subset distributed for this lab was filtered to "European-ancestry" participants upstream — but **upstream filters can leak**: a few individuals with recent admixture or genotyping artefacts often slip through. Always verify with a PC plot before regressing on a PGS, because a PGS trained on EUR predicts much less well outside the cluster and can drag your interaction estimates around.

```r
pcs <- paste0("pc", 1:10)

ggplot(d, aes(pc1, pc2)) +
  geom_point(alpha = 0.4, size = 0.6, colour = "steelblue") +
  labs(x = "PC1", y = "PC2",
       title = "Genetic PCs — HRS lab subset") +
  theme_minimal()
```

You should see one tight blob centred near 0 (the EUR cluster) plus a tail of points sitting away from the centroid.

### 1.5.2 Distribution of PC1 and PC2 — flag extreme values

Look at PC1 and PC2 marginally. Each has a long-tailed distribution: most people sit close to the centre, but a small fraction sit several standard deviations out. Those tail observations are very likely individuals whose genetic ancestry is unusual relative to the bulk of the sample (despite all of them self-reporting as non-Hispanic white).

```r
# Marginal histograms
ggplot(d, aes(pc1)) + geom_histogram(bins = 60, fill = "steelblue") +
  labs(x = "PC1", y = "Count") + theme_minimal()

ggplot(d, aes(pc2)) + geom_histogram(bins = 60, fill = "steelblue") +
  labs(x = "PC2", y = "Count") + theme_minimal()
```

A simple, transparent rule: flag anyone who is **more than 4 standard deviations from the mean on PC1 or PC2**.

```r
d[, pc_outlier := abs(scale(pc1)) > 4 | abs(scale(pc2)) > 4]
table(d$pc_outlier)
```

> **Checkpoint.** With this rule you should flag $\approx 210$ individuals ($\approx 2.5\%$). Tighter (|z|>5) gives almost no one; looser (|z|>3) flags too many. |z|>4 is a reasonable working threshold.

Re-plot PC1 vs PC2, highlighting the flagged points:

```r
ggplot(d, aes(pc1, pc2, colour = pc_outlier)) +
  geom_point(alpha = 0.6, size = 0.7) +
  scale_colour_manual(values = c(`FALSE` = "steelblue", `TRUE` = "red"),
                      labels = c(`FALSE` = "EUR-like", `TRUE` = "PC outlier")) +
  labs(x = "PC1", y = "PC2", colour = NULL,
       title = "Possible non-EUR individuals (|z| > 4 on PC1 or PC2)",
       subtitle = sprintf("%d flagged out of %d (%.1f%%)",
                          sum(d$pc_outlier), nrow(d),
                          100 * mean(d$pc_outlier))) +
  theme_minimal()
```

> **What this is and isn't.** The threshold is purely descriptive: it isolates points sitting far from the EUR centroid in this sample. It does *not* tell you what ancestry those points belong to — for that you'd need to project against a labelled reference panel like 1000 Genomes (Week 7). It also does not "fix" the EUR-trained PGS — the PGS still works less well in those individuals. Treat the flag as a **sensitivity control**: re-run the regressions on `d[pc_outlier == FALSE]` and check that the headline numbers ($\beta_{\text{PGS}} \approx 1.5$, $\beta_{GE} \approx 0.025$) are essentially unchanged.

---

## 2. Main effect: does the PGS predict BMI?

A baseline OLS without interaction:

```r
pcs <- paste0("pc", 1:10)
fmla_main <- as.formula(paste(
  "BMI_AV ~ pgs_bmi + birth_year + sex +",
  paste(pcs, collapse = " + ")
))

m_main <- lm(fmla_main, data = d)
coef(summary(m_main))["pgs_bmi", ]
```

> **Why no age control here?** Age-at-measurement and `birth_year` are mechanically tied in a panel like HRS: people born earlier are older at any given wave. Including both in the same regression produces strong collinearity and the cohort coefficient gets unstable. We keep `birth_year` (the cohort marker) and don't add a separate age control.

> **Checkpoint.** You should find $\beta_{\text{PGS}} \approx 1.5$ BMI units per SD of PGS, with $p < 10^{-100}$. A 1-SD higher PGS is associated with $\approx 1.5$ kg/m² higher BMI on average.

What does this *not* tell us? Whether the PGS effect is the same for everyone — across cohorts, sexes, environments. That is a G$\times$E question.

---

## 3. The interaction model

Walter's hypothesis: the slope of BMI on PGS depends on birth cohort. Centring birth year at 1944 makes the main effect of PGS interpretable as the slope *at* the Walter pivot:

```r
d[, by_c := birth_year - 1944]

fmla_gxe <- as.formula(paste(
  "BMI_AV ~ pgs_bmi * by_c + sex +",
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
              modx.values = c(-30, -10, 10, 30),
              modx.labels = c("1914 cohort", "1934 cohort",
                              "1954 cohort", "1974 cohort"),
              interval = TRUE, int.width = 0.95,
              x.label = "PGS-BMI (SD)", y.label = "Predicted BMI",
              legend.main = "Birth cohort")
```

You should see a fan: lines for older cohorts are *flatter*, lines for younger cohorts are *steeper*. Same PGS, more BMI in younger cohorts.

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
  "BMI_AV ~ pgs_bmi * by_c * sex +",
  paste(pcs, collapse = " + ")
))

m_3way <- lm(fmla_3way, data = d)
summary(m_3way)
anova(m_gxe, m_3way)
```

Interpret the three-way coefficient `pgs_bmi:by_c:sexfemale`. Does the obesogenic-environment effect on the BMI PGS look the same for men and women? Plot the predicted BMI by PGS, faceted by sex and split by birth-year tertile.

---

## 8. G$\times$E with education as the moderator

Birth year is one environmental story (the post-war obesogenic shift). A different class of E variables that the sociogenomics literature has explored heavily is **socio-economic position** — most commonly **education**. The intuition: in a high-education environment, people may have more information, money, and time to push back against an inherited tendency toward weight gain; the PGS slope should then be *attenuated* at higher education levels.

Three specifications. Run them in turn and compare the interaction terms to Walter's $\beta_{GE} \approx 0.025$ from § 3.

### 8.1 Continuous, z-scored years of schooling

```r
d[, edu_z := scale(raedyrs)]    # mean 0, SD 1

fmla_edu <- as.formula(paste(
  "BMI_AV ~ pgs_bmi * edu_z + birth_year + sex +",
  paste(pcs, collapse = " + ")
))
m_edu <- lm(fmla_edu, data = d)
coef(summary(m_edu))[c("pgs_bmi", "edu_z", "pgs_bmi:edu_z"), ]
```

> **Checkpoint.** Strong **main** effect of education on BMI ($\beta_{\text{edu}} \approx -0.43$, $p < 10^{-13}$): a one-SD higher education is associated with $\approx 0.4$ kg/m² lower BMI. But the **interaction** $\beta_{\text{PGS:edu}} \approx -0.04$, $p \approx 0.4$ — not distinguishable from zero. Education shifts the *level* of BMI but does not visibly moderate the PGS slope.

### 8.2 Three-level education (handle non-linearity)

Maybe the relationship isn't linear in years — perhaps only college-or-more matters. Recode:

```r
d[, edu3 := cut(raedyrs,
                breaks = c(-Inf, 12, 15, Inf),
                labels = c("HS_or_less", "some_coll", "coll_plus"))]
table(d$edu3)

fmla_e3 <- as.formula(paste(
  "BMI_AV ~ pgs_bmi * edu3 + birth_year + sex +",
  paste(pcs, collapse = " + ")
))
m_e3 <- lm(fmla_e3, data = d)
summary(m_e3)$coef[grep("pgs_bmi|edu3", rownames(summary(m_e3)$coef)), ]
```

> **Checkpoint.** College-or-more is associated with $\approx 1.1$ kg/m² lower BMI than HS-or-less, holding everything else fixed. But neither of the two `pgs_bmi:edu3*` interaction terms is significant ($p > 0.5$). Even the non-linear coding gives a null G$\times$Education.

### 8.3 Parental education (childhood SES)

The respondent's *own* years of schooling is shaped by childhood SES, by the same genes that drive the PGS, and by their adult choices. **Parental education** is a cleaner proxy for the *childhood* environment — closer to a true E variable in the rGE sense, because it cannot be a function of the respondent's own genotype.

```r
d_par <- d[!is.na(rameduc) & !is.na(rafeduc)]
d_par[, par_edu := scale((rameduc + rafeduc) / 2)]

fmla_pe <- as.formula(paste(
  "BMI_AV ~ pgs_bmi * par_edu + birth_year + sex +",
  paste(pcs, collapse = " + ")
))
m_pe <- lm(fmla_pe, data = d_par)
coef(summary(m_pe))[c("pgs_bmi", "par_edu", "pgs_bmi:par_edu"), ]
```

> **Checkpoint.** $n \approx 7{,}300$ after dropping rows missing parental education. The main effect is again negative ($\beta_{\text{par\_edu}} \approx -0.40$). The interaction $\beta_{\text{PGS:par\_edu}} \approx +0.07$, $p \approx 0.2$ — *positive* (i.e. PGS slope slightly steeper for kids of better-educated parents) and still not significant.

### 8.4 Visualise (specification 8.1) and discuss

```r
interact_plot(m_edu, pred = "pgs_bmi", modx = "edu_z",
              modx.values = c(-1, 0, 1),
              interval = TRUE, int.width = 0.95,
              x.label = "PGS-BMI (SD)", y.label = "Predicted BMI",
              legend.main = "Education (SD)")
```

You should see three nearly-parallel lines, just shifted vertically. The flat-fan pattern is what a non-significant G$\times$E looks like — contrast it with the wide fan you produced from `m_gxe` in § 4.1.

### 8.5 What does the contrast Walter (cohort) vs Education tell us?

A non-result is also a result. With the *same* PGS, in the *same* sample:

* **Birth year** moderates the PGS slope strongly and significantly ($\beta_{GE} \approx 0.025$, $p \approx 10^{-7}$).
* **Years of education**, **3-level education**, and **parental education** do *not* — every interaction is small and non-significant.

Interpretation: the **secular obesogenic shift** (cheap calories, sedentary jobs, motorisation) appears to "unlock" genetic predisposition for everyone in roughly the same way, regardless of where they sit in the education distribution. That contrasts with **education-attainment PGS**, where prior literature (Domingue et al., Tropf et al.) finds that environmental opportunities *do* gate genetic potential. The mechanism for BMI — exposure to ubiquitous environmental change — is harder to escape via SES than the mechanism for educational attainment.

> **Caveat on causal status.** Unlike `birth_year`, *neither* `raedyrs` nor `par_edu` is exogenous — they are themselves shaped by genes and environment. A G$\times$E coefficient on education is therefore *descriptive* (does the PGS slope differ across education levels?) and not *causal* (does education unlock the PGS?). Even if the interaction were significant, identifying it as causal would require an instrument (e.g. compulsory-schooling reforms, as in Davies et al. 2018).

---

## References

* Walter, S., Mejía-Guevara, I., Estrada, K., Liu, S. Y., & Glymour, M. M. (2016). *Association of a Genetic Risk Score with Body Mass Index Across Different Birth Cohorts.* **JAMA**, 316(1), 63–69.
* Herd, P., Freese, J., Sicinski, K., Domingue, B. W., Mullan Harris, K., Wei, C., & Hauser, R. M. (2019). *Genes, Gender Inequality, and Educational Attainment.* **American Sociological Review**, 84(6), 1069–1098.
* Domingue, B. W., Conley, D., Fletcher, J., & Boardman, J. D. (2016). *Cohort effects in the genetic influence on smoking.* **Behavior Genetics**, 46(1), 31–42.
* Tropf, F. C. (2019). *Applying Polygenic Scores: G$\times$E applications.* In Conley & Fletcher (eds), *The Genome Factor / textbook companion materials* — chapter and code on which this lab is built.
* Health and Retirement Study (HRS). Public-use dataset. Polygenic Score release 3 (PGS3), Ware et al.\ 2018.
