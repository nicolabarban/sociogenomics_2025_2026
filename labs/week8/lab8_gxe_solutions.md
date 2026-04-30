# Lab 8 (part B) — Solutions

These solutions report the numerical results expected from the anonymised HRS subset distributed with the lab ($N = 8{,}448$, birth years 1905–1980, BMI PGS standardised). Your numbers should match to within rounding.

---

## 2. Main effect

```
Estimate (pgs_bmi)  : ~ 1.50
Std. Error          : ~ 0.058
t value             : ~ 26
p value             : < 1e-100
```

A 1-SD higher PGS-BMI is associated with $\approx 1.5$ kg/m² higher mean BMI, controlling for sex, birth year, and 10 PCs. (We drop `Age_AV` here: it is mechanically tied to `birth_year` in HRS and including both produces collinearity that distorts the cohort coefficient.)

The $R^2$ of `m_main` is around $0.12$. PGS alone (with controls) explains $\approx 12\%$ of the variance in BMI in this older sample — typical for current-generation BMI PGS.

---

## 3. Interaction model

With `by_c = birth_year - 1944`:

```
                  Estimate   Std. Error   t value   p value
pgs_bmi             1.606      0.061       26.5    < 1e-100
by_c                0.097      0.005       21.1    < 1e-90
pgs_bmi:by_c        0.0250     0.0045       5.6    ~ 3e-08
```

* For someone born in 1944, the PGS slope is $\approx 1.60$.
* Each additional year of birth raises that slope by $\approx 0.025$ BMI units.
* Across the 75-year span 1905→1980, the implied total shift in slope is about $0.025 \times 75 \approx 1.9$ — i.e., the PGS slope nearly *doubles* from oldest to youngest cohort.
* The main effect of `by_c` ($\approx 0.097$) absorbs the secular obesity trend now that `Age_AV` is gone: BMI rises by about 1 kg/m² per decade of later birth.

`anova(m_main, m_gxe)`: $F \approx 31$, $p \approx 3 \times 10^{-8}$. The interaction is a real signal.

$\Delta R^2$ from adding the interaction: $\approx 0.0025$. **Tiny in variance, large in mechanism** — a classic feature of G$\times$E. The interaction barely moves $R^2$ but it changes who the coefficient applies to.

---

## 5. Walter design

```
group     beta    se     n
pre1944   1.20   0.066   5,262
post1944  1.98   0.108   3,186
```

* The PGS slope is **0.78 BMI units larger** in the post-1944 cohort.
* In relative terms, the PGS effect is **65% larger** in the post-war cohort than in the pre-war one.
* Both 95% CIs are well above zero and do not overlap each other (pre-1944 CI ≈ [1.07, 1.33]; post-1944 CI ≈ [1.77, 2.19]).

This is the central Walter et al.\ (2016) finding, replicated almost exactly with the same data and a slightly cleaner specification.

---

## 6. Discussion (one possible reading)

**1. rGE / selection.** Birth year is exogenous to the individual's genotype (you cannot choose when you were born). So *strict* rGE — genes selecting into environments — is not the worry here. The serious threat is **survival selection**: the pre-1944 group in HRS is observed only because they survived to enrolment (typically age 50+). If high-BMI-PGS individuals in earlier cohorts died at higher rates, the surviving pre-1944 sample is enriched for low-PGS people, and the smaller pre-1944 slope partly reflects mortality selection rather than a weaker biological effect of the PGS at the time. Sensitivity analysis: estimate the same model in waves where both cohorts are observable at the same age, and compare.

**2. Population vs individual.** No. The 1.98-per-SD slope is a *population statistic*: across people who differ in PGS by 1 SD, expected BMI differs by 1.98. It does **not** decompose any one person's BMI into "65% genetic". For a specific individual, genes and environment are non-separable inputs to a developmental process — both are 100% necessary. This is the Lewontin distinction at work (and the "two kinds of why" from class).

**3. Mechanism.** Plausible candidates between 1920 and 1980:
   * Cheap energy-dense food and the post-war shift in diet composition (refined carbs, sugars, vegetable oils);
   * Decline in physical activity tied to mechanisation, motorisation, sedentary work, suburbanisation.

   To distinguish, look for *exogenous variation* in one channel: for example, exploit the timing of price drops in particular food groups, or the introduction of fast-food chains in different US states. Then test whether $\beta_{GE}$ is larger in places/times where that channel changed most. The cleanest existing example uses fast-food density as $E$.

**4. Portability.** Two predictions, working in opposite directions:
   * The PGS itself was *trained* on a European-ancestry GWAS (GIANT 2015), so it predicts BMI **less well** in East Asian samples — a portability/LD-mismatch problem. So $\beta_{\text{PGS}}$ would be *attenuated* mechanically, before any G$\times$E story.
   * On top of that, the obesogenic-environment shift in East Asia happened later and faster (1980s onward in mainland China). So the *cohort* at which the PGS slope "switches on" should be displaced rightward by 30–40 years.
   You would not expect the same numbers; you would expect the same *shape*, shifted in time, and weaker in absolute terms.

---

## 7. Optional: 3-way (PGS × cohort × sex)

With `sex` as factor (reference category = `female`, so the contrast variable is `sexmale`):

```
                          Estimate   S.E.    t     p
pgs_bmi                     1.78    0.080  22.2  <1e-100
by_c                        0.076   0.019   4.1   ~5e-05
sexmale                     0.33    0.116   2.8   ~0.005
pgs_bmi:by_c                0.034   0.006   5.9   ~4e-09
pgs_bmi:sexmale            -0.38    0.116  -3.3   ~0.001
by_c:sexmale               -0.028   0.009  -3.0   ~0.003
pgs_bmi:by_c:sexmale       -0.020   0.009  -2.2   ~0.026
```

`anova(m_gxe, m_3way)`: $F \approx 7.2$, $p \approx 8 \times 10^{-5}$. The 3-way model fits significantly better.

**Reading the 3-way coefficient.** `pgs_bmi:by_c:sexmale` $\approx -0.020$ means the cohort-by-PGS amplification (the core Walter effect) is **smaller in men than in women**. For women, the BMI-PGS slope grows by $\approx 0.034$ per year of birth; for men it grows by only $\approx 0.034 - 0.020 = 0.014$. The post-war "obesogenic switch" did more to amplify women's genetic risk for high BMI than men's.

That is the opposite of the Herd et al.\ (2019) finding for **education**, where genetic potential was historically *more* gated for women and the gap closed over time. The contrast is instructive: each trait's G$\times$E story depends on which structural constraint is being relaxed, and for whom. There is no single "cohort unlocks genes" pattern — it has to be argued case by case.

---

## A note for the instructor

* The teaching subset is anonymised (HRS IDs dropped, PGS standardised, synthetic IDs assigned). It is shared with enrolled students under the HRS DUA via a private link — do **not** push it to the GitHub student repo.
* If you re-prepare the subset, regenerate it with `lab_prep/prepare_hrs_subset.R` (uses Felix Tropf's `hrs.csv` as the source).
* All the numbers above were computed against the subset built on 2026-04-27 with `set.seed(20260427)`.
