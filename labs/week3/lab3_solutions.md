# Lab 3 — Solutions

---

## Exercise 1 (Linkage Disequilibrium and Population Comparison)

### 1. What is the r-squared between rs994335 and rs2379903 in the full sample? How does it compare across EUR, AFR, and EAS?

**Full sample:**

```bash
plink --bfile hapmap3_qc \
      --ld rs994335 rs2379903 \
      --out ld_pair
grep "R-sq" ld_pair.log
```

Report the R-sq value from the log output.

**Per population:**

```bash
plink --bfile hapmap3_qc --keep samples_EUR.txt \
      --ld rs994335 rs2379903 --out ld_EUR
grep "R-sq" ld_EUR.log

plink --bfile hapmap3_qc --keep samples_AFR.txt \
      --ld rs994335 rs2379903 --out ld_AFR
grep "R-sq" ld_AFR.log

plink --bfile hapmap3_qc --keep samples_EAS.txt \
      --ld rs994335 rs2379903 --out ld_EAS
grep "R-sq" ld_EAS.log
```

**Expected:** The r-squared values will differ across populations. AFR typically shows lower LD than EUR or EAS for the same SNP pair, because African populations have larger effective population sizes and more historical recombination events that break down LD.

### 2. On chromosome 22, find the 5 SNP pairs with the highest r-squared in EUR. Are the same pairs in high LD in AFR?

```bash
# Top 5 pairs in EUR
sort -k7 -n -r ld_chr22_EUR.ld | head -5
```

Note the SNP names from columns 3 and 6. Then check if those same pairs appear in the AFR output:

```bash
# Example: check a specific pair (replace SNP names with your results)
sort -k7 -n -r ld_chr22_AFR.ld | head -5
```

Alternatively, extract the top EUR pairs and look them up in AFR:

```bash
# Get the top 5 EUR pairs as SNP_A-SNP_B
sort -k7 -n -r ld_chr22_EUR.ld | head -5 | awk '{print $3, $6}' > top_eur_pairs.txt

# Search for each pair in the AFR results
while read snpA snpB; do
    grep -w "$snpA" ld_chr22_AFR.ld | grep -w "$snpB" || echo "$snpA $snpB: not in AFR output (r2 < 0.05)"
done < top_eur_pairs.txt
```

**Expected:** Some SNP pairs that are in high LD in EUR may have much lower LD in AFR, or may not even appear in the AFR output (because their r-squared falls below the 0.05 threshold). This reflects the shorter LD blocks in African populations.

### 3. Does LD increase or decrease with distance? Which population shows the fastest LD decay and why?

LD **decreases** with physical distance. This is expected because recombination is more likely between distant SNPs, breaking down allelic associations over generations.

**AFR shows the fastest LD decay** (lowest mean r-squared at every distance bin). This is because:

- African populations have the **largest effective population size** among the three groups
- They have the **oldest population history**, meaning more generations of recombination have occurred
- Out-of-Africa populations (EUR, EAS) experienced **bottlenecks** during migration, which reduced diversity and created longer LD blocks
- The smaller effective population size after the bottleneck means fewer recombination events broke down LD

To verify this, compare the mean r-squared in the first distance bin (0-10 kb):

```bash
for POP in EUR AFR EAS; do
    echo -n "$POP 0-10kb mean r2: "
    awk 'NR>1 {
        dist = ($5 - $2); if (dist < 0) dist = -dist
        if (dist < 10000) { count++; sumr2 += $7 }
    } END { printf "%.4f (n=%d)\n", sumr2/count, count }' ld_chr22_${POP}.ld
done
```

---

## Exercise 2 (Detecting Related Individuals)

### 1. How many pairs have PI_HAT > 0.20? Is the top pair parent-offspring or siblings?

```bash
# Count pairs with PI_HAT > 0.20
awk 'NR>1 && $10 > 0.20' ibd_results.genome | wc -l

# Show the pair with highest PI_HAT, including Z values
awk 'NR>1 && $10 > 0.20' ibd_results.genome | sort -k10 -n -r | head -1
```

To determine the relationship type, look at the Z columns:

- If **Z0 ≈ 0, Z1 ≈ 1, Z2 ≈ 0** → **parent-offspring** (they always share exactly 1 allele IBD)
- If **Z0 ≈ 0.25, Z1 ≈ 0.50, Z2 ≈ 0.25** → **full siblings** (they can share 0, 1, or 2 alleles IBD)

Print Z values clearly:

```bash
awk 'NR>1 && $10 > 0.40 {print $1, $2, $3, $4, "Z0="$7, "Z1="$8, "Z2="$9, "PI_HAT="$10}' ibd_results.genome | sort -k8 -n -r
```

**Interpretation:** If Z1 is close to 1.0 and Z0 and Z2 are close to 0, the pair is parent-offspring. If Z0, Z1, Z2 are roughly 0.25, 0.50, 0.25, the pair is full siblings.

### 2. Using the KING output, find all pairs with kinship > 0.177. How many are there?

```bash
# From the between-family output
awk 'NR>1 && $NF > 0.177' king_results.kin0 | wc -l

# From the within-family output
awk 'NR>1 && $NF > 0.177' king_results.kin | wc -l
```

The total number of first-degree relative pairs is the sum from both files.

Note that KING kinship > 0.177 corresponds to approximately PI_HAT > 0.354, which captures parent-offspring and full sibling relationships.

### 3. If you identify two first-degree relatives, which one would you remove and why?

You should remove the individual with the **higher rate of missing genotypes** (higher `F_MISS` in the `.imiss` file). This is because:

- The individual with more complete data contributes more information to downstream analyses (association tests, PCA, etc.)
- Removing the lower-quality sample minimises data loss
- You can check this with:

```bash
# Look up the missing rate for each individual in the pair
grep "IID1" qc_missing.imiss
grep "IID2" qc_missing.imiss
```

Replace `IID1` and `IID2` with the actual individual IDs from the related pair. Remove the one with the larger `F_MISS` value.

If both individuals have similar missingness, other criteria can be used: remove the one with more extreme heterozygosity (F coefficient), or the one that is related to more other individuals in the sample.
