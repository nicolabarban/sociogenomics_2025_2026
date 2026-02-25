# Lab 3. Quality Control, Relatedness, and Linkage Disequilibrium

In this lab you will perform a full quality control (QC) pipeline on genetic data, investigate patterns of linkage disequilibrium (LD), and detect related individuals using PLINK and KING. By the end you will be able to run a production-quality QC workflow, compute and visualise LD between SNPs, and identify cryptic relatedness.

We continue using the **HapMap Phase III** dataset (`hapmap3`) from the course repository.

---

## 0. Getting started

Open [Google Cloud Shell](https://shell.cloud.google.com/) in your browser.

Update the course repository and pull the data files:

```bash
cd ~/sociogenomics_2025_2026
git pull
git lfs pull
cd $HOME
```

Make sure your project directories exist:

```bash
mkdir -p ~/Sociogenomics/Data ~/Sociogenomics/Results ~/Sociogenomics/Scripts
```

Copy the HapMap data (if not already present from Lab 2):

```bash
cp ~/sociogenomics_2025_2026/data/hapmap3.bed ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/hapmap3.bim ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/hapmap3.fam ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/BMI_pheno.txt ~/Sociogenomics/Data/
cp ~/sociogenomics_2025_2026/data/1kg_samples.txt ~/Sociogenomics/Data/
```

Move to your working data directory:

```bash
cd ~/Sociogenomics/Data
```

Make sure PLINK is available:

```bash
plink --version
```

If not, reinstall it:

```bash
bash ~/sociogenomics_2025_2026/scripts/setup_plink19.sh
source ~/.bashrc
```

---

## Part I. Quality Control Pipeline

A standard GWAS quality control pipeline removes low-quality SNPs and individuals before any analysis. The steps below follow best-practice recommendations (Anderson et al. 2010; Marees et al. 2018).

The typical QC order is:

1. **Per-individual missing rate** — remove samples with too much missing data
2. **Per-SNP missing rate** — remove SNPs with too much missing data
3. **Minor allele frequency (MAF)** — remove very rare variants
4. **Hardy-Weinberg equilibrium (HWE)** — remove SNPs with genotyping errors
5. **LD pruning** — reduce redundancy before relatedness and heterozygosity checks
6. **Heterozygosity / inbreeding** — remove outlier samples
7. **Relatedness** — identify and remove duplicates or close relatives

### Step 1 — Compute missing rates

```bash
plink --bfile hapmap3 \
      --missing \
      --out qc_missing
```

Inspect the output:

```bash
head qc_missing.imiss    # per-individual missingness
head qc_missing.lmiss    # per-SNP missingness
```

The key column is `F_MISS` — the fraction of missing genotypes.

Sort individuals by missingness (worst first):

```bash
sort -k6 -n -r qc_missing.imiss | head -20
```

Sort SNPs by missingness:

```bash
sort -k5 -n -r qc_missing.lmiss | head -20
```

### Step 2 — Compute allele frequencies

```bash
plink --bfile hapmap3 \
      --freq \
      --out qc_freq
```

Count how many SNPs fall in each frequency bin:

```bash
awk 'NR>1 {print $5}' qc_freq.frq | \
  awk '{
    if ($1 < 0.01) rare++
    else if ($1 < 0.05) lowfreq++
    else common++
  } END {
    print "Rare (<0.01):", rare
    print "Low-freq (0.01-0.05):", lowfreq
    print "Common (>=0.05):", common
  }'
```

### Step 3 — Compute Hardy-Weinberg statistics

```bash
plink --bfile hapmap3 \
      --hardy \
      --out qc_hwe
```

How many SNPs significantly deviate from HWE (p < 1e-6)?

```bash
awk 'NR>1 && $9 < 1e-6' qc_hwe.hwe | wc -l
```

### Step 4 — Apply all filters in one command

Now apply the standard thresholds to produce a clean dataset:

```bash
plink --bfile hapmap3 \
      --mind 0.05 \
      --geno 0.02 \
      --maf  0.01 \
      --hwe  1e-6 \
      --make-bed \
      --out hapmap3_qc
```

Check the PLINK log to see how many SNPs and individuals were removed at each step:

```bash
cat hapmap3_qc.log | grep -E "removed|remaining|pass"
```

**Question:** How many SNPs and individuals passed QC?

### Step 5 — LD pruning

Many downstream analyses (relatedness, PCA, heterozygosity) work better on a set of approximately independent SNPs. LD pruning removes one SNP from each pair with $r^2 > 0.2$.

```bash
plink --bfile hapmap3_qc \
      --indep-pairwise 50 5 0.2 \
      --out hapmap3_pruned
```

This creates two files:
- `hapmap3_pruned.prune.in` — SNPs to keep
- `hapmap3_pruned.prune.out` — SNPs removed (in high LD)

How many SNPs remain after pruning?

```bash
wc -l hapmap3_pruned.prune.in
wc -l hapmap3_pruned.prune.out
```

Extract the pruned SNP set:

```bash
plink --bfile hapmap3_qc \
      --extract hapmap3_pruned.prune.in \
      --make-bed \
      --out hapmap3_pruned_set
```

### Step 6 — Heterozygosity and inbreeding check

Individuals with unusually high or low heterozygosity may indicate sample contamination (high) or inbreeding / population stratification (low).

Compute the inbreeding coefficient $F$ on the pruned SNP set:

```bash
plink --bfile hapmap3_pruned_set \
      --het \
      --out qc_het
```

Inspect the output:

```bash
head qc_het.het
```

The columns are: FID, IID, O(HOM) observed homozygotes, E(HOM) expected homozygotes, N(NM) non-missing SNPs, F inbreeding coefficient.

Find outliers — flag individuals with $F < -0.15$ or $F > 0.15$:

```bash
awk 'NR>1 && ($6 < -0.15 || $6 > 0.15) {print $1, $2, $6}' qc_het.het
```

---

## Part II. Linkage Disequilibrium

**Linkage disequilibrium (LD)** is the non-random association of alleles at different loci. It arises from:
- **Physical linkage** (nearby SNPs tend to be co-inherited)
- **Population history** (founder effects, bottlenecks, admixture)
- **Selection** (variants near selected loci have inflated LD)

The most common measure is $r^2$ — the squared Pearson correlation between allele counts at two SNPs.

### Pairwise LD between specific SNPs

Calculate $r^2$ between two SNPs on chromosome 22:

```bash
plink --bfile hapmap3_qc \
      --ld rs4819391 rs9605927 \
      --out ld_pair
cat ld_pair.log | grep -A5 "LD"
```

### Compute LD across a chromosome

Compute pairwise $r^2$ for all SNP pairs on chromosome 22:

```bash
plink --bfile hapmap3_qc \
      --chr 22 \
      --r2 \
      --ld-window 100 \
      --ld-window-kb 1000 \
      --ld-window-r2 0.05 \
      --out ld_chr22
```

Options used:
- `--ld-window 100` — only consider pairs within 100 SNPs of each other
- `--ld-window-kb 1000` — only consider pairs within 1000 kb
- `--ld-window-r2 0.05` — only report pairs with $r^2 > 0.05$

Inspect the output:

```bash
head ld_chr22.ld
wc -l ld_chr22.ld
```

The columns are: CHR\_A, BP\_A, SNP\_A, CHR\_B, BP\_B, SNP\_B, R2.

Find the pair of SNPs with the highest $r^2$:

```bash
sort -k7 -n -r ld_chr22.ld | head -5
```

### LD as a function of physical distance

LD generally decreases with physical distance between SNPs.
To examine this, compute the mean $r^2$ in 10 kb distance bins:

```bash
awk 'NR>1 {
    dist = ($5 - $2)
    if (dist < 0) dist = -dist
    bin = int(dist/10000)  # 10kb bins
    count[bin]++
    sumr2[bin] += $7
} END {
    for (b in count) printf "%d\t%.4f\n", b*10, sumr2[b]/count[b]
}' ld_chr22.ld | sort -k1 -n | head -30
```

This shows mean $r^2$ by distance (in kb). You should see a clear **LD decay** — $r^2$ decreasing as SNPs get further apart.

### Compare LD across populations

LD patterns differ across populations due to different demographic histories (bottlenecks, admixture, effective population size). African populations typically have **shorter LD blocks** because of their larger effective population size and older history, while out-of-Africa populations (European, East Asian) tend to have **longer LD blocks** due to founder effects.

The course repository includes a sample information file from the **1000 Genomes Project** that maps individuals to populations and superpopulations.

First, copy it to your working directory:

```bash
cp ~/sociogenomics_2025_2026/data/1kg_samples.txt ~/Sociogenomics/Data/
```

Explore the file:

```bash
head -5 1kg_samples.txt
```

The columns are: Sample name, Sex, Biosample ID, Population code, Population name, Superpopulation code, Superpopulation name, Data collections.

Count individuals per superpopulation:

```bash
awk -F'\t' 'NR>1 {print $6, $7}' 1kg_samples.txt | sort | uniq -c | sort -n -r
```

You should see five groups: **AFR** (African), **EUR** (European), **EAS** (East Asian), **SAS** (South Asian), **AMR** (American/Admixed).

### Create population-specific sample lists

We need to create PLINK-compatible sample lists (two columns: FID IID) for each superpopulation. We match the sample names from `1kg_samples.txt` against the individuals in our dataset.

Create sample lists for three superpopulations:

```bash
awk -F'\t' 'NR>1 && $6 == "EUR" {print $1, $1}' 1kg_samples.txt > samples_EUR.txt
awk -F'\t' 'NR>1 && $6 == "AFR" {print $1, $1}' 1kg_samples.txt > samples_AFR.txt
awk -F'\t' 'NR>1 && $6 == "EAS" {print $1, $1}' 1kg_samples.txt > samples_EAS.txt
```

Check how many individuals from each group are in our dataset:

```bash
plink --bfile hapmap3_qc --keep samples_EUR.txt --make-just-fam --out check_EUR
wc -l check_EUR.fam

plink --bfile hapmap3_qc --keep samples_AFR.txt --make-just-fam --out check_AFR
wc -l check_AFR.fam

plink --bfile hapmap3_qc --keep samples_EAS.txt --make-just-fam --out check_EAS
wc -l check_EAS.fam
```

### Pairwise LD in different populations

Compare the $r^2$ between `rs4819391` and `rs9605927` across populations:

```bash
plink --bfile hapmap3_qc --keep samples_EUR.txt \
      --ld rs4819391 rs9605927 --out ld_EUR
grep "R-sq" ld_EUR.log

plink --bfile hapmap3_qc --keep samples_AFR.txt \
      --ld rs4819391 rs9605927 --out ld_AFR
grep "R-sq" ld_AFR.log

plink --bfile hapmap3_qc --keep samples_EAS.txt \
      --ld rs4819391 rs9605927 --out ld_EAS
grep "R-sq" ld_EAS.log
```

**Question:** Is the $r^2$ the same across all three populations? Which population shows the lowest LD?

### LD decay across populations

Compute LD on chromosome 22 separately for each population and compare the decay patterns:

```bash
for POP in EUR AFR EAS; do
    plink --bfile hapmap3_qc \
          --keep samples_${POP}.txt \
          --chr 22 \
          --r2 \
          --ld-window 100 \
          --ld-window-kb 1000 \
          --ld-window-r2 0.05 \
          --out ld_chr22_${POP}
done
```

Compute mean $r^2$ by distance bin for each population:

```bash
for POP in EUR AFR EAS; do
    echo "=== $POP ==="
    awk 'NR>1 {
        dist = ($5 - $2)
        if (dist < 0) dist = -dist
        bin = int(dist/10000)
        count[bin]++
        sumr2[bin] += $7
    } END {
        for (b in count) printf "%d kb\t%.4f\n", b*10, sumr2[b]/count[b]
    }' ld_chr22_${POP}.ld | sort -k1 -n | head -15
    echo ""
done
```

You should observe that:
- **AFR** has the fastest LD decay (lowest $r^2$ at every distance) — reflecting the large ancestral effective population size
- **EUR** and **EAS** have slower decay (higher $r^2$) — reflecting bottlenecks during the out-of-Africa migration
- All populations show the same general pattern: LD decreases with distance

### Exercise 1

1. What is the $r^2$ value between `rs4819391` and `rs9605927` in the full sample? How does it compare in EUR, AFR, and EAS?
2. On chromosome 22, find the 5 SNP pairs with the highest $r^2$ in the EUR population. Are the same pairs also in high LD in AFR?
3. Does LD generally increase or decrease with physical distance between SNPs? Which population shows the fastest LD decay and why?

---

## Part III. Detecting Related Individuals

Cryptic relatedness can inflate association statistics and bias heritability estimates. We use two approaches: PLINK's IBD estimation and the KING program.

### IBD estimation with PLINK (`--genome`)

PLINK estimates identity-by-descent (IBD) sharing between all pairs of individuals. The key output metric is **PI\_HAT** = proportion of genome shared IBD.

| PI\_HAT | Relationship |
|---------|-------------|
| ~1.00   | MZ twin / duplicate |
| ~0.50   | Parent-offspring or full siblings |
| ~0.25   | Half-siblings, avuncular, grandparent |
| ~0.125  | First cousins |
| < 0.10  | Effectively unrelated |

Run on the pruned dataset:

```bash
plink --bfile hapmap3_pruned_set \
      --genome \
      --out ibd_results
```

Inspect the output:

```bash
head ibd_results.genome
```

The key columns are: IID1, IID2, Z0 (prob. 0 alleles IBD), Z1 (prob. 1 allele IBD), Z2 (prob. 2 alleles IBD), PI\_HAT.

Find all pairs with PI\_HAT > 0.20 (second-degree relatives or closer):

```bash
awk 'NR>1 && $10 > 0.20 {print $1, $2, $3, $4, $10}' ibd_results.genome | sort -k5 -n -r
```

Count them:

```bash
awk 'NR>1 && $10 > 0.20' ibd_results.genome | wc -l
```

Find likely parent-offspring or full sibling pairs (PI\_HAT between 0.40 and 0.60):

```bash
awk 'NR>1 && $10 > 0.40 && $10 < 0.60 {print $1, $2, $3, $4, $10}' ibd_results.genome
```

You can distinguish parent-offspring from full siblings using the Z columns:
- **Parent-offspring:** Z0 ≈ 0, Z1 ≈ 1, Z2 ≈ 0 (always share exactly 1 allele IBD)
- **Full siblings:** Z0 ≈ 0.25, Z1 ≈ 0.50, Z2 ≈ 0.25 (share 0, 1, or 2 alleles)

Print the Z values for the closest pairs:

```bash
awk 'NR>1 && $10 > 0.40 {print $1, $2, $3, $4, "Z0="$7, "Z1="$8, "Z2="$9, "PI_HAT="$10}' ibd_results.genome | sort -k8 -n -r
```

**Note:** In a QC pipeline, you would remove one individual from each related pair. Typically, you keep the individual with lower missingness.

You can also run IBD with a minimum PI\_HAT threshold to speed up large datasets:

```bash
plink --bfile hapmap3_pruned_set \
      --genome \
      --min 0.125 \
      --out ibd_related_only
```

### KING — more robust kinship estimation

PLINK's `--genome` assumes homogeneous population structure. In samples with ancestry differences, PI\_HAT can be inflated. **KING** (Kinship-based INference for Gwas; Manichaikul et al. 2010) estimates kinship coefficients using a method robust to population structure.

> **Important:** KING requires genome-wide SNPs — do **not** use the LD-pruned set here. Use `hapmap3_qc` instead.

Install KING on Google Cloud Shell:

```bash
cd ~/Sociogenomics
wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar -xzf Linux-king.tar.gz
chmod +x king
```

Run KING kinship estimation (use the QC-filtered dataset, not the LD-pruned one):

```bash
./king -b ~/Sociogenomics/Data/hapmap3_qc.bed \
       --kinship \
       --prefix king_results
```

For large datasets, use `--related` instead, which is faster and focuses on close relatives:

```bash
./king -b ~/Sociogenomics/Data/hapmap3_qc.bed \
       --related \
       --degree 2 \
       --prefix king_related
```

The `--degree 2` flag reports up to 2nd-degree relatives (kinship > 0.0884).

Inspect the between-family output (pairs from different family IDs):

```bash
head king_results.kin0
```

The within-family output (pairs sharing a family ID in the FAM file):

```bash
head king_results.kin
```

Key columns in both files:
- `N_SNP` — number of SNPs used
- `HetHet` — proportion of SNPs where both individuals are heterozygous
- `IBS0` — proportion of SNPs with 0 alleles shared identical by state
- `Kinship` — estimated kinship coefficient

### Understanding KING kinship coefficients

The **kinship coefficient** measures the probability that two alleles sampled at random from two individuals are identical by descent. It is related to PI\_HAT but on a different scale:

$$\text{Kinship} \approx \frac{\text{PI\_HAT}}{2}$$

For example, parent-offspring pairs share half their genome IBD (PI\_HAT ≈ 0.50), so their kinship coefficient is ≈ 0.25.

| Kinship | Relationship | Equivalent PI\_HAT |
|---------|-------------|-------------------|
| > 0.354 | Duplicate / MZ twin | > 0.708 |
| 0.177 – 0.354 | 1st degree (parent-offspring, full siblings) | 0.354 – 0.708 |
| 0.0884 – 0.177 | 2nd degree (half-sibs, grandparent, avuncular) | 0.177 – 0.354 |
| 0.0442 – 0.0884 | 3rd degree (first cousins) | 0.0884 – 0.177 |
| < 0.0442 | Unrelated | < 0.0884 |

**Negative kinship values** are expected in KING output. They occur when two individuals come from different populations and share fewer alleles than average — this is not an error.

Find all related pairs (2nd degree or closer):

```bash
awk 'NR>1 && $NF > 0.0884 {print}' king_results.kin0 | sort -k8 -n -r | head -20
```

**Why use KING instead of PLINK?** PLINK's `--genome` estimates allele frequencies from the full sample. If the sample includes multiple ancestry groups, the frequency estimates are wrong for each group, which inflates PI\_HAT between unrelated individuals from the same population. KING avoids this by using a method that does not depend on allele frequency estimates — it compares genotypes directly between each pair.

For more details on KING: [https://www.kingrelatedness.com/manual.shtml](https://www.kingrelatedness.com/manual.shtml)

### Exercise 2

1. How many pairs of individuals in HapMap3 have a PI\_HAT > 0.20? Look at the Z0, Z1, Z2 columns of the pair with the highest PI\_HAT — is it a parent-offspring pair or a sibling pair?
2. Using the KING output, find all pairs with kinship > 0.177 (first-degree relatives). How many are there?
3. In a GWAS study design, if you identify two first-degree relatives in your sample, which one would you remove and why?

---

## Part IV. External Resources for LD Lookup

You do not always need to compute LD yourself. Several web tools let you look up LD for specific SNPs in reference populations:

### LDlink (NIH)

**[https://ldlink.nci.nih.gov/](https://ldlink.nci.nih.gov/)**

LDlink is a web-based tool backed by the 1000 Genomes Project data. It offers:

- **LDpair** — $r^2$ and $D'$ between two SNPs in a chosen population
- **LDproxy** — all SNPs in LD with a query SNP (useful after a GWAS hit)
- **LDpop** — compare LD across populations
- **LDmatrix** — pairwise LD matrix for a set of SNPs
- **SNPchip** — check if a SNP is on common genotyping arrays
- **RegulomeDB integration** — functional annotation of LD proxies

**Try it:** Look up `rs4819391` and `rs9605927` in the **JPT+CHB** population and compare with the PLINK result you computed above.

---

## Summary: Full QC Pipeline

Here is the complete recommended pipeline in a single block for reference:

```bash
# Step 1. Remove low-quality individuals and SNPs
plink --bfile hapmap3 \
      --mind 0.05 --geno 0.02 --maf 0.01 --hwe 1e-6 \
      --make-bed --out hapmap3_qc

# Step 2. LD pruning (for steps 3 and 4 only)
plink --bfile hapmap3_qc \
      --indep-pairwise 50 5 0.2 \
      --out hapmap3_pruned
plink --bfile hapmap3_qc \
      --extract hapmap3_pruned.prune.in \
      --make-bed --out hapmap3_pruned_set

# Step 3. Check heterozygosity on pruned set
plink --bfile hapmap3_pruned_set --het --out qc_het

# Step 4. Check relatedness on pruned set (PLINK)
plink --bfile hapmap3_pruned_set --genome --min 0.125 --out ibd_results

# Step 4b. Check relatedness with KING (use QC set, not pruned)
./king -b ~/Sociogenomics/Data/hapmap3_qc.bed --related --degree 2 --prefix king_related

# Step 5. (After manually checking outliers from steps 3-4)
# Remove outlier individuals and re-apply QC:
# plink --bfile hapmap3_qc --remove outliers.txt --make-bed --out hapmap3_final
```

> The final analysis dataset (`hapmap3_final`) uses the full set of QC-passing SNPs, **not** the LD-pruned set — unless PCA or heritability estimation specifically requires it.
