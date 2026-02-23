# Lab 3. Quality Control, Relatedness, and Linkage Disequilibrium

In this lab you will perform a full quality control (QC) pipeline on genetic data, detect related individuals using PLINK and KING, and investigate patterns of linkage disequilibrium (LD). By the end you will be able to run a production-quality QC workflow, identify cryptic relatedness, and compute and visualise LD between SNPs.

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

Find outliers (more than 3 SD from the mean F):

```bash
awk 'NR>1 {
    sum += $6; sumsq += $6^2; n++
  } END {
    mean = sum/n
    sd = sqrt(sumsq/n - mean^2)
    print "Mean F:", mean, "SD:", sd
    print "Thresholds:", mean - 3*sd, "to", mean + 3*sd
  }' qc_het.het

awk 'NR==1 {next}
  {sum += $6; n++} END {mean=sum/n; print mean}' qc_het.het > /tmp/het_mean.txt

awk 'NR==1 {next}
  {sumsq += ($6)^2; n++} END {mean=ENVIRON["m"]; print sqrt(sumsq/n - mean^2)}' \
  qc_het.het

```

A simpler approach — flag individuals with $F < -0.15$ or $F > 0.15$:

```bash
awk 'NR>1 && ($6 < -0.15 || $6 > 0.15) {print $1, $2, $6}' qc_het.het
```

### Exercise 1

1. How many SNPs and individuals remain in `hapmap3_qc` after the combined QC filter?
2. How many SNPs are removed by LD pruning? What fraction of post-QC SNPs does this represent?
3. Are there any individuals with an extreme inbreeding coefficient ($|F| > 0.15$) in the HapMap data? What might this indicate?
4. Why do we apply LD pruning *before* the heterozygosity and relatedness checks, rather than after?

---

## Part II. Detecting Related Individuals

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

Run on the pruned dataset (using all pairs):

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

Find likely duplicates or MZ twins (PI\_HAT > 0.90):

```bash
awk 'NR>1 && $10 > 0.90 {print $1, $2, $3, $4, $10}' ibd_results.genome
```

**Note:** In a QC pipeline, you would remove one individual from each related pair. Typically, you keep the individual with lower missingness.

You can also run IBD in a computationally efficient way using a minimum PI\_HAT threshold:

```bash
plink --bfile hapmap3_pruned_set \
      --genome \
      --min 0.125 \
      --out ibd_related_only
```

This only reports pairs above the threshold (faster for large datasets).

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

KING kinship coefficients correspond to:

| Kinship | Relationship |
|---------|-------------|
| > 0.354 | Duplicate / MZ twin |
| 0.177 – 0.354 | 1st degree (parent-offspring, full siblings) |
| 0.0884 – 0.177 | 2nd degree (half-sibs, grandparent, avuncular) |
| 0.0442 – 0.0884 | 3rd degree (first cousins) |
| < 0.0442 | Unrelated |

Find all related pairs (2nd degree or closer):

```bash
awk 'NR>1 && $NF > 0.0884 {print}' king_results.kin0 | sort -k8 -n -r | head -20
```

**Key difference from PLINK:** KING reports **kinship coefficients** (approximately PI\_HAT / 2), so a parent-offspring pair has kinship ≈ 0.25, not 0.5. Negative kinship values indicate unrelated individuals from different populations — this is expected and not an error.

For more details on KING commands and options: [https://www.kingrelatedness.com/manual.shtml](https://www.kingrelatedness.com/manual.shtml)

### GCTA — Genetic Relationship Matrix

**GCTA** (Genome-wide Complex Trait Analysis; Yang et al. 2011) is primarily used for heritability estimation and mixed-model association testing, but it also computes a **Genetic Relationship Matrix (GRM)** that serves as a precise measure of pairwise relatedness. The GRM is the central object in GCTA's linear mixed model (LMM) analyses.

The off-diagonal element $A_{ij}$ of the GRM is:

$$A_{ij} = \frac{1}{N} \sum_{k=1}^{N} \frac{(x_{ik} - 2p_k)(x_{jk} - 2p_k)}{2p_k(1-p_k)}$$

where $x_{ik}$ is the genotype count (0/1/2) for individual $i$ at SNP $k$, $p_k$ is the allele frequency of SNP $k$, and $N$ is the number of SNPs. For unrelated individuals $A_{ij} \approx 0$; for identical twins $A_{ij} \approx 1$.

> **Note:** Use the QC-filtered dataset (not LD-pruned) for GRM computation. GCTA recommends retaining as many SNPs as possible to get a precise estimate.

#### Install GCTA

```bash
cd ~/Sociogenomics
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.2-linux-kernel-3-x86_64.zip
unzip gcta-1.94.2-linux-kernel-3-x86_64.zip
ln -s ~/Sociogenomics/gcta-1.94.2-linux-kernel-3-x86_64/gcta64 ~/Sociogenomics/gcta64
```

Verify the installation:

```bash
~/Sociogenomics/gcta64 --help | head -20
```

#### Compute the GRM

```bash
cd ~/Sociogenomics/Data

~/Sociogenomics/gcta64 \
    --bfile hapmap3_qc \
    --make-grm \
    --out grm_hapmap
```

This produces three binary files:
- `grm_hapmap.grm.bin` — the GRM values (lower triangle, binary)
- `grm_hapmap.grm.N.bin` — number of SNPs used for each pair
- `grm_hapmap.grm.id` — individual IDs (FID, IID), one per line

Check how many individuals are in the GRM:

```bash
wc -l grm_hapmap.grm.id
head grm_hapmap.grm.id
```

#### Convert GRM to text and inspect

```bash
~/Sociogenomics/gcta64 \
    --grm grm_hapmap \
    --make-grm-gz \
    --out grm_hapmap_text
```

This creates a gzipped text file `grm_hapmap_text.grm.gz` with 4 columns: individual index 1, individual index 2, number of SNPs, relatedness value.

Inspect some values:

```bash
zcat grm_hapmap_text.grm.gz | head -20
```

The diagonal elements (where index 1 = index 2) represent each individual's estimated inbreeding coefficient plus 1:

```bash
zcat grm_hapmap_text.grm.gz | awk '$1 == $2 {print $3, $4}' | head -10
```

Off-diagonal elements for unrelated individuals should be close to 0:

```bash
zcat grm_hapmap_text.grm.gz | awk '$1 != $2 && $4 > 0.05 {print}' | wc -l
```

#### Identify and remove related individuals

GCTA can iteratively remove one individual from each related pair until no pairs exceed a given threshold:

```bash
~/Sociogenomics/gcta64 \
    --grm grm_hapmap \
    --grm-cutoff 0.05 \
    --make-grm \
    --out grm_unrelated
```

This outputs a new GRM (`grm_unrelated`) and a list of individuals that were **removed** (`grm_unrelated.grm.id` contains those retained). Check how many were removed:

```bash
wc -l grm_hapmap.grm.id
wc -l grm_unrelated.grm.id
```

To export the list of unrelated individuals for use with PLINK:

```bash
# Get retained IDs from GCTA
cat grm_unrelated.grm.id > unrelated_ids.txt

# Apply to PLINK dataset
plink --bfile hapmap3_qc \
      --keep unrelated_ids.txt \
      --make-bed \
      --out hapmap3_unrelated
```

#### Comparison of relatedness tools

| Feature | PLINK `--genome` | KING | GCTA GRM |
|---------|-----------------|------|----------|
| Metric | PI\_HAT (IBD proportion) | Kinship coefficient | Genetic relatedness ($A_{ij}$) |
| Input | LD-pruned recommended | QC-filtered (not pruned) | QC-filtered (not pruned) |
| Population structure | Sensitive | Robust | Method-of-moments |
| Primary use | Relatedness QC | Relatedness QC | Heritability / LMM / QC |
| Unrelated threshold | PI\_HAT < 0.10 | Kinship < 0.0442 | $A_{ij}$ < 0.05 |
| Speed (large $N$) | Slow (all pairs) | Fast | Moderate |

For a QC pipeline where you primarily want to remove close relatives before GWAS, **KING** is generally preferred because it is robust to population structure. **GCTA's GRM** is essential when you plan to run mixed-model analyses (GCTA-REML for heritability, or BOLT-LMM / SAIGE for association) since those methods require the GRM as input.

For more details on GCTA: [https://yanglab.westlake.edu.cn/software/gcta/](https://yanglab.westlake.edu.cn/software/gcta/)

### Exercise 2

1. How many pairs of individuals in HapMap3 have a PI\_HAT > 0.20? What is the likely family relationship for the pair with the highest PI\_HAT?
2. What is the theoretical PI\_HAT for:
   - Monozygotic (identical) twins?
   - Parent–offspring pairs?
   - Full siblings?
   - First cousins?
3. Why does KING perform better than PLINK's `--genome` in samples with mixed ancestry? (*Hint: think about what happens to allele sharing when two individuals are from different populations.*)
4. In a GWAS study design, if you identify two first-degree relatives in your sample, which one would you remove and why?
5. Look at the diagonal elements of the GCTA GRM (`$1 == $2`). What do values above 1.0 indicate? What about values substantially below 1.0?
6. What is the difference between the GRM off-diagonal element $A_{ij}$ and KING's kinship coefficient? How do they relate to PI\_HAT?

---

## Part III. Linkage Disequilibrium

**Linkage disequilibrium (LD)** is the non-random association of alleles at different loci. It arises from:
- **Physical linkage** (nearby SNPs tend to be co-inherited)
- **Population history** (founder effects, bottlenecks, admixture)
- **Selection** (variants near selected loci have inflated LD)

The most common measure is $r^2$ — the squared Pearson correlation between allele counts at two SNPs.

### Pairwise LD between specific SNPs

Calculate $r^2$ between two specific SNPs:

```bash
plink --bfile hapmap3_qc \
      --ld rs1048488 rs3115850 \
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
To examine this, let's compute the mean $r^2$ in distance bins:

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

### Compare LD between population groups

HapMap3 contains individuals from multiple ancestry groups. LD patterns differ across populations. Let's compare LD for the same SNP pair in different populations.

First, check what population labels are available (the 3rd column of the FAM file often contains population codes):

```bash
awk '{print $1}' hapmap3.fam | sort | uniq -c | sort -n -r | head -20
```

Extract European-ancestry individuals (CEU in HapMap):

```bash
awk '$1 == "CEU" {print $1, $2}' hapmap3.fam > samples_CEU.txt
wc -l samples_CEU.txt

plink --bfile hapmap3_qc \
      --keep samples_CEU.txt \
      --ld rs1048488 rs3115850 \
      --out ld_CEU
cat ld_CEU.log | grep -A5 "R-sq"
```

Extract Yoruba (YRI) individuals:

```bash
awk '$1 == "YRI" {print $1, $2}' hapmap3.fam > samples_YRI.txt
wc -l samples_YRI.txt

plink --bfile hapmap3_qc \
      --keep samples_YRI.txt \
      --ld rs1048488 rs3115850 \
      --out ld_YRI
cat ld_YRI.log | grep -A5 "R-sq"
```

**Question:** Are the LD values the same in CEU and YRI? Why might they differ?

### Exercise 3

1. What is the $r^2$ value between `rs1048488` and `rs3115850` in the full sample? Is it the same in CEU and YRI?
2. On chromosome 22, find the 5 SNP pairs with the highest $r^2$. What are their physical distances?
3. Does LD generally increase or decrease with physical distance between SNPs? Is this what you expected?
4. After LD pruning (from Part I), how many SNPs remained? Compare with the number before pruning. What window size and $r^2$ threshold did we use?

---

## Part IV. External Resources for LD Lookup

You do not always need to compute LD yourself. Several web tools let you look up LD for specific SNPs in reference populations:

### LDlink (NIH)

🔗 **[https://ldlink.nci.nih.gov/](https://ldlink.nci.nih.gov/)**

LDlink is a web-based tool backed by the 1000 Genomes Project data. It offers:

- **LDpair** — $r^2$ and $D'$ between two SNPs in a chosen population
- **LDproxy** — all SNPs in LD with a query SNP (useful after a GWAS hit)
- **LDpop** — compare LD across populations
- **LDmatrix** — pairwise LD matrix for a set of SNPs
- **SNPchip** — check if a SNP is on common genotyping arrays
- **RegulomeDB integration** — functional annotation of LD proxies

**Try it:** Look up `rs1048488` and `rs3115850` in the **CEU** population and compare with the PLINK result you computed above.

### Ensembl LD Calculator

🔗 **[https://www.ensembl.org/Homo_sapiens/Tools/LD](https://www.ensembl.org/Homo_sapiens/Tools/LD)**

Ensembl provides LD visualisation alongside genomic annotation (genes, regulatory elements, variants). Useful when you want to see LD in its genomic context.

### HaploReg

🔗 **[https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php](https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php)**

HaploReg annotates LD-linked variants with:
- eQTL data (gene expression)
- Regulatory element overlaps (enhancers, promoters)
- Motif disruptions
- Conservation scores

Useful after GWAS: find all SNPs in LD with your hit, and see which might be functionally relevant.

### SNAP (SNP Annotation and Proxy Search)

🔗 **[https://www.broadinstitute.org/snap/snap](https://www.broadinstitute.org/snap/snap)**

SNAP finds proxy SNPs in LD with a query SNP, given a population and $r^2$ threshold. Useful for identifying tagging SNPs.

### LDlink API (for scripts)

LDlink provides an API you can query from the command line. First register for a free API token at [https://ldlink.nci.nih.gov/?tab=apiaccess](https://ldlink.nci.nih.gov/?tab=apiaccess).

Then query from Google Cloud Shell:

```bash
TOKEN="your_token_here"   # replace with your token

curl -k "https://ldlink.nci.nih.gov/LDlinkRest/ldpair?var1=rs1048488&var2=rs3115850&pop=CEU&token=${TOKEN}"
```

This returns JSON with $r^2$, $D'$, haplotype frequencies, and a chi-squared test.

---

## Summary: Full QC Pipeline

Here is the complete recommended pipeline in a single block for reference:

```bash
# Step 1. Start from raw data
# Step 2. Remove low-quality individuals and SNPs
plink --bfile hapmap3 \
      --mind 0.05 --geno 0.02 --maf 0.01 --hwe 1e-6 \
      --make-bed --out hapmap3_qc

# Step 3. LD pruning (for steps 4 and 5 only)
plink --bfile hapmap3_qc \
      --indep-pairwise 50 5 0.2 \
      --out hapmap3_pruned
plink --bfile hapmap3_qc \
      --extract hapmap3_pruned.prune.in \
      --make-bed --out hapmap3_pruned_set

# Step 4. Check heterozygosity on pruned set
plink --bfile hapmap3_pruned_set --het --out qc_het

# Step 5. Check relatedness on pruned set (PLINK)
plink --bfile hapmap3_pruned_set --genome --min 0.125 --out ibd_results

# Step 5b. Check relatedness with KING (use QC set, not pruned)
./king -b ~/Sociogenomics/Data/hapmap3_qc.bed --related --degree 2 --prefix king_related

# Step 5c. Compute GRM with GCTA (use QC set, not pruned)
~/Sociogenomics/gcta64 --bfile hapmap3_qc --make-grm --out grm_hapmap
~/Sociogenomics/gcta64 --grm grm_hapmap --grm-cutoff 0.05 --make-grm --out grm_unrelated

# Step 6. (After manually checking outliers from steps 4-5)
# Remove outlier individuals and re-apply QC:
# plink --bfile hapmap3_qc --remove outliers.txt --make-bed --out hapmap3_final
```

> The final analysis dataset (`hapmap3_final`) uses the full set of QC-passing SNPs, **not** the LD-pruned set — unless PCA or heritability estimation specifically requires it.

---

## References

- Anderson, C.A. et al. (2010). Data quality control in genetic case-control association studies. *Nature Protocols*, 5, 1564–1573.
- Marees, A.T. et al. (2018). A tutorial on conducting genome-wide association studies: Quality control and statistical analysis. *International Journal of Methods in Psychiatric Research*, 27, e1608.
- Manichaikul, A. et al. (2010). Robust relationship inference in genome-wide association studies. *Bioinformatics*, 26(22), 2867–2873. [KING paper]
- Purcell, S. et al. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses. *American Journal of Human Genetics*, 81(3), 559–575.
- Yang, J. et al. (2011). GCTA: a tool for genome-wide complex trait analysis. *American Journal of Human Genetics*, 88(1), 76–82.
- GWAS Tutorial QC Reference: [https://cloufield.github.io/GWASTutorial/04_Data_QC/](https://cloufield.github.io/GWASTutorial/04_Data_QC/)
- GCTA documentation: [https://yanglab.westlake.edu.cn/software/gcta/](https://yanglab.westlake.edu.cn/software/gcta/)
- LDlink: [https://ldlink.nci.nih.gov/](https://ldlink.nci.nih.gov/)
