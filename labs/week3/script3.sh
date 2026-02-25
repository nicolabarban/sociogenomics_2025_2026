#!/usr/bin/env bash
# Lab 3 — Quality Control, Relatedness, and Linkage Disequilibrium
# Run from: ~/Sociogenomics/Data
set -euo pipefail

cd ~/Sociogenomics/Data

echo "============================================"
echo "Part I. Quality Control Pipeline"
echo "============================================"

# Step 1 — Compute missing rates
echo "--- Step 1: Computing missing rates ---"
plink --bfile hapmap3 \
      --missing \
      --out qc_missing

echo "Per-individual missingness (top 10 worst):"
sort -k6 -n -r qc_missing.imiss | head -10

echo ""
echo "Per-SNP missingness (top 10 worst):"
sort -k5 -n -r qc_missing.lmiss | head -10

# Step 2 — Compute allele frequencies
echo ""
echo "--- Step 2: Computing allele frequencies ---"
plink --bfile hapmap3 \
      --freq \
      --out qc_freq

echo "MAF distribution:"
awk 'NR>1 {print $5}' qc_freq.frq | \
  awk '{
    if ($1 < 0.01) rare++
    else if ($1 < 0.05) lowfreq++
    else common++
  } END {
    print "Rare (<0.01):", rare+0
    print "Low-freq (0.01-0.05):", lowfreq+0
    print "Common (>=0.05):", common+0
  }'

# Step 3 — Compute Hardy-Weinberg statistics
echo ""
echo "--- Step 3: Hardy-Weinberg equilibrium ---"
plink --bfile hapmap3 \
      --hardy \
      --out qc_hwe

echo "SNPs deviating from HWE (p < 1e-6):"
awk 'NR>1 && $9 < 1e-6' qc_hwe.hwe | wc -l

# Step 4 — Apply all QC filters
echo ""
echo "--- Step 4: Applying QC filters ---"
plink --bfile hapmap3 \
      --mind 0.05 \
      --geno 0.02 \
      --maf  0.01 \
      --hwe  1e-6 \
      --make-bed \
      --out hapmap3_qc

echo "QC summary:"
grep -E "removed|remaining|pass" hapmap3_qc.log || true

# Step 5 — LD pruning
echo ""
echo "--- Step 5: LD pruning ---"
plink --bfile hapmap3_qc \
      --indep-pairwise 50 5 0.2 \
      --out hapmap3_pruned

echo "SNPs kept after pruning:"
wc -l < hapmap3_pruned.prune.in
echo "SNPs removed by pruning:"
wc -l < hapmap3_pruned.prune.out

plink --bfile hapmap3_qc \
      --extract hapmap3_pruned.prune.in \
      --make-bed \
      --out hapmap3_pruned_set

# Step 6 — Heterozygosity check
echo ""
echo "--- Step 6: Heterozygosity / inbreeding ---"
plink --bfile hapmap3_pruned_set \
      --het \
      --out qc_het

echo "Outlier individuals (|F| > 0.15):"
awk 'NR>1 && ($6 < -0.15 || $6 > 0.15) {print $1, $2, $6}' qc_het.het || echo "(none)"


echo ""
echo "============================================"
echo "Part II. Linkage Disequilibrium"
echo "============================================"

# Pairwise LD between two SNPs
echo "--- Pairwise LD: rs994335 vs rs2379903 ---"
plink --bfile hapmap3_qc \
      --ld rs994335 rs2379903 \
      --out ld_pair
grep -A5 "LD" ld_pair.log || true

# LD across chromosome 22
echo ""
echo "--- LD across chromosome 22 ---"
plink --bfile hapmap3_qc \
      --chr 22 \
      --r2 \
      --ld-window 100 \
      --ld-window-kb 1000 \
      --ld-window-r2 0.05 \
      --out ld_chr22

echo "Number of SNP pairs with r2 > 0.05 on chr 22:"
wc -l < ld_chr22.ld

echo ""
echo "Top 5 SNP pairs by r2:"
sort -k7 -n -r ld_chr22.ld | head -5

echo ""
echo "LD decay (mean r2 by distance in kb):"
awk 'NR>1 {
    dist = ($5 - $2)
    if (dist < 0) dist = -dist
    bin = int(dist/10000)
    count[bin]++
    sumr2[bin] += $7
} END {
    for (b in count) printf "%d kb\t%.4f\n", b*10, sumr2[b]/count[b]
}' ld_chr22.ld | sort -k1 -n | head -30

# Population-specific LD comparison
echo ""
echo "--- Compare LD across populations (using 1kg_samples.txt) ---"
cp ~/sociogenomics_2025_2026/data/1kg_samples.txt ~/Sociogenomics/Data/ 2>/dev/null || true

echo "Superpopulation counts:"
awk -F'\t' 'NR>1 {print $6, $7}' 1kg_samples.txt | sort | uniq -c | sort -n -r

echo ""
echo "Creating population sample lists..."
awk -F'\t' 'NR>1 && $6 == "EUR" {print $1, $1}' 1kg_samples.txt > samples_EUR.txt
awk -F'\t' 'NR>1 && $6 == "AFR" {print $1, $1}' 1kg_samples.txt > samples_AFR.txt
awk -F'\t' 'NR>1 && $6 == "EAS" {print $1, $1}' 1kg_samples.txt > samples_EAS.txt

for POP in EUR AFR EAS; do
    n=$(plink --bfile hapmap3_qc --keep samples_${POP}.txt --make-just-fam --out check_${POP} 2>/dev/null && wc -l < check_${POP}.fam)
    echo "$POP: $n individuals in dataset"
done

echo ""
echo "--- Pairwise LD by population ---"
for POP in EUR AFR EAS; do
    plink --bfile hapmap3_qc --keep samples_${POP}.txt \
          --ld rs994335 rs2379903 --out ld_${POP} 2>/dev/null
    echo -n "$POP: "
    grep "R-sq" ld_${POP}.log || echo "(not available)"
done

echo ""
echo "--- LD decay by population (chr 22) ---"
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


echo ""
echo "============================================"
echo "Part III. Detecting Related Individuals"
echo "============================================"

# IBD estimation with PLINK
echo "--- IBD estimation (PLINK --genome) ---"
plink --bfile hapmap3_pruned_set \
      --genome \
      --out ibd_results

echo "Pairs with PI_HAT > 0.20:"
awk 'NR>1 && $10 > 0.20 {print $1, $2, $3, $4, "PI_HAT="$10}' ibd_results.genome | sort -k5 -n -r

echo ""
echo "Parent-offspring or sibling pairs (PI_HAT 0.40-0.60):"
awk 'NR>1 && $10 > 0.40 && $10 < 0.60 {print $1, $2, $3, $4, "Z0="$7, "Z1="$8, "Z2="$9, "PI_HAT="$10}' ibd_results.genome | sort -k8 -n -r

echo ""
echo "IBD with min threshold (faster for large datasets):"
plink --bfile hapmap3_pruned_set \
      --genome \
      --min 0.125 \
      --out ibd_related_only

echo "Related pairs (PI_HAT >= 0.125):"
awk 'NR>1' ibd_related_only.genome | wc -l

# KING kinship estimation
echo ""
echo "--- KING kinship estimation ---"
if [ ! -f ~/Sociogenomics/king ]; then
    echo "Installing KING..."
    cd ~/Sociogenomics
    wget -q https://www.kingrelatedness.com/Linux-king.tar.gz
    tar -xzf Linux-king.tar.gz
    chmod +x king
    cd ~/Sociogenomics/Data
fi

~/Sociogenomics/king -b ~/Sociogenomics/Data/hapmap3_qc.bed \
       --kinship \
       --prefix king_results

echo ""
echo "Between-family pairs (kin0):"
head king_results.kin0

echo ""
echo "Within-family pairs (kin):"
head king_results.kin

echo ""
echo "Related pairs (kinship > 0.0884, 2nd degree or closer):"
awk 'NR>1 && $NF > 0.0884 {print}' king_results.kin0 | sort -k8 -n -r | head -20

echo ""
echo "============================================"
echo "Lab 3 complete!"
echo "============================================"
