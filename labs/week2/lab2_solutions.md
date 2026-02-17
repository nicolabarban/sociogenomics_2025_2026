# Lab 2 — Solutions

## Exercise 1

**1. How many individuals? How many males and females?**

```
wc -l 1kg_hm3.fam
# 1092 individuals

awk '$5 == 1' 1kg_hm3.fam | wc -l
# Males: 525

awk '$5 == 2' 1kg_hm3.fam | wc -l
# Females: 567
```

**2. How many SNPs?**

```
wc -l 1kg_hm3.bim
# 851,065 SNPs
```

**3. How many SNPs on chromosome 1?**

```
awk '$1 == 1' 1kg_hm3.bim | wc -l
# 70,130 SNPs
```

**4. Find rs9930506**

```
grep rs9930506 1kg_hm3.bim
# 16	rs9930506	0	53830465	G	A
# Chromosome 16, alleles G and A
```

---

## Exercise 2

**1. MAF of rs9930506**

```
grep rs9930506 allele_freq.frq
```

**2. How many SNPs with MAF below 0.05?**

```
awk 'NR > 1 && $5 < 0.05' allele_freq.frq | wc -l
```

**3. Individual with highest missingness**

```
sort -k6 -n -r missing_report.imiss | head -2
```

The first line after the header shows the individual with the highest `F_MISS` value.

---

## Exercise 3

**1. Chromosome 6 SNPs with MAF > 0.10**

```
plink --bfile 1kg_hm3 --chr 6 --maf 0.10 --make-bed --out chr6_common
wc -l chr6_common.bim
```

**2. Extract three SNPs of your choice**

```
echo -e "rs1048488\nrs4970383\nrs9930506" > my_snps.txt
plink --bfile 1kg_hm3 --extract my_snps.txt --make-bed --out my_selected
wc -l my_selected.bim
# 3 SNPs
```

**3. Full QC pipeline**

```
plink --bfile 1kg_hm3 \
  --geno 0.05 \
  --mind 0.10 \
  --maf 0.01 \
  --hwe 1e-6 \
  --make-bed \
  --out 1kg_hm3_qc
```

Check the PLINK log output for the number of SNPs and individuals removed at each step and the final counts.
