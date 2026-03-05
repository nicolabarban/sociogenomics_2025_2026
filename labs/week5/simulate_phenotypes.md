# Simulating Phenotypes from Genotype Data

This document describes three approaches for simulating phenotypes from real genotype data,
suitable for GWAS teaching exercises where the true causal variants are known.

---

## Option 1 — GCTA `--simu-qt` (recommended)

Download GCTA: [yanglab.westlake.edu.cn/software/gcta](https://yanglab.westlake.edu.cn/software/gcta/)

**Step 1 — Pick causal SNPs** (randomly sample from your BIM file):

```bash
# Pick 20 random causal SNPs
awk 'NR>1 {print $2}' hapmap3_qc.bim | shuf | head -20 > causal_snps.txt
```

**Step 2 — Simulate quantitative phenotype:**

```bash
gcta64 --bfile hapmap3_qc \
       --simu-qt \
       --simu-causal-loci causal_snps.txt \
       --simu-hsq 0.5 \
       --simu-rep 1 \
       --out simulated_pheno
```

| Flag | Meaning |
|------|---------|
| `--simu-qt` | quantitative trait (use `--simu-cc N_cases N_controls` for binary) |
| `--simu-causal-loci` | file with causal SNP IDs |
| `--simu-hsq` | heritability (0–1); 0.5 = half variance explained by genetics |
| `--simu-rep` | number of independent phenotype replicates |

**Output:** `simulated_pheno.phen` (FID IID phenotype) — plug directly into PLINK `--pheno`.

---

## Option 2 — R (full control, no extra software)

This approach lets you set effect sizes, heritability, and noise explicitly.
Students can verify their GWAS hits against the known causal variants.

```r
library(BEDMatrix)  # install.packages("BEDMatrix")

G <- BEDMatrix("hapmap3_qc")   # individuals × SNPs matrix

# Pick 20 causal SNPs randomly
set.seed(42)
n_causal   <- 20
causal_idx <- sample(ncol(G), n_causal)
G_causal   <- as.matrix(G[, causal_idx])

# Standardise genotypes (mean 0, sd 1)
G_std <- scale(G_causal)

# Simulate effects from N(0,1)
beta <- rnorm(n_causal)

# Genetic component
g <- G_std %*% beta

# Scale to desired heritability h2 = 0.5
h2  <- 0.5
Vg  <- var(g)
Ve  <- Vg * (1 - h2) / h2
eps <- rnorm(nrow(G), sd = sqrt(Ve))

# Final phenotype
pheno <- g + eps

# Save in PLINK format (FID IID pheno)
fam <- read.table("hapmap3_qc.fam")
out <- data.frame(FID = fam$V1, IID = fam$V2, PHENO = pheno)
write.table(out, "simulated_pheno.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
```

---

## Option 3 — PLINK `--simulate` (synthetic data from scratch)

Generates both genotypes and phenotype without real input data.
Fast for demonstration but no real LD structure.

```bash
# Create a simulation parameter file
cat > sim.params << 'EOF'
10000 SNP_BLOCK 0.1 0.2 0.3 1
EOF

plink --simulate sim.params \
      --simulate-qt \
      --out synthetic_gwas
```

---

## Recommendation

Use **GCTA (Option 1)**:
- Works directly on existing `hapmap3_qc` data
- Preserves real LD structure
- Heritability is set to a known value — useful for teaching
- Output `.phen` file is directly compatible with PLINK `--pheno`
- Students can compare GWAS hits against the true causal SNP list

---

## References

- [GCTA documentation](https://yanglab.westlake.edu.cn/software/gcta/)
- [Biostars: simulate phenotype from genotype](https://www.biostars.org/p/428314/)
- [Biostars: simulate from 1000 Genomes](https://www.biostars.org/p/406828/)
- [PhenotypeSimulator — Bioinformatics 2018](https://academic.oup.com/bioinformatics/article/34/17/2951/4956348)
