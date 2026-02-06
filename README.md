# Sociogenomics 2025/2026 — Computer Labs

This repository contains the lab materials for the Sociogenomics course (2025/2026).

## Quick start (Cloud Shell)
```
git clone https://github.com/nicolabarban/sociogenomics_2025_2026.git
cd sociogenomics_2025_2026
bash scripts/setup_plink19.sh
```

Open a lab:
```
ls labs
```

## Quick start (Google Colab)
In a Colab notebook cell:
```
!git clone https://github.com/nicolabarban/sociogenomics_2025_2026.git
%cd sociogenomics_2025_2026
!bash scripts/setup_plink19.sh
```

## Data
Small datasets used in the labs are included in `data/`.
Large datasets are hosted on Dropbox with download instructions in `data/README.md`.

## Repository layout
- `labs/` lab instructions and scripts by week
- `data/` small datasets for exercises
- `scripts/` setup and download utilities
- `lectures/` lab lecture slides
- `papers/` required/optional readings
