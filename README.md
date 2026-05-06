# Mutation Dynamics of gyrA QRDR in Campylobacter spp.
## Genomic Insights from Three Decades of Global One Health Surveillance

**Authors:** Julien Alban Nguinkal, Evariste Bako, Yedomon Ange Bovys Zoclanclounon,
Muna Affara, Florian Gehre, Jürgen May

---

## Repository structure

```
campylobacter-gyra-fqr/
│
├── data/
│   ├── raw/                  # Accession lists from BV-BRC and NCBI GenBank
│   ├── processed/            # Deduplicated master dataset (post-QC, n = 5,963)
│   └── metadata/             # Host, country, year, isolation-source annotations
│
├── scripts/
│   ├── 01_data_collection/   # Download and deduplication
│   ├── 02_qc/                # BUSCO / CheckM wrappers
│   ├── 03_resistance_screening/  # abritAMR, AMRFinderPlus, QRDR extraction
│   ├── 04_mlst_cc/           # MLST and clonal complex assignment
│   ├── 05_phylogeny/         # Prokka → Panaroo → Gubbins → IQ-TREE pipeline
│   ├── 06_statistical_analysis/  # Prevalence, bootstrap CI, Fisher, BH-FDR (R)
│   └── 07_figures/           # All plotting scripts (R / Python)
│
├── results/
│   ├── qc/                   # Assembly quality summaries
│   ├── resistance/           # Per-genome gyrA mutation calls
│   ├── mlst/                 # ST and CC assignment table
│   ├── phylogeny/            # Tree files and annotation tables
│   ├── prevalence/           # Regional and host-stratified prevalence tables
│   └── statistical_tests/    # Test outputs (Fisher, permutation, BH-FDR)
│
├── figures/
│   ├── main/                 # Manuscript figures (Figures 1–4)
│   └── supplementary/        # Supplementary figures
│
├── supplementary/
│   ├── tables/               # Supplementary tables (xlsx / csv)
│   └── data/                 # Supplementary datasets
│
├── envs/                     # Conda environments / renv lockfile
├── .gitignore
├── LICENSE
└── README.md
```

## Reproducibility

All analyses were run under:
- **R** v4.5.1 (packages: ggplot2, dplyr, tidyr, scales, forcats, ape, ggtree)
- **Python** v3.12
- **Prokka** v1.14 | **Panaroo** v1.3 | **Gubbins** v3.4.3 | **IQ-TREE2** v2.2.0
- **mlst** v2.19.0 | **abritAMR** v1.0.19 | **AMRFinderPlus** v4.0.23

See `envs/` for conda environment files.

## Data availability

Raw genome assemblies are publicly available at NCBI GenBank and BV-BRC.
Accession numbers are listed in `data/raw/accessions.txt`.
The deduplicated master metadata table is provided in `data/processed/`.

## Citation

> Nguinkal JA, Bako E, et al. Mutation Dynamics of gyrA QRDR in *Campylobacter* spp.:
> Genomic Insights from Three Decades of Global One Health Surveillance. *[Journal]*, 2025.

## License

[MIT / CC-BY 4.0 — choose and update]
