# CITE-seq Malaria Immunology Pipeline

Analysis pipeline for 10x Genomics CITE-seq data from a controlled human malaria infection (CHMI) study. Characterises peripheral blood immune cell composition and transcriptional responses across three timepoints (pre-infection, acute, and convalescent) in patients with febrile (FEB) and chronic-tolerant (CHR) malaria phenotypes.

## Overview

| Script | Purpose |
|--------|---------|
| `config.R` | All file paths and tunable parameters — **edit this first** |
| `utils.R` | Shared helper functions (sourced by each script) |
| `01_integration.R` | Load, QC, doublet removal, SCT/DSB normalisation, RNA+ADT RPCA integration, WNN clustering |
| `02_subset_analysis.R` | Iterative WNN sub-clustering of T/NK, CD8, NK, γδ T, B, CD4/Treg lineages |
| `03_deg_mast.R` | MAST differential expression (CHR vs FEB) per celltype × timepoint |
| `04_visualisation.R` | UMAP SVGs, RNA/ADT dotplots, relative abundance plots, ComplexHeatmaps |
| `05_milo_da.R` | MiloR differential abundance testing |

## Data

Input: per-run Seurat objects (`.rds`) produced upstream by HTO demultiplexing. Ten runs (RO1–RO10) covering 10 patients across three timepoints:

- **C1** — Day 1 (pre-infection baseline)
- **C9** — Day 9 (acute infection)
- **C14** — Day 14 (convalescent)

Patient conditions:
- **FEB** — febrile malaria (n = 9)
- **CHR** — chronic tolerant (n = 1, patient P06)

## Requirements

```r
# CRAN
install.packages(c("Seurat", "tidyverse", "dplyr", "ggplot2",
                   "svglite", "patchwork", "gtable",
                   "Polychrome", "RColorBrewer", "openxlsx"))

# Bioconductor
BiocManager::install(c("DoubletFinder", "dsb", "STACAS",
                       "ComplexHeatmap", "circlize",
                       "MAST", "miloR", "scGate",
                       "SingleCellExperiment", "SummarizedExperiment"))
```

## Usage

1. Edit `config.R` — set `RDS_DIR` to your input data folder and adjust `OUTPUT_DIR` / `PLOT_DIR` as needed.
2. Run scripts in order:

```r
source("01_integration.R")
source("02_subset_analysis.R")
source("03_deg_mast.R")
source("04_visualisation.R")
source("05_milo_da.R")
```

Each script saves `.rds` checkpoints to `OUTPUT_DIR` so steps can be re-run independently.

## Methods summary

- **QC**: cells filtered on nFeature_RNA (200–6000) and mitochondrial reads (< 20 %)
- **Doublet removal**: DoubletFinder with pK optimised per sample (~5 % expected doublet rate)
- **RNA normalisation**: SCTransform (v2)
- **Protein normalisation**: DSB `ModelNegativeADTnorm` (background estimation from isotype controls)
- **Integration**: Seurat RPCA anchor-based integration (separate RNA and ADT)
- **Clustering**: Weighted Nearest Neighbour (WNN) combining RNA and protein modalities
- **Differential expression**: MAST with sex as a latent covariate (where both sexes present per subset)
- **Differential abundance**: MiloR neighbourhood testing per timepoint

## Output structure

```
output/
  integ_rna_annotated.rds
  tnk_all.rds
  cd8nk_subset.rds
  nk_focus.rds / gd_object.rds / cd8_cells.rds / b_cells.rds / CD4_and_treg.rds
  deg_results_mast.rds

plots/
  umap/        — per-object UMAP + legend SVGs
  dotplots/    — RNA and ADT dotplots + DEG CSVs
  freq/        — relative abundance line plots
  heatmaps/    — ComplexHeatmap SVGs
  milo/        — MiloR DA UMAP and beeswarm SVGs
```
