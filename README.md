# Disulfidptosis as a Prognostic Biomarker in Human Cancers
### A Systematic Review and Meta-Analysis with Bioinformatic Validation

> Code repository for the manuscript submitted for peer review.  
> All analyses were conducted in R ≥ 4.2.

---

## Repository structure

```
disulfidptosis-hcc/
│
├── README.md
│
├── 1_TCGA_Master/
│   └── Disulfidptosis_LIHC_Master_FIXED.R
│
├── 2_TCGA_EPIC_Patch/
│   └── Recorrer_Correcciones_EPIC.R
│
├── 3_GEO_Validation/
│   └── Fase4_Validacion_Externa_GEO.R
│
├── 4_Meta_Analysis/
│   └── Disulfidptosis_Metaanalysis.R
│
└── 5_QUIPS/
│   └── QUIPS_Disulfidptosis_plots.R
│
└── 6_Material_Suplementario/
    └── ...
```

---

## Script descriptions

| Script | Purpose | External data required |
|--------|---------|----------------------|
| `1_TCGA_Master` | Downloads TCGA-LIHC data, computes disulfidptosis scores, clustering, survival analysis (KM + Cox), tumor microenvironment deconvolution (EPIC), DEG, GSEA GO/KEGG, ORA, and supplementary Excel tables | TCGA-LIHC (auto-downloaded via TCGAbiolinks) |
| `2_TCGA_EPIC_Patch` | Re-runs only Blocks 5, 7, and 8 of the Master with EPIC corrections. **Run this after the Master if EPIC output needs to be regenerated** — does not re-download TCGA data | Output files from Script 1 |
| `3_GEO_Validation` | External validation in GSE14520 and GSE76427 (HCC cohorts); nomogram; immune checkpoint correlations | GSE14520, GSE76427 (auto-downloaded via GEOquery); output files from Script 1 |
| `4_Meta_Analysis` | Quantitative meta-analysis of 37 studies (HR pooled, subgroup analyses, leave-one-out, publication bias, forest plot, funnel plot) | None — all data embedded in the script |
| `5_QUIPS` | QUIPS risk-of-bias assessment figures for all 78 included studies | None — all data embedded in the script |

---

## Execution order

Scripts **must** be run in the following order:

```
1_TCGA_Master  →  2_TCGA_EPIC_Patch  →  3_GEO_Validation
```

Scripts 4 and 5 are **independent** and can be run at any time:

```
4_Meta_Analysis    (standalone)
5_QUIPS            (standalone)
```

> **Note on Script 2:** This patch is only needed if you need to regenerate the EPIC/TME figures or the supplementary Excel table after the Master has already run. It reads intermediate `.rds` files produced by Script 1 and overwrites only the affected outputs (Figures 07–09, Sheet S8 of the Excel file).

---

## Requirements

### R version
R ≥ 4.2.0 and ggplot2 ≥ 3.4.0 (required for `linewidth` parameter).

### CRAN packages
```r
install.packages(c(
  "meta", "metafor",
  "ggplot2", "ggpubr", "ggrepel", "patchwork", "ggtext", "cowplot",
  "dplyr", "tidyr", "stringr", "tibble", "forcats",
  "scales", "grid", "gridExtra", "RColorBrewer",
  "survival", "survminer",
  "cluster", "factoextra",
  "openxlsx", "data.table",
  "rms", "showtext", "remotes"
))
```

### Bioconductor packages
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "TCGAbiolinks", "SummarizedExperiment",
  "GEOquery", "Biobase", "limma",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot",
  "maftools"   # optional — only needed for TMB analysis
))
```

### GitHub packages
```r
# dmetar — optional, used in meta-analysis (fallback available without it)
remotes::install_github("MathiasHarrer/dmetar", upgrade = "never")

# IOBR — required for EPIC deconvolution (pinned version for reproducibility)
remotes::install_github("IOBR/IOBR", ref = "v0.99.9", upgrade = "never")
```

---

## Data availability

| Dataset | Source | Access |
|---------|--------|--------|
| TCGA-LIHC (RNA-seq + clinical) | GDC Data Portal | Public — auto-downloaded by Script 1 via `TCGAbiolinks` |
| GSE14520 | NCBI GEO | Public — auto-downloaded by Script 3 via `GEOquery` |
| GSE76427 | NCBI GEO | Public — auto-downloaded by Script 3 via `GEOquery` |
| Meta-analysis HR data (37 studies) | Embedded in Script 4 | See manuscript Table 1 for primary sources |
| QUIPS assessments (78 studies) | Embedded in Script 5 | See manuscript Supplementary Material S5 |

---

## Outputs

### Script 1 — TCGA Master
Saved to `resultados/`:
- `06_figuras/` — 14 figures (PNG 300 dpi + PDF)
- `07_tablas_suplementarias/Tablas_Suplementarias_Disulfidptosis_LIHC.xlsx` — 9-sheet Excel file
- `03_supervivencia/Datos_supervivencia.rds` — survival data (used by Scripts 2 and 3)
- `01_datos_procesados/TCGA_LIHC_log2fpkm_SYMBOL.rds` — normalized expression matrix (used by Scripts 2 and 3)

### Script 3 — GEO Validation
Saved to `resultados_fase4/`:
- `06_figuras_fase4/` — KM curves, nomogram, checkpoint heatmap

### Script 4 — Meta-analysis
Saved to `Disulfidptosis_Figures/`:
- 9 figures (PDF + PNG) + 2 CSV summary tables

### Script 5 — QUIPS
Saved to working directory:
- 4 individual figures + 1 combined panel (PDF)

---

## Known issues and notes

- **Script 1, Block 5:** A variable naming bug (`cibersort_res`) exists in the original Master script. **Script 2 contains the corrected version** and should be run after Script 1 to regenerate EPIC outputs.
- **`MY_LIB` path:** Lines referencing a local R library path (`C:/Users/...`) have been removed. If your R installation requires a custom library path, add `.libPaths("your/path")` at the top of the relevant scripts before running.
- **Internet access:** Scripts 1 and 3 require internet access on the first run to download TCGA and GEO data. Subsequent runs use cached `.rds` files.
- **RAM:** Script 1 requires ~16 GB RAM. Script 3 requires ~8 GB RAM.
- **TMB analysis (Script 3, Block 4.5):** Requires a MAF file for TCGA-LIHC. If not present, this block is skipped automatically with an informational message.

---

## License

Code released under the [MIT License](LICENSE).  
Data from TCGA and GEO are subject to their respective access policies.
