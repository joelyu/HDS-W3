# Healthcare Data Science Module 3 Assignment

## Immune Cell Infiltration in METABRIC

This is the submission GitHub repository for the Module 3 Assignment for the MSt program in Healthcare Data Science, exploring the infiltration of immune cells of tumour samples of the METABRIC breast cancer dataset using the marker gene scoring from Danaher et al. (2017). A live build of the rendered page is [here](https://joelyu.github.io/HDS-W3/HDS_03_YuChungYan_2602.html).

### Environment Setup

Requires [Quarto](https://quarto.org/docs/get-started/) and [Miniforge](https://github.com/conda-forge/miniforge) (mamba).

```bash
# Create environment
mamba env create -f environment.yml

# Activate
mamba activate immfiltration

# Install Bioconductor packages (bioconda dependency chains broken on osx-arm64)
Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install(c("ComplexHeatmap", "cBioPortalData"))'
```

Bioconductor packages are also auto-installed via `if(!require)` in the qmd setup chunk.

### Rendering

```bash
mamba activate immfiltration
quarto render HDS_03_YuChungYan_2602.qmd
```

The submission qmd sources scripts in order. No need to render intermediate files.

### Interactive Testing

Scripts don't include `library()` calls — the qmd setup chunk handles package loading. To run scripts in the console (e.g. Positron), load packages first:

```r
library(tidyverse)    # covers scripts 00-04
library(gtsummary)    # needed for script 02 (Table 1)
library(multcompView) # needed for script 03 (compact letter display)
library(survival)     # needed from script 05 onward
```

Scripts must run in order: each depends on objects produced by earlier scripts.

### File Structure

```
├── HDS_03_YuChungYan_2602.qmd       # Submission document — sources scripts/
├── scripts/
│   ├── _immune_markers.R             # Danaher gene lists (single source of truth)
│   ├── 00-data-cleaning.R            # Local-first data loading, API fallback
│   ├── 01-immune-cell-scoring.R      # 14 immune cell type scores (incl. CD4)
│   ├── 02-task1-exploration.R        # Cohort table, UpSet, gene coverage
│   ├── 03-q1-pam50-immune.R          # Q1: Heatmap, violins by PAM50 subtype
│   ├── 04-mutation-exploration.R     # Exploratory: mutation-immune associations
│   ├── 05-q2-survival.R              # Q2: Univariate Cox screening
│   ├── 06-q3-within-subtype.R        # Q3: Per-subtype Cox, interaction model
│   └── 07-extension.R                # Extension: immune clustering
├── data/processed/                   # Cleaned datasets (committed)
│   ├── clinical.csv                  # Patient + sample clinical merged
│   ├── mutations.maf                 # Stripped MAF for maftools
│   └── expression_immune_markers.csv # 60 marker genes, log2 intensity
├── references.bib
├── environment.yml                   # Mamba environment spec
├── nature.csl / vancouver.csl        # Citation styles
└── README.md
```

### Data

METABRIC dataset (Curtis et al. 2012) obtained from [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric) via the cBioPortal API. Cleaned datasets are committed in `data/processed/`. The data pipeline is local-first: `scripts/00-data-cleaning.R` uses the committed files if they exist, and only calls the cBioPortal API when they are missing.

### License

Code and analysis: CC BY 4.0. METABRIC data: ODbL v1.0. See `LICENSE.txt` for details.
