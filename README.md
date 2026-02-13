# Healthcare Data Science Module 3 Assignment

## Immune Cell Infiltration in METABRIC

This is the submission GitHub repository for the Module 3 Assignment for the MSt program in Healthcare Data Science, exploring the infiltration of immune cells of tumour samples of the METABRIC breast cancer dataset using the marker gene scoring from Danaher et al. (2017).

### Environment Setup

Requires [Quarto](https://quarto.org/docs/get-started/) and [Miniforge](https://github.com/conda-forge/miniforge) (mamba).

```bash
# Create environment
mamba env create -f environment.yml

# Activate
mamba activate immfiltration

# Install Bioconductor packages (bioconda dependency chains broken on osx-arm64)
Rscript -e 'if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install(c("ComplexHeatmap", "maftools", "cBioPortalData"))'
```

Bioconductor packages are also auto-installed via `if(!require)` in the qmd setup chunk.

### Rendering

```bash
mamba activate immfiltration
quarto render HDS_03_YuChungYan_2602.qmd
```

The submission qmd sources scripts in order. No need to render intermediate files.

### File Structure

```
├── HDS_03_YuChungYan_2602.qmd            # Final submission — sources scripts/
├── scripts/
│   ├── 00-data-cleaning.R                # Raw cBioPortal → data/processed/
│   ├── 01-danaher-scoring.R              # Compute 13 immune cell type scores
│   ├── 02-task1-exploration.R            # Cohort table, UpSet plot, bias check
│   ├── 03-q1-pam50-immune.R              # Heatmap, violins by PAM50
│   ├── 04-q2-mutations.R                 # Oncoplot, mutation-immune tests
│   ├── 05-q3-survival.R                  # Cox screening, KM curves, forest plot
│   └── 06-extension.R                    # Within-subtype immune stratification
├── data/processed/                       # Cleaned datasets (committed to git)
│   ├── clinical.csv                      # Patient + sample clinical merged
│   ├── mutations.maf                     # Stripped MAF for maftools
│   └── expression_danaher_markers.tsv    # 60 marker genes, log2 intensity
├── environment.yml                       # Mamba environment spec
├── README.md
└── .gitignore
```

### Data

METABRIC dataset (Curtis et al. 2012) obtained from [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric) via the cBioPortal API. Cleaned datasets are in `data/processed/` and committed to the repo. `scripts/00-data-cleaning.R` fetches fresh data from the API on render; if the API is unavailable, the committed `data/processed/` files are used as fallback.
