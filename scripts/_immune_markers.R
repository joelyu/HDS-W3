# =============================================================================
# _immune_markers.R
# Single source of truth: Danaher et al. 2017 marker gene list (Table 1)
#
# 60 genes across 14 cell types. CD4 T cells derived in 01-immune-cell-scoring.R
# as T-cell score minus CD8 score (not a separate gene set).
#
# Gene symbol verification:
#   All 60 symbols checked against the HGNC REST API (rest.genenames.org)
#   using /fetch/symbol, /search/prev_symbol, and /search/alias_symbol endpoints.
#   - KIAA0125 is a previous symbol for FAM30A (B-cells) — resolved to current name
#   - All other 59 symbols confirmed as current HGNC-approved names
#
# Platform coverage (Illumina HT-12 v3 microarray):
#   58/60 genes present in METABRIC expression data.
#   - TPSB2 (Mast cells, 1 of 5 markers) — absent from platform
#   - XCL2 (NK cells, 1 of 3 markers) — absent from platform
#   Both retained in this list so the coverage table (Table S1) reports 58/60.
#   Scoring proceeds with available markers (4/5 Mast, 2/3 NK).
# =============================================================================

immune_markers <- list(
  "B-cells"         = c("BLK", "CD19", "FCRL2", "MS4A1", "FAM30A",
                         "TNFRSF17", "TCL1A", "SPIB", "PNOC"),
  "CD45"            = c("PTPRC"),
  "Cytotoxic cells" = c("PRF1", "GZMA", "GZMB", "NKG7", "GZMH",
                         "KLRK1", "KLRB1", "KLRD1", "CTSW", "GNLY"),
  "DC"              = c("CCL13", "CD209", "HSD11B1"),
  "Exhausted CD8"   = c("LAG3", "CD244", "EOMES", "PTGER4"),
  "Macrophages"     = c("CD68", "CD84", "CD163", "MS4A4A"),
  "Mast cells"      = c("TPSAB1", "TPSB2", "CPA3", "MS4A2", "HDC"),
  "Neutrophils"     = c("FPR1", "SIGLEC5", "CSF3R", "FCAR", "FCGR3B",
                         "CEACAM3", "S100A12"),
  "NK CD56dim"      = c("KIR2DL3", "KIR3DL1", "KIR3DL2", "IL21R"),
  "NK cells"        = c("XCL1", "XCL2", "NCR1"),
  "T-cells"         = c("CD6", "CD3D", "CD3E", "SH2D1A", "TRAT1", "CD3G"),
  "Th1 cells"       = c("TBX21"),
  "Treg"            = c("FOXP3"),
  "CD8 T cells"     = c("CD8A", "CD8B")
)
