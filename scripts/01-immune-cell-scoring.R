# =============================================================================
# 01-immune-cell-scoring.R
# Compute immune cell type scores using Danaher et al. (2017) marker genes
#
# Method: mean log2 expression of cell-type-specific marker genes per patient
# Reference: Danaher P et al. J Immunother Cancer 5, 18 (2017)
#
# Input:  data/processed/expression_immune_markers.csv (from 00-data-cleaning.R)
# Output: data/processed/immune_scores.csv
#         Objects in memory: scores_df, scores_mat, expr_mat, immune_markers
#
# Expects: proc_dir, expr_file defined by parent qmd
# =============================================================================

# --- 14 cell types, 60 marker genes (Danaher et al. 2017, Table 1) ---
# CD4 derived as T-cell score minus CD8 score (r=0.65 vs flow cytometry, Table 2)
# Works with log2 intensity: subtraction ≈ log2(T-cell / CD8) ratio
# Gene alias: KIAA0125 → FAM30A (verified via HGNC)
# Platform-absent: TPSB2, XCL2 (not on Illumina HT-12 v3)
immune_markers <- list(
  "B-cells"         = c("BLK", "CD19", "FCRL2", "MS4A1", "FAM30A",
                         "TNFRSF17", "TCL1A", "SPIB", "PNOC"),
  "CD45"            = c("PTPRC"),
  "Cytotoxic cells" = c("PRF1", "GZMA", "GZMB", "NKG7", "GZMH",
                         "KLRK1", "KLRB1", "KLRD1", "CTSW", "GNLY"),
  "DC"              = c("CCL13", "CD209", "HSD11B1"),
  "Exhausted CD8"   = c("LAG3", "CD244", "EOMES", "PTGER4"),
  "Macrophages"     = c("CD68", "CD84", "CD163", "MS4A4A"),
  "Mast cells"      = c("TPSAB1", "CPA3", "MS4A2", "HDC"),
  "Neutrophils"     = c("S100A9", "S100A8", "CEACAM3", "SPI1", "FPR1",
                         "SIGLEC5", "CSF3R", "FCAR", "FCGR3B"),
  "NK CD56dim"      = c("KIR2DL3", "KIR3DL1", "KIR3DL2", "IL21R"),
  "NK cells"        = c("XCL1", "NCR1"),
  "T-cells"         = c("CD6", "CD3D", "CD3E", "SH2D1A", "TRAT1", "CD3G"),
  "Th1 cells"       = c("TBX21"),
  "Treg"            = c("FOXP3"),
  "CD8 T cells"     = c("CD8A", "CD8B")
)

# --- Load expression data (pre-processed marker gene subset, log2 intensity) ---
expr_raw <- read.csv(expr_file, check.names = FALSE)
gene_col <- colnames(expr_raw)[1]

# Handle KIAA0125/FAM30A alias — API may return either name
if ("KIAA0125" %in% expr_raw[[gene_col]] && !"FAM30A" %in% expr_raw[[gene_col]]) {
  expr_raw[[gene_col]][expr_raw[[gene_col]] == "KIAA0125"] <- "FAM30A"
}

expr_mat <- expr_raw %>%
  column_to_rownames(var = gene_col) %>%
  as.matrix()

# --- Verify marker gene coverage ---
all_markers <- unique(unlist(immune_markers))
found   <- all_markers[all_markers %in% rownames(expr_mat)]
missing <- all_markers[!all_markers %in% rownames(expr_mat)]

# --- Compute scores: mean log2 expression of available markers per patient ---
compute_cell_score <- function(expr_matrix, genes) {
  available <- genes[genes %in% rownames(expr_matrix)]
  if (length(available) == 0) return(rep(NA_real_, ncol(expr_matrix)))
  if (length(available) == 1) return(expr_matrix[available, ])
  colMeans(expr_matrix[available, , drop = FALSE], na.rm = TRUE)
}

scores_list <- lapply(immune_markers, compute_cell_score, expr_matrix = expr_mat)

# CD4: T-cell score minus CD8 score (Danaher et al. 2017, Table 1)
scores_list[["CD4 T cells"]] <- scores_list[["T-cells"]] - scores_list[["CD8 T cells"]]

scores_mat <- do.call(rbind, scores_list)
colnames(scores_mat) <- colnames(expr_mat)

# Transpose to patient-rows format
scores_df <- scores_mat %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("PATIENT_ID")

# Clean column names (spaces → underscores)
colnames(scores_df) <- gsub(" ", "_", colnames(scores_df))

# --- Save to CSV ---
write.csv(scores_df, file.path(proc_dir, "immune_scores.csv"), row.names = FALSE)

# --- Summary ------------------------------------------------------------------
message(sprintf(
  "Immune scoring: %d genes x %d patients | Coverage: %d/%d markers | %d scores (14 cell types + T-cells parent)%s",
  nrow(expr_mat), ncol(expr_mat),
  length(found), length(all_markers),
  nrow(scores_mat),
  if (length(missing) > 0) paste0(" | Missing: ", paste(missing, collapse = ", ")) else ""
))
