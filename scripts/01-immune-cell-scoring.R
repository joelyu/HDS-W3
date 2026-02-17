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
# Expects: proc_dir, expr_file (falls back to defaults if not set by parent qmd)
# =============================================================================

if (!exists("proc_dir"))  proc_dir  <- file.path("data", "processed")
if (!exists("expr_file")) expr_file <- file.path(proc_dir, "expression_immune_markers.csv")

# --- 14 cell types, 62 marker genes (single source of truth) -----------------
# CD4 derived below as T-cell score minus CD8 score (r=0.65 vs flow cytometry)
source("scripts/_immune_markers.R")

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

# Clean column names (spaces/hyphens → underscores for safe read.csv round-trips)
colnames(scores_df) <- gsub("[ -]", "_", colnames(scores_df))

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
