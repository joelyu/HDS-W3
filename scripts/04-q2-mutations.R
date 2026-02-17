# =============================================================================
# 04-q2-mutations.R
# Q2: Mutation Patterns and Immune Associations — COMPUTATION ONLY
#
# Produces objects: maf_obj, clin_for_maf, wilcox_results, sig_genes,
#                   effect_mat, fdr_mat, testable_genes,
#                   scores_q2, mut_mat, score_cols_q2
#
# Expects: scores_df from 01-immune-cell-scoring.R
#          clinical from 02-task1-exploration.R
# =============================================================================

if (!exists("proc_dir")) proc_dir <- file.path("data", "processed")
if (!exists("mut_file")) mut_file <- file.path(proc_dir, "mutations.maf")

# --- 1. Load MAF and verify Hugo_Symbol --------------------------------------
maf_df <- read.delim(mut_file, comment.char = "#")
stopifnot("Hugo_Symbol" %in% colnames(maf_df))

# Clinical annotation for maftools (needs Tumor_Sample_Barcode column)
clin_for_maf <- clinical %>%
  select(patient_id, pam50_subtype) %>%
  dplyr::rename(Tumor_Sample_Barcode = patient_id) %>%
  filter(!is.na(pam50_subtype) & pam50_subtype != "" &
         pam50_subtype != "NC" & pam50_subtype != "claudin-low")

maf_obj <- read.maf(maf = mut_file, clinicalData = clin_for_maf, verbose = FALSE)

# --- 2. Build mutation binary matrix ------------------------------------------
mut_binary <- maf_df %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  mutate(mutated = 1L) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = mutated, values_fill = 0L)

common_patients <- intersect(mut_binary$Tumor_Sample_Barcode, scores_df$PATIENT_ID)
mut_mat <- mut_binary %>% filter(Tumor_Sample_Barcode %in% common_patients)
scores_q2 <- scores_df %>% filter(PATIENT_ID %in% common_patients)
mut_mat <- mut_mat[match(scores_q2$PATIENT_ID, mut_mat$Tumor_Sample_Barcode), ]

mut_genes <- setdiff(colnames(mut_mat), "Tumor_Sample_Barcode")
score_cols_q2 <- setdiff(colnames(scores_q2), c("PATIENT_ID", "T_cells"))

# --- 3. Wilcoxon: mutated vs WT for each gene x cell type --------------------
gene_counts <- colSums(mut_mat[, mut_genes, drop = FALSE])
testable_genes <- names(gene_counts[gene_counts >= 10])

wilcox_results <- list()
for (gene in testable_genes) {
  mut_idx <- mut_mat[[gene]] == 1
  for (ct in score_cols_q2) {
    mut_vals <- scores_q2[[ct]][mut_idx]
    wt_vals  <- scores_q2[[ct]][!mut_idx]
    wt <- wilcox.test(mut_vals, wt_vals, exact = FALSE)
    wilcox_results[[length(wilcox_results) + 1]] <- data.frame(
      gene = gene, cell_type = ct,
      n_mut = sum(mut_idx), n_wt = sum(!mut_idx),
      p_value = wt$p.value,
      effect  = median(mut_vals, na.rm = TRUE) - median(wt_vals, na.rm = TRUE)
    )
  }
}
wilcox_results <- bind_rows(wilcox_results)
wilcox_results$fdr <- p.adjust(wilcox_results$p_value, method = "BH")

# --- 4. Effect size matrices for heatmap --------------------------------------
sig_genes <- wilcox_results %>%
  filter(fdr < 0.05) %>% pull(gene) %>% unique()

if (length(sig_genes) > 0) {
  # Show all significant genes if ≤30, otherwise top 30 by most significant pair
  max_genes <- 30
  top_sig <- wilcox_results %>%
    filter(gene %in% sig_genes) %>%
    group_by(gene) %>%
    summarise(min_fdr = min(fdr), .groups = "drop") %>%
    arrange(min_fdr)
  if (nrow(top_sig) > max_genes) top_sig <- head(top_sig, max_genes)
  plot_genes <- top_sig$gene

  effect_mat <- wilcox_results %>%
    filter(gene %in% plot_genes) %>%
    select(gene, cell_type, effect) %>%
    pivot_wider(names_from = cell_type, values_from = effect) %>%
    column_to_rownames("gene") %>% as.matrix()

  fdr_mat <- wilcox_results %>%
    filter(gene %in% plot_genes) %>%
    select(gene, cell_type, fdr) %>%
    pivot_wider(names_from = cell_type, values_from = fdr) %>%
    column_to_rownames("gene") %>% as.matrix()

  colnames(effect_mat) <- gsub("_", " ", colnames(effect_mat))
  colnames(fdr_mat)    <- gsub("_", " ", colnames(fdr_mat))
} else {
  effect_mat <- NULL
  fdr_mat <- NULL
}

# --- Summary ------------------------------------------------------------------
message(sprintf(
  "Q2: %d patients | %d testable genes (>=10 muts) | %d gene x cell pairs tested | %d FDR < 0.05 | %d genes with ≥1 sig association | showing %s",
  nrow(scores_q2), length(testable_genes),
  nrow(wilcox_results), sum(wilcox_results$fdr < 0.05),
  length(sig_genes),
  if (length(sig_genes) > 0) paste(length(plot_genes), "genes in heatmap") else "table fallback"
))
