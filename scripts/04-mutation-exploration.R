# =============================================================================
# 04-mutation-exploration.R
# Exploratory: Mutation Patterns and Immune Associations — COMPUTATION ONLY
#
# Produces objects: wilcox_results (with adj_effect, adj_p, adj_fdr columns),
#                   sig_genes, mut_mat
#
# Expects: scores_df from 01-immune-cell-scoring.R
#          clinical from 02-task1-exploration.R
# =============================================================================

if (!exists("proc_dir")) proc_dir <- file.path("data", "processed")
if (!exists("mut_file")) mut_file <- file.path(proc_dir, "mutations.maf")

# --- 1. Load MAF and verify Hugo_Symbol --------------------------------------
maf_df <- read.delim(mut_file, comment.char = "#")
stopifnot("Hugo_Symbol" %in% colnames(maf_df))

# --- 2. Build mutation binary matrix ------------------------------------------
mut_binary <- maf_df %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  mutate(mutated = 1L) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = mutated, values_fill = 0L)

common_patients <- intersect(mut_binary$Tumor_Sample_Barcode, scores_df$patient_id)
mut_mat <- mut_binary %>% filter(Tumor_Sample_Barcode %in% common_patients)
scores_q2 <- scores_df %>% filter(patient_id %in% common_patients)
mut_mat <- mut_mat[match(scores_q2$patient_id, mut_mat$Tumor_Sample_Barcode), ]

mut_genes <- setdiff(colnames(mut_mat), "Tumor_Sample_Barcode")
score_cols_q2 <- setdiff(colnames(scores_q2), c("patient_id", "T_cells"))

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

# --- 4. Identify significant genes for dot plot -------------------------------
sig_genes <- wilcox_results %>%
  filter(fdr < 0.05) %>% pull(gene) %>% unique()

# --- 5. PAM50-adjusted linear models -----------------------------------------
# Test whether mutation-immune associations are independent of intrinsic subtype
# Model: immune_score ~ mutation_status + PAM50_subtype
pam50_q2 <- clinical$pam50_subtype[match(scores_q2$patient_id, clinical$patient_id)]
valid_pam50 <- c("Basal-like", "HER2-Enriched", "Luminal A", "Luminal B", "Normal-like")
valid_idx <- pam50_q2 %in% valid_pam50
n_adj_valid <- sum(valid_idx)

adj_results <- list()
for (gene in testable_genes) {
  mut_status <- mut_mat[[gene]][valid_idx]
  if (length(unique(mut_status)) < 2) next
  pam50_status <- factor(pam50_q2[valid_idx])
  for (ct in score_cols_q2) {
    score_vals <- scores_q2[[ct]][valid_idx]
    df <- data.frame(score = score_vals, mutated = factor(mut_status), pam50 = pam50_status)
    fit <- tryCatch(lm(score ~ mutated + pam50, data = df), error = function(e) NULL)
    if (is.null(fit)) next
    coefs <- summary(fit)$coefficients
    if ("mutated1" %in% rownames(coefs)) {
      adj_results[[length(adj_results) + 1]] <- data.frame(
        gene = gene, cell_type = ct,
        adj_effect = coefs["mutated1", "Estimate"],
        adj_p = coefs["mutated1", "Pr(>|t|)"]
      )
    }
  }
}
adj_results <- bind_rows(adj_results)
adj_results$adj_fdr <- p.adjust(adj_results$adj_p, method = "BH")

# Merge adjusted results into wilcox_results
wilcox_results <- wilcox_results %>%
  left_join(adj_results %>% select(gene, cell_type, adj_effect, adj_p, adj_fdr),
            by = c("gene", "cell_type"))

# --- Summary ------------------------------------------------------------------
n_unadj_sig <- sum(wilcox_results$fdr < 0.05, na.rm = TRUE)
n_adj_sig   <- sum(wilcox_results$fdr < 0.05 & wilcox_results$adj_fdr < 0.05, na.rm = TRUE)
n_confounded <- n_unadj_sig - n_adj_sig
message(sprintf(
  "Exploratory: %d patients | %d testable genes (>=10 muts) | %d gene x cell pairs tested | %d FDR < 0.05 unadjusted | %d survive PAM50 adjustment | %d confounded | %d genes with ≥1 sig",
  nrow(scores_q2), length(testable_genes),
  nrow(wilcox_results), n_unadj_sig, n_adj_sig, n_confounded,
  length(sig_genes)
))
