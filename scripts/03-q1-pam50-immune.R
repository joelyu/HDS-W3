# =============================================================================
# 03-q1-pam50-immune.R
# Q1: Immune Cell Expression Across PAM50 Subtypes — COMPUTATION ONLY
#
# Produces objects: scores_clin, cell_types, heatmap_z, pam50_vals,
#                   scores_long, kw_results, pw_results, cld_df, pam50_colours
#
# Expects: scores_df, scores_mat, immune_markers from 01-immune-cell-scoring.R
#          clinical from 02-task1-exploration.R
# =============================================================================

# --- PAM50 colour palette (ColorBrewer Set2, used by qmd for Q1-Q3 plots) ----
pam50_colours <- setNames(
  RColorBrewer::brewer.pal(5, "Set2"),
  c("Basal-like", "HER2-Enriched", "Luminal A", "Luminal B", "Normal-like")
)

# --- 1. Merge scores with PAM50 -----------------------------------------------
scores_clin <- scores_df %>%
  inner_join(
    clinical %>% select(patient_id, pam50_subtype),
    by = "patient_id"
  ) %>%
  filter(!is.na(pam50_subtype) & pam50_subtype != "" &
    pam50_subtype != "NC" & pam50_subtype != "claudin-low")

# --- 2. Heatmap data: 14 cell types (drop T-cells parent) --------------------
cell_types <- setdiff(rownames(scores_mat), "T_cells")

pam50_map <- setNames(clinical$pam50_subtype, clinical$patient_id)
heatmap_patients <- intersect(
  colnames(scores_mat),
  names(pam50_map[!is.na(pam50_map) & pam50_map != "" &
    pam50_map != "NC" & pam50_map != "claudin-low"])
)
heatmap_mat <- scores_mat[cell_types, heatmap_patients]
heatmap_z <- t(scale(t(heatmap_mat)))
pam50_vals <- pam50_map[heatmap_patients]

# --- 3. Long format for violin plots ------------------------------------------
score_cols <- setdiff(colnames(scores_df), c("patient_id", "T_cells"))

scores_long <- scores_clin %>%
  select(patient_id, pam50_subtype, all_of(score_cols)) %>%
  pivot_longer(cols = all_of(score_cols), names_to = "cell_type", values_to = "score") %>%
  mutate(cell_type = gsub("_", " ", cell_type))

# --- 4. Kruskal-Wallis per cell type with BH FDR correction ------------------
kw_results <- scores_long %>%
  group_by(cell_type) %>%
  summarise(
    {
      kt <- kruskal.test(score ~ pam50_subtype)
      data.frame(kw_stat = kt$statistic, kw_p = kt$p.value)
    },
    .groups = "drop"
  ) %>%
  mutate(kw_fdr = p.adjust(kw_p, method = "BH"))

# --- 5. Pairwise Wilcoxon (for in-text reference) ----------------------------
pw_results <- scores_long %>%
  group_by(cell_type) %>%
  summarise(
    pw = list(pairwise.wilcox.test(score, pam50_subtype, p.adjust.method = "BH")),
    .groups = "drop"
  )

# --- 6. Compact letter display from pairwise Wilcoxon ------------------------
cld_df <- purrr::map(seq_len(nrow(pw_results)), function(i) {
  ct <- pw_results$cell_type[i]
  pmat <- pw_results$pw[[i]]$p.value

  # Replace hyphens in group names — multcompLetters splits on "-"
  all_groups <- unique(c(colnames(pmat), rownames(pmat)))
  safe_map <- setNames(gsub("-", ".", all_groups), all_groups)
  rev_map <- setNames(names(safe_map), safe_map)
  rownames(pmat) <- safe_map[rownames(pmat)]
  colnames(pmat) <- safe_map[colnames(pmat)]

  # Extract pairwise p-values as named vector
  pvals <- c()
  for (r in seq_len(nrow(pmat))) {
    for (cc in seq_len(ncol(pmat))) {
      if (!is.na(pmat[r, cc])) {
        pvals[paste(rownames(pmat)[r], colnames(pmat)[cc], sep = "-")] <- pmat[r, cc]
      }
    }
  }

  letters <- multcompView::multcompLetters(pvals, threshold = 0.05)$Letters

  data.frame(
    cell_type = ct,
    pam50_subtype = rev_map[names(letters)],
    cld_letter = unname(letters),
    stringsAsFactors = FALSE
  )
}) %>% dplyr::bind_rows()

# --- Summary ------------------------------------------------------------------
message(sprintf(
  "Q1: %d patients x %d cell types | %d/%d significant at FDR < 0.05",
  nrow(scores_clin), length(cell_types),
  sum(kw_results$kw_fdr < 0.05), nrow(kw_results)
))
