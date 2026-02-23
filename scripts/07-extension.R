# =============================================================================
# 07-extension.R
# Extension: Immune Clustering for Prognostic Improvement — COMPUTATION ONLY
#
# Question: Do multi-dimensional immune profiles predict survival better than
#           any single cell type?
#
# Produces objects: cluster_scaled, optimal_k, km_final, cluster_centroids,
#                   sig_km_final, sig_optimal_k, comparison_df
#
# Expects: surv_data from 05-q2-survival.R
#          scores_df from 01-immune-cell-scoring.R
#          adj_sig_cells from 06-q3-within-subtype.R
# =============================================================================

# --- 1. Prepare data for clustering ------------------------------------------
score_cols_ext <- setdiff(colnames(scores_df), c("patient_id", "T_cells"))
cluster_data <- surv_data[, score_cols_ext, drop = FALSE]
cluster_scaled <- scale(cluster_data)

# --- 2. Evaluate k = 2, 3, 4 via silhouette width ----------------------------
set.seed(42)
sil_results <- sapply(2:4, function(k) {
  km <- kmeans(cluster_scaled, centers = k, nstart = 25, iter.max = 100)
  ss <- cluster::silhouette(km$cluster, dist(cluster_scaled))
  mean(ss[, "sil_width"])
})
names(sil_results) <- 2:4
optimal_k <- as.integer(names(which.max(sil_results)))

# --- 3. Final clustering with optimal k --------------------------------------
set.seed(42)
km_final <- kmeans(cluster_scaled, centers = optimal_k, nstart = 25, iter.max = 100)
surv_data$immune_cluster <- factor(km_final$cluster)

# --- 4. Cluster centroids (mean z-scores per cell type) -----------------------
cluster_centroids <- km_final$centers # k x 14 matrix, column names preserved

# Label clusters by overall immune level (mean across all cell types)
cluster_means <- rowMeans(cluster_centroids)
cluster_rank <- rank(-cluster_means) # 1 = highest immune
cluster_labels <- paste0(
  "Cluster ", cluster_rank, " (n=",
  table(km_final$cluster)[as.character(seq_len(optimal_k))], ")"
)
levels(surv_data$immune_cluster) <- cluster_labels[as.integer(levels(surv_data$immune_cluster))]
rownames(cluster_centroids) <- cluster_labels

# --- 5. Sig-cells-only clustering ---------------------------------------------
sig_score_cols <- adj_sig_cells$cell_type
sig_cluster_scaled <- scale(surv_data[, sig_score_cols, drop = FALSE])

set.seed(42)
sig_sil_results <- sapply(2:4, function(k) {
  km <- kmeans(sig_cluster_scaled, centers = k, nstart = 25, iter.max = 100)
  ss <- cluster::silhouette(km$cluster, dist(sig_cluster_scaled))
  mean(ss[, "sil_width"])
})
names(sig_sil_results) <- 2:4
sig_optimal_k <- as.integer(names(which.max(sig_sil_results)))

set.seed(42)
sig_km_final <- kmeans(sig_cluster_scaled, centers = sig_optimal_k, nstart = 25, iter.max = 100)
surv_data$sig_cluster <- factor(sig_km_final$cluster)

# Label sig clusters by overall immune level
sig_cl_means <- rowMeans(sig_km_final$centers)
sig_cl_rank <- rank(-sig_cl_means)
sig_cl_labels <- paste0(
  "Sig Cluster ", sig_cl_rank, " (n=",
  table(sig_km_final$cluster)[as.character(seq_len(sig_optimal_k))], ")"
)
levels(surv_data$sig_cluster) <- sig_cl_labels[as.integer(levels(surv_data$sig_cluster))]

# --- 8. Compare: PAM50 baseline vs single cell types vs clustering -----------
cox_pam50_only <- coxph(Surv(os_months, os_event) ~ pam50_subtype, data = surv_data)

# All FDR-significant single cell types
single_results <- lapply(adj_sig_cells$cell_type, function(ct) {
  surv_data$score_z <- scale(surv_data[[ct]])[, 1]
  fit <- coxph(Surv(os_months, os_event) ~ score_z + pam50_subtype, data = surv_data)
  lr <- anova(cox_pam50_only, fit)
  data.frame(
    model = paste0("PAM50 + ", gsub("_", " ", ct)),
    df = 5L,
    concordance = summary(fit)$concordance[1],
    AIC = AIC(fit),
    lr_p = lr[["Pr(>|Chi|)"]][2],
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
surv_data$score_z <- NULL

# Clustering models
cox_cluster_all14 <- coxph(Surv(os_months, os_event) ~ immune_cluster + pam50_subtype, data = surv_data)
cox_cluster_sig <- coxph(Surv(os_months, os_event) ~ sig_cluster + pam50_subtype, data = surv_data)

cluster_results <- data.frame(
  model = c(
    paste0("PAM50 + Cluster (", nrow(adj_sig_cells), " sig, k=", sig_optimal_k, ")"),
    paste0("PAM50 + Cluster (all 14, k=", optimal_k, ")")
  ),
  df = c(4L + sig_optimal_k - 1L, 4L + optimal_k - 1L),
  concordance = c(
    summary(cox_cluster_sig)$concordance[1],
    summary(cox_cluster_all14)$concordance[1]
  ),
  AIC = c(AIC(cox_cluster_sig), AIC(cox_cluster_all14)),
  lr_p = c(
    anova(cox_pam50_only, cox_cluster_sig)[["Pr(>|Chi|)"]][2],
    anova(cox_pam50_only, cox_cluster_all14)[["Pr(>|Chi|)"]][2]
  ),
  stringsAsFactors = FALSE
)

# Baseline row
baseline_row <- data.frame(
  model = "PAM50 subtype only", df = 4L,
  concordance = summary(cox_pam50_only)$concordance[1],
  AIC = AIC(cox_pam50_only), lr_p = NA_real_,
  stringsAsFactors = FALSE
)

comparison_df <- bind_rows(baseline_row, single_results, cluster_results) %>%
  mutate(delta_AIC = AIC - AIC[1])

# --- Summary ------------------------------------------------------------------
message(sprintf(
  "Extension: %d sig single cell types all improve on PAM50 (best: %s, ΔAIC=%.1f) | Cluster (sig k=%d) ΔAIC=%.1f | Cluster (all14 k=%d) ΔAIC=%.1f",
  nrow(adj_sig_cells), gsub("_", " ", adj_sig_cells$cell_type[1]),
  min(comparison_df$delta_AIC, na.rm = TRUE),
  sig_optimal_k, cluster_results$AIC[1] - baseline_row$AIC,
  optimal_k, cluster_results$AIC[2] - baseline_row$AIC
))
