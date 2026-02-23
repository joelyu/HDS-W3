# =============================================================================
# 06-q3-within-subtype.R
# Q3: PAM50-Adjusted Immune Prognostic Analysis — COMPUTATION ONLY
#
# Question: How does accounting for PAM50 subtype change the prognostic
#           landscape of immune cell infiltration?
#
# Produces objects: adjusted_cox, adj_sig_cells, per_subtype_cox,
#                   interaction_tests
#
# Expects: surv_data, cox_results, sig_cells from 05-q2-survival.R
#          scores_df from 01-immune-cell-scoring.R
# =============================================================================

# --- 1. PAM50-adjusted Cox for all 14 cell types -----------------------------
cell_types_q3 <- setdiff(colnames(scores_df), c("patient_id", "T_cells"))

adjusted_cox <- lapply(cell_types_q3, function(ct) {
  surv_data$score_z <- scale(surv_data[[ct]])[, 1]
  fit <- coxph(Surv(os_months, os_event) ~ score_z + pam50_subtype, data = surv_data)
  s <- summary(fit)
  data.frame(
    cell_type = ct,
    hr = s$conf.int["score_z", 1],
    hr_lower = s$conf.int["score_z", 3],
    hr_upper = s$conf.int["score_z", 4],
    p_value = s$coefficients["score_z", 5],
    concordance = s$concordance[1],
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

surv_data$score_z <- NULL

# --- 2. Identify FDR-significant adjusted cell types --------------------------
adj_sig_cells <- adjusted_cox %>% filter(fdr < 0.05)

# --- 3. Per-subtype Cox for significant cell types ----------------------------
subtypes_q3 <- c("Basal-like", "HER2-Enriched", "Luminal A", "Luminal B", "Normal-like")
sig_ct_names <- adj_sig_cells$cell_type

# Pre-compute cohort-wide mean/SD for each sig cell type
cohort_stats <- data.frame(
  cell_type = sig_ct_names,
  ct_mean = sapply(sig_ct_names, function(ct) mean(surv_data[[ct]], na.rm = TRUE)),
  ct_sd = sapply(sig_ct_names, function(ct) sd(surv_data[[ct]], na.rm = TRUE)),
  stringsAsFactors = FALSE
)

per_subtype_cox <- lapply(sig_ct_names, function(ct) {
  ct_mean <- cohort_stats$ct_mean[cohort_stats$cell_type == ct]
  ct_sd <- cohort_stats$ct_sd[cohort_stats$cell_type == ct]
  lapply(subtypes_q3, function(st) {
    st_data <- surv_data %>% filter(pam50_subtype == st)
    st_data$score_z <- (st_data[[ct]] - ct_mean) / ct_sd
    fit <- coxph(Surv(os_months, os_event) ~ score_z, data = st_data)
    s <- summary(fit)
    data.frame(
      cell_type = ct, subtype = st,
      n = nrow(st_data), events = s$nevent,
      hr = s$conf.int[1, 1], hr_lower = s$conf.int[1, 3], hr_upper = s$conf.int[1, 4],
      p_value = s$coefficients[1, 5],
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
}) %>% bind_rows()

# --- 4. Interaction tests for significant cell types --------------------------
interaction_tests <- lapply(sig_ct_names, function(ct) {
  surv_data$score_z_int <- scale(surv_data[[ct]])[, 1]
  fit_main <- coxph(Surv(os_months, os_event) ~ score_z_int + pam50_subtype, data = surv_data)
  fit_int <- coxph(Surv(os_months, os_event) ~ score_z_int * pam50_subtype, data = surv_data)
  lr <- anova(fit_main, fit_int)
  data.frame(
    cell_type = ct,
    interaction_p = lr[["Pr(>|Chi|)"]][2],
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  mutate(interaction_fdr = p.adjust(interaction_p, method = "BH"))

surv_data$score_z_int <- NULL

# --- Summary ------------------------------------------------------------------
message(sprintf(
  "Q3: %d/%d cell types FDR < 0.05 after PAM50 adjustment (%s)",
  nrow(adj_sig_cells), length(cell_types_q3),
  paste(gsub("_", " ", adj_sig_cells$cell_type), collapse = ", ")
))
