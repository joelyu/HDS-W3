# =============================================================================
# 05-q2-survival.R
# Q2: Immune Scores and Clinical Outcomes — COMPUTATION ONLY
#
# Produces objects: surv_data, cox_results, sig_cells
#
# Expects: scores_df from 01-immune-cell-scoring.R
#          clinical from 02-task1-exploration.R
# =============================================================================

# --- 1. Prepare survival data ------------------------------------------------
surv_data <- scores_df %>%
  inner_join(
    clinical %>%
      select(patient_id, os_months, os_status, pam50_subtype) %>%
      filter(!is.na(os_months) & !is.na(os_status) & os_status != ""),
    by = "patient_id"
  ) %>%
  mutate(
    os_event = as.integer(sub(":.*", "", os_status))
  ) %>%
  filter(!is.na(pam50_subtype) & pam50_subtype != "" &
         pam50_subtype != "NC" & pam50_subtype != "claudin-low")

# --- 2. Univariate Cox screening — 14 cell types (per-SD HR) -----------------
cell_types_q3 <- setdiff(colnames(scores_df), c("patient_id", "T_cells"))

cox_results <- lapply(cell_types_q3, function(ct) {
  surv_data$score_z <- scale(surv_data[[ct]])[, 1]
  fit <- coxph(Surv(os_months, os_event) ~ score_z, data = surv_data)
  s <- summary(fit)
  data.frame(
    cell_type   = ct,
    hr          = s$conf.int[1, 1],
    hr_lower    = s$conf.int[1, 3],
    hr_upper    = s$conf.int[1, 4],
    p_value     = s$coefficients[1, 5],
    concordance = s$concordance[1],
    n           = s$n,
    events      = s$nevent
  )
}) %>% bind_rows() %>%
  mutate(fdr = p.adjust(p_value, method = "BH"))

# Clean up temporary column left by loop
surv_data$score_z <- NULL

# --- 3. Identify significant cell types for KM curves -------------------------
sig_cells <- cox_results %>%
  filter(fdr < 0.05) %>%
  arrange(p_value)

# --- Summary ------------------------------------------------------------------
message(sprintf(
  "Q2: %d patients, %d events | %d/%d cell types FDR < 0.05",
  cox_results$n[1], cox_results$events[1],
  sum(cox_results$fdr < 0.05), nrow(cox_results)
))
