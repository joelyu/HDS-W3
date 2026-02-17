# =============================================================================
# 02-task1-exploration.R
# Data Exploration — COMPUTATION ONLY
#
# Produces objects: clinical, membership, n_expr_pam50, n_expr_mut,
#                   tbl1, coverage_df
#
# Expects: scores_df, expr_mat, immune_markers from 01-immune-cell-scoring.R
# =============================================================================

if (!exists("proc_dir")) {
  proc_dir <- file.path("data", "processed")
}
if (!exists("clin_file")) {
  clin_file <- file.path(proc_dir, "clinical.csv")
}
if (!exists("mut_file")) {
  mut_file <- file.path(proc_dir, "mutations.maf")
}

# --- Load clinical data (already cleaned by 00-data-cleaning.R) ---
clinical <- read.csv(clin_file)
clin_id_col <- "patient_id"

# --- Expand abbreviated PAM50 labels to full names ---
pam50_remap <- c(
  "Basal" = "Basal-like", "Her2" = "HER2-Enriched",
  "LumA" = "Luminal A", "LumB" = "Luminal B", "Normal" = "Normal-like"
)
clinical$pam50_subtype <- ifelse(
  clinical$pam50_subtype %in% names(pam50_remap),
  pam50_remap[clinical$pam50_subtype],
  clinical$pam50_subtype
)

# --- Get patient IDs per data type ---
expr_patients <- colnames(expr_mat) # from 01-immune-cell-scoring.R

maf <- read.delim(mut_file, comment.char = "#")
mut_patients <- unique(maf$Tumor_Sample_Barcode)

# PAM50 set for UpSet: includes claudin-low (raw data availability view).
# Claudin-low is excluded from downstream analysis (Table 1, Q1-Q3) — see filters below.
pam50_patients <- clinical[[clin_id_col]][
  !is.na(clinical$pam50_subtype) & clinical$pam50_subtype != ""
]

# --- Data overlap ---
all_patients <- unique(c(
  as.character(clinical[[clin_id_col]]),
  expr_patients,
  as.character(mut_patients)
))

membership <- data.frame(
  Patient = all_patients,
  Clinical = as.integer(
    all_patients %in% as.character(clinical[[clin_id_col]])
  ),
  Expression = as.integer(all_patients %in% expr_patients),
  Mutations = as.integer(all_patients %in% as.character(mut_patients)),
  PAM50 = as.integer(all_patients %in% as.character(pam50_patients)),
  stringsAsFactors = FALSE
)

n_expr_pam50 <- sum(membership$Expression == 1 & membership$PAM50 == 1)
n_expr_mut <- sum(membership$Expression == 1 & membership$Mutations == 1)

# --- Table 1: cohort characteristics by PAM50 ---
clinical$has_expression <- as.character(clinical[[clin_id_col]]) %in%
  expr_patients

table_data <- clinical %>%
  filter(
    !is.na(pam50_subtype) &
      pam50_subtype != "" &
      pam50_subtype != "NC" &
      pam50_subtype != "claudin-low"
  ) %>%
  mutate(
    deceased = factor(
      ifelse(sub(":.*", "", os_status) == "1", "Yes", "No"),
      levels = c("No", "Yes")
    )
  )

table_vars <- c(
  "age_at_diagnosis",
  "er_status",
  "grade",
  "cellularity",
  "os_months",
  "deceased"
)
avail_vars <- intersect(table_vars, colnames(table_data))

tbl1 <- table_data %>%
  select(all_of(c("pam50_subtype", avail_vars))) %>%
  tbl_summary(
    by = pam50_subtype,
    missing = "ifany",
    label = list(
      age_at_diagnosis ~ "Age at diagnosis",
      er_status ~ "ER status",
      grade ~ "Tumour Grade",
      cellularity ~ "Cellularity",
      os_months ~ "Follow-up (months)",
      deceased ~ "Deceased"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    )
  ) %>%
  add_overall() %>%
  add_p() %>%
  bold_labels()

# --- Danaher marker gene coverage ---
coverage_df <- tibble(
  Cell_type = rep(names(immune_markers), sapply(immune_markers, length)),
  Gene = unlist(immune_markers)
) %>%
  mutate(In_METABRIC = Gene %in% rownames(expr_mat))

# --- Summary ------------------------------------------------------------------
message(sprintf(
  "Cohort: %d clinical | %d expression | %d mutations | %d PAM50 | Working: Q1/Q3 %d, Q2 %d",
  nrow(clinical),
  length(expr_patients),
  length(mut_patients),
  length(pam50_patients),
  n_expr_pam50,
  n_expr_mut
))
