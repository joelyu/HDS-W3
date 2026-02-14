# =============================================================================
# 02-task1-exploration.R
# Task 1: Data Exploration — cohort table, data completeness, missing data bias
#
# Expects: proc_dir, clin_file, mut_file defined by parent qmd
#          scores_df, expr_mat, immune_markers from 01-immune-cell-scoring.R
# Produces: Table 1, UpSet plot, bias check table, gene coverage table
# =============================================================================

# --- Load clinical data (already cleaned by 00-data-cleaning.R) ---
clinical <- read.csv(clin_file)
clin_id_col <- "patient_id"

# --- Get patient IDs per data type ---
expr_patients <- colnames(expr_mat)  # from 01-immune-cell-scoring.R

maf <- read.delim(mut_file, comment.char = "#")
mut_patients <- unique(maf$Tumor_Sample_Barcode)

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
  Patient    = all_patients,
  Clinical   = as.integer(all_patients %in% as.character(clinical[[clin_id_col]])),
  Expression = as.integer(all_patients %in% expr_patients),
  Mutations  = as.integer(all_patients %in% as.character(mut_patients)),
  PAM50      = as.integer(all_patients %in% as.character(pam50_patients)),
  stringsAsFactors = FALSE
)

n_expr_pam50 <- sum(membership$Expression == 1 & membership$PAM50 == 1)
n_expr_mut   <- sum(membership$Expression == 1 & membership$Mutations == 1)

# --- UpSet plot: data completeness ---
upset(
  membership[, -1],
  sets = c("Clinical", "Expression", "Mutations", "PAM50"),
  order.by = "freq",
  main.bar.color = "#4A708B",
  sets.bar.color = "#4A708B",
  text.scale = c(1.5, 1.3, 1.2, 1.2, 1.5, 1.2),
  mb.ratio = c(0.6, 0.4)
)

# --- Table 1: cohort characteristics by PAM50 ---
clinical$has_expression <- as.character(clinical[[clin_id_col]]) %in% expr_patients

table_data <- clinical %>%
  filter(!is.na(pam50_subtype) & pam50_subtype != "")

table_vars <- c("age_at_diagnosis", "menopausal_state", "er_status",
                "her2_status", "grade", "tumor_stage", "has_expression")
avail_vars <- intersect(table_vars, colnames(table_data))

tbl1 <- table_data %>%
  select(all_of(c("pam50_subtype", avail_vars))) %>%
  tbl_summary(
    by = pam50_subtype,
    missing = "ifany",
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    )
  ) %>%
  add_overall() %>%
  add_p() %>%
  bold_labels()
tbl1

# --- Missing data bias check ---
clinical$expr_group <- ifelse(
  clinical$has_expression,
  "Has expression", "No expression"
)

bias_vars <- c("age_at_diagnosis", "er_status", "her2_status",
               "grade", "tumor_stage", "pam50_subtype")
bias_avail <- intersect(bias_vars, colnames(clinical))

tbl_bias <- clinical %>%
  select(all_of(c("expr_group", bias_avail))) %>%
  tbl_summary(
    by = expr_group,
    missing = "ifany",
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    )
  ) %>%
  add_p() %>%
  bold_labels()
tbl_bias

# --- Danaher marker gene coverage ---
coverage_df <- tibble(
  Cell_type = rep(names(immune_markers), sapply(immune_markers, length)),
  Gene = unlist(immune_markers)
) %>%
  mutate(In_METABRIC = Gene %in% rownames(expr_mat))

coverage_df %>%
  group_by(Cell_type) %>%
  summarise(
    Markers = n(),
    Found = sum(In_METABRIC),
    Missing = paste(Gene[!In_METABRIC], collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(Missing = ifelse(Missing == "", "\u2014", Missing)) %>%
  knitr::kable(caption = "Immune marker gene coverage in METABRIC (Danaher et al. 2017, after HGNC alias mapping)")

# --- Summary ------------------------------------------------------------------
message(sprintf(
  "Cohort: %d clinical | %d expression | %d mutations | %d PAM50 | Working: Q1/Q3 %d, Q2 %d",
  nrow(clinical), length(expr_patients), length(mut_patients), length(pam50_patients),
  n_expr_pam50, n_expr_mut
))
