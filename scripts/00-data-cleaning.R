# =============================================================================
# 00-data-cleaning.R
# Fetch METABRIC data from cBioPortal API → data/processed/
#
# Strategy:
#   1. Try cBioPortal API (pulls only what we need — no bulk download)
#   2. If API fails, use existing data/processed/ files (committed to git)
#
# Output:  data/processed/clinical.csv
#          data/processed/mutations.maf
#          data/processed/expression_immune_markers.csv   (log2 intensity)
#
# Expects: proc_dir defined by parent qmd
# =============================================================================

if (!exists("proc_dir")) proc_dir <- file.path("data", "processed")
dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)

# --- Danaher marker genes (single source of truth in _immune_markers.R) ------
source("scripts/_immune_markers.R")
danaher_genes <- unique(c(unlist(immune_markers), "KIAA0125"))  # KIAA0125 = API fallback alias for FAM30A

STUDY_ID <- "brca_metabric"
.log <- character()

# --- Skip API if processed data already exists --------------------------------
proc_files <- c("clinical.csv", "mutations.maf", "expression_immune_markers.csv")
all_exist  <- all(file.exists(file.path(proc_dir, proc_files)))

if (all_exist) {
  .log <- c(.log, sprintf("Source: existing files in %s (delete to force re-download)", proc_dir))
  api_success <- TRUE
} else {

# --- Try cBioPortal API ------------------------------------------------------
api_success <- tryCatch({
  if (!require(cBioPortalData, quietly = TRUE)) {
    BiocManager::install("cBioPortalData", ask = FALSE, update = FALSE)
    library(cBioPortalData)
  }

  cbio <- cBioPortal()
  .log <- c(.log, "Source: cBioPortal API")

  # ---- 1. Clinical data -----------------------------------------------------
  clin_raw <- clinicalData(cbio, studyId = STUDY_ID)

  clinical_clean <- clin_raw %>%
    transmute(
      patient_id       = patientId,
      age_at_diagnosis = as.numeric(AGE_AT_DIAGNOSIS),
      sex              = SEX,
      menopausal_state = INFERRED_MENOPAUSAL_STATE,
      er_status        = ER_IHC,
      her2_status      = HER2_SNP6,
      pr_status        = if ("PR_STATUS" %in% names(.)) PR_STATUS else NA_character_,
      grade            = if ("GRADE" %in% names(.)) GRADE else NA_character_,
      tumor_stage      = if ("TUMOR_STAGE" %in% names(.)) as.character(TUMOR_STAGE) else NA_character_,
      tumor_size       = if ("TUMOR_SIZE" %in% names(.)) as.numeric(TUMOR_SIZE) else NA_real_,
      lymph_nodes_pos  = as.numeric(LYMPH_NODES_EXAMINED_POSITIVE),
      npi              = as.numeric(NPI),
      cellularity      = CELLULARITY,
      chemotherapy     = CHEMOTHERAPY,
      hormone_therapy  = HORMONE_THERAPY,
      radiotherapy     = RADIO_THERAPY,
      breast_surgery   = BREAST_SURGERY,
      histological_subtype = HISTOLOGICAL_SUBTYPE,
      # NOTE: CLAUDIN_SUBTYPE is cBioPortal's "PAM50 + Claudin-low subtype" field —
      # a combined classification whose provenance is undocumented. We retain the
      # raw values here; claudin-low is excluded in downstream analysis scripts.
      pam50_subtype    = CLAUDIN_SUBTYPE,
      intclust_subtype = INTCLUST,
      cohort           = COHORT,
      os_months        = as.numeric(OS_MONTHS),
      os_status        = OS_STATUS,
      rfs_months       = as.numeric(RFS_MONTHS),
      rfs_status       = RFS_STATUS,
      tmb_nonsynonymous = if ("TMB_NONSYNONYMOUS" %in% names(.)) as.numeric(TMB_NONSYNONYMOUS) else NA_real_
    ) %>%
    mutate(across(where(is.character), ~ gsub("Positve", "Positive", .x)))

  write.csv(clinical_clean, file.path(proc_dir, "clinical.csv"), row.names = FALSE)
  .log <- c(.log, sprintf("  clinical.csv: %d patients x %d cols", nrow(clinical_clean), ncol(clinical_clean)))

  # ---- 2. Mutation data -----------------------------------------------------
  # mutationData() requires entrezGeneIds — use REST API /fetch endpoint
  mut_resp <- httr::POST(
    sprintf("https://www.cbioportal.org/api/molecular-profiles/%s_mutations/mutations/fetch",
            STUDY_ID),
    body = jsonlite::toJSON(
      list(sampleListId = paste0(STUDY_ID, "_all")),
      auto_unbox = TRUE
    ),
    httr::content_type_json(),
    httr::add_headers(accept = "application/json")
  )
  httr::stop_for_status(mut_resp, "fetch mutation data from cBioPortal")
  mut_raw <- jsonlite::fromJSON(
    httr::content(mut_resp, as = "text", encoding = "UTF-8"),
    flatten = TRUE
  )

  # API returns entrezGeneId but not Hugo symbols — look them up separately
  entrez_col <- grep("entrezGeneId", colnames(mut_raw), value = TRUE)[1]
  if (is.na(entrez_col)) stop("No entrezGeneId column found in mutation data")

  gene_ids <- unique(as.integer(mut_raw[[entrez_col]]))
  gene_resp <- httr::POST(
    "https://www.cbioportal.org/api/genes/fetch",
    body = jsonlite::toJSON(gene_ids),
    httr::content_type_json(),
    httr::add_headers(accept = "application/json")
  )
  httr::stop_for_status(gene_resp, "fetch gene symbols from cBioPortal")
  gene_map <- jsonlite::fromJSON(
    httr::content(gene_resp, as = "text", encoding = "UTF-8")
  )
  mut_raw$Hugo_Symbol <- gene_map$hugoGeneSymbol[
    match(as.integer(mut_raw[[entrez_col]]), gene_map$entrezGeneId)
  ]

  # Map REST API column names to standard MAF names for maftools
  maf_rename <- c(
    chr = "Chromosome",
    startPosition = "Start_Position", endPosition = "End_Position",
    referenceAllele = "Reference_Allele", variantAllele = "Tumor_Seq_Allele2",
    mutationType = "Variant_Classification", variantType = "Variant_Type",
    sampleId = "Tumor_Sample_Barcode", proteinChange = "HGVSp_Short"
  )

  avail_cols <- intersect(names(maf_rename), colnames(mut_raw))
  maf_clean <- mut_raw[, c("Hugo_Symbol", avail_cols), drop = FALSE]
  for (old_name in avail_cols) {
    colnames(maf_clean)[colnames(maf_clean) == old_name] <- maf_rename[old_name]
  }

  # Verify Hugo_Symbol is populated
  if (all(is.na(maf_clean$Hugo_Symbol))) {
    stop("Hugo_Symbol lookup failed — all values are NA")
  }

  write.table(maf_clean, file.path(proc_dir, "mutations.maf"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  .log <- c(.log, sprintf("  mutations.maf: %d mutations x %d cols", nrow(maf_clean), ncol(maf_clean)))

  # ---- 3. Expression: log2 intensity for Danaher markers --------------------
  expr_raw <- getDataByGenes(
    cbio,
    studyId    = STUDY_ID,
    genes      = danaher_genes,
    by         = "hugoGeneSymbol",
    molecularProfileIds = paste0(STUDY_ID, "_mrna")
  )

  df <- expr_raw[[paste0(STUDY_ID, "_mrna")]]
  if (is.null(df) || nrow(df) == 0) stop("No expression data returned from API")

  genes_returned <- unique(df$hugoGeneSymbol)
  genes_missing  <- setdiff(danaher_genes, genes_returned)
  .log <- c(.log, sprintf("  expression_immune_markers.csv: %d genes x %d patients",
                           length(genes_returned), length(unique(df$sampleId))))
  if (length(genes_missing) > 0) {
    .log <- c(.log, sprintf("  Missing genes: %s", paste(genes_missing, collapse = ", ")))
  }

  # Pivot to genes-as-rows, patients-as-cols format
  expr_wide <- df %>%
    select(hugoGeneSymbol, sampleId, value) %>%
    pivot_wider(names_from = sampleId, values_from = value,
                values_fn = ~ mean(.x, na.rm = TRUE))

  write.csv(expr_wide, file.path(proc_dir, "expression_immune_markers.csv"),
            row.names = FALSE)

  TRUE
}, error = function(e) {
  .log <<- c(.log, sprintf("  API failed: %s", conditionMessage(e)),
                           "  Falling back to existing processed data")
  FALSE
})

} # end if/else all_exist

# ---- Verify processed data exists -------------------------------------------
for (f in proc_files) {
  fp <- file.path(proc_dir, f)
  if (file.exists(fp)) {
    .log <- c(.log, sprintf("  ok %s (%s)", f, format(file.size(fp), big.mark = ",")))
  } else {
    .log <- c(.log, sprintf("  MISSING %s", f))
    warning(sprintf("MISSING %s — downstream scripts may fail.", f))
  }
}

# --- Summary ------------------------------------------------------------------
message(paste(.log, collapse = "\n"))
