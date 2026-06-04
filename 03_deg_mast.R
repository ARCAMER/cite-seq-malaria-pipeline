# =============================================================================
# 03_deg_mast.R
# Differential expression (CHR vs FEB) per celltype x timepoint
# using MAST with sex as a latent variable
# =============================================================================

source("config.R")
source("utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

install_bioc_if_missing("MAST")
suppressPackageStartupMessages(library(MAST))

# ---- Load objects ----

objects_to_run <- list(
  gd_object    = readRDS(file.path(OUTPUT_DIR, "gd_object.rds")),
  b_cells      = readRDS(file.path(OUTPUT_DIR, "b_cells.rds")),
  cd8_cells    = readRDS(file.path(OUTPUT_DIR, "cd8_cells.rds")),
  nk_focus     = readRDS(file.path(OUTPUT_DIR, "nk_focus.rds")),
  CD4_and_treg = readRDS(file.path(OUTPUT_DIR, "CD4_and_treg.rds")),
  mono_obj     = readRDS(file.path(OUTPUT_DIR, "mono_obj.rds"))
)

# ---- Helper: add sex annotation from patient ID ----

add_sex_from_patient <- function(seu, sex_map, patient_col = "patient") {
  pat       <- as.character(seu@meta.data[[patient_col]])
  seu$sex   <- unname(sex_map[pat])
  seu$sex   <- factor(seu$sex, levels = c("male", "female"))
  seu
}

# ---- Helper: basic metadata validation ----

assert_vars_present <- function(seu, vars) {
  miss <- setdiff(vars, colnames(seu@meta.data))
  if (length(miss)) stop("Missing metadata columns: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

# ---- Run MAST per object x celltype x timepoint ----

col_celltype  <- "celltype"
col_day       <- "day"
col_condition <- "condition"
col_patient   <- "patient"

deg_results_all <- list()

for (nm in names(objects_to_run)) {
  seu <- objects_to_run[[nm]]
  seu <- add_sex_from_patient(seu, SEX_MAP, patient_col = col_patient)

  message("\n====================")
  message("OBJECT: ", nm)
  print(table(seu$sex, useNA = "ifany"))

  assay_used <- if ("RNA" %in% names(seu@assays)) "RNA" else DefaultAssay(seu)
  subtypes   <- sort(unique(as.character(seu@meta.data[[col_celltype]])))
  days       <- sort(unique(as.character(seu@meta.data[[col_day]])))
  fm_formals <- names(formals(Seurat::FindMarkers))

  obj_res <- list()

  for (subtype in subtypes) {
    for (d in days) {
      sub_obj <- subset_cells_by_celltype_day(seu, subtype, d,
                                               col_celltype = col_celltype,
                                               col_day      = col_day)
      if (is.null(sub_obj) || ncol(sub_obj) == 0) next

      assert_vars_present(sub_obj, c(col_condition, "sex"))

      cond_vec   <- as.character(sub_obj@meta.data[[col_condition]])
      cond_count <- table(cond_vec)

      if (!(all(c("CHR", "FEB") %in% names(cond_count)) &&
            all(cond_count[c("CHR", "FEB")] >= MAST_MIN_CELLS))) next

      Idents(sub_obj)       <- col_condition
      DefaultAssay(sub_obj) <- assay_used
      if (is.factor(sub_obj$sex)) sub_obj$sex <- droplevels(sub_obj$sex)

      sex_levels <- unique(na.omit(as.character(sub_obj$sex)))

      args <- list(
        object          = sub_obj,
        ident.1         = "CHR",
        ident.2         = "FEB",
        assay           = assay_used,
        test.use        = "MAST",
        min.pct         = MAST_MIN_PCT,
        logfc.threshold = MAST_LOGFC_THRESH
      )
      if (length(sex_levels) > 1) args$latent.vars <- "sex"
      if ("layer" %in% fm_formals) args$layer <- "data"

      message(nm, " | ", subtype, " | ", d,
              " | CHR=", cond_count["CHR"], " FEB=", cond_count["FEB"],
              if (length(sex_levels) > 1) " | +sex" else " | no-sex")

      markers <- tryCatch(
        do.call(Seurat::FindMarkers, args),
        error = function(e) {
          message("FAIL ", nm, " | ", subtype, " | ", d, " : ", conditionMessage(e))
          NULL
        }
      )

      if (!is.null(markers) && nrow(markers) > 0) {
        key          <- paste(subtype, d, sep = "__")
        obj_res[[key]] <- markers
      }
    }
  }

  deg_results_all[[nm]] <- obj_res
  message("Stored DEG tables for ", nm, ": ", length(obj_res))
}

message("\nDONE. Objects processed: ", paste(names(deg_results_all), collapse = ", "))

saveRDS(deg_results_all, file.path(OUTPUT_DIR, "deg_results_mast.rds"))
