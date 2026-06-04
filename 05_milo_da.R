# =============================================================================
# 05_milo_da.R
# Differential abundance testing with MiloR (CHR vs FEB, per timepoint)
# Outputs: DA result CSVs, UMAP DA plots, beeswarm plots (all as SVGs)
# =============================================================================

source("config.R")
source("utils.R")

install_bioc_if_missing(c("miloR", "SingleCellExperiment", "SummarizedExperiment"))
install_cran_if_missing(c("svglite", "tidyr"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(svglite)
  library(miloR)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

objects_to_run <- list(
  gd_object    = readRDS(file.path(OUTPUT_DIR, "gd_object.rds")),
  b_cells      = readRDS(file.path(OUTPUT_DIR, "b_cells.rds")),
  cd8_cells    = readRDS(file.path(OUTPUT_DIR, "cd8_cells.rds")),
  nk_focus     = readRDS(file.path(OUTPUT_DIR, "nk_focus.rds")),
  CD4_and_treg = readRDS(file.path(OUTPUT_DIR, "CD4_and_treg.rds")),
  mono_obj     = readRDS(file.path(OUTPUT_DIR, "mono_obj.rds")),
  integ_rna    = readRDS(file.path(OUTPUT_DIR, "integ_rna_annotated.rds"))
)


# ---- Neighbourhood labelling helpers ----

label_nhoods_by_celltype <- function(milo, celltypes, min_prop = 0.5) {
  stopifnot(length(celltypes) == ncol(milo))
  nh_idx   <- miloR::nhoodIndex(milo)
  nhood_ids <- if (!is.null(names(nh_idx)) && any(nzchar(names(nh_idx))))
    as.character(names(nh_idx))
  else
    as.character(seq_along(nh_idx))

  lab <- vapply(nh_idx, function(ii) {
    ct <- celltypes[ii]
    ct <- ct[!is.na(ct) & nzchar(ct)]
    if (!length(ct)) return(NA_character_)
    tab  <- sort(table(ct), decreasing = TRUE)
    prop <- as.numeric(tab[1]) / sum(tab)
    if (!is.null(min_prop) && prop < min_prop) return("Mixed")
    names(tab)[1]
  }, character(1))

  prop_top <- vapply(nh_idx, function(ii) {
    ct <- celltypes[ii]
    ct <- ct[!is.na(ct) & nzchar(ct)]
    if (!length(ct)) return(NA_real_)
    tab <- sort(table(ct), decreasing = TRUE)
    as.numeric(tab[1]) / sum(tab)
  }, numeric(1))

  data.frame(Nhood    = nhood_ids,
             celltype = unname(lab),
             top_prop = unname(prop_top),
             stringsAsFactors = FALSE)
}

summarise_da_by_celltype <- function(da_df, nhood_ct,
                                      fdr_thresh = 0.05) {
  fdr_col <- intersect(c("SpatialFDR","FDR","adj.P.Val","p_val_adj"),
                        colnames(da_df))[1]
  lfc_col <- intersect(c("logFC","logFC.condition","logFC_condition",
                          "logFC_conditionCHR","logFC_conditionFEB"),
                        colnames(da_df))[1]
  if (is.na(fdr_col)) stop("No FDR column in DA results: ", paste(colnames(da_df), collapse=", "))
  if (is.na(lfc_col)) stop("No logFC column in DA results: ", paste(colnames(da_df), collapse=", "))

  da_df$Nhood <- as.character(da_df$Nhood %||% rownames(da_df))
  nhood_ct$Nhood <- as.character(nhood_ct$Nhood)

  out     <- dplyr::left_join(da_df, nhood_ct, by = "Nhood") %>%
    dplyr::mutate(sig = .data[[fdr_col]] < fdr_thresh)
  out_sig <- out %>% dplyr::filter(sig)

  sum_df <- if (!nrow(out_sig)) {
    data.frame(celltype = character(0), n_sig = integer(0),
               n_pos = integer(0), n_neg = integer(0))
  } else {
    out_sig %>%
      dplyr::group_by(celltype) %>%
      dplyr::summarise(
        n_sig = dplyr::n(),
        n_pos = sum(.data[[lfc_col]] > 0, na.rm = TRUE),
        n_neg = sum(.data[[lfc_col]] < 0, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(n_sig))
  }
  list(per_nhood = out, summary = sum_df, fdr_col = fdr_col, lfc_col = lfc_col)
}


# ---- Per-day Milo run ----

run_milo_one_day <- function(obj, obj_name, DAY,
                              svg_w = 6.5, svg_h = 5.5,
                              beeswarm_w = 6.5, beeswarm_h = 5.5,
                              condition_levels   = c("FEB", "CHR"),
                              fdr_thresh         = 0.05,
                              nhood_min_prop     = 0.5) {
  outd <- file.path(PLOT_DIR, "milo", obj_name)
  dir.create(outd, showWarnings = FALSE, recursive = TRUE)

  DefaultAssay(obj) <- "RNA"
  sce      <- Seurat::as.SingleCellExperiment(obj, assay = "RNA")
  wnn_umap <- Seurat::Embeddings(obj, "wnn.umap")[colnames(sce), , drop = FALSE]
  SingleCellExperiment::reducedDim(sce, "UMAP") <- wnn_umap

  sce_day <- sce[, sce$day == DAY, drop = FALSE]
  if (ncol(sce_day) == 0) {
    warning(obj_name, ": 0 cells for day=", DAY, " — skipping.")
    return(NULL)
  }

  cd <- SummarizedExperiment::colData(sce_day)
  cd$condition <- factor(cd$condition, levels = condition_levels)
  SummarizedExperiment::colData(sce_day) <- cd

  milo <- miloR::Milo(sce_day)
  milo <- miloR::buildGraph(milo, k = 30, d = 30, reduced.dim = "PCA")
  milo <- miloR::makeNhoods(milo, prop = 0.1, k = 30, d = 30,
                              refined = TRUE, reduced_dims = "PCA")
  milo <- miloR::countCells(
    milo,
    meta.data = as.data.frame(SummarizedExperiment::colData(milo)[, c("Sample","condition")]),
    sample    = "Sample"
  )

  design_df <- SummarizedExperiment::colData(milo) %>%
    as.data.frame() %>%
    dplyr::distinct(Sample, condition) %>%
    tidyr::drop_na(Sample, condition)

  if (length(unique(design_df$condition)) < 2) {
    warning(obj_name, " day=", DAY, ": only one condition — skipping Milo.")
    return(NULL)
  }

  rownames(design_df) <- design_df$Sample
  design_df$Sample    <- NULL
  design <- model.matrix(~ condition, data = design_df)

  da   <- miloR::testNhoods(milo, design = design, design.df = design_df,
                              reduced.dim = "UMAP")
  da$Nhood <- as.character(da$Nhood)
  milo <- miloR::buildNhoodGraph(milo)

  csv_path <- file.path(outd, paste0(obj_name, "_milo_DA_", DAY, "_CHR_vs_FEB.csv"))
  write.csv(da, csv_path, row.names = FALSE)

  ct_vec        <- setNames(as.character(obj@meta.data[colnames(obj), "celltype"]),
                             colnames(obj))
  ct_vec        <- ct_vec[colnames(sce_day)]
  nhood_ct      <- label_nhoods_by_celltype(milo, ct_vec, min_prop = nhood_min_prop)
  da_sum        <- summarise_da_by_celltype(da, nhood_ct, fdr_thresh = fdr_thresh)

  sum_path <- file.path(outd, paste0(obj_name, "_milo_DA_", DAY, "_celltype_summary.csv"))
  write.csv(da_sum$summary, sum_path, row.names = FALSE)

  # UMAP DA plot
  p_milo <- miloR::plotNhoodGraphDA(
    milo, milo_res = da, alpha = fdr_thresh,
    title = paste0(obj_name, " — CHR vs FEB — ", DAY, " (WNN UMAP)")
  )
  svg_path <- file.path(outd, paste0(obj_name, "_milo_DA_umap_", DAY, ".svg"))
  svglite::svglite(svg_path, width = svg_w, height = svg_h)
  print(p_milo)
  grDevices::dev.off()
  message("Saved: ", svg_path)

  # Beeswarm plot
  da_annotated <- da %>% dplyr::left_join(nhood_ct, by = "Nhood")
  p_beeswarm   <- miloR::plotDAbeeswarm(da_annotated, group.by = "celltype",
                                         alpha = fdr_thresh) +
    ggplot2::labs(title = paste0(obj_name, " — CHR vs FEB — ", DAY)) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(size = 10, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8)
    )
  bee_path <- file.path(outd, paste0(obj_name, "_milo_beeswarm_", DAY, ".svg"))
  svglite::svglite(bee_path, width = beeswarm_w, height = beeswarm_h)
  print(p_beeswarm)
  grDevices::dev.off()
  message("Saved: ", bee_path)

  invisible(list(day          = DAY,
                 milo         = milo,
                 da           = da,
                 plot_umap    = p_milo,
                 plot_bees    = p_beeswarm,
                 nhood_ct     = nhood_ct,
                 da_summary   = da_sum$summary))
}

run_milo_all_days <- function(obj, obj_name, days = NULL,
                               condition_levels = c("CHR", "FEB"), ...) {
  if (is.null(days)) days <- sort(unique(as.character(obj@meta.data$day)))
  res        <- lapply(days, function(D) run_milo_one_day(obj, obj_name, D,
                                                           condition_levels = condition_levels, ...))
  names(res) <- days
  res
}

# ---- Run ----

milo_results <- lapply(names(objects_to_run), function(nm) {
  run_milo_all_days(objects_to_run[[nm]], nm, condition_levels = c("CHR", "FEB"))
})
names(milo_results) <- names(objects_to_run)

message("\nDONE. Milo completed for: ", paste(names(milo_results), collapse = ", "))
