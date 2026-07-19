# =============================================================================
# 05_milo_da.R
# Differential abundance testing with MiloR (CHR vs FEB, per timepoint)
#
# Neighbourhoods are built on a combined RNA + ADT PCA embedding (the same
# dimensions WNN was built with, via pc_cutoff) and on the day-subsetted WNN
# graph. Neighbourhood k / proportion auto-scale to the number of cells.
# The DA test is run at multiple FDR thresholds.
#
# Outputs: DA result CSVs, per-celltype summaries, UMAP DA plots and beeswarm
# plots (all as SVGs), one output folder per object x FDR threshold.
# =============================================================================

source("config.R")
source("utils.R")

install_bioc_if_missing(c("miloR", "SingleCellExperiment", "SummarizedExperiment",
                          "BiocParallel"))
install_cran_if_missing(c("svglite", "tidyr", "igraph", "ggbeeswarm", "ggplot2"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggbeeswarm)
  library(svglite)
  library(igraph)
  library(methods)
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


# ---- Auto-scale neighbourhood parameters to object size ----

auto_milo_params <- function(n_cells) {
  nhood_k <- dplyr::case_when(
    n_cells <  1000 ~ 15,
    n_cells <  5000 ~ 25,
    n_cells < 10000 ~ 35,
    n_cells < 30000 ~ 45,
    n_cells < 60000 ~ 55,
    TRUE            ~ 60
  )
  nhood_prop <- round(min(0.35, max(0.05, 300 / n_cells)), 3)
  list(nhood_k = as.integer(nhood_k), nhood_prop = nhood_prop)
}


# ---- Neighbourhood labelling helpers ----

label_nhoods_by_celltype <- function(milo, celltypes, min_prop = 0.5) {
  nh_mat <- miloR::nhoods(milo)
  k      <- ncol(nh_mat)
  if (k == 0) stop("No neighbourhoods found in milo object.")

  lab <- vapply(seq_len(k), function(j) {
    ii <- which(nh_mat[, j] > 0)
    ct <- celltypes[ii]
    ct <- ct[!is.na(ct) & nzchar(ct)]
    if (!length(ct)) return(NA_character_)
    tab  <- sort(table(ct), decreasing = TRUE)
    prop <- as.numeric(tab[1]) / sum(tab)
    if (!is.null(min_prop) && prop < min_prop) return("Mixed")
    names(tab)[1]
  }, character(1))

  prop_top <- vapply(seq_len(k), function(j) {
    ii <- which(nh_mat[, j] > 0)
    ct <- celltypes[ii]
    ct <- ct[!is.na(ct) & nzchar(ct)]
    if (!length(ct)) return(NA_real_)
    tab <- sort(table(ct), decreasing = TRUE)
    as.numeric(tab[1]) / sum(tab)
  }, numeric(1))

  data.frame(Nhood    = as.character(seq_len(k)),
             celltype = unname(lab),
             top_prop = unname(prop_top),
             stringsAsFactors = FALSE)
}

summarise_da_by_celltype <- function(da_df, nhood_ct, fdr_thresh = 0.05) {
  fdr_col <- intersect(c("SpatialFDR", "FDR", "adj.P.Val", "p_val_adj"),
                       colnames(da_df))[1]
  lfc_col <- intersect(c("logFC", "logFC.condition", "logFC_condition",
                         "logFC_conditionCHR", "logFC_conditionFEB"),
                       colnames(da_df))[1]
  if (is.na(fdr_col)) stop("No FDR column in DA results: ",   paste(colnames(da_df), collapse = ", "))
  if (is.na(lfc_col)) stop("No logFC column in DA results: ", paste(colnames(da_df), collapse = ", "))

  da_df$Nhood    <- as.character(da_df$Nhood %||% rownames(da_df))
  nhood_ct$Nhood <- as.character(nhood_ct$Nhood)

  out     <- dplyr::left_join(da_df, nhood_ct, by = "Nhood") %>%
    dplyr::mutate(sig = .data[[fdr_col]] < fdr_thresh)
  out_sig <- dplyr::filter(out, sig)

  sum_df <- if (!nrow(out_sig)) {
    data.frame(celltype = character(0), n_sig = integer(0),
               n_pos = integer(0), n_neg = integer(0),
               stringsAsFactors = FALSE)
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
                             condition_levels = c("FEB", "CHR"),
                             fdr_thresh       = 0.05,
                             nhood_min_prop   = 0.5,
                             svg_w = 7, svg_h = 6,
                             beeswarm_w = 7, beeswarm_h = 6) {
  outd <- file.path(PLOT_DIR, "milo", obj_name)
  dir.create(outd, showWarnings = FALSE, recursive = TRUE)
  message("\n  --- ", obj_name, " | day = ", DAY, " ---")

  required_cols <- c("day", "condition", "Sample", "celltype")
  missing_cols  <- setdiff(required_cols, colnames(obj@meta.data))
  if (length(missing_cols))
    stop("Missing metadata columns in ", obj_name, ": ",
         paste(missing_cols, collapse = ", "))

  # Auto-detect PCA dims using pc_cutoff (mirrors run_wnn in utils.R)
  rna_dims_use <- pc_cutoff(obj[["pca"]]@stdev)
  adt_dims_use <- pc_cutoff(obj[["adt.pca"]]@stdev)
  d_combined   <- rna_dims_use + adt_dims_use
  message("    RNA PCA dims: 1:", rna_dims_use,
          " | ADT PCA dims: 1:", adt_dims_use,
          " | combined: ", d_combined, "D")

  # Convert to SCE, store WNN UMAP
  DefaultAssay(obj) <- "RNA"
  sce      <- Seurat::as.SingleCellExperiment(obj, assay = "RNA")
  wnn_umap <- Seurat::Embeddings(obj, "wnn.umap")[colnames(sce), , drop = FALSE]
  SingleCellExperiment::reducedDim(sce, "UMAP") <- wnn_umap

  # Subset to day
  sce_day <- sce[, sce$day == DAY, drop = FALSE]
  if (ncol(sce_day) == 0) {
    warning(obj_name, ": 0 cells for day=", DAY, " — skipping.")
    return(NULL)
  }
  message("    Cells in day ", DAY, ": ", ncol(sce_day))

  # Auto-scale neighbourhood parameters to cell count
  auto       <- auto_milo_params(ncol(sce_day))
  nhood_k    <- auto$nhood_k
  nhood_prop <- auto$nhood_prop
  message("    Auto params — nhood_k: ", nhood_k, " | nhood_prop: ", nhood_prop)

  cells_day <- colnames(sce_day)

  # Build combined embedding (same dims WNN was built with)
  rna_embed      <- Seurat::Embeddings(obj, "pca")[cells_day, seq_len(rna_dims_use), drop = FALSE]
  adt_embed      <- Seurat::Embeddings(obj, "adt.pca")[cells_day, seq_len(adt_dims_use), drop = FALSE]
  combined_embed <- cbind(rna_embed, adt_embed)
  colnames(combined_embed) <- paste0("dim_", seq_len(d_combined))
  SingleCellExperiment::reducedDim(sce_day, "combined_pca") <- combined_embed

  # Fix condition factor
  cd           <- SummarizedExperiment::colData(sce_day)
  cd$condition <- factor(as.character(cd$condition), levels = condition_levels)
  SummarizedExperiment::colData(sce_day) <- cd

  # Subset WNN graph to day cells
  wnn_graph_name <- grep("wknn", names(obj@graphs), value = TRUE)[1]
  if (is.na(wnn_graph_name))
    stop("No wknn graph in ", obj_name,
         ". Available: ", paste(names(obj@graphs), collapse = ", "))
  message("    WNN graph: ", wnn_graph_name)

  wnn_sub <- obj@graphs[[wnn_graph_name]][cells_day, cells_day]
  wnn_ig  <- igraph::graph_from_adjacency_matrix(
    methods::as(wnn_sub, "dgCMatrix"),
    weighted = TRUE, mode = "undirected"
  )

  # Build Milo object
  milo               <- miloR::Milo(sce_day)
  miloR::graph(milo) <- wnn_ig

  milo <- miloR::makeNhoods(milo, prop = nhood_prop, k = nhood_k, d = d_combined,
                            refined = TRUE, reduced_dims = "combined_pca")
  message("    Neighbourhoods: ", ncol(miloR::nhoods(milo)))

  milo <- miloR::countCells(
    milo,
    meta.data = as.data.frame(
      SummarizedExperiment::colData(milo))[, c("Sample", "condition"), drop = FALSE],
    sample    = "Sample"
  )

  # Design matrix
  design_df <- as.data.frame(SummarizedExperiment::colData(milo)) %>%
    dplyr::distinct(Sample, condition) %>%
    tidyr::drop_na(Sample, condition)

  if (length(unique(design_df$condition)) < 2) {
    warning(obj_name, " day=", DAY, ": fewer than 2 conditions — skipping.")
    return(NULL)
  }

  rownames(design_df) <- design_df$Sample
  design_df$Sample    <- NULL
  design <- model.matrix(~ condition, data = design_df)

  # Spatial-FDR nhood distances + DA test
  milo <- miloR::calcNhoodDistance(milo, d = d_combined, reduced.dim = "combined_pca")
  da   <- miloR::testNhoods(milo, design = design, design.df = design_df,
                            reduced.dim   = "combined_pca",
                            fdr.weighting = "neighbour-distance")
  da$Nhood <- as.character(seq_len(nrow(da)))
  milo     <- miloR::buildNhoodGraph(milo)

  # Save DA CSV
  csv_path <- file.path(outd, paste0(obj_name, "_milo_DA_", DAY, "_CHR_vs_FEB.csv"))
  write.csv(da, csv_path, row.names = FALSE)
  message("    Saved: ", basename(csv_path))

  # Label neighbourhoods by celltype and summarise
  ct_vec        <- as.character(obj@meta.data[cells_day, "celltype"])
  names(ct_vec) <- cells_day

  nhood_ct <- label_nhoods_by_celltype(milo, ct_vec, min_prop = nhood_min_prop)
  da_sum   <- summarise_da_by_celltype(da, nhood_ct, fdr_thresh = fdr_thresh)

  sum_path <- file.path(outd, paste0(obj_name, "_milo_DA_", DAY, "_celltype_summary.csv"))
  write.csv(da_sum$summary, sum_path, row.names = FALSE)

  n_sig <- sum(da[[da_sum$fdr_col]] < fdr_thresh, na.rm = TRUE)
  message("    Saved: ", basename(sum_path), "  (FDR col: ", da_sum$fdr_col, ")")
  message("    Significant nhoods: ", n_sig, " / ", nrow(da))

  # UMAP DA plot (neighbourhood centroids on WNN UMAP)
  umap_coords <- SingleCellExperiment::reducedDim(milo, "UMAP")
  nh_mat      <- miloR::nhoods(milo)
  nhood_pos   <- t(vapply(seq_len(ncol(nh_mat)), function(j) {
    ii <- which(nh_mat[, j] > 0)
    colMeans(umap_coords[ii, , drop = FALSE])
  }, numeric(2)))

  plot_df <- data.frame(
    UMAP1      = nhood_pos[, 1],
    UMAP2      = nhood_pos[, 2],
    logFC      = da$logFC,
    SpatialFDR = da$SpatialFDR,
    sig        = da$SpatialFDR < fdr_thresh & !is.na(da$SpatialFDR)
  )

  p_umap <- ggplot2::ggplot(plot_df,
      ggplot2::aes(x = UMAP1, y = UMAP2, fill = logFC, size = sig, alpha = sig)) +
    ggplot2::geom_point(shape = 21, stroke = 0.2, color = "black") +
    ggplot2::scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                                  midpoint = 0, na.value = "grey80", name = "log FC") +
    ggplot2::scale_size_manual(values = c("TRUE" = 2.5, "FALSE" = 1), guide = "none") +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.3), guide = "none") +
    ggplot2::labs(title = paste0(obj_name, " — CHR vs FEB — day ", DAY)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5))

  umap_path <- file.path(outd, paste0(obj_name, "_milo_DA_umap_", DAY, ".svg"))
  save_svg_plot(p_umap, umap_path, width = svg_w, height = svg_h)

  # Beeswarm plot
  da_annotated <- da %>%
    dplyr::left_join(nhood_ct, by = "Nhood") %>%
    dplyr::mutate(sig = SpatialFDR < fdr_thresh & !is.na(SpatialFDR))

  da_nonsig <- dplyr::filter(da_annotated, !sig)
  da_sig    <- dplyr::filter(da_annotated,  sig)

  p_beeswarm <- ggplot2::ggplot(da_annotated,
      ggplot2::aes(x = logFC, y = celltype)) +
    ggbeeswarm::geom_quasirandom(data = da_nonsig, size = 0.4,
                                 color = "grey80", groupOnX = FALSE) +
    ggbeeswarm::geom_quasirandom(data = da_sig,
                                 ggplot2::aes(color = logFC),
                                 size = 0.8, groupOnX = FALSE) +
    ggplot2::scale_color_gradient2(low = "red", mid = "lightgrey", high = "blue",
                                   midpoint = 0, name = "log FC") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::labs(
      title    = paste0(obj_name, " — CHR vs FEB — day ", DAY),
      subtitle = paste0(n_sig, " / ", nrow(da), " nhoods significant"),
      x        = "log fold change",
      y        = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 10, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8),
      axis.text.y   = ggplot2::element_text(size = 8)
    )

  bee_path <- file.path(outd, paste0(obj_name, "_milo_DA_beeswarm_", DAY, ".svg"))
  save_svg_plot(p_beeswarm, bee_path, width = beeswarm_w, height = beeswarm_h)

  invisible(list(day                 = DAY,
                 milo                = milo,
                 da                  = da,
                 nhood_celltype      = nhood_ct,
                 da_celltype_summary = da_sum$summary,
                 plot_umap           = p_umap,
                 plot_beeswarm       = p_beeswarm))
}

run_milo_all_days <- function(obj, obj_name, days = NULL,
                              condition_levels = c("FEB", "CHR"),
                              fdr_thresh = 0.05, ...) {
  if (is.null(days)) days <- sort(unique(as.character(obj@meta.data$day)))
  message("\n========== ", obj_name, " | days: ",
          paste(days, collapse = ", "), " ==========")

  res <- lapply(days, function(D) {
    tryCatch(
      run_milo_one_day(obj, obj_name, D,
                       condition_levels = condition_levels,
                       fdr_thresh       = fdr_thresh, ...),
      error = function(e) {
        message("  ERROR in ", obj_name, " day=", D, ": ", conditionMessage(e))
        NULL
      }
    )
  })
  names(res) <- days
  invisible(res)
}


# ---- Run (all objects x all FDR thresholds) ----

DAYS             <- c("C1", "C9", "C14")
CONDITION_LEVELS <- c("FEB", "CHR")
FDR_THRESHOLDS   <- c(0.05, 0.1)

milo_results <- list()

for (fdr in FDR_THRESHOLDS) {
  fdr_tag <- paste0("fdr", gsub("\\.", "p", as.character(fdr)))
  message("\n########## FDR threshold: ", fdr, " (", fdr_tag, ") ##########")

  for (nm in names(objects_to_run)) {
    milo_results[[fdr_tag]][[nm]] <- run_milo_all_days(
      obj              = objects_to_run[[nm]],
      obj_name         = paste0(nm, "_", fdr_tag),
      days             = DAYS,
      condition_levels = CONDITION_LEVELS,
      fdr_thresh       = fdr
    )
  }
}

message("\nDONE. Milo completed for: ",
        paste(names(objects_to_run), collapse = ", "),
        " at FDR ", paste(FDR_THRESHOLDS, collapse = ", "))
