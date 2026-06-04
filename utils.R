# =============================================================================
# utils.R — Shared helper functions used across pipeline scripts
# Source this after config.R: source("utils.R")
# =============================================================================

# ---- Null-coalescing operator ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- String helpers ----
canon_ct      <- function(x)    trimws(as.character(x))
uniq_keep_order <- function(x)  x[!duplicated(x)]
present_in    <- function(f, v) f[match(f, v, nomatch = 0L) > 0L]

# ---- PC cutoff heuristic ----
# Returns the elbow point in a scree plot: whichever comes first —
# (a) cumulative variance > cumu_thr AND per-PC variance < pct_thr, or
# (b) the last PC where the variance drop exceeds drop_thr.
pc_cutoff <- function(stdev, cumu_thr = 90, pct_thr = 5, drop_thr = 0.1) {
  pct  <- stdev / sum(stdev) * 100
  cumu <- cumsum(pct)
  co1  <- which(cumu > cumu_thr & pct < pct_thr)[1]
  drops <- pct[-length(pct)] - pct[-1]
  co2   <- if (any(drops > drop_thr)) tail(which(drops > drop_thr), 1) + 1L else length(stdev)
  min(co1, co2, na.rm = TRUE)
}

# ---- Preprocessing helpers ----

#' Standard RNA PCA (scales on 'integrated' assay)
prep_rna_pca <- function(obj, hvf_layer = "data", hvgs = NULL,
                         scale_vars = NULL, n_pcs = 50, pca_name = "pca") {
  if (is.null(hvgs)) {
    DefaultAssay(obj) <- "RNA"
    obj  <- FindVariableFeatures(obj, layer = hvf_layer)
    hvgs <- VariableFeatures(obj)
  }
  DefaultAssay(obj) <- "integrated"
  obj <- ScaleData(obj, features = hvgs, vars.to.regress = scale_vars)
  obj <- RunPCA(obj, reduction.name = pca_name, features = hvgs, npcs = n_pcs)
  obj
}

#' ADT PCA (scales on 'integratedADT' assay, drops isotype controls)
prep_adt_pca <- function(obj, nfeatures, scale_vars = NULL,
                         n_pcs = 50, pca_name = "adt.pca") {
  DefaultAssay(obj) <- "ADT"
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  DefaultAssay(obj) <- "integratedADT"
  adt_hvfs <- setdiff(VariableFeatures(obj), ISOTYPES)
  obj <- ScaleData(obj, features = adt_hvfs, do.scale = FALSE, do.center = TRUE,
                   vars.to.regress = scale_vars)
  obj <- RunPCA(obj, features = adt_hvfs, reduction.name = pca_name, npcs = n_pcs)
  obj
}

#' WNN graph, clustering, and UMAP in one call
run_wnn <- function(obj, rna_red = "pca", adt_red = "adt.pca",
                    resolution, algo = 3, umap_name = "wnn.umap") {
  rna_pcs <- pc_cutoff(obj[[rna_red]]@stdev)
  adt_pcs <- pc_cutoff(obj[[adt_red]]@stdev)
  obj <- FindMultiModalNeighbors(
    obj,
    reduction.list = list(rna_red, adt_red),
    dims.list      = list(seq_len(rna_pcs), seq_len(adt_pcs))
  )
  obj <- FindClusters(obj, graph.name = "wsnn", resolution = resolution, algorithm = algo)
  obj <- RunUMAP(obj, nn.name = "weighted.nn",
                 reduction.name = umap_name, reduction.key = "WNN_", return.model = TRUE)
  obj
}

# ---- Marker discovery ----

#' Run FindAllMarkers for RNA and ADT, returning top hits per cluster
find_markers_and_top <- function(obj,
                                 rna_only_pos = TRUE, adt_only_pos = TRUE,
                                 rna_n = 20, adt_n = 10,
                                 rna_abs = FALSE, adt_abs = FALSE) {
  Idents(obj) <- "seurat_clusters"

  DefaultAssay(obj) <- "RNA"
  rna_markers <- FindAllMarkers(obj, only.pos = rna_only_pos)

  DefaultAssay(obj) <- "ADT"
  adt_markers <- FindAllMarkers(obj, only.pos = adt_only_pos)

  top_fun <- function(df, n, use_abs) {
    df %>%
      group_by(cluster) %>%
      filter(p_val_adj < 0.05) %>%
      { if (use_abs) slice_max(., abs(avg_log2FC), n = n)
        else         slice_max(., avg_log2FC,       n = n) }
  }

  list(
    rna     = rna_markers,
    adt     = adt_markers,
    top_rna = top_fun(rna_markers, rna_n, rna_abs),
    top_adt = top_fun(adt_markers, adt_n, adt_abs)
  )
}

# ---- Subset helpers ----

#' Subset by celltype x day without using subset= expressions (avoids scoping bugs)
subset_cells_by_celltype_day <- function(seu, subtype, day,
                                         col_celltype = "celltype",
                                         col_day = "day") {
  md   <- seu@meta.data
  keep <- (as.character(md[[col_celltype]]) == subtype) &
    (as.character(md[[col_day]]) == day)
  cells_use <- rownames(md)[which(keep)]
  if (!length(cells_use)) return(NULL)
  subset(seu, cells = cells_use)
}

#' Subset by cluster identity using cell barcodes (no subset= expression)
subset_by_clusters_cells <- function(seu, cluster_col, clusters_to_keep) {
  if (is.null(clusters_to_keep)) return(seu)
  md   <- seu@meta.data
  stopifnot(cluster_col %in% colnames(md))
  keep <- canon_ct(md[[cluster_col]]) %in% canon_ct(clusters_to_keep)
  cells_use <- rownames(md)[which(keep)]
  if (!length(cells_use)) stop("0 cells left after filtering clusters_to_keep.")
  subset(seu, cells = cells_use)
}

# ---- Celltype factor ordering ----

#' Regex-based lineage + differentiation-stage ordering for celltypes
first_match_rank <- function(x, regex_vec) {
  x2 <- tolower(x)
  r  <- rep(length(regex_vec) + 1L, length(x2))
  for (i in seq_along(regex_vec)) {
    hit      <- grepl(regex_vec[[i]], x2, perl = TRUE)
    r[hit & r == (length(regex_vec) + 1L)] <- i
  }
  r
}

reorder_celltype_by_name <- function(obj) {
  labs <- sort(unique(as.character(obj$celltype)))

  lineage_names <- c("cd4","treg","cd8","nk","gd","mait","mono","dc","b","plasma","other")
  lineage_patterns <- c(
    cd4    = "\\bcd4\\b|helper",
    treg   = "treg|foxp3",
    cd8    = "\\bcd8\\b|ctl|cytotoxic",
    nk     = "\\bnk\\b|nkg7|gnly",
    gd     = "gamma.?delta|\\bgd\\b|trdc|trgc",
    mait   = "mait",
    mono   = "mono|cd14|fcgr3a|lyz|macroph",
    dc     = "\\bdc\\b|dendritic|fcerg1a|cst3",
    b      = "\\bb\\b|cd79a|ms4a1|igh",
    plasma = "plasma|plasmablast|jchain|mzb1|xbp1",
    other  = ".*"
  )

  stage_names <- c("naive","cm","memory","effector","temra","exhausted","cycling","ifn","stress","other")
  stage_patterns <- c(
    naive     = "naive|\\btn\\b",
    cm        = "tcm|central[ _-]?memory",
    memory    = "memory|\\btem\\b|effector[ _-]?memory|trm|gzm(k)?\\+",
    effector  = "effector|teff|cytotoxic|ctl|gzm|prf1|fgfbp2",
    temra     = "temra",
    exhausted = "exhaust|\\btex\\b|pdcd1|lag3|tigit|havcr2|ctla4|entpd1",
    cycling   = "cycle|cycling|mki67|top2a|stmn1",
    ifn       = "\\bifn\\b|interferon|isg",
    stress    = "stress|heat.?shock|\\bhsp\\b|^mt-|mitochond",
    other     = ".*"
  )

  lin_rank   <- first_match_rank(labs, lineage_patterns[lineage_names])
  stage_rank <- first_match_rank(labs, stage_patterns[stage_names])
  levs <- labs[order(lin_rank, stage_rank, labs)]

  obj$celltype <- factor(canon_ct(as.character(obj$celltype)), levels = levs, ordered = TRUE)
  Idents(obj)  <- obj$celltype
  obj
}

# ---- Palette helpers ----

#' Relative luminance of a hex colour (sRGB → linear → Y)
hex_luminance <- function(hex) {
  hex <- gsub("^#", "", hex)
  rgb <- grDevices::col2rgb(paste0("#", hex)) / 255
  lin <- ifelse(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  as.numeric(0.2126 * lin[1, ] + 0.7152 * lin[2, ] + 0.0722 * lin[3, ])
}

#' k Glasbey colours with near-white filtered out
distinct_palette_no_white <- function(levs, lum_max = 0.90) {
  levs       <- as.character(levs)
  k          <- length(levs)
  extra      <- max(50L, ceiling(k * 2.0))
  candidates <- Polychrome::glasbey.colors(k + extra)
  L          <- hex_luminance(candidates)
  keep       <- candidates[L <= lum_max]
  if (length(keep) < k) keep <- candidates[L <= min(0.95, lum_max + 0.05)]
  keep <- keep[!toupper(keep) %in% c("#FFFFFF", "#FFFFFE", "#FFFFFD")]
  if (length(keep) < k) stop("Not enough non-white colours. Try lum_max = 0.95 or reduce k.")
  cols       <- keep[seq_len(k)]
  names(cols) <- levs
  cols
}

get_desired_levels <- function(obj, celltype_col = "celltype") {
  stopifnot(celltype_col %in% colnames(obj@meta.data))
  ct_raw <- obj@meta.data[[celltype_col]]
  ct_chr <- canon_ct(as.character(ct_raw))
  ct_chr <- ct_chr[nzchar(ct_chr)]
  levs <- if (is.factor(ct_raw)) {
    canon_ct(levels(ct_raw))
  } else if (is.factor(Seurat::Idents(obj))) {
    canon_ct(levels(Seurat::Idents(obj)))
  } else {
    unique(ct_chr)
  }
  levs  <- levs[nzchar(levs)]
  miss  <- setdiff(unique(ct_chr), levs)
  if (length(miss)) levs <- c(levs, miss)
  levs
}

hard_assign_celltype_palette <- function(obj, celltype_col = "celltype",
                                         force = TRUE, lum_max = 0.90) {
  stopifnot(celltype_col %in% colnames(obj@meta.data))
  levs     <- get_desired_levels(obj, celltype_col)
  ct_chr   <- canon_ct(as.character(obj@meta.data[[celltype_col]]))
  obj@meta.data[[celltype_col]] <- factor(ct_chr, levels = levs)
  Seurat::Idents(obj) <- obj@meta.data[[celltype_col]]

  if (!force) {
    pal_old <- obj@misc$celltype_palette %||% NULL
    if (!is.null(pal_old) &&
        setequal(canon_ct(names(pal_old)), levs) &&
        length(pal_old) == length(levs)) {
      names(pal_old) <- canon_ct(names(pal_old))
      obj@misc$celltype_palette <- pal_old[levs]
      return(obj)
    }
  }

  pal_new <- distinct_palette_no_white(levs, lum_max = lum_max)
  obj@misc$celltype_palette <- pal_new
  obj
}

get_celltype_palette <- function(obj, celltype_col = "celltype", lum_max = 0.90) {
  obj  <- hard_assign_celltype_palette(obj, celltype_col, force = TRUE, lum_max = lum_max)
  levs <- levels(Seurat::Idents(obj))
  obj@misc$celltype_palette[levs]
}

get_stored_palette_strict <- function(obj, celltype_col = "celltype") {
  levs <- get_desired_levels(obj, celltype_col)
  pal  <- obj@misc$celltype_palette %||% NULL
  if (is.null(pal) || !length(pal))
    stop("obj@misc$celltype_palette is missing. Run hard_assign_celltype_palette() first.")
  names(pal) <- canon_ct(names(pal))
  missing    <- setdiff(levs, names(pal))
  if (length(missing))
    stop("Palette missing levels: ", paste(missing, collapse = ", "))
  pal[levs]
}

# ---- Save helpers ----

save_svg_plot <- function(p, path, width = 8, height = 6) {
  svglite::svglite(path, width = width, height = height)
  print(p)
  grDevices::dev.off()
  message("Saved: ", path)
}

save_svg_grob <- function(g, path, width = 8, height = 5) {
  svglite::svglite(path, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(g)
  grDevices::dev.off()
  message("Saved: ", path)
}

save_heatmap_svg <- function(ht, path, outer_pad_mm = 1.5, w_cm = 19.05, h_cm = 15.23) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  svglite::svglite(path, width = w_cm / 2.54, height = h_cm / 2.54)
  ComplexHeatmap::draw(ht, padding = grid::unit(rep(outer_pad_mm, 4), "mm"))
  grDevices::dev.off()
  message("Saved: ", path)
}

save_legends_svg <- function(legend_grob, path, w_cm = 10.29, h_cm = 2.12) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  svglite::svglite(path, width = w_cm / 2.54, height = h_cm / 2.54)
  grid::grid.newpage()
  ComplexHeatmap::draw(legend_grob)
  grDevices::dev.off()
  message("Saved: ", path)
}

# ---- Package installation helpers ----

install_cran_if_missing <- function(pkgs, repos = "https://cloud.r-project.org") {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) install.packages(missing_pkgs, repos = repos)
  invisible(TRUE)
}

install_bioc_if_missing <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing_pkgs)) return(invisible(TRUE))
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE)
  invisible(TRUE)
}

# ---- Assay QC ----

#' Suffix duplicate rownames in an assay to prevent downstream failures
fix_assay_duplicate_rownames <- function(obj, assay_name) {
  assay <- obj[[assay_name]]
  rn    <- rownames(assay)
  if (!any(duplicated(rn))) return(obj)
  tab      <- table(rn)
  dups     <- names(tab[tab > 1])
  new_names <- rn
  for (sym in dups) {
    idx <- which(rn == sym)
    for (k in seq(2, length(idx))) new_names[idx[k]] <- paste0(sym, "__dup", k)
  }
  if (length(assay@counts))   rownames(assay@counts)     <- new_names
  if (length(assay@data))     rownames(assay@data)       <- new_names
  if (nrow(assay@scale.data)) rownames(assay@scale.data) <- new_names
  obj[[assay_name]] <- assay
  obj
}

# ---- Grid layout helper ----

#' Choose a balanced ncol x nrow grid for k panels
choose_patchwork_grid <- function(k, max_cols = 4, min_cols = 1) {
  k <- as.integer(k)
  if (k <= 2) return(list(ncol = k, nrow = 1))
  if (k <= 4) return(list(ncol = 2, nrow = ceiling(k / 2)))
  ncol <- max(min_cols, min(max_cols, ceiling(sqrt(k))))
  nrow <- ceiling(k / ncol)
  last_row <- k - (nrow - 1L) * ncol
  if (last_row == 1L && ncol > 2L) {
    ncol2    <- ncol - 1L
    nrow2    <- ceiling(k / ncol2)
    last_row2 <- k - (nrow2 - 1L) * ncol2
    if (last_row2 > last_row) { ncol <- ncol2; nrow <- nrow2 }
  }
  list(ncol = ncol, nrow = nrow)
}
