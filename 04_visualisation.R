# =============================================================================
# 04_visualisation.R
# Produces publication-quality SVG figures:
#   A. UMAP + standalone legend per lineage object
#   B. Ordered RNA + ADT dotplots
#   C. Relative cell-type abundance (line plots per patient group x timepoint)
#   D. ComplexHeatmap (gene expression z-score, cluster x day x condition)
# =============================================================================

source("config.R")
source("utils.R")

install_cran_if_missing(c("Seurat", "dplyr", "ggplot2", "svglite",
                           "gtable", "Polychrome", "patchwork",
                           "RColorBrewer", "openxlsx"))
install_bioc_if_missing(c("ComplexHeatmap", "circlize"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(svglite)
  library(patchwork)
  library(gtable)
  library(Polychrome)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

# ---- Load objects ----

objects_to_run <- list(
  gd_object    = readRDS(file.path(OUTPUT_DIR, "gd_object.rds")),
  b_cells      = readRDS(file.path(OUTPUT_DIR, "b_cells.rds")),
  cd8_cells    = readRDS(file.path(OUTPUT_DIR, "cd8_cells.rds")),
  nk_focus     = readRDS(file.path(OUTPUT_DIR, "nk_focus.rds")),
  CD4_and_treg = readRDS(file.path(OUTPUT_DIR, "CD4_and_treg.rds")),
  mono_obj     = readRDS(file.path(OUTPUT_DIR, "mono_obj.rds"))
)

deg_results_all <- readRDS(file.path(OUTPUT_DIR, "deg_results_mast.rds"))


# ========== A. UMAP SVGs ==========

make_umap_svg <- function(obj, obj_name,
                           reduction       = "wnn.umap",
                           pt_size         = 1,
                           umap_w_cm       = 8.76,
                           umap_h_cm       = 8.83,
                           legend_nrow     = 3,
                           legend_key_in   = 0.22,
                           legend_text_size  = 10,
                           legend_title_size = 10,
                           legend_w_cm     = 8.76,
                           legend_h_cm     = 3.45,
                           celltype_col    = "celltype",
                           lum_max         = 0.90) {
  outd <- file.path(PLOT_DIR, "umap", obj_name)
  dir.create(outd, showWarnings = FALSE, recursive = TRUE)

  obj <- hard_assign_celltype_palette(obj, celltype_col = celltype_col,
                                       force = TRUE, lum_max = lum_max)
  pal <- get_celltype_palette(obj, celltype_col = celltype_col, lum_max = lum_max)

  p_umap <- Seurat::DimPlot(obj, reduction = reduction, group.by = celltype_col,
                             label = FALSE, cols = pal,
                             repel = TRUE, raster = FALSE, pt.size = pt_size) +
    Seurat::NoLegend() +
    theme(plot.margin = margin(0, 0, 0, 0))

  save_svg_plot(p_umap,
                file.path(outd, paste0(obj_name, "_UMAP.svg")),
                width  = umap_w_cm / 2.54,
                height = umap_h_cm / 2.54)

  p_leg <- Seurat::DimPlot(obj, reduction = reduction, group.by = celltype_col,
                            label = FALSE, cols = pal,
                            repel = TRUE, raster = FALSE, pt.size = pt_size) +
    guides(colour = guide_legend(nrow = legend_nrow, byrow = TRUE,
                                  override.aes = list(size = 3))) +
    theme_void() +
    theme(
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.title      = element_text(size = legend_title_size),
      legend.text       = element_text(size = legend_text_size),
      legend.key.width  = unit(legend_key_in, "in"),
      legend.key.height = unit(legend_key_in, "in"),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin       = margin(0, 0, 0, 0)
    )

  legend_grob <- gtable::gtable_filter(ggplotGrob(p_leg), "guide-box", trim = TRUE)
  save_svg_grob(legend_grob,
                file.path(outd, paste0(obj_name, "_UMAP_LEGEND.svg")),
                width  = legend_w_cm / 2.54,
                height = legend_h_cm / 2.54)

  assign(obj_name, obj, envir = .GlobalEnv)
  invisible(obj)
}

lapply(names(objects_to_run),
       function(nm) make_umap_svg(objects_to_run[[nm]], nm, pt_size = 0.3))


# ========== B. RNA + ADT dotplots ==========

# Per-object manually curated significant genes to add to dotplots
sig_map <- list(
  gd_object = list(RNA = c("STAT1", "RNF213", "PIM1", "KLRC1")),
  nk_focus  = list(RNA = c("DUSP1", "ZFP36", "FKBP5", "STAT1", "GBP5",
                             "PARP9", "IRF1", "PDIA3", "NFKBIA", "FTH1")),
  cd8_cells = list(RNA = c("KLRC2", "LILRB1", "TYROBP", "TNFRSF9",
                             "IRF4", "PDCD1", "HAVCR2", "STAT1")),
  b_cells   = list(RNA = c("HLA-DRB5", "PRDM1", "XBP1", "MZB1",
                             "IL4R", "BLNK", "ZFP36", "FKBP5"))
)

run_deg_and_dotplots_svg <- function(seu, obj_name,
                                      out_dir  = file.path(PLOT_DIR, "dotplots"),
                                      sig_map  = list(),
                                      top_n    = 10) {
  stopifnot("celltype" %in% colnames(seu@meta.data))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  Idents(seu) <- "celltype"
  ct_order    <- names(sort(table(Idents(seu)), decreasing = TRUE))
  seu$celltype <- factor(seu$celltype, levels = ct_order)
  Idents(seu) <- "celltype"

  # RNA DEGs
  DefaultAssay(seu) <- "RNA"
  rna_deg <- FindAllMarkers(seu, assay = "RNA", only.pos = TRUE,
                             logfc.threshold = 0.25, min.pct = 0.20,
                             return.thresh = 0.05, test.use = "wilcox") %>%
    group_by(cluster) %>%
    distinct(gene, .keep_all = TRUE) %>%
    ungroup()

  write.csv(rna_deg,
            file.path(out_dir, paste0("DEG_RNA__", obj_name, ".csv")),
            row.names = FALSE)

  rna_top <- rna_deg %>%
    group_by(cluster) %>%
    arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()

  # ADT DEGs
  DefaultAssay(seu) <- "ADT"
  adt_deg <- FindAllMarkers(seu, assay = "ADT", only.pos = TRUE,
                             logfc.threshold = 0.25, min.pct = 0.20,
                             return.thresh = 0.05, test.use = "wilcox") %>%
    group_by(cluster) %>%
    distinct(gene, .keep_all = TRUE) %>%
    ungroup()

  write.csv(adt_deg,
            file.path(out_dir, paste0("DEG_ADT__", obj_name, ".csv")),
            row.names = FALSE)

  adt_top <- adt_deg %>%
    group_by(cluster) %>%
    arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()

  # Build ordered feature lists
  grp_levels <- levels(seu$celltype)
  sm         <- sig_map[[obj_name]] %||% list()

  rna_feats <- character(0)
  for (g in grp_levels) {
    top_g <- rna_top %>% filter(cluster == g) %>% pull(gene)
    sig_g <- intersect(sm$RNA %||% character(0), rownames(seu[["RNA"]]))
    rna_feats <- c(rna_feats, top_g, setdiff(sig_g, top_g))
  }

  adt_feats <- character(0)
  for (g in grp_levels) {
    top_g <- adt_top %>% filter(cluster == g) %>% pull(gene)
    sig_g <- intersect(sm$ADT %||% character(0), rownames(seu[["ADT"]]))
    adt_feats <- c(adt_feats, top_g, setdiff(sig_g, top_g))
  }

  seu          <- fix_assay_duplicate_rownames(seu, "RNA")
  seu          <- fix_assay_duplicate_rownames(seu, "ADT")
  rna_markers  <- uniq_keep_order(present_in(uniq_keep_order(rna_feats), rownames(seu[["RNA"]])))
  adt_markers  <- uniq_keep_order(present_in(uniq_keep_order(adt_feats), rownames(seu[["ADT"]])))

  dot_theme <- theme_minimal() +
    theme(
      axis.text.x      = element_text(angle = 90, hjust = 1, size = 9),
      axis.text.y      = element_text(size = 12),
      axis.title       = element_blank(),
      plot.title       = element_text(hjust = 0.5, size = 13),
      panel.grid.major = element_line(color = "darkgray", linewidth = 0.5),
      panel.grid.minor = element_line(color = "darkgray", linewidth = 0.25)
    )

  DefaultAssay(seu) <- "RNA"
  p_rna <- DotPlot(seu, features = rna_markers, group.by = "celltype", assay = "RNA") +
    scale_color_gradientn(colors = c("blue", "lightgrey", "red")) +
    labs(title = paste0("RNA Expression — ", obj_name)) +
    dot_theme

  DefaultAssay(seu) <- "ADT"
  p_adt <- DotPlot(seu, features = adt_markers, group.by = "celltype",
                   assay = "ADT", scale = TRUE, cluster.idents = FALSE) +
    scale_color_gradientn(colors = c("blue", "lightgrey", "red")) +
    scale_x_discrete(labels = function(x) gsub("^Hu(?:MsRt|Ms)?\\.?", "", x, perl = TRUE)) +
    labs(title = paste0("ADT Expression — ", obj_name)) +
    dot_theme

  rna_svg <- file.path(out_dir, paste0("DOT_RNA__", obj_name, ".svg"))
  adt_svg <- file.path(out_dir, paste0("DOT_ADT__", obj_name, ".svg"))
  ggsave(rna_svg, p_rna, width = 10.91, height = 6.06, units = "in")
  ggsave(adt_svg, p_adt, width = 10.91, height = 6.06, units = "in")
  message("Saved SVGs for: ", obj_name)

  list(object = seu, rna_dotplot = p_rna, adt_dotplot = p_adt)
}

dotplot_results <- list()
for (nm in names(objects_to_run)) {
  res               <- run_deg_and_dotplots_svg(objects_to_run[[nm]], nm,
                                                 sig_map = sig_map)
  dotplot_results[[nm]] <- res
  objects_to_run[[nm]]  <- res$object
}


# ========== C. Relative abundance plots ==========

make_freq_svg <- function(obj_name, freq_df,
                           out_dir        = file.path(PLOT_DIR, "freq"),
                           w_cm           = 10.29,
                           h_cm           = 10.16,
                           day_levels     = c("C1", "C9", "C14"),
                           max_cols       = 4,
                           celltypes_keep = NULL,
                           celltype_levels = NULL) {
  stopifnot(all(c("day", "pct_total", "condition", "celltype", "n") %in% colnames(freq_df)))

  df <- freq_df %>%
    mutate(day      = factor(as.character(day), levels = day_levels),
           celltype = as.character(celltype)) %>%
    filter(!is.na(day))

  if (!is.null(celltypes_keep)) df <- df %>% filter(celltype %in% celltypes_keep)

  ct_order <- celltype_levels %||%
    (df %>%
       group_by(celltype) %>%
       summarise(total = sum(n, na.rm = TRUE), .groups = "drop") %>%
       arrange(desc(total)) %>%
       pull(celltype))

  df$celltype <- factor(df$celltype, levels = intersect(ct_order, unique(df$celltype)))

  plots <- split(df, df$celltype) |>
    lapply(function(dd) {
      ggplot(dd, aes(day, pct_total, colour = condition, group = condition)) +
        stat_summary(fun = mean, geom = "line",     linewidth = 1) +
        stat_summary(fun = mean, geom = "point",    size = 2) +
        stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25, linewidth = 0.6) +
        scale_colour_manual(values = COND_COLS) +
        scale_x_discrete(drop = FALSE) +
        labs(y = "% of all cells", x = NULL, title = unique(dd$celltype)) +
        theme_minimal(base_size = 9) +
        theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
              legend.position = "none",
              plot.margin     = margin(1, 1, 1, 1))
    })

  k    <- length(plots)
  if (k == 0) stop(obj_name, ": 0 celltypes to plot.")
  grid <- choose_patchwork_grid(k, max_cols = max_cols)

  p <- patchwork::wrap_plots(plots, ncol = grid$ncol) +
    patchwork::plot_annotation(
      title = obj_name,
      theme = theme(plot.title = element_text(face = "bold", size = 10, hjust = 0))
    )

  outd     <- file.path(out_dir, obj_name)
  dir.create(outd, showWarnings = FALSE, recursive = TRUE)
  svg_path <- file.path(outd, paste0(obj_name, "_freq.svg"))

  svglite::svglite(svg_path, width = w_cm / 2.54, height = h_cm / 2.54)
  print(p)
  grDevices::dev.off()

  message("Panels=", k, " | ncol=", grid$ncol, " | nrow=", grid$nrow)
  message("Saved: ", svg_path)
  invisible(TRUE)
}

# freq_df_list should be a named list of data frames with columns:
# day, pct_total, condition, celltype, n
# Build it from objects_to_run, e.g.:
# freq_df_list <- lapply(objects_to_run, compute_freq_df)  # user-defined


# ========== D. ComplexHeatmap ==========

# Gene selection helpers

select_top_deg_in_order <- function(deg_list, top_per_group,
                                     adj_p_cutoff, log2fc_cutoff,
                                     unique_genes = TRUE) {
  genes <- character(0)
  seen  <- new.env(parent = emptyenv())
  for (nm in names(deg_list)) {
    df <- deg_list[[nm]]
    if (is.null(df) || !nrow(df)) next
    if (!all(c("p_val_adj", "avg_log2FC") %in% colnames(df))) next
    df <- df[df$p_val_adj < adj_p_cutoff & abs(df$avg_log2FC) > log2fc_cutoff, , drop = FALSE]
    if (!nrow(df)) next
    df        <- df[order(df$p_val_adj), , drop = FALSE]
    top_genes <- rownames(head(df, top_per_group))
    top_genes <- top_genes[!is.na(top_genes) & nzchar(top_genes)]
    if (unique_genes) {
      top_genes <- top_genes[!vapply(top_genes, exists, logical(1),
                                      envir = seen, inherits = FALSE)]
      for (g in top_genes) assign(g, TRUE, envir = seen)
    }
    genes <- c(genes, top_genes)
  }
  unique(genes)
}

get_sig_genes_from_deg <- function(deg_list, adj_p_cutoff, log2fc_cutoff) {
  sig <- character(0)
  for (nm in names(deg_list)) {
    df <- deg_list[[nm]]
    if (is.null(df) || !nrow(df)) next
    if (!all(c("p_val_adj", "avg_log2FC") %in% colnames(df))) next
    df  <- df[df$p_val_adj < adj_p_cutoff & abs(df$avg_log2FC) > log2fc_cutoff, , drop = FALSE]
    sig <- c(sig, rownames(df))
  }
  unique(sig[!is.na(sig) & nzchar(sig)])
}

# Manual gene lists for text annotation on heatmaps
manual_de_genes_map <- list(
  mono_obj = c(
    "TLR7","TLR8","ALPK1","JAK2","STAT1","STAT2","IRF1","IRF2","IRF7",
    "IFI44L","IFIT1","IFIT2","IFIT3","IFITM3","ISG15","OAS1","OAS2","OAS3",
    "GBP1","GBP2","GBP4","GBP5","LY6E","EPSTI1","NLRC5","TAP1","PSME2",
    "HLA-DQA1","HLA-DRB5","LAP3","FCGR1A","FCGR1B","VAMP5","PLAAT4","PLAC8",
    "CALHM6","PSTPIP2","DYNLL1","MYOF","SERPING1","APOL3","APOL6","MT2A",
    "WARS","GIMAP4","C5AR1","FCAR","PHACTR1","DNMBP","ARL4C","ADGRE2",
    "CXCL8","CSF3R","ECE1","TREM1","LY86","TNFRSF21","TOB1","CEBPD","NR4A1",
    "MAFB","ZFHX3","ZNF331","PPARG","G0S2","CHKA","PDE4D","RNF152","VNN1",
    "B3GNT5","XYLT1","EEPD1","MBNL2","ZFAS1","LDLRAD3","TMCC3","CDA","THBS1"
  ),
  nk_focus = c(
    "DUSP1","ZFP36","FKBP5","STAT1","GBP5","PARP9","IRF1","PDIA3","NFKBIA",
    "FTH1","RAB8B","ITGA4","RPL35A","RAPGEF6","FYB1","TRIO","ITGAX","CDK13",
    "CRIP1","PSME1","CCDC91","UBL5","RPL35","DHRS7","S100A4"
  ),
  gd_object = c("STAT1","RNF213","PIM1","KLRC1"),
  cd8_cells = c(
    "KLRC2","LILRB1","TYROBP","TNFRSF9","IRF4","PDCD1","HAVCR2","XCL1",
    "XCR1","CCL4L2","PROK2","IL7","IKZF2","HHEX","ZBTB16","STAT1","IFI44L",
    "IFI27","GBP1","GBP4","EPSTI1","ETV7","ALPK1","SOCS2","SOCS3","EGLN3",
    "MTHFD1L","ABCG2","PBX1","RAI2","SAMD4A","ZNF667","ZNF793","ZNF175",
    "ZNF135","SEMA3A","PLXNA4","SDK1","PCDH9","PTK7","DAB1","SERPINE2",
    "MXRA8","HAPLN3"
  ),
  b_cells = c("HLA-DRB5","HLA-DQA2","PRDM1","XBP1","MZB1","IL4R","BLNK",
              "ZFP36","FKBP5","TXNIP","SIMC1")
)

plot_heatmap <- function(seurat_obj, deg_results, obj_name,
                          cluster_col      = "celltype",
                          condition_col    = "condition",
                          day_col          = "day",
                          day_levels       = c("C1", "C9", "C14"),
                          clusters_to_keep = NULL,
                          adj_p_cutoff     = HEAT_ADJ_P,
                          log2fc_cutoff    = HEAT_LOGFC,
                          top_per_group    = HEAT_TOP_N,
                          top_n_genes      = HEAT_MAX_GENES,
                          unique_genes     = TRUE,
                          assay            = "RNA",
                          show_gene_names  = TRUE,
                          gene_fontsize    = 9,
                          anno_bar_mm      = 1.4,
                          anno_gap_mm      = 0.6,
                          palette_full     = NULL,
                          id_levels_full   = NULL) {
  message("[", obj_name, "] Building heatmap")

  stopifnot(all(c(day_col, cluster_col, condition_col) %in% colnames(seurat_obj@meta.data)))
  if (is.null(palette_full))   stop("palette_full required.")
  if (is.null(id_levels_full)) stop("id_levels_full required.")
  names(palette_full) <- canon_ct(names(palette_full))
  id_levels_full      <- canon_ct(id_levels_full)
  if (!setequal(id_levels_full, names(palette_full)))
    stop("id_levels_full and names(palette_full) differ.")

  seurat_obj@meta.data[[cluster_col]] <- canon_ct(seurat_obj@meta.data[[cluster_col]])
  seurat_obj <- subset_by_clusters_cells(seurat_obj, cluster_col, clusters_to_keep)
  seurat_obj@meta.data[[condition_col]] <- toupper(canon_ct(seurat_obj@meta.data[[condition_col]]))
  seurat_obj@meta.data[[day_col]] <- factor(as.character(seurat_obj@meta.data[[day_col]]),
                                             levels = day_levels)

  present_genes <- rownames(seurat_obj[[assay]])

  # Combine top DEG genes with manually curated genes
  deg_genes   <- intersect(
    select_top_deg_in_order(deg_results, top_per_group, adj_p_cutoff,
                             log2fc_cutoff, unique_genes),
    present_genes
  )
  sig_universe <- get_sig_genes_from_deg(deg_results, adj_p_cutoff, log2fc_cutoff)
  text_genes   <- intersect(
    intersect(manual_de_genes_map[[obj_name]] %||% character(0), present_genes),
    sig_universe
  )

  all_genes <- if (is.null(top_n_genes)) {
    unique(c(text_genes, deg_genes))
  } else {
    remaining  <- max(0L, top_n_genes - length(text_genes))
    unique(c(text_genes, head(setdiff(deg_genes, text_genes), remaining)))
  }

  if (length(all_genes) < 2) stop("[", obj_name, "] Too few genes after filtering.")
  message("[", obj_name, "] ", length(all_genes), " genes selected.")

  Idents(seurat_obj) <- factor(seurat_obj@meta.data[[cluster_col]], levels = id_levels_full)
  seurat_obj$group_label <- paste(
    canon_ct(as.character(Idents(seurat_obj))),
    seurat_obj@meta.data[[condition_col]],
    as.character(seurat_obj@meta.data[[day_col]]),
    sep = "|"
  )
  keep_cells  <- !is.na(seurat_obj$group_label) & nzchar(seurat_obj$group_label)
  seurat_obj  <- subset(seurat_obj, cells = colnames(seurat_obj)[keep_cells])

  has_layer <- "layer" %in% names(formals(Seurat::AverageExpression))
  avg_args  <- list(object   = seurat_obj, features = all_genes,
                    assays   = assay,      group.by = "group_label")
  if (has_layer) avg_args$layer <- "data" else avg_args$slot <- "data"

  expr_mat <- do.call(Seurat::AverageExpression, avg_args)[[assay]][all_genes, , drop = FALSE]
  expr_z   <- t(scale(t(expr_mat)))
  expr_z   <- expr_z[rowSums(!is.na(expr_z)) > 0, , drop = FALSE]

  parts    <- strsplit(colnames(expr_z), "\\|")
  annot_df <- data.frame(
    cluster   = canon_ct(sapply(parts, `[`, 1)),
    condition = canon_ct(sapply(parts, `[`, 2)),
    day       = sapply(parts, `[`, 3),
    row.names = colnames(expr_z),
    stringsAsFactors = FALSE
  )
  annot_df$cluster   <- factor(annot_df$cluster,            levels = id_levels_full)
  annot_df$condition <- factor(toupper(annot_df$condition), levels = c("CHR", "FEB"))
  annot_df$day       <- factor(annot_df$day,                levels = day_levels)

  ord      <- with(annot_df, order(cluster, day, condition))
  expr_z   <- expr_z[, ord, drop = FALSE]
  annot_df <- annot_df[ord, , drop = FALSE]

  clusters_present <- levels(droplevels(annot_df$cluster))
  cluster_colors   <- palette_full[clusters_present]
  if (any(is.na(cluster_colors)))
    stop("[", obj_name, "] palette missing: ", paste(clusters_present[is.na(cluster_colors)], collapse = ", "))
  names(cluster_colors) <- clusters_present

  day_colors       <- setNames(brewer.pal(3, "Set2")[1:3], day_levels)
  condition_colors <- c(CHR = "steelblue", FEB = "firebrick")

  ha <- HeatmapAnnotation(
    Cluster   = droplevels(annot_df$cluster),
    Day       = annot_df$day,
    Condition = annot_df$condition,
    col       = list(Cluster   = cluster_colors,
                     Day       = day_colors,
                     Condition = condition_colors),
    annotation_height = unit(rep(anno_bar_mm, 3), "mm"),
    gap               = unit(anno_gap_mm, "mm"),
    show_legend       = FALSE
  )

  lab_gp         <- gpar(fontsize = gene_fontsize, lineheight = 1.0)
  row_name_width <- max_text_width(rownames(expr_z), gp = lab_gp) + unit(0.5, "mm")
  col_fun        <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

  ht <- Heatmap(
    expr_z,
    name                = "Z-score",
    col                 = col_fun,
    top_annotation      = ha,
    show_row_names      = show_gene_names,
    row_names_gp        = lab_gp,
    row_names_max_width = row_name_width,
    show_column_names   = FALSE,
    cluster_rows        = TRUE,
    cluster_columns     = FALSE,
    show_heatmap_legend = FALSE,
    rect_gp             = gpar(col = "grey90", lwd = 0.3)
  )

  lgd_z    <- Legend(title = "Z-score", col_fun = col_fun,
                     at = c(-4, 0, 4), direction = "horizontal")
  lgd_day  <- Legend(title = "Day",  labels = names(day_colors),
                     legend_gp = gpar(fill = day_colors),  direction = "horizontal")
  lgd_cond <- Legend(title = "Condition", labels = names(condition_colors),
                     legend_gp = gpar(fill = condition_colors), direction = "horizontal")
  legends  <- packLegend(lgd_z, lgd_day, lgd_cond,
                          direction = "horizontal", gap = unit(3, "mm"))

  invisible(list(heatmap = ht, legends = legends))
}

# Run heatmaps for all objects
heatmap_objs <- list()
exclude_map  <- list(
  nk_focus = "CD3+ CD8+ T"   # contaminant cluster — exclude from heatmap
)

for (nm in names(objects_to_run)) {
  if (is.null(deg_results_all[[nm]]) || length(deg_results_all[[nm]]) == 0) {
    message("Skipping ", nm, ": no DEG results.")
    next
  }

  seu          <- objects_to_run[[nm]]
  pal_full     <- get_stored_palette_strict(seu)
  id_levels    <- names(pal_full)
  cl_keep      <- setdiff(id_levels, exclude_map[[nm]] %||% character(0))

  objs <- plot_heatmap(
    seurat_obj       = seu,
    deg_results      = deg_results_all[[nm]],
    obj_name         = nm,
    clusters_to_keep = cl_keep,
    palette_full     = pal_full,
    id_levels_full   = id_levels
  )

  obj_out <- file.path(PLOT_DIR, "heatmaps", nm)
  dir.create(obj_out, showWarnings = FALSE, recursive = TRUE)
  save_heatmap_svg(objs$heatmap,
                   file.path(obj_out, paste0(nm, "_heatmap.svg")))
  save_legends_svg(objs$legends,
                   file.path(obj_out, paste0(nm, "_heatmap_legends.svg")))
  heatmap_objs[[nm]] <- objs
}

message("\nDONE. Heatmaps built for: ", paste(names(heatmap_objs), collapse = ", "))
