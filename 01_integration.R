# =============================================================================
# 01_integration.R
# Malaria CITE-seq integration pipeline
#
# Steps:
#   1. Load per-run RDS objects (post-HTO demultiplexing)
#   2. Merge time-points per patient (Day 1 / 9 / 14)
#   3. QC filtering (nFeature, % mitochondrial reads)
#   4. Doublet removal with DoubletFinder
#   5. RNA normalisation with SCTransform
#   6. ADT normalisation with DSB ModelNegativeADTnorm
#   7. RNA integration (RPCA anchor-based)
#   8. ADT integration (RPCA anchor-based)
#   9. WNN graph, UMAP, Louvain clustering
#  10. Cluster marker discovery and manual annotation
# =============================================================================

source("config.R")
source("utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(DoubletFinder)
  library(dsb)
  library(STACAS)
  library(dplyr)
})


# ---- 1. Load RDS files ----

rds_files <- list(
  RO1  = "FilterMatrix.RO-1_QC_HTOseperate.rds",
  RO2  = "FilterMatrix.RO-2_QC_HTOseperate.rds",
  RO3  = "FilterMatrix.RO-3_QC_HTOseperate.rds",
  RO4  = "FilterMatrix.RO-4_QC_HTOseperate.rds",
  RO5  = "FilterMatrix.RO-5_QC_HTOseperate.rds",
  RO6  = "FilterMatrix.RO-6_QC_HTOseperate.rds",
  RO7  = "RO7_QC_HTOseperate.rds",
  RO8  = "RO8_QC_HTOseperate.rds",
  RO9  = "RO9_QC_HTOseperate.rds",
  RO10 = "RO10_QC_HTOseperate.rds"
)

runs <- lapply(rds_files, function(f) readRDS(file.path(RDS_DIR, f)))
list2env(runs, envir = .GlobalEnv)


# ---- 2. Build per-patient objects (merge Day 1 / 9 / 14) ----

merge_patient <- function(d1, d9, d14, pid) {
  d1$day  <- "C1"
  d9$day  <- "C9"
  d14$day <- "C14"
  merged  <- merge(d1, y = list(d9, d14),
                   add.cell.ids = c("C1", "C9", "C14"),
                   merge.data   = FALSE)
  merged$Patient <- pid
  merged$patient <- paste(pid, merged$day, sep = "_")
  merged
}

multis <- list(
  "P01" = merge_patient(subset(RO5,  patient == "P01"),
                        subset(RO4,  patient == "P01"),
                        subset(RO6,  patient == "P01"), "P01"),
  "P02" = merge_patient(subset(RO7,  patient == "P02"),
                        subset(RO9,  patient == "P02"),
                        subset(RO10, patient == "P02"), "P02"),
  "P03" = merge_patient(subset(RO10, patient == "P03"),
                        subset(RO8,  patient == "P03"),
                        subset(RO7,  patient == "P03"), "P03"),
  "P04" = merge_patient(subset(RO6,  patient == "P04"),
                        subset(RO5,  patient == "P04"),
                        subset(RO4,  patient == "P04"), "P04"),
  "P05" = merge_patient(subset(RO4,  patient == "P05"),
                        subset(RO6,  patient == "P05"),
                        subset(RO5,  patient == "P05"), "P05"),
  "P06" = merge_patient(subset(RO8,  patient == "P06"),
                        subset(RO10, patient == "P06"),
                        subset(RO9,  patient == "P06"), "P06"),
  "P07" = merge_patient(subset(RO1,  patient == "P07"),
                        subset(RO3,  patient == "P07"),
                        subset(RO2,  patient == "P07"), "P07"),
  "P08" = merge_patient(subset(RO3,  patient == "P08"),
                        subset(RO2,  patient == "P08"),
                        subset(RO1,  patient == "P08"), "P08"),
  "P09" = merge_patient(subset(RO2,  patient == "P09"),
                        subset(RO1,  patient == "P09"),
                        subset(RO3,  patient == "P09"), "P09"),
  "P10" = merge_patient(subset(RO9,  patient == "P10"),
                        subset(RO7,  patient == "P10"),
                        subset(RO8,  patient == "P10"), "P10")
)

# P06 is the chronic-tolerant (CHR) patient; all others are febrile (FEB)
multis[["P06"]]$condition <- "CHR"
for (nm in setdiff(names(multis), "P06")) multis[[nm]]$condition <- "FEB"


# ---- 3. QC filtering ----

qc_filter <- function(s) {
  s <- JoinLayers(s)
  s[["percent.mt"]] <- PercentageFeatureSet(s, "^MT-")
  subset(s,
         nFeature_RNA > 200  &
           nFeature_RNA < 6000 &
           percent.mt   < 20)
}

multis <- lapply(multis, qc_filter)


# ---- 4. Doublet removal (DoubletFinder) ----

mark_doublets <- function(s) {
  DefaultAssay(s) <- "RNA"
  s <- NormalizeData(s, verbose = FALSE) |>
    FindVariableFeatures(verbose = FALSE) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(npcs = N_PCS_RNA, verbose = FALSE)

  sweep.list  <- paramSweep(s, PCs = seq_len(N_PCS_RNA))
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn       <- find.pK(sweep.stats)
  optimal_pk  <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  nExp      <- round(DOUBLET_RATE * ncol(s))
  before_md <- colnames(s@meta.data)

  s <- doubletFinder(s, PCs = seq_len(N_PCS_RNA), pN = 0.25, pK = optimal_pk,
                     nExp = nExp, reuse.pANN = FALSE, sct = FALSE)

  new_md  <- setdiff(colnames(s@meta.data), before_md)
  df_cols <- grep("^DF\\.classifications", new_md, value = TRUE)
  if (!length(df_cols))
    stop("DoubletFinder did not produce a DF.classifications column.")

  singlets <- rownames(s@meta.data)[s@meta.data[[df_cols[1]]] == "Singlet"]
  subset(s, cells = singlets)
}

multis <- lapply(multis, mark_doublets)

saveRDS(multis, file.path(OUTPUT_DIR, "after_doublet_removal.rds"))


# ---- 5. RNA normalisation (SCTransform) ----

multis <- lapply(multis, function(s) {
  SCTransform(s, assay = "RNA", new.assay.name = "SCT", verbose = TRUE)
})


# ---- 6. ADT normalisation (DSB ModelNegativeADTnorm) ----

# Remove HTO rows from the ADT assay before normalisation
hto_feats <- grep("^HTO", rownames(multis[[1]][["ADT"]]), value = TRUE)

multis <- lapply(multis, function(s) {
  if (!"ADT" %in% names(s@assays)) return(s)
  counts     <- GetAssayData(s, assay = "ADT", layer = "counts")
  keep_feats <- setdiff(rownames(counts), hto_feats)
  s[["ADT"]] <- CreateAssayObject(counts = counts[keep_feats, , drop = FALSE])
  s
})

raw_counts_list <- lapply(multis, function(s) {
  if ("ADT" %in% names(s@assays)) GetAssayData(s, assay = "ADT", layer = "counts")
  else NULL
})

dsb_list <- lapply(raw_counts_list, function(mat) {
  if (is.null(mat)) return(NULL)
  dsb::ModelNegativeADTnorm(
    cell_protein_matrix      = as.matrix(mat),
    denoise.counts           = TRUE,
    use.isotype.control      = TRUE,
    isotype.control.name.vec = ISOTYPES,
    define.pseudocount       = FALSE,
    return.stats             = FALSE
  )
})
names(dsb_list) <- names(raw_counts_list)

multis <- mapply(function(s, raw, dsb_out) {
  if (is.null(dsb_out)) return(s)
  s[["ADT"]] <- CreateAssayObject(counts = raw)
  s <- SetAssayData(s, assay = "ADT", slot = "data", new.data = dsb_out)
  s
}, multis, raw_counts_list, dsb_list, SIMPLIFY = FALSE)


# ---- 7. RNA integration (RPCA) ----

rna_features <- SelectIntegrationFeatures(multis, nfeatures = N_HVG_RNA)
multis       <- PrepSCTIntegration(multis, anchor.features = rna_features)

multis <- lapply(multis, function(s) {
  RunPCA(s, assay = "SCT", features = rna_features,
         reduction.name = "pca", npcs = N_PCS_RNA, verbose = FALSE)
})

rna_anchors <- FindIntegrationAnchors(
  object.list          = multis,
  normalization.method = "SCT",
  anchor.features      = rna_features,
  reduction            = "rpca",
  dims                 = seq_len(N_PCS_RNA),
  k.filter             = 200
)

integ_rna <- IntegrateData(rna_anchors, normalization.method = "SCT")
saveRDS(integ_rna, file.path(OUTPUT_DIR, "integ_rna.rds"))


# ---- 8. ADT integration (RPCA) ----

adt_list <- lapply(multis, function(s) {
  DefaultAssay(s) <- "ADT"
  VariableFeatures(s) <- head(rownames(s[["ADT"]]), N_HVG_ADT)
  s
})

adt_features <- setdiff(
  SelectIntegrationFeatures(adt_list, nfeatures = 200),
  ISOTYPES
)

adt_list <- lapply(adt_list, function(s) {
  DefaultAssay(s) <- "ADT"
  s <- ScaleData(s, assay = "ADT", features = adt_features, do.scale = FALSE)
  s <- RunPCA(s, assay = "ADT", features = adt_features,
              reduction.name = "pca", npcs = N_PCS_ADT, verbose = FALSE)
  s
})

adt_anchors <- FindIntegrationAnchors(
  adt_list,
  assay                = rep("ADT", length(adt_list)),
  anchor.features      = adt_features,
  normalization.method = "LogNormalize",
  reduction            = "rpca",
  dims                 = 1:30,
  scale                = FALSE
)

integ_adt <- IntegrateData(adt_anchors, normalization.method = "LogNormalize",
                           new.assay.name = "integratedADT")
DefaultAssay(integ_adt) <- "integratedADT"
integ_adt <- ScaleData(integ_adt, do.scale = FALSE) |>
  RunPCA(npcs = N_PCS_ADT, reduction.name = "adt.pca")

saveRDS(integ_adt, file.path(OUTPUT_DIR, "integ_adt.rds"))

# Transfer ADT assays into the RNA-integrated object
integ_rna[["ADT"]]           <- integ_adt[["ADT"]]
integ_rna[["integratedADT"]] <- integ_adt[["integratedADT"]]
integ_rna[["adt.pca"]]       <- integ_adt[["adt.pca"]]


# ---- 9. WNN graph, UMAP and clustering ----

DefaultAssay(integ_rna) <- "integrated"
integ_rna <- RunPCA(integ_rna, verbose = FALSE)

rna_pcs <- pc_cutoff(integ_rna[["pca"]]@stdev)
adt_pcs <- pc_cutoff(integ_rna[["adt.pca"]]@stdev)
message("RNA PCs used: ", rna_pcs, "  |  ADT PCs used: ", adt_pcs)

integ_rna <- FindMultiModalNeighbors(
  integ_rna,
  reduction.list = list("pca", "adt.pca"),
  dims.list      = list(seq_len(rna_pcs), seq_len(adt_pcs))
)
integ_rna <- FindClusters(integ_rna, graph.name = "wsnn", resolution = WNN_RES_GLOBAL)
integ_rna <- RunUMAP(integ_rna, nn.name = "weighted.nn",
                     reduction.name = "wnn.umap", reduction.key = "WNN_",
                     return.model = TRUE)

DimPlot(integ_rna, reduction = "wnn.umap", label = TRUE)


# ---- 10. Cluster markers and manual annotation ----

# 10A. RNA markers from SCT residuals
Idents(integ_rna)       <- "seurat_clusters"
DefaultAssay(integ_rna) <- "RNA"
integ_rna <- PrepSCTFindMarkers(integ_rna, assay = "RNA")

rna_markers <- FindAllMarkers(
  integ_rna,
  assay           = "SCT",
  only.pos        = TRUE,
  test.use        = "wilcox",
  min.pct         = 0.50,
  logfc.threshold = 0.50
)

# 10B. ADT markers
DefaultAssay(integ_rna) <- "ADT"
adt_markers <- FindAllMarkers(
  integ_rna,
  assay           = "ADT",
  only.pos        = TRUE,
  test.use        = "wilcox",
  min.pct         = 0.40,
  logfc.threshold = 0.40
)

# Top 20 RNA and top 10 ADT markers per cluster
top_rna <- rna_markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(avg_log2FC, n = 20) %>%
  arrange(cluster, desc(avg_log2FC))

top_adt <- adt_markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(avg_log2FC, n = 10) %>%
  arrange(cluster, desc(avg_log2FC))

print(top_rna, n = Inf)
print(top_adt, n = Inf)

# 10C. Manual cell-type annotation
cluster_ids <- c(
  `0`  = "CD14+ Monocytes",
  `1`  = "CD4 T EM",
  `2`  = "KLRC2+ NK",
  `3`  = "Naïve CD4 T",
  `4`  = "Memory CD8 T",
  `5`  = "Naïve CD8 T",
  `6`  = "Memory B",
  `7`  = "CD16+ Monocytes",
  `8`  = "KLRC1+ NK",
  `9`  = "Naïve B",
  `10` = "Vδ2 T",
  `11` = "Naïve CD4 T",
  `12` = "MAIT/γδ T",
  `13` = "Regulatory T",
  `14` = "Erythroid",
  `15` = "Cycling T",
  `16` = "GZMB+ CD8 T",
  `17` = "Stressed cells",
  `18` = "Plasma cells",
  `19` = "Transitional B",
  `20` = "Platelets",
  `21` = "Progenitors",
  `22` = "pDC"
)

integ_rna <- RenameIdents(integ_rna, cluster_ids)
integ_rna$celltype <- Idents(integ_rna)

saveRDS(integ_rna, file.path(OUTPUT_DIR, "integ_rna_annotated.rds"))
message("01_integration.R complete.")
