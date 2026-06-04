# =============================================================================
# 02_subset_analysis.R
# Iterative WNN sub-clustering for each major immune lineage
#
# Lineages processed:
#   A. T/NK (TNK_ALL) — all T and NK clusters from the global object
#   B. CD8 / NK / MAIT  (cd8nk_subset)
#   C. NK focus          (nk_focus)
#   D. γδ T / MAIT       (gd_object)
#   E. B cells           (b_cells)
#   F. CD8 T cells       (cd8_cells)
#   G. CD4 T + Tregs     (CD4_and_treg)
# =============================================================================

source("config.R")
source("utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scGate)
})

integ_rna <- readRDS(file.path(OUTPUT_DIR, "integ_rna_annotated.rds"))


# ---- A. T / NK subset ----

tnk_all <- subset(integ_rna,
                  subset = seurat_clusters %in% c(1, 2, 3, 4, 5, 8, 9, 11, 12, 14, 18))

DefaultAssay(tnk_all) <- "RNA"
tnk_all <- prep_rna_pca(tnk_all)
tnk_all <- prep_adt_pca(tnk_all, nfeatures = 140)
tnk_all <- run_wnn(tnk_all, resolution = WNN_RES_TNK)

DimPlot(tnk_all, reduction = "wnn.umap", label = TRUE) + NoLegend()

# Cluster annotation
tnk_ids <- c(
  "0"  = "CD4_Tmem/Th",
  "1"  = "CD4_Tnaive/TCM",
  "2"  = "NK_CD56dim_cytotoxic",
  "3"  = "CD8_GZMK_TEM",
  "4"  = "CD8_Tnaive",
  "5"  = "NK_CD56bright_cytokine_resp",
  "6"  = "NK_KIRpos_CD57pos",
  "7"  = "CD4_IL7R_Tmem",
  "8"  = "CD8_CX3CR1_TEMRA",
  "9"  = "MAIT",
  "10" = "gammaDelta_T_Vd2",
  "11" = "Treg",
  "12" = "CD4_CCR6_Th17like",
  "13" = "InnateLike_PLZF_T",
  "14" = "CD8_Tnaive/TCM_2",
  "15" = "NK_CD56bright_activated",
  "16" = "CD4_Activated_PD1hi",
  "17" = "CD8_Tnaive/TCM_3"
)
tnk_all <- RenameIdents(tnk_all, tnk_ids)

lineage_map_tnk <- c(
  "CD4_Tmem/Th"           = "CD4",
  "CD4_Tnaive/TCM"        = "CD4",
  "CD4_IL7R_Tmem"         = "CD4",
  "CD4_CCR6_Th17like"     = "CD4",
  "CD4_Activated_PD1hi"   = "CD4",
  "Treg"                  = "Treg",
  "CD8_GZMK_TEM"          = "CD8",
  "CD8_Tnaive"            = "CD8",
  "CD8_CX3CR1_TEMRA"      = "CD8",
  "CD8_Tnaive/TCM_2"      = "CD8",
  "CD8_Tnaive/TCM_3"      = "CD8",
  "NK_CD56dim_cytotoxic"  = "NK",
  "NK_KIRpos_CD57pos"     = "NK",
  "NK_CD56bright_cytokine_resp" = "NK",
  "NK_CD56bright_activated"     = "NK",
  "MAIT"                  = "MAIT/gdT",
  "gammaDelta_T_Vd2"      = "MAIT/gdT",
  "InnateLike_PLZF_T"     = "MAIT/gdT"
)
tnk_all$lineage <- unname(lineage_map_tnk[as.character(Idents(tnk_all))])

DimPlot(tnk_all, group.by = "lineage", label = TRUE)
VlnPlot(tnk_all,
        features  = c("Hu.CD8", "Hu.CD4-RPA.T4", "Hu.CD3-UCHT1",
                      "Hu.TCR.Va7.2", "Hu.TCR.Vd2"),
        pt.size   = 0,
        group.by  = "seurat_clusters")

saveRDS(tnk_all, file.path(OUTPUT_DIR, "tnk_all.rds"))


# ---- B. CD8 / NK / MAIT subset ----

cd8nk_subset <- subset(tnk_all, subset = lineage %in% c("CD8", "NK", "MAIT/gdT"))
cd8nk_subset <- prep_rna_pca(cd8nk_subset)
cd8nk_subset <- prep_adt_pca(cd8nk_subset, nfeatures = 134)
cd8nk_subset <- run_wnn(cd8nk_subset, resolution = WNN_RES_CD8NK)

DimPlot(cd8nk_subset, reduction = "wnn.umap", label = TRUE)

mk_cd8nk <- find_markers_and_top(cd8nk_subset, rna_n = 20, adt_n = 10)
print(mk_cd8nk$top_rna, n = Inf)
print(mk_cd8nk$top_adt, n = Inf)

VlnPlot(cd8nk_subset, features = "Hu.CD3-UCHT1", pt.size = 0)
FeaturePlot(cd8nk_subset, features = "Hu.TCR.Vd2", reduction = "wnn.umap")

Idents(cd8nk_subset) <- "seurat_clusters"
hard_map_cd8nk <- c(
  `0`  = "NK_CD56dim_KIR+_CD16hi",
  `1`  = "CD8_Tcm_CCR7+",
  `2`  = "CD8_EM_cytotoxic",
  `3`  = "CD8_GZMK+_Tem",
  `4`  = "NK_CD56bright_activated",
  `5`  = "CD8_KIR+_innate-like",
  `6`  = "CD8_activated_HLA-DR+",
  `7`  = "gamma-delta_T_Vd2",
  `8`  = "CD4_Tcm_contaminant",
  `9`  = "CD8_Tcm_MYC+",
  `10` = "NK_CD11c+",
  `11` = "MAIT_activated",
  `12` = "MAIT_canonical",
  `13` = "CD8_Tcm_LEF1+",
  `14` = "CD8_Tcm_CD62Lhi",
  `15` = "NK_CD56bright_IL12RB2+"
)
cd8nk_subset <- RenameIdents(cd8nk_subset, hard_map_cd8nk)
cd8nk_subset$celltype <- as.character(Idents(cd8nk_subset))

lineage_map_cd8nk <- c(
  "CD8_Tcm_CCR7+"           = "CD8",
  "CD8_Tcm_MYC+"            = "CD8",
  "CD8_Tcm_LEF1+"           = "CD8",
  "CD8_Tcm_CD62Lhi"         = "CD8",
  "CD8_EM_cytotoxic"        = "CD8",
  "CD8_GZMK+_Tem"           = "CD8",
  "CD8_KIR+_innate-like"    = "CD8",
  "CD8_activated_HLA-DR+"   = "CD8",
  "NK_CD56dim_KIR+_CD16hi"  = "NK",
  "NK_CD56bright_activated" = "NK",
  "NK_CD56bright_IL12RB2+"  = "NK",
  "NK_CD11c+"               = "NK",
  "gamma-delta_T_Vd2"       = "γδ_T",
  "MAIT_activated"          = "MAIT",
  "MAIT_canonical"          = "MAIT",
  "CD4_Tcm_contaminant"     = "CD4_T"
)
cd8nk_subset$lineage <- unname(lineage_map_cd8nk[cd8nk_subset$celltype])
cd8nk_subset$lineage <- factor(cd8nk_subset$lineage,
                                levels = c("CD8", "NK", "MAIT", "γδ_T", "CD4_T"))

Idents(cd8nk_subset) <- "lineage"
DimPlot(cd8nk_subset, reduction = "wnn.umap", label = TRUE, repel = TRUE) +
  ggtitle("Refined Lineages: CD8, NK, MAIT, γδ T, CD4") + NoLegend()

saveRDS(cd8nk_subset, file.path(OUTPUT_DIR, "cd8nk_subset.rds"))


# ---- C. NK focus ----

nk_focus <- subset(cd8nk_subset, subset = lineage == "NK")
nk_focus <- prep_rna_pca(nk_focus)
nk_focus <- prep_adt_pca(nk_focus, nfeatures = 134)
nk_focus <- run_wnn(nk_focus, resolution = WNN_RES_NK)

DimPlot(nk_focus, reduction = "wnn.umap", label = TRUE, split.by = "condition")

mk_nk <- find_markers_and_top(nk_focus, rna_n = 20, adt_n = 10)
print(mk_nk$top_rna, n = Inf)
print(mk_nk$top_adt, n = Inf)

FeaturePlot(nk_focus, features = c("Hu.CD3-UCHT1", "Hu.TCR.AB"), pt.size = 0)

nk_labels <- c(
  `0` = "CD57+ KIR+ NK",
  `1` = "CD161+ PLZF+ NK",
  `2` = "HLA-DR+ NK",
  `3` = "CD11c+ TOX+ NK",
  `4` = "CD3+ CD8+ T",         # contaminant
  `5` = "CD56+ GZMK+ NK",
  `6` = "CD161+ PLZF+ NK",
  `7` = "TNF+ IFNG+ NK"
)
nk_focus <- RenameIdents(nk_focus, nk_labels)
nk_focus$celltype <- Idents(nk_focus)
nk_focus <- subset(nk_focus, subset = celltype == "CD3+ CD8+ T", invert = TRUE)

DimPlot(nk_focus, reduction = "wnn.umap", label = TRUE) + NoLegend()
saveRDS(nk_focus, file.path(OUTPUT_DIR, "nk_focus.rds"))


# ---- D. γδ T / MAIT ----

gd_object <- subset(cd8nk_subset, subset = lineage == "γδ_T")
gd_object <- prep_rna_pca(gd_object, scale_vars = "patient")
gd_object <- prep_adt_pca(gd_object, nfeatures = 100, scale_vars = "patient")
gd_object <- run_wnn(gd_object, resolution = WNN_RES_GD)

DimPlot(gd_object, reduction = "wnn.umap", label = TRUE)

mk_gd <- find_markers_and_top(gd_object, rna_n = 30, adt_n = 20,
                               rna_abs = TRUE, adt_abs = TRUE)
print(mk_gd$top_rna, n = Inf)
print(mk_gd$top_adt, n = Inf)

VlnPlot(gd_object,
        features = c("Hu.CD3-UCHT1", "Hu.TCR.AB", "Hu.TCR.Va7.2", "Hu.TCR.Vd2"))

Idents(gd_object) <- "seurat_clusters"

# NOTE: cluster 0 is CD314+ and cluster 1 is CD226+ (marker annotations are
# labelled in the reverse order from what the cluster names imply; retained
# here to match the manuscript figure labels).
gd_ids <- c(
  "0" = "CD226+ GD T",
  "1" = "CD314+ GD T",
  "2" = "CD127+ CD45RO+ GD T"
)
gd_object <- RenameIdents(gd_object, gd_ids)
gd_object$celltype <- Idents(gd_object)

DimPlot(gd_object, reduction = "wnn.umap", label = TRUE, repel = TRUE,
        split.by = "condition") + NoLegend()
table(gd_object$celltype, gd_object$condition, gd_object$day)

VlnPlot(gd_object, features = "HAVCR2", assay = "RNA",
        pt.size = 0, split.by = "condition")

saveRDS(gd_object, file.path(OUTPUT_DIR, "gd_object.rds"))


# ---- E. B cells ----

b_cells <- subset(integ_rna, subset = seurat_clusters %in% c(6, 10, 13, 20))
b_cells <- prep_rna_pca(b_cells)
b_cells <- prep_adt_pca(b_cells, nfeatures = 134)
b_cells <- run_wnn(b_cells, resolution = 0.15)

# Remove contaminating non-B clusters identified after initial clustering
b_cells <- subset(b_cells, subset = seurat_clusters %in% c(3, 4, 6), invert = TRUE)
b_cells <- prep_rna_pca(b_cells)
b_cells <- prep_adt_pca(b_cells, nfeatures = 134)
b_cells <- run_wnn(b_cells, resolution = WNN_RES_B)

DimPlot(b_cells, reduction = "wnn.umap", label = TRUE)

bcell_ids <- c(
  `0` = "CD27+ Memory B",
  `1` = "CD23+ Transitional B",
  `2` = "CD11c+ Memory B",
  `3` = "SLAMF7+ Memory B",
  `4` = "Naive B",
  `5` = "CD69+ B"
)
b_cells <- RenameIdents(b_cells, bcell_ids)
b_cells$celltype <- Idents(b_cells)

saveRDS(b_cells, file.path(OUTPUT_DIR, "b_cells.rds"))


# ---- F. CD8 T cells ----

install_bioc_if_missing("scGate")

cd8_cells <- subset(cd8nk_subset, subset = lineage == "CD8")
cd8_cells <- prep_rna_pca(cd8_cells)
cd8_cells <- prep_adt_pca(cd8_cells, nfeatures = 140)
cd8_cells <- run_wnn(cd8_cells, resolution = WNN_RES_CD8)

DefaultAssay(cd8_cells) <- "RNA"
scGate_models_DB <- scGate::get_scGateDB()
cd8_cells <- scGate::scGate(cd8_cells,
                             model = scGate_models_DB$human$generic$Tcell.alphabeta)
DimPlot(cd8_cells, group.by = "is.pure", label = TRUE, repel = TRUE)
cd8_cells <- subset(cd8_cells, subset = is.pure == "Pure")

cd8_cells <- prep_rna_pca(cd8_cells)
cd8_cells <- prep_adt_pca(cd8_cells, nfeatures = 140)
cd8_cells <- run_wnn(cd8_cells, resolution = WNN_RES_CD8)

cd8_names <- c(
  "0" = "CD8 Naive",
  "1" = "CD8 Effector",
  "2" = "CD8 Effector-Memory",
  "3" = "CD8 Innate-like",
  "4" = "CD8 TEMRA",
  "5" = "CD8 Naive",
  "6" = "CD8 Naive"
)
cd8_cells$celltype <- unname(cd8_names[as.character(cd8_cells$seurat_clusters)])
Idents(cd8_cells) <- "celltype"

DimPlot(cd8_cells, reduction = "wnn.umap", label = TRUE) + NoLegend()
saveRDS(cd8_cells, file.path(OUTPUT_DIR, "cd8_cells.rds"))


# ---- G. CD4 T cells and Tregs ----

CD4_and_treg <- subset(tnk_all, subset = lineage %in% c("Treg", "CD4"))
CD4_and_treg <- prep_rna_pca(CD4_and_treg)
CD4_and_treg <- prep_adt_pca(CD4_and_treg, nfeatures = 134)
CD4_and_treg <- run_wnn(CD4_and_treg, resolution = WNN_RES_CD4)

DimPlot(CD4_and_treg, reduction = "wnn.umap", label = TRUE)

cd4_names <- c(
  "0"  = "CD4 Naive",
  "1"  = "CD4 Naive-early LEF1hi",
  "2"  = "CD4 CCR6+ CXCR5+",
  "3"  = "CD4 IFNG-AS1+ PD-1+",
  "4"  = "CD4 Th2/CCR4+ Effector",
  "5"  = "CD4 Th1/Th17-activated",
  "6"  = "CD4 Th17 CCR6+ IL18R1+",
  "7"  = "CD4 Naive LEF1/BACH2hi",
  "8"  = "Treg Activated",
  "9"  = "Treg Resting",
  "10" = "Cytotoxic/Innate-like"
)
CD4_and_treg$celltype <- unname(cd4_names[as.character(CD4_and_treg$seurat_clusters)])
Idents(CD4_and_treg) <- "celltype"

saveRDS(CD4_and_treg, file.path(OUTPUT_DIR, "CD4_and_treg.rds"))
message("02_subset_analysis.R complete.")
