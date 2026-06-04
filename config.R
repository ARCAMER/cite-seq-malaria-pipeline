# =============================================================================
# config.R — Set all paths and global parameters here before running
# =============================================================================

# ---- Input / output directories ----
RDS_DIR       <- "/datastore/tdo/Honours/data"          # RDS files from HTO demultiplexing
OUTPUT_DIR    <- "output"                               # local outputs (CSVs, RDS checkpoints)
PLOT_DIR      <- "plots"                                # SVG / PNG outputs

# Create output directories if missing
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PLOT_DIR, "umap"),     showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PLOT_DIR, "freq"),     showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PLOT_DIR, "heatmaps"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PLOT_DIR, "milo"),     showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PLOT_DIR, "dotplots"), showWarnings = FALSE, recursive = TRUE)

# ---- Analysis parameters ----
N_PCS_RNA  <- 30          # PCs for RNA integration / DoubletFinder
N_PCS_ADT  <- 50          # PCs for ADT integration
N_HVG_RNA  <- 3000        # highly variable genes for RNA integration
N_HVG_ADT  <- 134         # variable features for ADT integration
DOUBLET_RATE <- 0.05      # expected doublet rate (5 %)

# WNN clustering resolutions
WNN_RES_GLOBAL <- 0.2
WNN_RES_TNK    <- 0.3
WNN_RES_CD8NK  <- 0.35
WNN_RES_NK     <- 0.3
WNN_RES_GD     <- 0.1
WNN_RES_B      <- 0.2
WNN_RES_CD8    <- 0.2
WNN_RES_CD4    <- 0.3

# MAST DEG parameters
MAST_MIN_CELLS   <- 10    # minimum cells per condition to run MAST
MAST_MIN_PCT     <- 0.3
MAST_LOGFC_THRESH <- 0.5

# Heatmap / dotplot gene selection
HEAT_ADJ_P     <- 0.05
HEAT_LOGFC     <- 0.5
HEAT_TOP_N     <- 5
HEAT_MAX_GENES <- 60

# ---- Isotype controls ----
ISOTYPES <- c(
  "Isotype-MOPC.173", "Isotype-MOPC.21",
  "Isotype-MPC.11",   "Isotype-RTK2071",
  "Isotype-RTK2758",  "Isotype-RTK4530"
)

# ---- Patient sex map ----
SEX_MAP <- c(
  "P01" = "female",
  "P02" = "male",
  "P03" = "male",
  "P04" = "male",
  "P05" = "male",
  "P06" = "male",
  "P07" = "female",
  "P08" = "female",
  "P09" = "male",
  "P10" = "male"
)

# ---- Condition colours (for frequency and beeswarm plots) ----
COND_COLS <- c(FEB = "firebrick", CHR = "steelblue")

# ---- Misc ----
options(future.globals.maxSize = 100 * 1024^3)   # 100 GB for parallel processing
