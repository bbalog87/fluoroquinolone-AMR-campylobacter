# ============================================================
# 01_cc_landscape_deduplicated.R
# Clonal complex (CC) landscape after genome deduplication
#
# PURPOSE
#   1) Attach clonal complex (CC) assignments to each genome using ST → CC mapping
#   2) Quantify CC diversity and size distribution per species
#   3) Identify “analysis-ready” CCs using minimum-data thresholds
#
# WHY THIS STAGE MATTERS
#   - Prevents overinterpretation of rare lineages (small CCs)
#   - Provides denominators and transparency for downstream phylogeny/subsampling
#   - Produces Supplementary tables and a figure for CC size distribution
#
# INPUTS (expected under data/)
#   - cleaned_campy_data_deduplicated.xlsx
#       One row per genome (Isolate accession unique), includes:
#       Isolate, Species, ST, and GyrA-Mutation (0/1)
#   - ST_data_with_CC.xlsx
#       Contains ST and clonal complex assignment (may have repeated ST entries)
#
# OUTPUTS
#   Tables:
#     - tables/cc_landscape_after_dedup.xlsx
#     - mlst_cc/eligible_ccs_after_dedup.xlsx
#     - mlst_cc/eligible_ccs_after_dedup.tsv
#   Figure:
#     - figures/cc_size_distribution_log10.png
#
# SUGGESTED THRESHOLDS (can be justified in Methods)
#   - Keep CCs with >= 20 genomes AND >= 5 gyrA mutants per species
#     Rationale: avoids unstable estimates from sparse CCs
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(writexl)
  library(ggplot2)
})

# -------------------------------
# 0) PATH SETUP
# -------------------------------

proj_root <- "~/PROJECT-2025/Campylobacter-Paper/gyA_analysis/REVISION-1/campy-rev"

data_dir   <- file.path(proj_root, "data")
tables_dir <- file.path(proj_root, "tables")
fig_dir    <- file.path(proj_root, "figures")
mlst_dir   <- file.path(proj_root, "mlst_cc")

dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(mlst_dir,   showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# 1) LOAD DEDUPLICATED MASTER DATA
# -------------------------------

campy <- read_excel(file.path(data_dir, "cleaned_campy_data_deduplicated.xlsx"))

# Sanity check: confirm the key columns exist
required_cols <- c("Isolate", "Species", "ST", "GyrA-Mutation")
missing_cols <- required_cols[!required_cols %in% names(campy)]
if (length(missing_cols) > 0) {
  stop("Missing required columns in cleaned_campy_data_deduplicated.xlsx: ",
       paste(missing_cols, collapse = ", "))
}

# IMPORTANT:
# The column name "GyrA-Mutation" contains a hyphen, which requires backticks in R.
# To avoid repeating backticks throughout the pipeline, we create a standard column:
campy <- campy %>%
  mutate(
    ST = as.character(ST),
    gyrA_mut = as.integer(`GyrA-Mutation`)   # enforce 0/1
  )

# Optional: confirm gyrA_mut is only 0/1/NA
if (any(!is.na(campy$gyrA_mut) & !campy$gyrA_mut %in% c(0L, 1L))) {
  stop("gyrA_mut contains values other than 0/1. Check 'GyrA-Mutation' coding.")
}

# -------------------------------
# 2) LOAD AND CLEAN ST → CC TABLE
# -------------------------------
# Your ST_data_with_CC.xlsx has repeated ST entries (normal if it was genome-level).
# We collapse it to a clean lookup: 1 row per ST.
#
# NOTE: The warnings you saw earlier were due to Excel parsing of some columns.
# We explicitly read only the two needed columns and force them to character.

st_cc_raw <- read_excel(file.path(data_dir, "ST_data_with_CC.xlsx"))

# Validate expected column names exist
if (!all(c("ST", "ClonalComplex_Assigned") %in% names(st_cc_raw))) {
  stop("ST_data_with_CC.xlsx must contain columns: ST and ClonalComplex_Assigned.")
}

st_cc <- st_cc_raw %>%
  transmute(
    ST = as.character(ST),
    CC = as.character(ClonalComplex_Assigned)
  ) %>%
  filter(!is.na(ST), ST != "") %>%
  distinct()

# Optional integrity check:
# If an ST maps to multiple CCs, that is biologically or curation-ambiguous.
# We flag it explicitly so it can be discussed or resolved.
st_to_multi_cc <- st_cc %>%
  group_by(ST) %>%
  summarise(n_cc = n_distinct(CC), .groups = "drop") %>%
  filter(n_cc > 1)

if (nrow(st_to_multi_cc) > 0) {
  message("WARNING: Some STs map to multiple CCs (ambiguous mapping).")
  message("These STs will still merge, but you should inspect them for curation issues.")
  # Save a diagnostic table
  write_xlsx(
    list(ambiguous_ST_to_CC = st_to_multi_cc),
    file.path(mlst_dir, "diagnostic_ambiguous_ST_to_CC.xlsx")
  )
}

# -------------------------------
# 3) MERGE CC INTO MASTER DATA
# -------------------------------
# LEFT JOIN keeps all genomes even if CC is missing for some STs.
# This is important to avoid silently dropping genomes.

campy_cc <- campy %>%
  left_join(st_cc, by = "ST")

# Sanity: join should not change number of genomes
stopifnot(nrow(campy_cc) == nrow(campy))

# -------------------------------
# 4) SUMMARISE CC LANDSCAPE
# -------------------------------
# We restrict to rows with CC assigned for CC-specific analyses.
# Genomes lacking CC can still be used elsewhere (e.g., species-level analyses),
# but not for CC-stratified inference.

cc_summary <- campy_cc %>%
  filter(!is.na(CC), CC != "") %>%
  group_by(Species, CC) %>%
  summarise(
    N_genomes      = n(),
    N_gyrA_mut     = sum(gyrA_mut == 1, na.rm = TRUE),
    Prop_gyrA_mut  = round(N_gyrA_mut / N_genomes, 3),
    .groups = "drop"
  ) %>%
  arrange(Species, desc(N_genomes))

# -------------------------------
# 5) FILTER “ELIGIBLE” CCs FOR ROBUST DOWNSTREAM ANALYSES
# -------------------------------
# Suggested thresholds:
#   >= 20 genomes per CC per species
#   >= 5 gyrA mutants per CC per species
#
# You can tune these thresholds later, but keep them fixed once decided,
# and report them explicitly in Methods.

min_genomes <- 20
min_mutants <- 5

cc_eligible <- cc_summary %>%
  filter(N_genomes >= min_genomes, N_gyrA_mut >= min_mutants) %>%
  arrange(Species, desc(N_genomes))

# -------------------------------
# 6) SAVE OUTPUT TABLES
# -------------------------------

write_xlsx(
  list(
    cc_landscape_all = cc_summary,
    cc_eligible      = cc_eligible
  ),
  file.path(tables_dir, "cc_landscape_after_dedup.xlsx")
)

write_xlsx(
  cc_eligible,
  file.path(mlst_dir, "eligible_ccs_after_dedup.xlsx")
)

write.table(
  cc_eligible,
  file.path(mlst_dir, "eligible_ccs_after_dedup.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# -------------------------------
# 7) FIGURE: CC SIZE DISTRIBUTION
# -------------------------------
# We plot CC size distribution (all species pooled, but CC rows are species-stratified).
# Log10 x-axis because CC sizes are typically heavy-tailed.

p <- ggplot(cc_summary, aes(x = N_genomes)) +
  geom_histogram(bins = 40) +
  scale_x_log10() +
  labs(
    x = "Genomes per CC (log10 scale)",
    y = "Number of CCs",
    title = "Clonal complex size distribution after genome deduplication"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(fig_dir, "cc_size_distribution_log10.png"),
  plot = p,
  width = 7,
  height = 5,
  dpi = 300
)

# -------------------------------
# 8) CONSOLE REPORT (FOR REVISION LOGGING)
# -------------------------------

message("CC landscape summary complete.")
message("Total genomes (deduplicated): ", nrow(campy))
message("Genomes with CC assigned: ", sum(!is.na(campy_cc$CC) & campy_cc$CC != ""))

message("Total CC rows (Species × CC): ", nrow(cc_summary))
message("Unique CC overall (ignoring species): ", n_distinct(cc_summary$CC))

message("Eligible CC rows (Species × CC): ", nrow(cc_eligible))
message("Eligible CC rows by species:")
print(table(cc_eligible$Species))

# End of script
