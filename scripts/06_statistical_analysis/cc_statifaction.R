# ============================================================
# 04_cc_stratification.R
# Purpose:
#   Identify clonal complexes (CCs) with sufficient data
#   for downstream evolutionary and phylogenetic analyses
#
# Inputs:
#   - cleaned_campy_data.xlsx
#   - ST_data_with_CC.xlsx
#
# Outputs:
#   - cc_summary_all.tsv
#   - cc_summary_eligible.tsv
#
# Author: Julien A. Nguinkal
# ============================================================

# ============================================================
# PROJECT PATH SETUP
# ============================================================

proj_root <- "~/PROJECT-2025/Campylobacter-Paper/gyA_analysis/REVISION-1/campy-rev"

data_dir     <- file.path(proj_root, "data")
scripts_dir  <- file.path(proj_root, "scripts")
tables_dir   <- file.path(proj_root, "tables")
figures_dir  <- file.path(proj_root, "figures")
models_dir   <- file.path(proj_root, "models")
phylo_dir    <- file.path(proj_root, "phylogeny")
gyrA_dir     <- file.path(proj_root, "gyrA")
mlst_dir     <- file.path(proj_root, "mlst_cc")
bias_dir     <- file.path(proj_root, "sampling_bias")
meta_dir     <- file.path(proj_root, "metadata")
qc_dir       <- file.path(proj_root, "qc")
docs_dir     <- file.path(proj_root, "docs")
revlog_dir   <- file.path(proj_root, "revision_log")

# Create folders if missing (only for output folders, not raw data or scripts)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(phylo_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(bias_dir, showWarnings = FALSE, recursive = TRUE)


## ============================================================
# 04_cc_stratification.R
# Identify clonal complexes with sufficient data
# ============================================================
#### Package loading
library(readxl)
library(dplyr)
library(writexl)

# -------------------------------
# PATH SETUP
# -------------------------------

proj_root <- "~/PROJECT-2025/Campylobacter-Paper/gyA_analysis/REVISION-1/campy-rev"

data_dir   <- file.path(proj_root, "data")
tables_dir <- file.path(proj_root, "tables")
mlst_dir   <- file.path(proj_root, "mlst_cc")

dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# LOAD DATA
# -------------------------------

campy <- read_excel(file.path(data_dir, "cleaned_campy_data.xlsx"))

st_cc <- read_excel(file.path(data_dir, "ST_data_with_CC.xlsx")) %>%
  select(ST, ClonalComplex_Assigned) %>%
  distinct() %>%
  rename(CC = ClonalComplex_Assigned)

# -------------------------------
# MERGE CC INTO MASTER TABLE
# -------------------------------

campy_cc <- campy %>%
  left_join(st_cc, by = "ST") %>%
  distinct(Isolate, .keep_all = TRUE)

# -------------------------------
# CC SUMMARY
# -------------------------------

cc_summary <- campy_cc %>%
  filter(!is.na(CC)) %>%
  group_by(Species, CC) %>%
  summarise(
    N_genomes  = n(),
    N_gyrA_mut = sum(GyrA_Mutation == 1, na.rm = TRUE),
    Prop_gyrA  = round(N_gyrA_mut / N_genomes, 3),
    .groups = "drop"
  )

# -------------------------------
# ELIGIBILITY FILTER
# -------------------------------

cc_eligible <- cc_summary %>%
  filter(
    N_genomes >= 20,
    N_gyrA_mut >= 5
  )

# -------------------------------
# SAVE OUTPUTS
# -------------------------------

write_xlsx(
  list(
    all_ccs      = cc_summary,
    eligible_ccs = cc_eligible
  ),
  file.path(tables_dir, "cc_stratification_summary.xlsx")
)

write.table(
  cc_eligible,
  file.path(mlst_dir, "eligible_ccs.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# -------------------------------
# REPORT
# -------------------------------

message("CC stratification complete.")
message("Eligible CC counts by species:")
print(table(cc_eligible$Species))



nrow(campy)
nrow(campy_cc)


campy %>%
  count(Isolate) %>%
  filter(n > 1) %>%
  arrange(desc(n)) %>%
  head(20)

campy %>%
  left_join(st_cc, by = "ST") %>%
  count(Isolate) %>%
  filter(n > 1) %>%
  arrange(desc(n)) %>%
  head(20)


st_cc_raw <- read_excel(file.path(data_dir, "ST_data_with_CC.xlsx"))

st_cc_raw %>%
  count(ST) %>%
  filter(n > 1) %>%
  arrange(desc(n))






















# ============================================================
# 00_deduplicate_master_genomes.R
# Collapse duplicated genome accessions and log removals
# ============================================================

library(readxl)
library(dplyr)
library(writexl)

# -------------------------------
# PATH SETUP
# -------------------------------

proj_root <- "~/PROJECT-2025/Campylobacter-Paper/gyA_analysis/REVISION-1/campy-rev"

data_dir   <- file.path(proj_root, "data")
revlog_dir <- file.path(proj_root, "revision_log")

dir.create(revlog_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# LOAD MASTER DATA
# -------------------------------

campy_raw <- read_excel(file.path(data_dir, "cleaned_campy_data.xlsx"))

# -------------------------------
# IDENTIFY DUPLICATES
# -------------------------------

dup_table <- campy_raw %>%
  count(Isolate, name = "n_rows") %>%
  filter(n_rows > 1)

# Extract all duplicated rows
dup_rows <- campy_raw %>%
  semi_join(dup_table, by = "Isolate") %>%
  arrange(Isolate)

# -------------------------------
# DEDUPLICATE (KEEP FIRST)
# -------------------------------

campy_dedup <- campy_raw %>%
  arrange(Isolate) %>%
  distinct(Isolate, .keep_all = TRUE)

# -------------------------------
# SAVE LOGS
# -------------------------------

write_xlsx(
  list(
    duplicated_accessions = dup_table,
    duplicated_rows       = dup_rows
  ),
  file.path(revlog_dir, "duplicate_genomes_audit.xlsx")
)

write_xlsx(
  campy_dedup,
  file.path(data_dir, "cleaned_campy_data_deduplicated.xlsx")
)

# -------------------------------
# REPORT
# -------------------------------

message("Initial genomes: ", nrow(campy_raw))
message("After deduplication: ", nrow(campy_dedup))
message("Removed duplicates: ", nrow(campy_raw) - nrow(campy_dedup))

