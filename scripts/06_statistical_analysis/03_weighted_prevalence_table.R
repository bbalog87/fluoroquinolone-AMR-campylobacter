# ============================================================
# 03_weighted_prevalence_table.R
# Final publication-quality weighted prevalence table
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(writexl)
})

# -------------------------------
# PATH SETUP
# -------------------------------

proj_root <- "~/PROJECT-2025/Campylobacter-Paper/gyA_analysis/REVISION-1/campy-rev"
data_dir  <- file.path(proj_root, "data")
tables_dir <- file.path(proj_root, "tables")

dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# LOAD REAL DATA
# -------------------------------

campy <- read_excel(
  file.path(data_dir, "cleaned_campy_data_deduplicated.xlsx")
) %>%
  mutate(
    gyrA_mut = as.integer(`GyrA-Mutation`)
  )

# Sanity check
if(any(!campy$gyrA_mut %in% c(0,1,NA))){
  stop("gyrA mutation column contains unexpected values.")
}

# ============================================================
# STEP 1 — RAW REGIONAL PREVALENCE
# ============================================================

raw_region <- campy %>%
  group_by(Region) %>%
  summarise(
    Genomes = n(),
    Mutant = sum(gyrA_mut == 1, na.rm = TRUE),
    Raw_Prevalence = Mutant / Genomes,
    .groups = "drop"
  )

# ============================================================
# STEP 2 — COUNTRY-BALANCED PREVALENCE
# ============================================================

country_level <- campy %>%
  group_by(Region, Country) %>%
  summarise(
    n_country = n(),
    mut_country = sum(gyrA_mut == 1, na.rm = TRUE),
    prev_country = mut_country / n_country,
    .groups = "drop"
  )

country_balanced <- country_level %>%
  group_by(Region) %>%
  summarise(
    Country_Balanced = mean(prev_country),
    n_countries = n(),
    .groups = "drop"
  )

# ============================================================
# STEP 3 — COUNTRY-YEAR BALANCED
# ============================================================

country_year_level <- campy %>%
  group_by(Region, Country, Year) %>%
  summarise(
    n_cy = n(),
    mut_cy = sum(gyrA_mut == 1, na.rm = TRUE),
    prev_cy = mut_cy / n_cy,
    .groups = "drop"
  )

country_year_balanced <- country_year_level %>%
  group_by(Region) %>%
  summarise(
    Country_Year_Balanced = mean(prev_cy),
    .groups = "drop"
  )

# ============================================================
# STEP 4 — BUILD FINAL TABLE
# ============================================================

final_table <- raw_region %>%
  left_join(country_balanced, by = "Region") %>%
  left_join(country_year_balanced, by = "Region") %>%
  mutate(
    Raw_Prevalence_pct = round(Raw_Prevalence * 100, 1),
    Country_Balanced_pct = round(Country_Balanced * 100, 1),
    Country_Year_Balanced_pct = round(Country_Year_Balanced * 100, 1),
    Raw_Label = paste0(Raw_Prevalence_pct, "% (", Mutant, "/", Genomes, ")")
  ) %>%
  select(
    Region,
    Genomes,
    Mutant,
    Raw_Label,
    Country_Balanced_pct,
    Country_Year_Balanced_pct
  ) %>%
  arrange(desc(Genomes))

# ============================================================
# SAVE TABLE
# ============================================================

write_xlsx(
  final_table,
  file.path(tables_dir, "Table2_weighted_prevalence_real_data.xlsx")
)

write.csv(
  final_table,
  file.path(tables_dir, "Table2_weighted_prevalence_real_data.csv"),
  row.names = FALSE
)

print(final_table)

message("Final weighted prevalence table generated using real dataset.")
