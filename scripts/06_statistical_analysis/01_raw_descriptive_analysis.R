# ============================================================
# 01_raw_descriptive_analysis.R
# RAW descriptive analysis of gyrA mutations
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(writexl)
  library(scales)
})

# -------------------------------
# PATH SETUP
# -------------------------------

proj_root <- "~/PROJECT-2025/Campylobacter-Paper/gyA_analysis/REVISION-1/campy-rev"

data_dir   <- file.path(proj_root, "data")
tables_dir <- file.path(proj_root, "tables")
fig_dir    <- file.path(proj_root, "figures", "raw_descriptive")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# LOAD DATA
# -------------------------------

campy <- read_excel(
  file.path(data_dir, "cleaned_campy_data_deduplicated.xlsx")
)

# Create clean mutation variable
campy <- campy %>%
  mutate(
    gyrA_mut = as.integer(`GyrA-Mutation`)
  )

# Sanity check
if(any(!campy$gyrA_mut %in% c(0,1,NA))){
  stop("gyrA mutation column contains unexpected values.")
}

message("Total genomes analysed: ", nrow(campy))

# ============================================================
# 1. OVERALL PREVALENCE
# ============================================================

overall <- campy %>%
  summarise(
    Total = n(),
    Mutant = sum(gyrA_mut == 1, na.rm = TRUE),
    Prevalence = Mutant / Total
  )

print(overall)

# ============================================================
# 2. SPECIES-LEVEL PREVALENCE
# ============================================================

species_prev <- campy %>%
  group_by(Species) %>%
  summarise(
    Total = n(),
    Mutant = sum(gyrA_mut == 1, na.rm = TRUE),
    Prevalence = Mutant / Total,
    .groups = "drop"
  ) %>%
  arrange(desc(Prevalence))

# ============================================================
# 3. CONTINENT-LEVEL PREVALENCE
# ============================================================

continent_prev <- campy %>%
  group_by(Region) %>%
  summarise(
    Total = n(),
    Mutant = sum(gyrA_mut == 1, na.rm = TRUE),
    Prevalence = Mutant / Total,
    .groups = "drop"
  ) %>%
  arrange(desc(Prevalence))

# ============================================================
# 4. CC-LEVEL PREVALENCE
# ============================================================

cc_prev <- campy %>%
  filter(!is.na(ClonalComplex_Assigned)) %>%
  group_by(Species, ClonalComplex_Assigned) %>%
  summarise(
    Total = n(),
    Mutant = sum(gyrA_mut == 1, na.rm = TRUE),
    Prevalence = Mutant / Total,
    .groups = "drop"
  ) %>%
  arrange(desc(Total))

# ============================================================
# 5. TEMPORAL TREND (with denominators)
# ============================================================

year_prev <- campy %>%
  group_by(Year) %>%
  summarise(
    Total = n(),
    Mutant = sum(gyrA_mut == 1, na.rm = TRUE),
    Prevalence = Mutant / Total,
    .groups = "drop"
  ) %>%
  filter(Total >= 20)   # avoid unstable early sparse years

# ============================================================
# SAVE TABLES
# ============================================================

write_xlsx(
  list(
    overall = overall,
    species = species_prev,
    continent = continent_prev,
    cc = cc_prev,
    yearly = year_prev
  ),
  file.path(tables_dir, "raw_descriptive_summary.xlsx")
)

# ============================================================
# FIGURE 1 — MAIN DESCRIPTIVE MULTI-PANEL
# ============================================================

# A) Continent barplot
p1 <- continent_prev %>%
  ggplot(aes(x = reorder(Region, Prevalence), y = Prevalence)) +
  geom_col(fill = "#2C7BB6") +
  scale_y_continuous(labels = percent_format()) +
  coord_flip() +
  labs(
    x = NULL,
    y = "gyrA mutant prevalence (%)",
    title = "Fluoroquinolone-associated gyrA mutations by continent"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(fig_dir, "Figure1A_continent_prevalence.png"),
       p1, width = 7, height = 5, dpi = 300)

# B) Temporal trend
p2 <- ggplot(year_prev, aes(x = Year, y = Prevalence)) +
  geom_line(color = "#D7191C", size = 1.2) +
  geom_point(size = 1.8) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Year",
    y = "gyrA mutant prevalence (%)",
    title = "Temporal dynamics of gyrA mutations"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(fig_dir, "Figure1B_temporal_trend.png"),
       p2, width = 8, height = 4.5, dpi = 300)

message("Raw descriptive analysis complete.")
