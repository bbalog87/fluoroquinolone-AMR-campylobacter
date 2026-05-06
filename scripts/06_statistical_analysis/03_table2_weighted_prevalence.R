# ============================================================
# 03_weighted_prevalence_table.R
# Final publication-quality prevalence table
#   - Raw prevalence (n/N)
#   - Country-balanced prevalence (mean across countries, within region)
#   - Country-year balanced prevalence (mean across country-year strata, within region)
#
# Reviewer-proof structure (with denominators):
#   Region
#   Genomes (N)
#   Mutant (n)
#   Raw prevalence % (n/N)
#   Countries (k)
#   Country-balanced prevalence % (k)
#   Country-year strata (m)
#   Country-year balanced prevalence % (m)
#
# INPUT
#   data/cleaned_campy_data_deduplicated.xlsx
#
# OUTPUT
#   tables/Table2_weighted_prevalence_real_data.xlsx
#   tables/Table2_weighted_prevalence_real_data.csv
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(writexl)
})

# -------------------------------
# 0) PATH SETUP
# -------------------------------
proj_root  <- "~/PROJECT-2025/Campylobacter-Paper/gyA_analysis/REVISION-1/campy-rev"
data_dir   <- file.path(proj_root, "data")
tables_dir <- file.path(proj_root, "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

infile <- file.path(data_dir, "cleaned_campy_data_deduplicated.xlsx")

# -------------------------------
# 1) LOAD DATA + CLEAN VARIABLES
# -------------------------------
campy <- read_excel(infile) %>%
  mutate(
    gyrA_mut = as.integer(`GyrA-Mutation`),
    Year     = suppressWarnings(as.integer(Year)),
    Country  = as.character(Country),
    Region   = as.character(Region)
  )

# Required columns sanity check
required_cols <- c("Region", "Country", "Year", "gyrA_mut")
missing_cols <- required_cols[!required_cols %in% names(campy)]
if (length(missing_cols) > 0) {
  stop("Missing required columns in input file: ", paste(missing_cols, collapse = ", "))
}

# Mutation sanity check
if (any(!is.na(campy$gyrA_mut) & !campy$gyrA_mut %in% c(0L, 1L))) {
  stop("gyrA_mut contains values other than 0/1/NA. Check 'GyrA-Mutation' coding.")
}

# Remove records missing key stratification fields
campy <- campy %>%
  filter(!is.na(Region), Region != "",
         !is.na(Country), Country != "",
         !is.na(Year))

message("Total genomes analysed (deduplicated): ", nrow(campy))

# -------------------------------
# 2) RAW PREVALENCE (REGION + GLOBAL)
# -------------------------------
raw_region <- campy %>%
  group_by(Region) %>%
  summarise(
    `Genomes (N)` = n(),
    `Mutant (n)`  = sum(gyrA_mut == 1, na.rm = TRUE),
    Raw_Prevalence = `Mutant (n)` / `Genomes (N)`,
    .groups = "drop"
  )

raw_global <- campy %>%
  summarise(
    `Genomes (N)` = n(),
    `Mutant (n)`  = sum(gyrA_mut == 1, na.rm = TRUE),
    Raw_Prevalence = `Mutant (n)` / `Genomes (N)`
  ) %>%
  mutate(Region = "Global")

# -------------------------------
# 3) COUNTRY-BALANCED PREVALENCE (REGION + GLOBAL)
#   Mean prevalence across countries within region.
#   Each country contributes equally (not proportional to #genomes).
# -------------------------------
country_level <- campy %>%
  group_by(Region, Country) %>%
  summarise(
    n_country   = n(),
    mut_country = sum(gyrA_mut == 1, na.rm = TRUE),
    prev_country = mut_country / n_country,
    .groups = "drop"
  )

country_balanced_region <- country_level %>%
  group_by(Region) %>%
  summarise(
    `Countries (k)` = n_distinct(Country),
    Country_Balanced = mean(prev_country, na.rm = TRUE),
    .groups = "drop"
  )

country_balanced_global <- country_level %>%
  summarise(
    `Countries (k)` = n_distinct(Country),
    Country_Balanced = mean(prev_country, na.rm = TRUE)
  ) %>%
  mutate(Region = "Global")

# -------------------------------
# 4) COUNTRY-YEAR BALANCED PREVALENCE (REGION + GLOBAL)
#   Mean prevalence across (Country × Year) strata within region.
#   Each country-year contributes equally.
# -------------------------------
country_year_level <- campy %>%
  group_by(Region, Country, Year) %>%
  summarise(
    n_cy   = n(),
    mut_cy = sum(gyrA_mut == 1, na.rm = TRUE),
    prev_cy = mut_cy / n_cy,
    .groups = "drop"
  )

country_year_balanced_region <- country_year_level %>%
  group_by(Region) %>%
  summarise(
    `Country-year strata (m)` = n(),
    Country_Year_Balanced = mean(prev_cy, na.rm = TRUE),
    .groups = "drop"
  )

country_year_balanced_global <- country_year_level %>%
  summarise(
    `Country-year strata (m)` = n(),
    Country_Year_Balanced = mean(prev_cy, na.rm = TRUE)
  ) %>%
  mutate(Region = "Global")

# -------------------------------
# 5) BUILD FINAL TABLE (FORMATTED)
# -------------------------------
final_table2 <- bind_rows(raw_region, raw_global) %>%
  left_join(bind_rows(country_balanced_region, country_balanced_global), by = "Region") %>%
  left_join(bind_rows(country_year_balanced_region, country_year_balanced_global), by = "Region") %>%
  mutate(
    `Raw prevalence % (n/N)` = paste0(
      round(Raw_Prevalence * 100, 1),
      "% (", `Mutant (n)`, "/", `Genomes (N)`, ")"
    ),
    `Country-balanced prevalence %` = paste0(
      round(Country_Balanced * 100, 1),
      "% (k=", `Countries (k)`, ")"
    ),
    `Country-year balanced prevalence %` = paste0(
      round(Country_Year_Balanced * 100, 1),
      "% (m=", `Country-year strata (m)`, ")"
    )
  ) %>%
  select(
    Region,
    `Genomes (N)`,
    `Mutant (n)`,
    `Raw prevalence % (n/N)`,
    `Countries (k)`,
    `Country-balanced prevalence %`,
    `Country-year strata (m)`,
    `Country-year balanced prevalence %`
  ) %>%
  arrange(desc(`Genomes (N)`))

# -------------------------------
# 6) SAVE
# -------------------------------
out_xlsx <- file.path(tables_dir, "Table2_weighted_prevalence_real_data.xlsx")
out_csv  <- file.path(tables_dir, "Table2_weighted_prevalence_real_data.csv")

write_xlsx(final_table2, out_xlsx)
write.csv(final_table2, out_csv, row.names = FALSE)

# -------------------------------
# 7) CONSOLE OUTPUT
# -------------------------------
print(final_table2)
message("Saved: ", out_xlsx)
message("Saved: ", out_csv)


























# ============================================================
# 03_weighted_prevalence_table.R
# Correct handling of missing Year
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(writexl)
})

proj_root  <- "~/PROJECT-2025/Campylobacter-Paper/gyA_analysis/REVISION-1/campy-rev"
data_dir   <- file.path(proj_root, "data")
tables_dir <- file.path(proj_root, "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

campy0 <- read_excel(file.path(data_dir, "cleaned_campy_data_deduplicated.xlsx")) %>%
  mutate(
    gyrA_mut = as.integer(`GyrA-Mutation`),
    Year     = suppressWarnings(as.integer(Year)),
    Country  = as.character(Country),
    Region   = as.character(Region)
  )

# Keep only what is needed for Region and Country based summaries
campy_RC <- campy0 %>%
  filter(!is.na(Region), Region != "",
         !is.na(Country), Country != "")

# Keep only what is needed for Country-Year balancing
campy_RCY <- campy0 %>%
  filter(!is.na(Region), Region != "",
         !is.na(Country), Country != "",
         !is.na(Year))

message("Rows for Raw and Country-balanced: ", nrow(campy_RC))
message("Rows for Country-Year balanced: ", nrow(campy_RCY))

# -------------------------------
# RAW PREVALENCE (Region + Global)
# -------------------------------
raw_region <- campy_RC %>%
  group_by(Region) %>%
  summarise(
    `Genomes (N)` = n(),
    `Mutant (n)`  = sum(gyrA_mut == 1, na.rm = TRUE),
    Raw_Prevalence = `Mutant (n)` / `Genomes (N)`,
    .groups = "drop"
  )

raw_global <- campy_RC %>%
  summarise(
    `Genomes (N)` = n(),
    `Mutant (n)`  = sum(gyrA_mut == 1, na.rm = TRUE),
    Raw_Prevalence = `Mutant (n)` / `Genomes (N)`
  ) %>%
  mutate(Region = "Global")

# -------------------------------
# COUNTRY-BALANCED (Region + Global)
# -------------------------------
country_level <- campy_RC %>%
  group_by(Region, Country) %>%
  summarise(
    n_country   = n(),
    mut_country = sum(gyrA_mut == 1, na.rm = TRUE),
    prev_country = mut_country / n_country,
    .groups = "drop"
  )

country_balanced_region <- country_level %>%
  group_by(Region) %>%
  summarise(
    `Countries (k)` = n_distinct(Country),
    Country_Balanced = mean(prev_country, na.rm = TRUE),
    .groups = "drop"
  )

country_balanced_global <- country_level %>%
  summarise(
    `Countries (k)` = n_distinct(Country),
    Country_Balanced = mean(prev_country, na.rm = TRUE)
  ) %>%
  mutate(Region = "Global")

# -------------------------------
# COUNTRY-YEAR BALANCED (Region + Global)
# -------------------------------
country_year_level <- campy_RCY %>%
  group_by(Region, Country, Year) %>%
  summarise(
    n_cy   = n(),
    mut_cy = sum(gyrA_mut == 1, na.rm = TRUE),
    prev_cy = mut_cy / n_cy,
    .groups = "drop"
  )

country_year_balanced_region <- country_year_level %>%
  group_by(Region) %>%
  summarise(
    `Country-year strata (m)` = n(),
    Country_Year_Balanced = mean(prev_cy, na.rm = TRUE),
    .groups = "drop"
  )

country_year_balanced_global <- country_year_level %>%
  summarise(
    `Country-year strata (m)` = n(),
    Country_Year_Balanced = mean(prev_cy, na.rm = TRUE)
  ) %>%
  mutate(Region = "Global")

# -------------------------------
# FINAL TABLE
# -------------------------------
final_table2 <- bind_rows(raw_region, raw_global) %>%
  left_join(bind_rows(country_balanced_region, country_balanced_global), by = "Region") %>%
  left_join(bind_rows(country_year_balanced_region, country_year_balanced_global), by = "Region") %>%
  mutate(
    `Raw prevalence % (n/N)` = paste0(
      round(Raw_Prevalence * 100, 1),
      "% (", `Mutant (n)`, "/", `Genomes (N)`, ")"
    ),
    `Country-balanced prevalence %` = paste0(
      round(Country_Balanced * 100, 1),
      "% (k=", `Countries (k)`, ")"
    ),
    `Country-year balanced prevalence %` = ifelse(
      is.na(`Country-year strata (m)`),
      NA_character_,
      paste0(round(Country_Year_Balanced * 100, 1),
             "% (m=", `Country-year strata (m)`, ")")
    )
  ) %>%
  select(
    Region,
    `Genomes (N)`,
    `Mutant (n)`,
    `Raw prevalence % (n/N)`,
    `Countries (k)`,
    `Country-balanced prevalence %`,
    `Country-year strata (m)`,
    `Country-year balanced prevalence %`
  ) %>%
  arrange(desc(`Genomes (N)`))

out_xlsx <- file.path(tables_dir, "Table2_weighted_prevalence_real_data2.xlsx")
out_csv  <- file.path(tables_dir, "Table2_weighted_prevalence_real_data2.csv")

write_xlsx(final_table2, out_xlsx)
write.csv(final_table2, out_csv, row.names = FALSE)

print(final_table2)
message("Saved: ", out_xlsx)
message("Saved: ", out_csv)




# ============================================================
# Species composition by region (donut plot)
# Based on real deduplicated dataset
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

library(patchwork)
library(ggplot2)
library(dplyr)

# -------------------------------
# PATH
# -------------------------------

proj_root <- "~/PROJECT-2025/Campylobacter-Paper/gyA_analysis/REVISION-1/campy-rev"
data_dir  <- file.path(proj_root, "data")
fig_dir   <- file.path(proj_root, "figures")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# LOAD REAL DATA
# -------------------------------

campy <- read_excel(
  file.path(data_dir, "cleaned_campy_data_deduplicated.xlsx")
)

# -------------------------------
# CLEAN & CLASSIFY SPECIES
# -------------------------------

campy_clean <- campy %>%
  filter(!is.na(Region), !is.na(Species)) %>%
  mutate(
    SpeciesGroup = case_when(
      Species == "Campylobacter jejuni" ~ "C. jejuni",
      Species == "Campylobacter coli" ~ "C. coli",
      TRUE ~ "Other Campylobacter"
    )
  )

# -------------------------------
# SUMMARISE PER REGION
# -------------------------------

species_summary <- campy_clean %>%
  group_by(Region, SpeciesGroup) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Region) %>%
  mutate(
    Total = sum(Count),
    Percent = Count / Total * 100
  ) %>%
  arrange(Region, desc(Percent)) %>%
  group_by(Region) %>%
  mutate(
    ymax = cumsum(Percent),
    ymin = lag(ymax, default = 0),
    label_pos = (ymax + ymin) / 2,
    Label = paste0(round(Percent, 1), "%")
  )

# -------------------------------
# CUSTOM COLORS (journal-safe)
# -------------------------------

species_colors <- c(
  "C. jejuni" = "#1b9e77",
  "C. coli"   = "#d95f02",
  "Other Campylobacter" = "#7570b3"
)

# -------------------------------
# DONUT PLOT
# -------------------------------

p_species <- ggplot(species_summary) +
  
  geom_rect(aes(
    ymin = ymin,
    ymax = ymax,
    xmin = 3,
    xmax = 4,
    fill = SpeciesGroup
  ),
  color = "white",
  linewidth = 0.5) +
  
  coord_polar(theta = "y") +
  facet_wrap(~Region, ncol = 3) +
  xlim(c(0, 4)) +
  
  # Center total n
  geom_text(
    data = species_summary %>% distinct(Region, Total),
    aes(x = 0, y = 0, label = paste0("n = ", Total)),
    size = 5,
    fontface = "bold"
  ) +
  
  # Percent labels
 # geom_text(
   # aes(x = 3.5, y = label_pos, label = Label),
   # size = 3.8,
   # fontface = "bold"
  #) +
  
  scale_fill_manual(values = species_colors) +
  
  theme_void(base_size = 16) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.margin = margin(10,10,10,10)
  )

# Save high-resolution figure
ggsave(
  file.path(fig_dir, "Figure_species_distribution_by_region.png"),
  p_species,
  width = 12,
  height = 8,
  dpi = 600
)

p_species




# --------------------------------------------------
# Ensure plots exist:
# p_species  -> Donut
# pA         -> Regional prevalence
# pB         -> Temporal trend
# --------------------------------------------------

final_figure <- 
  (p_species | pA) /
  pB +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        size = 18,
        face = "bold"
      )
    )
  )

final_figure





# ============================================================
# FINAL PUBLICATION-QUALITY MULTI-PANEL FIGURE
# Layout:
#   A = Species distribution donut
#   B = Regional gyrA prevalence
#   C = Temporal trend
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

# ------------------------------------------------------------
# HIGH-READABILITY THEME
# ------------------------------------------------------------

pub_theme <- theme_minimal(base_size = 24) +
  theme(
    axis.text = element_text(size = 20, face = "bold", color = "black"),
    axis.title = element_text(size = 22, face = "bold"),
    strip.text = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 24, face = "bold"),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------
# APPLY THEME TO PANELS
# ------------------------------------------------------------

p_species <- p_species + pub_theme
pA        <- pA + pub_theme
pB        <- pB +
  pub_theme +
  theme(
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 20, face = "bold")
  )

# ------------------------------------------------------------
# BUILD MULTI-PANEL LAYOUT
# ------------------------------------------------------------

final_figure <- 
  (p_species | pA) /
  pB +
  plot_layout(heights = c(1, 1.3)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        size = 30,
        face = "bold"
      )
    )
  )

# ------------------------------------------------------------
# SAVE VECTOR PUBLICATION PDF
# ------------------------------------------------------------

ggsave(
  filename = "Figure1_Final_Publication.pdf",
  plot = final_figure,
  width = 190/25.4,   # 190 mm width (full journal page)
  height = 240/25.4,  # tall layout
  units = "in",
  device = cairo_pdf
)

# Also save high-res PNG if needed
ggsave(
  filename = "Figure1_Final_Publication.png",
  plot = final_figure,
  width = 190/25.4,
  height = 240/25.4,
  units = "in",
  dpi = 600
)

# Print to viewer
final_figure



###############################################


suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(scales)
})

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------

df <- read_excel("CC_resitance.xlsx")

df_clean <- df %>%
  rename(
    CC = ClonalComplex,
    HostGroup = Host.Group
  ) %>%
  mutate(
    CC = ifelse(is.na(CC) | CC == "", "Orphan", CC)
  ) %>%
  filter(!is.na(HostGroup))

# Keep only the 3 host categories
df_clean <- df_clean %>%
  filter(HostGroup %in% c("Human", "Avian", "Nonhuman Mammal"))

# ------------------------------------------------------------
# IDENTIFY TOP 10 TRUE CC (excluding Orphan)
# ------------------------------------------------------------

top10_cc <- df_clean %>%
  filter(CC != "Orphan") %>%
  count(CC, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(CC)

# ------------------------------------------------------------
# MERGE RARE + ORPHAN INTO "Other"
# ------------------------------------------------------------

df_grouped <- df_clean %>%
  mutate(
    CC_grouped = ifelse(CC %in% top10_cc, CC, "Other")
  )

# ------------------------------------------------------------
# SUMMARISE HOST PROPORTIONS
# ------------------------------------------------------------

cc_summary <- df_grouped %>%
  group_by(CC_grouped, HostGroup) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CC_grouped) %>%
  mutate(
    Total = sum(n),
    Proportion = n / Total
  ) %>%
  ungroup()

# ------------------------------------------------------------
# ORDER BARS BY HUMAN COMPOSITION (DESCENDING)
# ------------------------------------------------------------

human_order <- cc_summary %>%
  filter(HostGroup == "Human") %>%
  arrange(desc(Proportion)) %>%
  pull(CC_grouped)

cc_summary$CC_grouped <- factor(
  cc_summary$CC_grouped,
  levels = human_order
)

# ------------------------------------------------------------
# NATURE MICROBIOLOGY STYLE PALETTE
# ------------------------------------------------------------

nm_palette <- c(
  "Human" = "#F58518",
  "Avian" = "#4C78A8",
  "Nonhuman Mammal" = "#54A24B"
)


nm_palette2 <- c(
  "Human" = "#1b9e77",
  "Avian" = "#d95f02",
  "Nonhuman Mammal" = "#7570b3"
)


# ------------------------------------------------------------
# PLOT
# ------------------------------------------------------------

p_cc <- ggplot(cc_summary,
               aes(x = CC_grouped,
                   y = Proportion,
                   fill = HostGroup)) +
  geom_col(width = 0.7,
           color = "white",
           linewidth = 0.3) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0,1),
    expand = c(0,0)
  ) +
  scale_fill_manual(values = nm_palette2) +
  labs(
    x = NULL,
    y = "Host composition within clonal complex (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank()
  )

p_cc




# ------------------------------------------------------------
# EXPORT – Double-column Nature size
# ------------------------------------------------------------

ggsave(
  "Figure_CC_Top10_Host_NatureStyle2.pdf",
  p_cc,
  width = 183/25.4,
  height = 110/25.4,
  units = "in",
  device = cairo_pdf
)

ggsave(
  "Figure_CC_Top10_Host_NatureStyle.png",
  p_cc,
  width = 183/25.4,
  height = 110/25.4,
  units = "in",
  dpi = 600
)









