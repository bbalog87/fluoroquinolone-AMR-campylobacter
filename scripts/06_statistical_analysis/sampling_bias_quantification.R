# ===============================
# 02_sampling_bias_quantification.R
# ===============================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(writexl)
  library(ineq)
})

# -------------------------------
# PATH SETUP
# -------------------------------

proj_root <- "~/PROJECT-2025/Campylobacter-Paper/gyA_analysis/REVISION-1/campy-rev"

data_dir   <- file.path(proj_root, "data")
fig_dir    <- file.path(proj_root, "figures", "sampling_bias")
tables_dir <- file.path(proj_root, "tables")
bias_dir   <- file.path(proj_root, "sampling_bias")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bias_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# LOAD DATA
# -------------------------------

campy <- read_excel(
  file.path(data_dir, "cleaned_campy_data_deduplicated.xlsx")
)

message("Total genomes analysed: ", nrow(campy))

# -------------------------------
# COUNTRY × YEAR MATRIX
# -------------------------------

country_year <- campy %>%
  count(Country, Year, name = "N_genomes") %>%
  arrange(desc(N_genomes)) %>% filter(Year>1990)

write.table(
  country_year,
  file.path(bias_dir, "country_year_counts.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# -------------------------------
# COUNTRY TOTALS + GINI
# -------------------------------

country_totals <- campy %>%
  count(Country, name = "N_total") %>%
  arrange(desc(N_total))

gini_country <- ineq(country_totals$N_total, type = "Gini")

# -------------------------------
# REGION TOTALS
# -------------------------------

region_totals <- campy %>%
  count(Region, name = "N_total") %>%
  arrange(desc(N_total))

# -------------------------------
# YEAR TOTALS
# -------------------------------

year_totals <- campy %>%
  count(Year, name = "N_total") %>%
  arrange(Year)

# -------------------------------
# SAVE TABLES
# -------------------------------

write_xlsx(
  list(
    country_totals = country_totals,
    region_totals  = region_totals,
    year_totals    = year_totals,
    country_year   = country_year,
    gini_country   = tibble(Gini = gini_country)
  ),
  file.path(tables_dir, "sampling_bias_summary.xlsx")
)

# -------------------------------
# FIGURES
# -------------------------------

p_country <- country_totals %>%
  slice_max(N_total, n = 25) %>%
  ggplot(aes(x = reorder(Country, N_total), y = N_total)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top 25 countries by genome count",
    x = "Country",
    y = "Number of genomes"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  file.path(fig_dir, "top25_countries_counts.png"),
  p_country,
  width = 7,
  height = 8,
  dpi = 300
)

p_region <- region_totals %>%
  ggplot(aes(x = reorder(Region, N_total), y = N_total)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Genome sampling by region",
    x = "Region",
    y = "Genomes"
  ) +
  theme_minimal()

ggsave(
  file.path(fig_dir, "region_counts.png"),
  p_region,
  width = 6,
  height = 5,
  dpi = 300
)

p_year <- year_totals %>%
  ggplot(aes(x = Year, y = N_total)) +
  geom_line() +
  geom_point(size = 1) +
  labs(
    title = "Genome sampling intensity over time",
    x = "Year",
    y = "Genomes"
  ) +
  theme_minimal()

ggsave(
  file.path(fig_dir, "year_counts.png"),
  p_year,
  width = 7,
  height = 4,
  dpi = 300
)

# -------------------------------
# LORENZ CURVE
# -------------------------------

png(
  file.path(fig_dir, "lorenz_country_sampling.png"),
  width = 700,
  height = 600
)

Lc(country_totals$N_total,
   main = "Lorenz curve of genome sampling by country")

dev.off()

# -------------------------------
# CONSOLE SUMMARY
# -------------------------------

message("----- Sampling bias summary -----")
message("Countries represented: ", n_distinct(campy$Country))
message("Regions represented:   ", n_distinct(campy$Region))
message("Years represented:     ", n_distinct(campy$Year))

message("Top 5 countries:")
print(head(country_totals, 5))

message("Gini coefficient (country imbalance): ", round(gini_country, 3))

message("Sampling bias analysis complete.")


















# -------------------------------
# PUBLICATION COLORS + THEME
# -------------------------------

cb_pal <- c(
  blue   = "#0072B2",
  orange = "#E69F00",
  green  = "#009E73",
  red    = "#D55E00",
  purple = "#CC79A7",
  sky    = "#56B4E9",
  grey   = "#666666"
)

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      legend.position = "top",
      panel.grid.major = element_blank()
    )
}

save_pub <- function(plot, file, width = 180, height = 120) {
  ggsave(
    file,
    plot = plot,
    width = width,
    height = height,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
}


p_country <- country_totals %>%
  slice_max(N_total, n = 25) %>%
  mutate(Country = reorder(Country, N_total)) %>%
  ggplot(aes(x = Country, y = N_total)) +
  geom_col(fill = cb_pal$blue) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Number of genomes",
    title = "Sampling intensity by country (top 25)"
  ) +
  theme_pub()

save_pub(
  p_country,
  file.path(fig_dir, "Fig_sampling_top25_countries.png"),
  width = 180,
  height = 140
)


# ===============================
# 02_sampling_bias_quantification.R
# ===============================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(writexl)
  library(ineq)
})

# -------------------------------
# PATH SETUP
# -------------------------------

proj_root <- "~/PROJECT-2025/Campylobacter-Paper/gyA_analysis/REVISION-1/campy-rev"

data_dir   <- file.path(proj_root, "data")
fig_dir    <- file.path(proj_root, "figures", "sampling_bias")
tables_dir <- file.path(proj_root, "tables")
bias_dir   <- file.path(proj_root, "sampling_bias")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bias_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# LOAD DATA
# -------------------------------

campy <- read_excel(
  file.path(data_dir, "cleaned_campy_data_deduplicated.xlsx")
)

message("Total genomes analysed: ", nrow(campy))

# -------------------------------
# COUNTRY × YEAR MATRIX
# -------------------------------

country_year <- campy %>%
  count(Country, Year, name = "N_genomes") %>%
  arrange(desc(N_genomes)) %>% 
  filter(Year > 1990)

write.table(
  country_year,
  file.path(bias_dir, "country_year_counts.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# -------------------------------
# COUNTRY TOTALS + GINI
# -------------------------------

country_totals <- campy %>%
  count(Country, name = "N_total") %>%
  arrange(desc(N_total))

gini_country <- ineq(country_totals$N_total, type = "Gini")

# -------------------------------
# REGION TOTALS
# -------------------------------

region_totals <- campy %>%
  count(Region, name = "N_total") %>%
  arrange(desc(N_total))

# -------------------------------
# YEAR TOTALS
# -------------------------------

year_totals <- campy %>%
  count(Year, name = "N_total") %>%
  arrange(Year)

# -------------------------------
# SAVE TABLES
# -------------------------------

write_xlsx(
  list(
    country_totals = country_totals,
    region_totals  = region_totals,
    year_totals    = year_totals,
    country_year   = country_year,
    gini_country   = tibble(Gini = gini_country)
  ),
  file.path(tables_dir, "sampling_bias_summary.xlsx")
)

# -------------------------------
# PUBLICATION COLORS + THEME
# -------------------------------

cb_pal <- c(
  blue   = "#0072B2",
  orange = "#E69F00",
  green  = "#009E73",
  red    = "#D55E00",
  purple = "#CC79A7",
  sky    = "#56B4E9",
  grey   = "#666666"
)

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      legend.position = "top",
      panel.grid.major = element_blank()
    )
}

save_pub <- function(plot, file, width = 180, height = 120) {
  ggsave(
    file,
    plot = plot,
    width = width,
    height = height,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
}

# -------------------------------
# FIGURES
# -------------------------------

p_country <- country_totals %>%
  slice_max(N_total, n = 25) %>%
  mutate(Country = reorder(Country, N_total)) %>%
  ggplot(aes(x = Country, y = N_total)) +
  geom_col(fill = cb_pal["blue"]) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Number of genomes",
    title = "Sampling intensity by country (top 25)"
  ) +
  theme_pub()

save_pub(
  p_country,
  file.path(fig_dir, "Fig_sampling_top25_countries.png"),
  width = 180,
  height = 140
)

p_region <- region_totals %>%
  ggplot(aes(x = reorder(Region, N_total), y = N_total)) +
  geom_col(fill = cb_pal["blue"]) +
  coord_flip() +
  labs(
    title = "Genome sampling by region",
    x = "Region",
    y = "Genomes"
  ) +
  theme_pub()

save_pub(
  p_region,
  file.path(fig_dir, "Fig_region_counts.png"),
  width = 180,
  height = 120
)

p_year <- year_totals %>%
  ggplot(aes(x = Year, y = N_total)) +
  geom_line(color = cb_pal["blue"]) +
  geom_point(size = 1.5, color = cb_pal["blue"]) +
  labs(
    title = "Genome sampling intensity over time",
    x = "Year",
    y = "Genomes"
  ) +
  theme_pub()

save_pub(
  p_year,
  file.path(fig_dir, "Fig_year_counts.png"),
  width = 180,
  height = 120
)

# -------------------------------
# LORENZ CURVE
# -------------------------------

png(
  file.path(fig_dir, "lorenz_country_sampling.png"),
  width = 700,
  height = 600
)

Lc(country_totals$N_total,
   main = "Lorenz curve of genome sampling by country")

dev.off()

# -------------------------------
# CONSOLE SUMMARY
# -------------------------------

message("----- Sampling bias summary -----")
message("Countries represented: ", n_distinct(campy$Country))
message("Regions represented:   ", n_distinct(campy$Region))
message("Years represented:     ", n_distinct(campy$Year))

message("Top 5 countries:")
print(head(country_totals, 5))

message("Gini coefficient (country imbalance): ", round(gini_country, 3))

message("Sampling bias analysis complete.")














