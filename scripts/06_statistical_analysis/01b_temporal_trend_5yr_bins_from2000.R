# ============================================================
# 01b_temporal_trend_5yr_bins_from2000.R
# 5-year binned prevalence of gyrA mutations (>=2000)
# Produces: tables + publication-quality figure
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(writexl)
  library(stringr)
})




suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(stringr)
})

# -------------------------------
# LOAD DATA
# -------------------------------

proj_root <- "~/PROJECT-2025/Campylobacter-Paper/gyA_analysis/REVISION-1/campy-rev"
data_dir  <- file.path(proj_root, "data")
fig_dir   <- file.path(proj_root, "figures")

campy <- read_excel(file.path(data_dir, "cleaned_campy_data_deduplicated.xlsx"))

campy <- campy %>%
  mutate(
    gyrA_mut = as.integer(`GyrA-Mutation`)
  )

# ============================================================
# PANEL A — CONTINENT PREVALENCE
# ============================================================

continent_prev <- campy %>%
  group_by(Region) %>%
  summarise(
    Total = n(),
    Mutant = sum(gyrA_mut == 1, na.rm = TRUE),
    Prevalence = Mutant / Total,
    .groups = "drop"
  ) %>%
  arrange(desc(Prevalence)) %>%
  mutate(
    Label = paste0(round(Prevalence*100,1), "%\n(", Mutant, "/", Total, ")")
  )

pA <- ggplot(continent_prev,
             aes(x = reorder(Region, Prevalence),
                 y = Prevalence)) +
  geom_col(fill = "#66c2a5", width = 0.7) +
  geom_text(aes(label = Label),
            hjust = -0.1,
            size = 4,
            fontface = "bold") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Genomic prevalence of gyrA QRDR substitutions (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

pA


# ============================================================
# PANEL B — TEMPORAL TREND (5-YEAR BINS SINCE 2000)
# ============================================================

temporal_data <- campy %>%
  filter(Year >= 2000) %>%
  mutate(
    Interval = paste0(floor(Year/5)*5, "-", floor(Year/5)*5 + 4)
  ) %>%
  group_by(Interval) %>%
  summarise(
    Total = n(),
    Mutant = sum(gyrA_mut == 1, na.rm = TRUE),
    Prevalence = Mutant / Total,
    .groups = "drop"
  ) %>%
  filter(Total >= 20) %>%
  mutate(
    IntervalMid = as.numeric(str_extract(Interval, "^\\d{4}")) + 2.5,
    Label = paste0(round(Prevalence*100,1), "%")
  )

model <- lm(Prevalence ~ IntervalMid, data = temporal_data)

slope <- round(coef(model)[2]*100, 2)
r2    <- round(summary(model)$r.squared, 3)
pval  <- signif(summary(model)$coefficients[2,4], 3)

pB <- ggplot(temporal_data,
             aes(x = IntervalMid,
                 y = Prevalence)) +
  geom_line(color = "#B22222", linewidth = 1.5) +
  geom_point(size = 3, color = "#B22222") +
  geom_smooth(method = "lm", linetype = "dashed",
              color = "black", se = TRUE) +
  scale_x_continuous(
    breaks = temporal_data$IntervalMid,
    labels = temporal_data$Interval
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  annotate("text",
           x = min(temporal_data$IntervalMid),
           y = max(temporal_data$Prevalence),
           hjust = 0,
           label = paste0("Slope = ", slope,
                          "% per 5-year interval\nR² = ", r2,
                          "\np = ", pval),
           size = 4) +
  labs(
    x = NULL,
    y = "Genomic prevalence (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

# ============================================================
# PANEL C — SPECIES PREVALENCE
# ============================================================

species_prev <- campy %>%
  group_by(Species) %>%
  summarise(
    Total = n(),
    Mutant = sum(gyrA_mut == 1, na.rm = TRUE),
    Prevalence = Mutant / Total,
    .groups = "drop"
  ) %>%
  arrange(desc(Prevalence)) %>%
  mutate(
    Label = paste0(round(Prevalence*100,1), "%\n(", Mutant, "/", Total, ")")
  )

pC <- ggplot(species_prev,
             aes(x = reorder(Species, Prevalence),
                 y = Prevalence)) +
  geom_col(fill = "#00A087FF", width = 0.7) +
  geom_text(aes(label = Label),
            hjust = -0.1,
            size = 4,
            fontface = "bold") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Genomic prevalence (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

# ============================================================
# COMBINE PANELS
# ============================================================

final_plot <- (pA | pB) / pC +
  plot_annotation(
    tag_levels = "A",
    title = "Global genomic prevalence and temporal dynamics of fluoroquinolone-associated gyrA QRDR substitutions"
  )

ggsave(
  file.path(fig_dir, "Figure1_Main_NatureStyle.png"),
  final_plot,
  width = 14,
  height = 10,
  dpi = 600
)

print(final_plot)




















# ============================================================
# PANEL B — TEMPORAL TREND WITH MULTIPLE TREND LINES
# ============================================================
temporal_data <- campy %>%
  filter(Year >= 2000) %>%
  mutate(
    Interval = paste0(floor(Year/5)*5, "-", floor(Year/5)*5 + 4)
  ) %>%
  group_by(Interval) %>%
  summarise(
    Total = n(),
    Mutant = sum(gyrA_mut == 1, na.rm = TRUE),
    Prevalence = Mutant / Total,
    .groups = "drop"
  ) %>%
  filter(Total >= 20) %>%
  mutate(
    IntervalMid = as.numeric(str_extract(Interval, "^\\d{4}")) + 2.5,
    Label = paste0(round(Prevalence*100, 1), "%")
  )

# Linear model statistics
model <- lm(Prevalence ~ IntervalMid, data = temporal_data)
slope <- round(coef(model)[2] * 100, 2)
r2    <- round(summary(model)$r.squared, 3)
pval  <- signif(summary(model)$coefficients[2, 4], 3)

# LOESS model for comparison
loess_model <- loess(Prevalence ~ IntervalMid, data = temporal_data, span = 0.75)
loess_r2 <- round(cor(temporal_data$Prevalence, 
                      predict(loess_model))^2, 3)





# -------------------------------
# 4) PUBLICATION-QUALITY FIGURE (WITH LEGEND + n/N LABELS)
# -------------------------------
# ============================================================
# PANEL B — TEMPORAL TREND WITH LEGEND
# ============================================================
library(ggrepel)
library(dplyr)
library(ggplot2)
library(scales)

# Build labels with numerator/denominator
temporal_data <- temporal_data %>%
  mutate(
    Label = paste0(
      round(Prevalence * 100, 1), "%\n(",
      Mutant, "/", Total, ")"
    )
  )

pB <- ggplot(temporal_data, aes(x = IntervalMid, y = Prevalence)) +
  
  # Observed points
  geom_point(aes(color = "Observed data"),
             size = 4.2) +
  
  # Labels
  geom_text_repel(
    aes(label = Label),
    size = 4.0,                # bigger text
    fontface = "bold",         # bold labels
    lineheight = 0.9,
    color = "black",
    nudge_y = 0.02,
    segment.color = "gray60",
    segment.size = 0.35,
    min.segment.length = 0,
    box.padding = 0.6,
    point.padding = 0.35,
    max.overlaps = Inf
  ) +
  
  # LOESS
  geom_smooth(aes(color = "LOESS smoothing", linetype = "LOESS smoothing"),
              method = "loess",
              span = 0.75,
              se = FALSE,
              linewidth = 1.4) +
  
  # Linear
  geom_smooth(aes(color = "Linear regression",
                  linetype = "Linear regression",
                  fill = "Linear regression"),
              method = "lm",
              linewidth = 1.4,
              alpha = 0.18) +
  
  scale_color_manual(
    name = NULL,
    values = c(
      "Observed data" = "#fc8d62",
      "Linear regression" = "#fc8d62",
      "LOESS smoothing" = "#66c2a5"
    ),
    breaks = c("Observed data", "LOESS smoothing", "Linear regression")
  ) +
  
  scale_linetype_manual(
    name = NULL,
    values = c(
      "Linear regression" = "dashed",
      "LOESS smoothing" = "solid"
    ),
    breaks = c("LOESS smoothing", "Linear regression")
  ) +
  
  scale_fill_manual(values = c("Linear regression" = "#cdcdcd"), guide = "none") +
  
  # Model annotation. Make it bold and bigger.
  annotate("text",
           x = min(temporal_data$IntervalMid),
           y = max(temporal_data$Prevalence) * 1.10,
           hjust = 0,
           vjust = 1,
           label = paste0("Slope = ", slope,
                          "% per 5-year interval\nR² = ", r2,
                          "\np = ", pval),
           size = 5.2,
           fontface = "bold",
           color = "black") +
  
  scale_x_continuous(breaks = temporal_data$IntervalMid,
                     labels = temporal_data$Interval) +
  
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.05, 0.25))) +
  
  labs(
    x = NULL,
    y = "gyrA mutant prevalence (%)"
  ) +
  
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(color = "black", face = "bold", size = 15),
    axis.text.y = element_text(color = "black", face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 18),
    legend.position = "bottom",
    legend.text = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 20, 10, 10)
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        shape = c(16, NA, NA),
        linetype = c("blank", "solid", "dashed"),
        linewidth = c(0, 1.4, 1.4)
      )
    ),
    linetype = "none"
  )

pB
