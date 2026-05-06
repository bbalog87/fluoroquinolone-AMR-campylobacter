###############################################
# GLOBAL + CC-AWARE SAMPLING WITH JEJUNI CAP 60
# + WT-ONLY EXPANSION TO ~110 GENOMES
# + GEOGRAPHIC + HOST STRATIFICATION
###############################################

library(tidyverse)
library(janitor)
library(readxl)

###############################################
# 1. Load raw datasets
###############################################

dat_raw <- read_excel("masterFile_campylobacter.xlsx") %>% clean_names()
cc_raw  <- read_excel("CC_resitance.xlsx") %>% clean_names()

###############################################
# 2. Deduplicate by isolate
###############################################

dat    <- dat_raw %>% arrange(isolate) %>% distinct(isolate, .keep_all = TRUE)
cc_tab <- cc_raw  %>% arrange(isolate) %>% distinct(isolate, .keep_all = TRUE)

###############################################
# 3. Build species_group
###############################################

cc_tab <- cc_tab %>%
  mutate(
    species_clean = str_to_lower(species),
    species_group = case_when(
      str_detect(species_clean, "campylobacter jejuni")        ~ "campylobacter_jejuni",
      str_detect(species_clean, "campylobacter coli")          ~ "campylobacter_coli",
      str_detect(species_clean, "campylobacter lari")          ~ "campylobacter_lari",
      str_detect(species_clean, "campylobacter upsaliensis")   ~ "campylobacter_upsaliensis",
      str_detect(species_clean, "campylobacter insulaenigrae") ~ "campylobacter_insulaenigrae",
      str_detect(species_clean, "campylobacter sp\\. jh")      ~ "campylobacter_sp_JH",
      str_detect(species_clean, "campylobacter sp\\. cfsan")   ~ "campylobacter_sp_CFSAN",
      str_detect(species_clean, "campylobacter sp\\.")         ~ "campylobacter_sp_other",
      TRUE                                                     ~ "campylobacter_unknown"
    )
  )

###############################################
# 4. Determine top 8 species_group
###############################################

top8_species <- cc_tab %>%
  count(species_group) %>%
  arrange(desc(n)) %>%
  slice_head(n = 8) %>%
  pull(species_group)

###############################################
# 5. Merge species_group + CC into master dataset
###############################################

# Remove species_group from dat only if it exists
if ("species_group" %in% names(dat)) {
  dat <- dat %>% select(-species_group)
}

dat <- dat %>%
  left_join(
    cc_tab %>% select(isolate, species_group, clonal_complex),
    by = "isolate"
  ) %>%
  filter(species_group %in% top8_species)

###############################################
# 6. Assign CC only for jejuni/coli
###############################################

dat <- dat %>%
  mutate(
    cc = case_when(
      species_group %in% c("campylobacter_jejuni", "campylobacter_coli") &
        !is.na(clonal_complex) ~ clonal_complex,
      species_group %in% c("campylobacter_jejuni", "campylobacter_coli") &
        is.na(clonal_complex)  ~ "unassigned",
      TRUE ~ NA_character_
    )
  )

###############################################
# 7. QRDR classification
###############################################

dat <- dat %>%
  mutate(
    qrdr_class = case_when(
      is.na(quinolone) | quinolone == "" ~ "WT",
      str_detect(quinolone, "T86I|T86V|T86A|T86K|D90N|D90Y|P104S") ~ "QRDR_mutant",
      TRUE ~ "WT"
    )
  )

###############################################
# 8. Define strata
###############################################

dat <- dat %>% mutate(stratum = species_group)

###############################################
# 9. Sampling functions (MUTANTS FIXED)
###############################################

# Non-jejuni species: WT expanded, mutants fixed at 4
sample_species <- function(df, n_wt = 10, n_mut = 4) {
  df_unique     <- df %>% distinct(isolate, species_group, .keep_all = TRUE)
  wt_available  <- df_unique %>% filter(qrdr_class == "WT")
  mut_available <- df_unique %>% filter(qrdr_class == "QRDR_mutant")
  
  wt  <- wt_available  %>% slice_sample(n = min(n_wt,  nrow(wt_available)))
  mut <- mut_available %>% slice_sample(n = min(n_mut, nrow(mut_available)))
  
  bind_rows(wt, mut)
}

# Jejuni/coli CC-aware: WT expanded, mutants fixed at 3
sample_species_cc <- function(df, n_wt = 6, n_mut = 3) {
  df_unique <- df %>% distinct(isolate, species_group, .keep_all = TRUE)
  
  df_unique %>%
    mutate(cc = if_else(is.na(cc), "unassigned", cc)) %>%
    group_by(cc) %>%
    group_modify(~ {
      wt_available  <- .x %>% filter(qrdr_class == "WT")
      mut_available <- .x %>% filter(qrdr_class == "QRDR_mutant")
      
      wt  <- wt_available  %>% slice_sample(n = min(n_wt,  nrow(wt_available)))
      mut <- mut_available %>% slice_sample(n = min(n_mut, nrow(mut_available)))
      
      bind_rows(wt, mut)
    }) %>%
    ungroup()
}

###############################################
# 10. Apply sampling across strata
###############################################

subset_list <- dat %>%
  group_by(stratum) %>%
  group_split() %>%
  map(~ {
    strat <- unique(.x$stratum)
    if (strat %in% c("campylobacter_jejuni", "campylobacter_coli")) {
      sample_species_cc(.x)
    } else {
      sample_species(.x)
    }
  })

subset_global <- bind_rows(subset_list)

###############################################
# 11. WT-ONLY EXPANSION TO ~110 GENOMES
###############################################

target_total <- 110
current_total <- nrow(subset_global)
needed <- target_total - current_total

if (needed > 0) {
  
  wt_pool <- dat %>%
    filter(qrdr_class == "WT") %>%
    anti_join(subset_global, by = "isolate")
  
  # STEP A: geographic stratification
  geo_wt <- wt_pool %>%
    filter(!is.na(region)) %>%
    group_by(region) %>%
    slice_sample(n = 1, replace = FALSE) %>%
    ungroup()
  
  # STEP B: host stratification
  host_wt <- wt_pool %>%
    filter(!is.na(host_group)) %>%
    anti_join(geo_wt, by = "isolate") %>%
    group_by(host_group) %>%
    slice_sample(n = 1, replace = FALSE) %>%
    ungroup()
  
  stratified_wt <- bind_rows(geo_wt, host_wt) %>%
    distinct(isolate, .keep_all = TRUE)
  
  remaining_needed <- needed - nrow(stratified_wt)
  
  filler_wt <- wt_pool %>%
    anti_join(stratified_wt, by = "isolate") %>%
    slice_sample(n = min(remaining_needed, nrow(.)))
  
  extra_wt <- bind_rows(stratified_wt, filler_wt) %>%
    distinct(isolate, .keep_all = TRUE)
  
  subset_global <- bind_rows(subset_global, extra_wt) %>%
    distinct(isolate, .keep_all = TRUE)
}

###############################################
# 12. Apply JEJUNI CAP (60 genomes)
###############################################

jejuni_cap <- 60

subset_global <- subset_global %>%
  group_by(species_group) %>%
  group_modify(~ {
    if (.y$species_group == "campylobacter_jejuni") {
      .x %>% slice_sample(n = min(jejuni_cap, nrow(.x)))
    } else {
      .x
    }
  }) %>%
  ungroup()

###############################################
# 13. Save global subset
###############################################

write_csv(subset_global, "subset_global_balanced_core_phylogeny_110.csv")

###############################################
# 14. Jejuni/coli zoom tree (no cap)
###############################################

subset_jejuni_coli <- subset_global %>%
  filter(species_group %in% c("campylobacter_jejuni", "campylobacter_coli"))

write_csv(subset_jejuni_coli, "subset_jejuni_coli_zoom_tree_110.csv")

###############################################
# 15. Sanity checks
###############################################

print(subset_global %>% count(species_group))
print(subset_global %>% count(qrdr_class))
print(subset_jejuni_coli %>% count(cc, qrdr_class))
print(nrow(subset_global))





















###############################################
# 16. Comprehensive dataset statistics
###############################################

library(vegan)

cat("\n=============================\n")
cat("SPECIES DIVERSITY\n")
cat("=============================\n")

species_stats <- subset_global %>% count(species_group)
print(species_stats)

cat("\nSpecies richness:", nrow(species_stats), "\n")
cat("Shannon diversity:", diversity(species_stats$n, index = "shannon"), "\n")
cat("Simpson diversity:", diversity(species_stats$n, index = "simpson"), "\n")


cat("\n=============================\n")
cat("QRDR / WT BALANCE\n")
cat("=============================\n")

qrdr_stats <- subset_global %>% count(qrdr_class)
print(qrdr_stats)

cat("\nWT proportion:", round(qrdr_stats$n[qrdr_stats$qrdr_class=="WT"] /
                                sum(qrdr_stats$n), 3), "\n")


cat("\n=============================\n")
cat("JEJUNI + COLI CC DIVERSITY\n")
cat("=============================\n")

cc_stats <- subset_jejuni_coli %>% count(cc)
print(cc_stats)

cat("\nNumber of CCs:", nrow(cc_stats), "\n")
cat("Shannon CC diversity:", diversity(cc_stats$n, index = "shannon"), "\n")


cat("\n=============================\n")
cat("GEOGRAPHIC DIVERSITY\n")
cat("=============================\n")

geo_stats <- subset_global %>% count(region)
print(geo_stats)

cat("\nNumber of regions:", nrow(geo_stats), "\n")


cat("\n=============================\n")
cat("HOST DIVERSITY\n")
cat("=============================\n")

host_stats <- subset_global %>% count(host_group)
print(host_stats)

cat("\nNumber of host groups:", nrow(host_stats), "\n")


cat("\n=============================\n")
cat("GENOME COMPLETENESS\n")
cat("=============================\n")

comp_stats <- subset_global %>% summarise(
  mean_completeness = mean(check_m_completeness, na.rm = TRUE),
  median_completeness = median(check_m_completeness, na.rm = TRUE),
  min_completeness = min(check_m_completeness, na.rm = TRUE),
  max_completeness = max(check_m_completeness, na.rm = TRUE)
)

print(comp_stats)


cat("\n=============================\n")
cat("SUMMARY\n")
cat("=============================\n")

cat("Total genomes:", nrow(subset_global), "\n")
cat("Total WT:", qrdr_stats$n[qrdr_stats$qrdr_class=="WT"], "\n")
cat("Total mutants:", qrdr_stats$n[qrdr_stats$qrdr_class=="QRDR_mutant"], "\n")
cat("Species richness:", nrow(species_stats), "\n")
cat("CC richness (jejuni+coli):", nrow(cc_stats), "\n")
cat("Geographic regions:", nrow(geo_stats), "\n")
cat("Host groups:", nrow(host_stats), "\n")

