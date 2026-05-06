############################################################
# Robust pipeline: select ~120 C. jejuni + C. coli genomes
# stratified by clonal complex × host group × region
#
# FIX LOG (v2):
#  - Hard-coded correct column names confirmed from actual data:
#      masterFile  -> Species, HostGroup, Region, Country
#      CC_resitance -> ClonalComplex, Host.Group, Geographic.Group, Species
#  - Removed fragile pick_col() helper; columns are selected explicitly
#  - Fixed species normalization (source has proper-case "Campylobacter jejuni")
#  - Fixed slice_sample() crash: replaced n=first(n_per_group) with
#    group_modify() + explicit min(nwant, nrow(.x))
#  - Fixed step-8 crash: target_per_sp cast to integer before slice_sample()
#  - Added data-quality checks & informative messages throughout
############################################################

library(tidyverse)
library(janitor)
library(readxl)
library(purrr)

# ── 1. Load & deduplicate ────────────────────────────────────────────────────

dat_raw <- read_excel("masterFile_campylobacter.xlsx") %>% clean_names()
cc_raw  <- read_excel("CC_resitance.xlsx")             %>% clean_names()

# clean_names() converts:
#   masterFile  : Species -> species, HostGroup -> host_group,
#                 Region -> region, Country -> country
#   CC_resitance: ClonalComplex -> clonal_complex,
#                 Host.Group -> host_group, Geographic.Group -> geographic_group,
#                 Species -> species

# Quick column audit — uncomment if you need to re-check:
# message("masterFile cols: ", paste(names(dat_raw), collapse = ", "))
# message("CC cols:         ", paste(names(cc_raw),  collapse = ", "))

dat    <- dat_raw %>% arrange(isolate) %>% distinct(isolate, .keep_all = TRUE)
cc_tab <- cc_raw  %>% arrange(isolate) %>% distinct(isolate, .keep_all = TRUE)

message("masterFile isolates (deduped): ", nrow(dat))
message("CC_resistance isolates (deduped): ", nrow(cc_tab))

# ── 2. Build clean metadata ──────────────────────────────────────────────────
#
# Priority for host  : host_group (masterFile) -> host_group (CC file)
# Priority for region: geographic_group (CC file) -> region (masterFile)
# Clonal complex comes ONLY from CC file

meta_clean <- dat %>%
  select(
    isolate,
    species,           # "Campylobacter jejuni" / "Campylobacter coli" etc.
    host_group,        # e.g. "Human", "Avian", "Nonhuman Mammal"
    region,            # broad geographic region from masterFile
    country            # fallback geography
  ) %>%
  left_join(
    cc_tab %>% select(isolate, clonal_complex, geographic_group),
    by = "isolate"
  ) %>%
  mutate(
    # Normalise species to snake_case canonical names
    species = case_when(
      str_detect(tolower(species), "jejuni") ~ "campylobacter_jejuni",
      str_detect(tolower(species), "coli")   ~ "campylobacter_coli",
      TRUE ~ tolower(str_replace_all(str_trim(species), "\\s+", "_"))
    ),
    
    # Prefer geographic_group from CC file; fall back to region from masterFile
    region = case_when(
      !is.na(geographic_group) & geographic_group != "" ~ geographic_group,
      !is.na(region)           & region != ""           ~ region,
      !is.na(country)          & country != ""          ~ country,
      TRUE ~ "Unknown"
    ),
    
    # Fill missing host_group
    host_group = if_else(
      is.na(host_group) | host_group == "" | host_group == "NA",
      "Unknown",
      as.character(host_group)
    ),
    
    clonal_complex = if_else(
      is.na(clonal_complex) | clonal_complex == "" | clonal_complex == "NA",
      "Unassigned",
      as.character(clonal_complex)
    )
  ) %>%
  select(-geographic_group, -country) %>%   # housekeeping
  distinct(isolate, .keep_all = TRUE)

# ── 3. Filter to C. jejuni + C. coli ────────────────────────────────────────

meta_jj_cc <- meta_clean %>%
  filter(species %in% c("campylobacter_jejuni", "campylobacter_coli"))

message("C. jejuni rows: ", sum(meta_jj_cc$species == "campylobacter_jejuni"))
message("C. coli rows:   ", sum(meta_jj_cc$species == "campylobacter_coli"))

if (nrow(meta_jj_cc) == 0) {
  stop("No C. jejuni or C. coli rows after filtering. ",
       "Check species column. Distinct values found:\n",
       paste(unique(meta_clean$species), collapse = "\n"))
}

# ── 4. Identify major CCs (≥5 genomes) ──────────────────────────────────────

cc_counts <- meta_jj_cc %>%
  filter(clonal_complex != "Unassigned") %>%
  count(species, clonal_complex, name = "n", sort = TRUE)

major_cc <- cc_counts %>%
  filter(n >= 5) %>%
  pull(clonal_complex) %>%
  unique()

message("Major CCs (≥5 genomes): ", length(major_cc),
        " -> ", paste(head(major_cc, 10), collapse = ", "),
        if (length(major_cc) > 10) "..." else "")

# ── 5. Stratified sampling: species × CC × host_group × region ──────────────
#
#   Major CC strata -> up to 4 isolates each
#   Minor CC strata -> up to 2 isolates each
#
# FIX: use group_modify() so n can reference per-group data;
#      min(nwant, nrow(.x)) prevents errors when a stratum is smaller than nwant

set.seed(123)

meta_jj_cc_annot <- meta_jj_cc %>%
  mutate(
    cc_class   = if_else(clonal_complex %in% major_cc, "major", "minor"),
    n_per_group = if_else(cc_class == "major", 4L, 2L)
  )

meta_selected <- meta_jj_cc_annot %>%
  group_by(species, clonal_complex, host_group, region) %>%
  group_modify(function(.x, .k) {
    nwant <- as.integer(.x$n_per_group[1])
    take  <- min(nwant, nrow(.x))
    slice_sample(.x, n = take)
  }) %>%
  ungroup()

message("After stratified sampling: ", nrow(meta_selected), " isolates")

# ── 6. Balance to ~120 total (≈60 jejuni + 60 coli) ─────────────────────────
#
# FIX: cast target_per_sp to integer — slice_sample() requires integer n

target_total  <- 120L
target_per_sp <- as.integer(target_total / 2L)

meta_balanced <- meta_selected %>%
  group_by(species) %>%
  group_modify(function(.x, .k) {
    slice_sample(.x, n = min(nrow(.x), target_per_sp))
  }) %>%
  ungroup()

n_jej  <- sum(meta_balanced$species == "campylobacter_jejuni")
n_col  <- sum(meta_balanced$species == "campylobacter_coli")
message("Final: ", nrow(meta_balanced), " isolates  (",
        n_jej, " jejuni + ", n_col, " coli)")

# ── 7. Export ────────────────────────────────────────────────────────────────

write_lines(meta_balanced$isolate,
            "jejuni_coli_selected_isolates.txt")

write_lines(paste0("genomes/", meta_balanced$isolate, ".fna"),
            "jejuni_coli_selected_fastas.txt")

write_csv(meta_balanced %>% select(-cc_class, -n_per_group),
          "jejuni_coli_selected_metadata.csv")

message("\nExported files:")
message("  jejuni_coli_selected_isolates.txt")
message("  jejuni_coli_selected_fastas.txt")
message("  jejuni_coli_selected_metadata.csv")

# ── 8. Summary tables ────────────────────────────────────────────────────────

cat("\n── CC distribution ──────────────────────────────\n")
meta_balanced %>%
  count(species, clonal_complex, name = "n") %>%
  arrange(species, desc(n)) %>%
  print(n = 200)

cat("\n── Host group distribution ──────────────────────\n")
meta_balanced %>%
  count(species, host_group, name = "n") %>%
  arrange(species, desc(n)) %>%
  print(n = 200)

cat("\n── Region distribution ──────────────────────────\n")
meta_balanced %>%
  count(species, region, name = "n") %>%
  arrange(species, desc(n)) %>%
  print(n = 200)

















# Load selected isolate IDs
selected_ids <- readLines("jejuni_coli_selected_isolates.txt")

# Directory containing all genomes
all_genomes_dir <-"/home/nguinkal/PROJECT-2025/Campylobacter-Paper/gyA_analysis/all_genomes"   # <-- adjust if needed

# Output directory
outdir <- "selected_genomes"
dir.create(outdir, showWarnings = FALSE)

# List all .fna files
all_fna <- list.files(all_genomes_dir, pattern = "\\.fna$", full.names = TRUE)

# Function to check if filename contains the isolate ID
matches_isolate <- function(file, isolate) {
  grepl(isolate, basename(file), fixed = TRUE)
}

# Copy genomes
for (iso in selected_ids) {
  hits <- all_fna[sapply(all_fna, matches_isolate, isolate = iso)]
  
  if (length(hits) == 0) {
    message("WARNING: No genome found for isolate: ", iso)
  } else if (length(hits) > 1) {
    message("WARNING: Multiple matches for isolate: ", iso, " → copying first")
    file.copy(hits[1], outdir)
  } else {
    file.copy(hits, outdir)
  }
}

message("Done. Copied genomes are in: ", outdir)
