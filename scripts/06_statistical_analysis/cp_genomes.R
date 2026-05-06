library(tidyverse)
library(fs)

subset <- read_csv("subset_global_balanced_core_phylogeny_110.csv")

# Adjust this if your column name is different:
selected_ids <- subset$assembly_accession

src_dir <- "/home/nguinkal/PROJECT-2025/Campylobacter-Paper/gyA_analysis/all_genomes"
dst_dir <- "/home/nguinkal/PROJECT-2025/Campylobacter-Paper/gyA_analysis/REVISION-1/phylogeny_110/selected_genomes"

dir_create(dst_dir)

for (id in selected_ids) {
  # Try .fna, .fa, .fasta
  candidates <- c(
    file.path(src_dir, paste0(id, ".fna")),
    file.path(src_dir, paste0(id, ".fa")),
    file.path(src_dir, paste0(id, ".fasta"))
  )
  
  existing <- candidates[file_exists(candidates)]
  
  if (length(existing) == 1) {
    file_copy(existing, file.path(dst_dir, basename(existing)), overwrite = TRUE)
  } else {
    message("Missing genome for: ", id)
  }
}
 


#### Phylo genetric trees


library(ape)
library(ggtree)
library(phangorn)


tree <- read.tree("iqtree/campy.treefile")
rooted <- midpoint(tree)
write.tree(rooted, file = "iqtree/campy_rooted.treefile")

ggtree(rooted) + geom_tiplab(size = 2)




library(readxl)
library(tidyverse)
library(ape)

raw <- read_excel("subset_global_balanced_core_phylogeny_110.xlsx")

meta <- raw %>%
  transmute(
    taxon = assembly_accession,   # or genome_name, depending on tree labels
    species = stratum,
    clonal_complex = clonal_complex,
    host_group = host_group,
    country = country,
    region = region,
    year = year,
    qrdr_class = qrdr_class,
    quinolone = quinolone,
    cmeB_status = florfenicol_quinolone,  # or create a clean column
    st = st,
    isolation_source = isolation_source
  )

write.csv(meta, "metadata_tree.txt", row.names = F)
write.csv(meta2, "metadata_tree_clade.txt", row.names = F)




library(ape)
library(tidyverse)

# Load tree
tree <- read.tree("iqtree/campy_rooted.treefile")

# Patristic distance matrix
D <- cophenetic(tree)

# Hierarchical clustering
hc <- hclust(as.dist(D), method = "average")

# Choose number of clades
groups <- cutree(hc, k = 12)

# Convert to tibble
clade_df <- tibble(
  taxon = names(groups),
  clade = paste0("Clade_", groups)
)


meta2 <- meta %>% left_join(clade_df, by = "taxon")







###############################################################
# GLOBAL CAMPYLOBACTER TREE — HOST-SHARING + CC + REGION
###############################################################
###############################################################
# CAMPYLOBACTER CC-CENTRIC ULTRAMETRIC CIRCULAR PHYLOGENY
# WITH REGION + ISOLATION SOURCE RINGS
# Julien — fully documented, end-to-end pipeline
###############################################################

library(ape)
library(tidyverse)
library(ggtree)
library(ggforce)      # for geom_arcbar (continuous arcs)
library(ggtreeExtra)  # for geom_fruit (annotation rings)
library(ggnewscale)

###############################################################
# 1. Load tree + metadata
###############################################################

tree <- read.tree("iqtree/campy_rooted.treefile")
#meta <- read_tsv("metadata.tsv")

###############################################################
# 2. Convert tree to ultrametric and save it
###############################################################

tree_ultra <- chronos(tree, lambda = 1)
write.tree(tree_ultra, file = "campy_ultrametric.tree")

###############################################################
# 3. Harmonize metadata (species, CC, region, source)
###############################################################

meta2 <- meta %>%
  mutate(
    clonal_complex = str_trim(as.character(clonal_complex)),
    
    species = case_when(
      species == "campylobacter_jejuni"        ~ "C. jejuni",
      species == "campylobacter_coli"          ~ "C. coli",
      species == "campylobacter_lari"          ~ "C. lari",
      species == "campylobacter_insulaenigrae" ~ "C. insulaenigrae",
      TRUE                                     ~ "Campylobacter sp."
    ),
    
    region = if_else(is.na(region), "Unknown", region),
    
    # Isolation source (host_group or more detailed variable)
    source = case_when(
      host_group == "Human" ~ "Human",
      host_group == "Avian" ~ "Avian",
      host_group == "Nonhuman Mammal" ~ "Mammal",
      is.na(host_group) ~ "Unknown",
      TRUE ~ host_group
    )
  )

###############################################################
# 4. Define CC ecological categories
###############################################################

cc_ecology <- meta2 %>%
  group_by(clonal_complex) %>%
  summarise(
    n_human  = sum(source == "Human", na.rm = TRUE),
    n_animal = sum(source != "Human" & source != "Unknown"),
    n_total  = n()
  ) %>%
  mutate(
    cc_category = case_when(
      n_human > 0 & n_animal == 0 ~ "Human-Only",
      n_animal > 0 & n_human == 0 ~ "Animal-Only",
      n_animal > 0 & n_human > 0  ~ "Shared",
      TRUE                        ~ "Orphan"
    )
  )

meta2 <- meta2 %>% left_join(cc_ecology, by = "clonal_complex")

###############################################################
# 5. Build circular ultrametric tree and extract tip positions
###############################################################

#p <- ggtree(tree_ultra, layout = "circular")
p <- ggtree(tree_ultra, layout = "circular") %<+% meta2


tip_df <- p$data %>%
  filter(isTip) %>%
  select(label, angle, y) %>%
  left_join(meta2 %>% select(taxon, clonal_complex, cc_category, region, source),
            by = c("label" = "taxon"))

###############################################################
# 6. Compute CC arc boundaries
###############################################################

cc_arcs <- tip_df %>%
  arrange(angle) %>%
  group_by(clonal_complex, cc_category) %>%
  summarise(
    start = min(angle),
    end   = max(angle),
    mid   = (min(angle) + max(angle)) / 2,
    .groups = "drop"
  )

###############################################################
# 7. Define colour palettes
###############################################################

cc_cat_cols <- c(
  "Human-Only"  = "#C44E52",
  "Animal-Only" = "#55A868",
  "Shared"      = "#4C72B0",
  "Orphan"      = "#DDDDDD"
)

region_cols <- c(
  "Europe"        = "#4C72B0",
  "Asia"          = "#55A868",
  "North America" = "#C44E52",
  "South America" = "#8172B2",
  "Oceania"       = "#64B5CD",
  "Unknown"       = "#DDDDDD"
)

source_cols <- c(
  "Human"   = "#C44E52",
  "Avian"   = "#55A868",
  "Mammal"  = "#4C72B0",
  "Unknown" = "#DDDDDD"
)

###############################################################
# 8. Add REGION and SOURCE rings inside the tree
###############################################################

# Ring 1: region
p_r1 <- p +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(x = 1, y = label, fill = region),
    width = 0.05, offset = 0.02
  ) +
  scale_fill_manual(values = region_cols)

# Ring 2: isolation source
p_r2 <- p_r1 +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(x = 2, y = label, fill = source),
    width = 0.05, offset = 0.02
  ) +
  scale_fill_manual(values = source_cols)

###############################################################
# 9. Add CC arcs around the tree
###############################################################

p_arcs <- p_r2 +
  new_scale_fill() +
  geom_arcbar(
    data = cc_arcs,
    aes(
      x0 = 0, y0 = 0,
      r0 = max(p$data$y) + 0.5,
      r  = max(p$data$y) + 1.2,
      start = start * pi/180,
      end   = end   * pi/180,
      fill  = cc_category
    ),
    color = NA,
    alpha = 0.9
  ) +
  scale_fill_manual(values = cc_cat_cols)

###############################################################
# 10. Add CC labels outside the circle
###############################################################

p_labels <- p_arcs +
  geom_text(
    data = cc_arcs,
    aes(
      x = max(p$data$y) + 1.4,
      y = mid,
      label = clonal_complex,
      angle = mid - 90
    ),
    size = 2.5,
    hjust = 0
  )

###############################################################
# 11. Final styling
###############################################################

p_final <- p_labels +
  theme(
    legend.position = "right",
    text = element_text(family = "Helvetica", size = 8),
    plot.background = element_rect(fill = "white", colour = NA)
  )

###############################################################
# 12. Save final figure
###############################################################

ggsave("campy_ultrametric_CC_region_source.pdf", p_final,
       width = 12, height = 12, units = "in", useDingbats = FALSE)

ggsave("campy_ultrametric_CC_region_source.png", p_final,
       width = 12, height = 12, units = "in", dpi = 600)
