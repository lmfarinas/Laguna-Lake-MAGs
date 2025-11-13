#Load packages:
library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)
library(xlsx)
library(ggplot2)
library(viridis)
library(UpSetR)

#### Summarize bin quality stats
# List all bin stat files
mags_stats_files <- list.files(path = "MAGs_stats", pattern = "\\.tabular$", full.names = TRUE)

# Import all .tsv files and combine into one data frame
mags_stats <- lapply(mags_stats_files, function(file) {
  sample_name <- str_split(basename(file), "_", simplify = TRUE)[1] # Extract sample name (everything before first "_")
  df <- read.delim(file, stringsAsFactors = FALSE) # Read file
  colnames(df)[1] <- "bin_number" # Rename 1st column to bin_number
  
  # Add sample_name as first column
  df <- df %>%
    mutate(sample_name = sample_name,
           bin_name = paste0(sample_name, ".", bin_number)) %>%
    select(sample_name, bin_number, bin_name, everything())
  
  return(df)
}) %>% bind_rows()

# Export MAGs quality containing genome, completeness, and contamination for dRep
mags_quality <- mags_stats %>%
  select(bin_name, completeness, contamination) %>%
  rename(genome = bin_name) %>%
  mutate(genome = paste0(genome, ".fasta"))

write.table(mags_quality,
            file = "mags_quality.csv",
            sep = ",",
            quote = FALSE,
            row.names = FALSE)

#### Summarize gtdb-tk files
# List all gtdb tsv files in the folder
gtdb_tsv_files <- list.files(path = "gtdbtk_summary", pattern = "\\.tsv$", full.names = TRUE)

# Import all .tsv files and combine into one data frame

gtdb_tsv <- lapply(gtdb_tsv_files, function(file) {
  sample_name <- str_split(basename(file), "_", simplify = TRUE)[1] # Extract sample name (everything before first "_")
  df <- read.delim(file, stringsAsFactors = FALSE) # Read file
  colnames(df)[1] <- "bin_number" # Rename 1st column to bin_number
  
  # Add sample_name as first column
  df <- df %>%
    mutate(sample_name = sample_name,
           bin_name = paste0(sample_name, ".", bin_number)) %>%
    select(sample_name, bin_number, bin_name, everything())
  
  return(df)
}) %>% bind_rows()


#### Combine Info
# Import pertinent data tables from drep_out
genomeInformation <- read.table(file = "drep_out/data_tables/genomeInformation.csv",
                                sep=",", header=TRUE)
widb <- read.table(file = "drep_out/data_tables/Widb.csv",
                   sep=",", header=TRUE)
cdb <- read.table(file = "drep_out/data_tables/Cdb.csv",
                   sep=",", header=TRUE)
sdb <- read.table(file = "drep_out/data_tables/Sdb.csv",
                  sep=",", header=TRUE)

# Combine information from metaWrap-checkM (mags_stats),
# dRep (genomeInformation), dRep (sdb), dRep (cdb), dRep (widb)
# and gtdb-tk (gtdb_tsv)
mags_info <- mags_stats %>%
  select(-lineage, -N50, -size) %>%  # remove info also found from gtdb-tk and drep
  # Add other MAGs statistics from dRep (genomeInformation)
  # N50 values from checkm and drep were different. Since drep used its internally calculated N50 for scoring, I chose to report it
  left_join(
    genomeInformation %>%
      mutate(bin_name = sub("\\.fasta$", "", genome)) %>%
      select(bin_name, length, N50, centrality),
    by = "bin_name") %>%
  # Add scoring from dRep (sdb)
  left_join(
    sdb %>%
      mutate(bin_name = sub("\\.fasta$", "", genome)) %>%
      select(bin_name, score),
    by = "bin_name") %>%
  # Add clustering info from dRep (cdb)
  left_join(
    cdb %>%
      mutate(bin_name = sub("\\.fasta$", "", genome)) %>%
      select(bin_name, secondary_cluster) %>% 
      rename(cluster = secondary_cluster),
    by = "bin_name") %>%
  # Identify SRGs through dRep(widb)
  mutate(is_SRG = bin_name %in% sub("\\.fasta$", "", widb$genome)) %>%
  mutate(SRGs = if_else(is_SRG, "★", "")) %>%
  select(-is_SRG) %>%
  # Add info from gtdb-tk
  left_join(
    gtdb_tsv %>% select(-sample_name, -bin_number, -other_related_references.genome_id.species_name.radius.ANI.AF.),
    by = "bin_name") %>%
  # Arrange by cluster and score
  separate(cluster, into = c("primary", "secondary"), sep = "_", convert = TRUE) %>%
  arrange(primary, secondary, desc(score)) %>%
  group_by(primary, secondary) %>%
  ungroup() %>%
  unite(cluster, primary, secondary, sep = "_", remove = TRUE)

#### create excel table - info on MAGs

# Create a new workbook
mags_info_wb <- createWorkbook()

# Add first sheet
sheet1 <- createSheet(mags_info_wb, sheetName = "sheet 1")
addDataFrame(as.data.frame(mags_info), sheet1, row.names = FALSE, col.names = TRUE)

# Save workbook
saveWorkbook(mags_info_wb, "mags_info.xlsx")

#### Plot
# get species representative genomes (SRGs) information
SRG_info <- genomeInformation %>%
  # keep only needed columns from genomeInformation
  select(genome, completeness, contamination, length, N50, centrality) %>%
  # keep only genomes that appear in widb (SRGs)
  semi_join(widb, by = "genome") %>%
  # attach SRG metadata from widb
  left_join(
    widb %>% select(genome, score, cluster, cluster_members),
    by = "genome"
  )

## generate Contamination vs Completeness for SRGs
gtdb_classification <- gtdb_tsv %>%
  select(bin_name, classification) %>%
  mutate(bin_name = paste0(bin_name, ".fasta")) %>%
  rename(genome = bin_name)

SRG_info_joined <- SRG_info %>%
  left_join(gtdb_classification, by = "genome") %>%
  separate(classification, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";", convert = FALSE)

# Define the quality regions as a data frame
quality_regions <- data.frame(
  xmin  = c(50, 90),
  xmax  = c(90, 100),
  ymin  = c(0, 0),
  ymax  = c(10, 5),
  quality = c("Medium-quality SRGs", "High-quality SRGs")
)

ggplot(SRG_info_joined, aes(x = completeness, y = contamination, color = phylum)) +
  # Background quality regions (with legend and actual transparency)
  geom_rect(
    data = quality_regions,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = quality),
    inherit.aes = FALSE,
    color = NA, alpha = 0.1
  ) +
  # Threshold lines
  geom_vline(xintercept = c(50, 90), linetype = "dashed", color = "grey35") +
  geom_hline(yintercept = c(10, 5),   linetype = "dashed", color = "grey35") +
  # Points
  geom_point(size = 3, alpha = 0.9) +
  labs(
    x = "Completeness (%)",
    y = "Contamination (%)",
    color = "Phylum",
    fill  = "Quality thresholds"
  ) +
  scale_x_continuous(limits = c(50, 100)) +
  # Colors for fills and points
  scale_color_viridis_d(option = "turbo") +
  scale_fill_manual(
    values = c("Medium-quality SRGs" = "orange3", "High-quality SRGs" = "blue3")
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 13),
    panel.grid.minor = element_blank()
  )

## generate Completeness and Contamination plot for MAGs and SRGs
# combine info for All MAGs and SRGs
SRGs_MAGs_combined <- bind_rows(
  genomeInformation %>% mutate(dataset = "MAGs") %>% select(genome, completeness, contamination, dataset),
  SRG_info %>% mutate(dataset = "SRGs") %>% select(genome, completeness, contamination, dataset)
)

# Plot for completeness
ggplot(SRGs_MAGs_combined, aes(x = completeness, fill = dataset)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 10, color = "black") +
  scale_fill_manual(values = c("MAGs" = "#a6cee3", "SRGs" = "#1f78b4")) +
  labs(
    x = "Completeness (%)",
    y = "Count",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Plot for contamination
ggplot(SRGs_MAGs_combined, aes(x = contamination, fill = dataset)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 10, color = "black") +
  scale_fill_manual(values = c("MAGs" = "#a6cee3", "SRGs" = "#1f78b4")) +
  labs(
    x = "Contamination (%)",
    y = "Count",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

## upset plot for Shared MAGs
# make wide boolean matrix (samples = columns, clusters = rows)
cluster_matrix <- mags_info %>%
  select(sample_name, cluster) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = sample_name, values_from = present, values_fill = 0)

# sample columns (all except cluster)
samples <- setdiff(colnames(cluster_matrix), "cluster")

UpSetR::upset(
  as.data.frame(cluster_matrix),
  sets = samples,
  nintersects = NA,   # how many intersections to show (top 30)
  order.by = "freq",   # order by most common intersections
  main.bar.color = "steelblue",
  sets.bar.color = "gray40",
  point.size = 3,
  line.size = 1.2,
  text.scale = 1.3,
  mainbar.y.label = "Number of shared MAGs",  # top bars
  sets.x.label = "MAGs per Sample",       # left bars
)
