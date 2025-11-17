#Load packages:
library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)
library(xlsx)
library(ggplot2)
library(viridis)
library(UpSetR)
library(jsonlite)
library(gt)

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

# Export putative high quality MAGs for bakta (predict rRNAs and tRNAs)
putative_highQ_mags <- mags_stats %>%
  filter(completeness > 90, contamination < 5) %>%
  mutate(bin_name = paste0(bin_name, ".fasta")) %>%
  pull(bin_name)

writeLines(putative_highQ_mags, "putative_highQ_mags.txt")

#### Get info from bakta json files (run after getting json files from bakta)
# Folder containing Bakta JSON files
json_dir <- "bakta_out/json_output"

# 20 tRNAs
standard_AA <- c(
  "Ala", "Arg", "Asn", "Asp", "Cys",
  "Gln", "Glu", "Gly", "His", "Ile",
  "Leu", "Lys", "Met", "Phe", "Pro",
  "Ser", "Thr", "Trp", "Tyr", "Val"
)

# Function to parse one Bakta JSON
parse_bakta <- function(json_file) {
  j <- fromJSON(json_file)
  
  # rRNA detection
  rrna_products <-as_tibble(j$features) %>%
    filter(type == "rRNA") %>%
    pull(product) %>%
    na.omit() %>%
    unique()
  
  has_5s  <- ifelse(any(str_detect(rrna_products, "5S ribosomal RNA")),"Yes", "No")
  has_16s <- ifelse(any(str_detect(rrna_products, "16S ribosomal RNA")),"Yes", "No")
  has_23s <- ifelse(any(str_detect(rrna_products, "23S ribosomal RNA")),"Yes", "No")
  
  # tRNA detection
  trna_AA <-as_tibble(j$features) %>%
    filter(type == "tRNA") %>%
    pull(amino_acid) %>%
    na.omit() %>%
    unique()
  
  # Validate against 20 standard tRNAs
  missing <- setdiff(standard_AA, trna_AA)
  
  # Create table
  tibble(
    bin_name = tools::file_path_sans_ext(basename(json_file)),
    rRNA_5S  = has_5s,
    rRNA_16S = has_16s,
    rRNA_23S = has_23s,
    num_unique_tRNAs = 20-length(missing),
    missing_tRNAs = paste(missing, collapse = ", ")
  )
}
# Apply to all JSON files
json_files <- list.files(json_dir, pattern = "\\.json$", full.names = TRUE)
trna_rrna_info <- map_dfr(json_files, parse_bakta)

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
# Import information from seqkit stats output
seqkit_bin_stats <- read.table(file="metawrap_50_10_bins/bins_seqkit_stats.txt", header = TRUE, sep = "\t")

# Import pertinent data tables from drep_out
genomeInformation <- read.table(file = "drep_out/data_tables/genomeInformation.csv",
                                sep=",", header=TRUE)
widb <- read.table(file = "drep_out/data_tables/Widb.csv",
                   sep=",", header=TRUE)
cdb <- read.table(file = "drep_out/data_tables/Cdb.csv",
                   sep=",", header=TRUE)
sdb <- read.table(file = "drep_out/data_tables/Sdb.csv",
                  sep=",", header=TRUE)

# Combine information from metaWrap-checkM (mags_stats), seqkit (seqkit_bin_stats)
# dRep (genomeInformation), dRep (sdb), dRep (cdb), dRep (widb)
# and gtdb-tk (gtdb_tsv)
mags_info <- mags_stats %>%
  select(-lineage, -N50, -size) %>%  # remove info also found from gtdb-tk and drep
  rename(MetaWRAP_binner = binner) %>%
  # N50 values from checkm were different from seqkit and drep. Since drep used its internally calculated N50 for scoring, I chose to report it.
  # Moreover, drep's N50 is similar with seqkit, which calculates N50 including the N's of the scaffolds.
  # Add info on rrnas and trnas from bakta
  left_join(trna_rrna_info, by = "bin_name") %>%
  # Add MIMAG Quality Classification
  mutate(
    MIMAG_qual = case_when(
      completeness > 90 & contamination < 5 &
        rRNA_5S == "Yes"  & rRNA_16S == "Yes" & rRNA_23S == "Yes" &
        num_unique_tRNAs >= 18 ~ "High",
      completeness > 90 & contamination < 5 ~ "nHigh",
      completeness >= 50 & contamination < 10 ~ "Medium",
      TRUE ~ ""
    )) %>%
  # Add statistics from seqkit
  left_join(
    seqkit_bin_stats %>%
      mutate(bin_name = sub("\\.fasta$", "", file)) %>%
      rename(num_scaffolds = num_seqs, total_length = sum_len, min_length = min_len,
             avg_length = avg_len, max_length = max_len, scaffold_N50 = N50,
             scaffold_L50 = N50_num, total_N = sum_n) %>%
      select(bin_name, num_scaffolds, total_length, min_length, avg_length, max_length,
             scaffold_N50, scaffold_L50, total_N),
    by = "bin_name") %>%
  # Add other MAGs statistics from dRep (genomeInformation)
  left_join(
    genomeInformation %>%
      mutate(bin_name = sub("\\.fasta$", "", genome)) %>%
      select(bin_name, centrality),
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
  # Add classification from gtdb-tk
  left_join(
    gtdb_tsv %>% select(bin_name, classification),
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

#### Plots
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

## upset plot for Shared SRGs
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
  sets.x.label = "Shared SRGs",       # left bars
)


data <- data.frame(
  Category = c("A", "B", "C"),
  Value1 = c(10, 25, 15),
  Value2 = c(5, 12, 8)
)

## Create GT table for MAG counts per quality
quality_counts <- mags_info %>%
  # count ALL MAGs per quality
  count(MIMAG_qual, name = "MAGs") %>%
  # join SRG counts
  left_join(
    mags_info %>%
      filter(SRGs == "★") %>%
      count(MIMAG_qual, name = "SRGs"),
    by = "MIMAG_qual"
  ) %>%
  # replace NA with 0 for SRGs
  mutate(SRGs = replace_na(SRGs, 0)) %>%
  # reorder factor levels
  mutate(MIMAG_qual = factor(MIMAG_qual, levels = c("High", "nHigh", "Medium"))) %>%
  arrange(MIMAG_qual) %>%
  # add totals row
  bind_rows(
    summarise(.,
              MIMAG_qual = "Total",
              MAGs = sum(MAGs),
              SRGs = sum(SRGs)
    )
  )

quality_counts_gt <- quality_counts %>%
  gt() %>%
  cols_label(MIMAG_qual = "") 
quality_counts_gt
