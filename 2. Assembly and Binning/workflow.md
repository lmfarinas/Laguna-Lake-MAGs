## Assembly
Galaxy - Metaspades

## Binning

Galaxy - Metawrap

downloaded quality information from galaxy into the `MAGs_stats/` folder:

- [sampleName]_metawrap_50_10_bins.tabular

Using R, quality information of the MAGs (>=50 completeness, <=10 contamination) was combined into a single dataframe.
```r
#Load packages:
library(dplyr)
library(stringr)

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
```
The script above produces the `mags_quality.csv` file, which will be used by dRep for genome quality information and to skip CheckM.

## Dereplication via dRep
`dRep version 3.6.1`

```bash
# Run dRep dereplication
# Output folder: drep_out
# Genomes from metawrap_50_10_bins/ALL/
# Genome quality info from mags_quality.csv
# Completeness ≥ 50%, contamination ≤ 10%
# Primary clustering ANI = 0.90, secondary = 0.99
# Using 6 threads

dRep dereplicate drep_output \
  -g metawrap_50_10_bins/ALL/*.fasta \
  --genomeInfo mags_quality.csv \
  -comp 50 -con 10 \
  -pa 0.90 -sa 0.95 \
  -p 6
```

output of drep in `drep_out/` folder.

## Taxonomic Classification of MAGs
gtdbtk via galaxy

**NOTE: rerun gtdbtk now using the dereplicated fasta**

gtdbtk classification files were downloaded into `gtdbtk_summary/` folder:

- [sampleName]_gtdbtk.bac120.summary.tsv

The files were tabulated on R
```r
#Load packages:
library(dplyr)
library(stringr)

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
           bin_name = paste0(sample_name, "_", bin_number)) %>%
    select(sample_name, bin_number, bin_name, everything())
  
  return(df)
}) %>% bind_rows()
```

#### Plots
Classification of MAGs: https://pmc.ncbi.nlm.nih.gov/articles/PMC12114606/#B40-microorganisms-13-00985

- High-quality MAGs: >90% completeness, <5% contamination.
- Medium-quality MAGs: >50% completeness, <10% contamination.
- Low-quality MAGs: <50% completeness, <10% contamination.

````r
#Load packages:
library(dplyr)

# Import pertinent data tables
genomeInformation <- read.table(file = "drep_out/data_tables/genomeInformation.csv",
                                sep=",", header=TRUE)

widb <- read.table(file = "drep_out/data_tables/Widb.csv",
                   sep=",", header=TRUE)

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

# generate QC plot (Contam vs Comple)
# note: categorize by phyla

ggplot(SRG_info, aes(x = completeness, y = contamination, color = cluster)) +
  # Low-quality region (<50% completeness, <10% contamination)
  geom_rect(aes(xmin = 0, xmax = 50, ymin = 0, ymax = 10),
            fill = "red4", alpha = 0.12, inherit.aes = FALSE) +
  # Medium-quality region (50–90% completeness, <10% contamination)
  geom_rect(aes(xmin = 50, xmax = 90, ymin = 0, ymax = 10),
            fill = "yellow4", alpha = 0.12, inherit.aes = FALSE) +
  # High-quality region (>90% completeness, <5% contamination)
  geom_rect(aes(xmin = 90, xmax = 100, ymin = 0, ymax = 5),
            fill = "green4", alpha = 0.12, inherit.aes = FALSE) +
  geom_point(size = 3, alpha = 0.9) +
  # Reference threshold lines
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 90, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "grey40") +
  labs(
    x = "Completeness (%)",
    y = "Contamination (%)",
    color = "Cluster",
    title = "Genome QC Plot (CheckM Quality Reference Ranges)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 13),
    panel.grid.minor = element_blank()
  )

```

