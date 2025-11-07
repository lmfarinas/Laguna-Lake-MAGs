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

