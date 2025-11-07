## Assembly
Galaxy - Metaspades

## Binning

Galaxy - Metawrap

downloaded quality information from galaxy into the `MAGs_stats/` folder:

- [sampleName]_metawrap_50_10_bins.tabular

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
## Dereplication via dRep
`dRep version 3.6.1`

```bash
dRep dereplicate drep_out \    # output folder: drep_out
  -g metawrap_50_10_bins/ALL/*.fasta \    # location of bins
  --genomeInfo mags_quality.csv \    # contains quality information of mags (completeness, contamination)
  -p 6    # 6 threads used
```
## Taxonomic Classification of MAGs
gtdbtk via galaxy

