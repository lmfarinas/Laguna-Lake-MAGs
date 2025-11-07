# mOTUs profiling
Microbial taxonomic profiling was performed using mOTUs v3.1.0, a tool that estimates the abundance of known and unknown microbial species based on single-copy phylogenetic marker genes ([Ruscheweyh et al., 2022](https://doi.org/10.1186/s40168-022-01410-z)). Reads were first mapped to the mOTUs marker gene database, and then the mapped reads were assigned to known microbial taxa. The abundance of each taxon, down to the genus level, was determined using the default counting algorithm. The relative abundance of each taxon was then calculated as the number of reads mapped to the taxon divided by the total number of mapped reads.

### Script: `motus_profile.sh`
- loops through the samples' clean reads in the `clean_reads/` folder
- runs ``bwa`` through ``motus profile`` against the mOTUs marker gene database and generates *.bam* files in the `taxo_profile/` folder.
- produces **relative abundance** data at phylum level in .txt format. `taxo_profile/phylum/profile_phylum.txt` contains the merged data for all samples.

```bash
#!/bin/bash

# List of sample names
samples=(AUGLS01 AUGLS05 AUGLS16 SEPLS01 SEPLS05 SEPLS16) 

# create the output folder
mkdir -p taxo_profile

# Loop through the clean reads for each sample
for sample in "${samples[@]}"; do
    echo "Processing $sample..."
    motus profile \
        -f "clean_reads/${sample}_1_cleaned.fastq.gz" \
        -r "clean_reads/${sample}_2_cleaned.fastq.gz" \
        -I "taxo_profile/${sample}.bam" \ # save BWA result as .bam file
        -t 8 \
	-k phylum \ # set taxonomic level
	-n "${sample}" \
	> "taxo_profile/phylum/${sample}_profile_phylum.txt"
done

# merge profiles per sample into a single text file
motus merge -d "taxo_profile/phylum" -o "taxo_profile/phylum/profile_phylum.txt"
```

### Create stacked barplot: `taxo_profile.R`
```r
# Load necessary library
library(readr)
library(tidyr)
library(dplyr)

# set '#' as indicator of the metadata lines in the motus output
# Clean the header of the motus output
{
  # Read file line by line
  lines <- readLines("taxo_profile/phylum/profile_phylum.txt")
  # Find the header line (starts with "#consensus_taxonomy")
  header_line <- grep("^#consensus_taxonomy", lines, value = TRUE)
  
  # remove leading '#'
  header <- gsub("^#", "", header_line)
  
  # Read the actual data, skipping metadata lines
  profile_phylum <- read_tsv(
    "taxo_profile/phylum/profile_phylum.txt",
    comment = "#",
    col_names = strsplit(header, "\t")[[1]]
  )
}

# Clean the dataframe
motus_abundance <- profile_phylum %>%
  filter(!if_all(AUGLS01:SEPLS16, ~ .x == 0))%>%  # remove 
  mutate(consensus_taxonomy = if_else(consensus_taxonomy == "Bacteria phylum incertae sedis",
    "Other Bacteria", consensus_taxonomy))  # rename Bacteria phylum incertae sedis to Other bacteria

```

```r
##### Plot for Unassigned Taxa ##### 
# Import CSV file
TaxoProfile_UA <- as.data.frame(read.csv(
  '/Users/elmashu/Documents/UP/Thesis/Data Analysis/R Files/Taxonomic Profile/TaxoProfile_phyla_unassigned.csv'))

#convert data to long format
TaxoProfile_UA_long <- pivot_longer(TaxoProfile_UA, cols = c("Unassigned", "Assigned"), 
                                 names_to = "Category", values_to = "Proportion")

#Set priority levels (what appears on top)
TaxoProfile_UA_long$Category <- factor(TaxoProfile_UA_long$Category, levels = c("Unassigned", "Assigned"))

#Convert decimal to percent
TaxoProfile_UA_long$Proportion <- TaxoProfile_UA_long$Proportion * 100

#Plot
ggplot(TaxoProfile_UA_long, aes(x = Sample, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Assigned" = "#80CBC4", "Unassigned" = "#D3D3D3")) +
  labs(y = "Relative Abundance (%)", x = NULL) +
  theme_minimal(base_family = "sans", base_size = 12) +
  theme(legend.position = "right", 
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        panel.background = element_blank())

##### Plot for Assigned Taxa w/ no unassigned ##### 
# Import CSV file
TaxoProfile_phyla <- as.data.frame(read.csv(
  '/Users/elmashu/Documents/UP/Thesis/Data Analysis/R Files/Taxonomic Profile/TaxoProfile_phyla_no_unassigned.csv'))

#convert data to long format
TaxoProfile_phyla_long <- pivot_longer(TaxoProfile_phyla,
                                    cols = c("Actinobacteria", "Proteobacteria", "Cyanobacteria", "Bacteroidetes", "Other.Bacteria", "Verrucomicrobia", "Planctomycetes", "Acidobacteria", "Gemmatimonadetes", "Chloroflexi", "Fusobacteria"), 
                                    names_to = "Category", values_to = "Proportion")

#Convert decimal to percent
TaxoProfile_phyla_long$Proportion <- TaxoProfile_phyla_long$Proportion * 100

#Set priority levels (what appears on top)
TaxoProfile_phyla_long$Category <- factor(TaxoProfile_phyla_long$Category,
                                               levels = rev(c("Actinobacteria", "Proteobacteria", "Cyanobacteria",
                                                          "Bacteroidetes", "Other.Bacteria", "Verrucomicrobia",
                                                          "Planctomycetes", "Acidobacteria", "Gemmatimonadetes",
                                                          "Chloroflexi", "Fusobacteria")))

# Plot
ggplot(TaxoProfile_phyla_long, aes(x = Sample, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Actinobacteria" = "#E69F00",
                               "Proteobacteria" = "#56B4E9",
                               "Cyanobacteria" = "#009E73",
                               "Bacteroidetes" = "#F0E442",
                               "Other.Bacteria" = "#0072B2",
                               "Verrucomicrobia" = "#D55E00",
                               "Planctomycetes" = "#CC79A7",
                               "Acidobacteria" = "#8B4513",
                               "Gemmatimonadetes" = "#800080",
                               "Chloroflexi" = "#CCFF33",
                               "Fusobacteria" = "#FF6666")) +
  labs(y = "Relative Abundance (%)", x = NULL, fill = "Phylum") +
  theme_minimal(base_family = "sans", base_size = 12) +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        panel.background = element_blank())
```
