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
	> "taxo_profile/${sample}_profile_phylum.txt"
done

# merge profiles per sample into a single text file
motus merge -d taxo_profile -o "taxo_profile/profile_phylum.txt"
