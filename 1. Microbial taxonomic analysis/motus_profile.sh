#!/bin/bash

# List of sample names
samples=(AUGLS01 AUGLS05 AUGLS16 SEPLS01 SEPLS05 SEPLS16) 

# Make sure the output folder exists
mkdir -p taxo_profile

# Loop through each sample
for sample in "${samples[@]}"; do
    echo "Processing $sample..."
    motus profile \
        -f "clean_reads/${sample}_1_cleaned.fastq.gz" \
        -r "clean_reads/${sample}_2_cleaned.fastq.gz" \
        -I "taxo_profile/${sample}.bam" \
        -t 6 
done
