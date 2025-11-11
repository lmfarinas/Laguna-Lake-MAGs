#!/bin/bash

# Set directories
input_dir="scaffolds"
output_stats="$input_dir/scaffolds_stats.txt"
filtered_dir="scaffolds_min1000"
filtered_stats="$filtered_dir/scaffolds_min1000_stats.txt"

# Create filtered output directory if it does not exist
mkdir -p "$filtered_dir"

# Initialize output stats files with headers
echo -e "sample\tfile\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len" > "$output_stats"
echo -e "sample\tfile\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len" > "$filtered_stats"

# Loop through all .fasta files
for fasta_file in "$input_dir"/*.fasta; do
    filename=$(basename "$fasta_file")
    sample="${filename%.*}"

    # Get stats for the original file
    stats=$(seqkit stats "$fasta_file" | tail -n +2)
    echo -e "$sample\t$stats" >> "$output_stats"

    # Filter sequences ≥1000 bp
    filtered_file="$filtered_dir/${sample}_min1000.fasta"
    seqkit seq -m 1000 "$fasta_file" -o "$filtered_file"

    # Get stats for the filtered file
    filtered_stats_line=$(seqkit stats "$filtered_file" | tail -n +2)
    echo -e "$sample\t$filtered_stats_line" >> "$filtered_stats"
done