# mOTUs profiling
Microbial taxonomic profiling was performed using mOTUs v3.1.0, a tool that estimates the abundance of known and unknown microbial species based on single-copy phylogenetic marker genes ([Ruscheweyh et al., 2022](https://doi.org/10.1186/s40168-022-01410-z)). Reads were first mapped to the mOTUs marker gene database, and then the mapped reads were assigned to known microbial taxa. The abundance of each taxon, down to the genus level, was determined using the default counting algorithm. The relative abundance of each taxon was then calculated as the number of reads mapped to the taxon divided by the total number of mapped reads.

### Script: `motus_profile.sh`
- loops through the samples' clean reads in the `clean_reads/` folder
- runs in  ``bwa`` against the mOTUs marker gene database and generates *.bam* files in the `taxo_profile/` folder.
