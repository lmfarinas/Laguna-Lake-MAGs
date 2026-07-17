# Laguna-Lake-MAGs

## 0. Pre-processing
Pre-processing of raw reads using **Fastp v0.24.0** to trim adapter sequences, remove poly-G and poly-X tails, perform paired-end base correction, and filter reads based on quality (Phred >15)
see: https://github.com/lmfarinas/laguna-lake-ARG-Analysis/blob/main/ARG-Quantification-Workflow.md

## 1. Microbial taxonomic analysis
see:  [mOTUs_workflow.md](1.%20Microbial%20taxonomic%20analysis/mOTUs_workflow.md)

uploaded data (folders) in Zenodo:10.5281/zenodo.21405147
1. clean_reads
2. taxo_profile

## 2. Assembly and Binning
see: [MAGs_workflow.md](https://github.com/lmfarinas/Laguna-Lake-MAGs/blob/main/2.%20Assembly%20and%20Binning/MAGs_workflow.md)

uploaded data (folders) in Zenodo:10.5281/zenodo.21405147
1. scaffolds
2. scaffolds_min1000
3. metawrap_50_10_bins
4. bakta_out
5. gtdbtk_summary
6. drep_out
7. amrfinder_out
8. coverage
