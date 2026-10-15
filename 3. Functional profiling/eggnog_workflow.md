# KEGG Profiling

## Set up of eggnog-mapper3
```bash
mamba create -n eggnog3_env \
    -c conda-forge -c bioconda \
    python=3.11 pip diamond mmseqs2 prodigal hmmer git \
    -y
```
