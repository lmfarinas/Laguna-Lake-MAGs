# KEGG Profiling

## Set up of eggnog-mapper3

emapper-3.0.0-beta6 / eggNOG DB version: 7.0.0 / Diamond version found: diamond version 2.2.8 / MMseqs2 version found: 18.8cc5c
```bash
# Set up a conda environment
mamba create -n eggnog3_env \
    -c conda-forge -c bioconda \
    python=3.11 pip diamond mmseqs2 prodigal hmmer git \
    -y

conda activate eggnog3_env

# Install eggnog-mapper
mkdir -p ~/software
cd software

git clone https://github.com/eggnogdb/eggnog-mapper.git
cd eggnog-mapper

python -m pip install .

# Download database
mkdir -p ~/databases/eggnog7

# there seems to be a bug in the download_eggnog_data.py of eggNOG-mapper v3.0.0-beta6
# It is trying to resolve a default database directory inside the Conda environment:
# /miniconda3/envs/eggnog3_env/lib/python3.11/site-packages/data

# resolved by creating the expected directory:
mkdir -p \
  /home/lance.farinas/miniconda3/envs/eggnog3_env/lib/python3.11/site-packages/data

# then download
download_eggnog_data.py \
    -y \
    --data_dir ~/databases/eggnog7
```
Create a symlink to point from the default to the actual directory of the database
```bash
DB_DEFAULT="$CONDA_PREFIX/lib/python3.11/site-packages/data"
DB_ACTUAL="/scratch2/lance.farinas/databases/eggnog7"

# Remove the empty default directory
rmdir "$DB_DEFAULT"

# Link it to your downloaded eggNOG v7 database directory
ln -s "$DB_ACTUAL" "$DB_DEFAULT"

# Verify the symlink and database path
ls -ld "$DB_DEFAULT"
readlink -f "$DB_DEFAULT/eggnog.db"
```
