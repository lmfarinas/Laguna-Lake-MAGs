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


```bash
cat > run_MAG_kegg_srun.sh <<'EOF'
#!/bin/bash
set -uo pipefail

ROOT="$(pwd)"
MAG_DIR="$ROOT/metawrap_50_10_bins/ALL_Samples"
MANIFEST="$ROOT/mag_list.txt"
OUT_DIR="$ROOT/kegg_MAG_annotations"
LOG_DIR="$ROOT/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: mag_list.txt not found."
    exit 1
fi

N=$(wc -l < "$MANIFEST" | tr -d ' ')
echo "MAGs in manifest: $N"

if [[ "$N" -eq 0 ]]; then
    echo "ERROR: mag_list.txt is empty."
    exit 1
fi

if [[ "${1:-}" != "--inside" ]]; then
    srun \
        --partition=batch \
        --qos=batch_default \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=8 \
        --mem=32G \
        --time=24:00:00 \
        bash "$ROOT/run_MAG_kegg_srun.sh" --inside \
        2>&1 | tee -a "$LOG_DIR/kegg_batch.log"
    exit "${PIPESTATUS[0]}"
fi

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate eggnog2_env

echo "Host: $(hostname)"
echo "Started: $(date)"
emapper.py --version

completed=0
skipped=0
failed=0

while IFS= read -r MAG || [[ -n "$MAG" ]]; do
    [[ -z "$MAG" ]] && continue

    NAME="$(basename "$MAG" .fasta)"
    ANNOT="$OUT_DIR/${NAME}.emapper.annotations"
    DONE="$OUT_DIR/${NAME}.done"

    if [[ -f "$DONE" && -s "$ANNOT" ]]; then
        echo "SKIP: $NAME (already completed)"
        ((skipped+=1))
        continue
    fi

    if [[ ! -f "$MAG" ]]; then
        echo "ERROR: MAG file not found: $MAG"
        ((failed+=1))
        continue
    fi

    echo "[$(date)] Starting $NAME"

    if emapper.py \
        -i "$MAG" \
        --itype genome \
        --genepred prodigal \
        --output "$NAME" \
        --output_dir "$OUT_DIR" \
        --cpu "$SLURM_CPUS_PER_TASK"
    then
        if [[ -s "$ANNOT" ]]; then
            touch "$DONE"
            echo "SUCCESS: $NAME"
            ((completed+=1))
        else
            echo "ERROR: Annotation file missing or empty: $NAME"
            ((failed+=1))
        fi
    else
        echo "ERROR: eggNOG-mapper failed: $NAME"
        ((failed+=1))
    fi

    echo "Progress: completed=$completed skipped=$skipped failed=$failed"
done < "$MANIFEST"

echo "Finished: $(date)"
echo "Newly completed: $completed"
echo "Previously completed: $skipped"
echo "Failed: $failed"

[[ "$failed" -eq 0 ]]
EOF

chmod +x run_MAG_kegg_srun.sh
tmux new -s kegg_MAGs
./run_MAG_kegg_srun.sh
```
