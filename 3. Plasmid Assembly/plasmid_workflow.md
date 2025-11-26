Perform plasmid assembly via metaplasmidspades (Galaxy Version 3.15.5+galaxy2) using default parameters.

Contigs were downloaded into `metaplasmidspades_out/contigs`. Filename format: [samplename]_metaplasmidspades_contigs.fasta

Contigs were then renamed, appending the sample name at the beginning of their fasta headers. Then, all contigs were combined into a single fasta file, `all_metaplasmidspades_contigs.fasta`.

```bash
for f in *.fasta; do
    sample="${f%%_*}"
    awk -v p="${sample}_" '
        /^>/ { sub(/^>/, ">" p) }
        { print }
    ' "$f" >> all_metaplasmidspades_contigs.fasta
done
```
Moreover, the putative plasmids were dereplicated using MMseqs2 Version: 18.8cc5c:

```bash
mmseqs easy-cluster all_metaplasmidspades_contigs.fasta plasmids99 tmp \
    --min-seq-id 0.99 \
    -c 0.9 \
    --cov-mode 1

mmseqs easy-cluster all_metaplasmidspades_contigs.fasta plasmids95 tmp \
    --min-seq-id 0.95 \
    -c 0.9 \
    --cov-mode 1
```


to do: https://github.com/ablab/plasmidVerify, https://github.com/ablab/viralVerify


use viralVerify bioconda v 1.1
```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
```

