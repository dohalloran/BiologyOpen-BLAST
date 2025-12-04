#!/usr/bin/env bash
set -euo pipefail

# Run from the BLAST directory:
#   cd /path/to/BLAST
#   ./run_blast_and_filter.sh

BLASTP="./ncbi-blast-2.10.0+/bin/blastp"
DB="./db/ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa"
QUERY="./MajorityVote_targets.fasta"
BLAST_OUT="./MajorityVote_targets_vs_ancylostoma.tsv"

echo "Running BLASTP..."
"$BLASTP" \
  -query "$QUERY" \
  -db "$DB" \
  -evalue 1e-2 \
  -max_target_seqs 10 \
  -num_threads 4 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -out "$BLAST_OUT"

echo "BLAST complete, results in $BLAST_OUT"

echo "Filtering BLAST hits at e <= 1e-4 and mapping back to DrugBank IDs..."
python3 filter_blast_and_map_drugbank.py \
  --blast "$BLAST_OUT" \
  --fasta "$QUERY" \
  --targets_csv "./DrugBank_Targets.csv" \
  --evalue_cutoff 1e-4 \
  --out_fasta "./MajorityVote_targets_filtered.fasta" \
  --out_csv "./MajorityVote_targets_filtered_mapping.csv"

echo "Done."
echo "  Filtered FASTA: MajorityVote_targets_filtered.fasta"
echo "  Mapping CSV   : MajorityVote_targets_filtered_mapping.csv"
