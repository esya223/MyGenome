#!/bin/bash
#
# Script: exclude_contigs.sh
# Description: Excludes contigs from a FASTA file based on a Contamination.txt file.
# Usage: ./exclude_contigs.sh <input.fasta> <Contamination.txt>
# Output: excluded.fasta
#

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input.fasta> <Contamination.txt>"
  exit 1
fi

input_fasta="$1"
contamination_txt="$2"
output_fasta="excluded.fasta"

# Extract the sequence names to exclude from the Contamination.txt file
exclude_names=$(grep "^Pd8825_contig" "$contamination_txt" | awk '{print $1}' | uniq)

# Use seqkit to filter out the sequences
singularity run --app seqkit2900 /share/singularity/images/ccs/conda/amd-conda19-rocky9.sinf \
  seqkit seq -w 0 -n -s -v -i "$input_fasta" \
  | singularity run --app seqkit2900 /share/singularity/images/ccs/conda/amd-conda19-rocky9.sinf \
    seqkit grep -v -f <(echo "$exclude_names") \
  > "$output_fasta"

echo "Excluded sequences and saved to $output_fasta"
