#!/bin/bash

FASTA_FILE="$1"
CONTAM_FILE="$2"
OUTFILE="cleaned_${FASTA_FILE}"

if [[ -z "$FASTA_FILE" || -z "$CONTAM_FILE" ]]; then
  echo "❌ Error: Missing input files."
  echo "Usage: $0 <fasta_file> <contamination_file>"
  exit 1
fi

if [[ ! -f "$FASTA_FILE" ]]; then
  echo "❌ Error: FASTA file '$FASTA_FILE' not found."
  exit 1
fi

if [[ ! -f "$CONTAM_FILE" ]]; then
  echo "❌ Error: Contamination report '$CONTAM_FILE' not found."
  exit 1
fi

echo "✅ Inputs found. Processing..."

awk '/^Pd8825_contig/ {
  split($3, coords, "\\.\\.");
  print $1 "\t" coords[1] "\t" coords[2];
}' "$CONTAM_FILE" > adaptor_regions.tsv

awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq = seq $0} END {print seq}' "$FASTA_FILE" > temp.fasta

> "$OUTFILE"
total_written=0
total_skipped=0

trim_Ns() {
  local seq="$1"
  # Remove leading Ns
  seq="${seq#"${seq%%[!N]*}"}"
  # Remove trailing Ns
  seq="${seq%"${seq##*[!N]}"}"
  echo "$seq"
}

while read -r header sequence; do
  contig=$(echo "$header" | cut -d' ' -f1 | sed 's/>//')
  awk -v contig="$contig" '$1 == contig' adaptor_regions.tsv > contig_regions.tsv
  seq_length=${#sequence}

  if [[ ! -s contig_regions.tsv ]]; then
    # No contamination, trim Ns and keep if long enough
    sequence=$(trim_Ns "$sequence")
    if (( ${#sequence} >= 200 )); then
      echo "$header" >> "$OUTFILE"
      echo "$sequence" >> "$OUTFILE"
      ((total_written++))
    else
      echo "⚠️ Skipped $contig (length after trimming Ns < 200 bp)"
      ((total_skipped++))
    fi
    continue
  fi

  cut_points=()
  while read -r name start end; do
    if (( start <= 10 )); then
      sequence="${sequence:end}"
    elif (( end >= seq_length - 10 )); then
      sequence="${sequence:0:start-1}"
    else
      cut_points+=("$start $end")
    fi
  done < contig_regions.tsv

  if [[ ${#cut_points[@]} -gt 0 ]]; then
    start_pos=1
    seg_num=1
    for cp in "${cut_points[@]}"; do
      s=$(echo "$cp" | cut -d' ' -f1)
      e=$(echo "$cp" | cut -d' ' -f2)
      part="${sequence:start_pos-1:s-start_pos}"
      part=$(trim_Ns "$part")
      if (( ${#part} >= 200 )); then
        echo ">${contig}_split${seg_num}" >> "$OUTFILE"
        echo "$part" >> "$OUTFILE"
        ((total_written++))
      else
        echo "⚠️ Skipped ${contig}_split${seg_num} (after trimming Ns < 200 bp)"
        ((total_skipped++))
      fi
      start_pos=$((e + 1))
      ((seg_num++))
    done
    part="${sequence:start_pos-1}"
    part=$(trim_Ns "$part")
    if (( ${#part} >= 200 )); then
      echo ">${contig}_split${seg_num}" >> "$OUTFILE"
      echo "$part" >> "$OUTFILE"
      ((total_written++))
    else
      echo "⚠️ Skipped ${contig}_split${seg_num} (final segment < 200 bp)"
      ((total_skipped++))
    fi
  else
    sequence=$(trim_Ns "$sequence")
    if (( ${#sequence} >= 200 )); then
      echo "$header" >> "$OUTFILE"
      echo "$sequence" >> "$OUTFILE"
      ((total_written++))
    else
      echo "⚠️ Skipped $contig (after trim + Ns < 200 bp)"
      ((total_skipped++))
    fi
  fi

done < <(paste - - < temp.fasta)

rm -f temp.fasta contig_regions.tsv adaptor_regions.tsv

echo "✅ Cleaning complete."
echo "📁 Output FASTA: $OUTFILE"
echo "🧾 Contigs written: $total_written"
echo "🚫 Contigs skipped: $total_skipped"

