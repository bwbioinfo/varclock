#!/usr/bin/env bash

###############################################################################
# tag_bams
#
# Description:
#   Extracts reads supporting specific variants from a TSV file,
#   tags them with sample information, and outputs sorted and indexed BAMs.
#
# Usage:
#   ./tag_bams --vcf VARIANTS.vcf --tsv SUPPORT.tsv --bam INPUT.bam --output OUTDIR
#
# Author: Nicholas Geoffrion
# Created: 2025-07-25
# Version: 1.1
###############################################################################

set -euo pipefail

usage() {
  echo "Usage: $0 --vcf VCF --tsv TSV --bam BAM --output OUTDIR"
  exit 1
}

VCF=""
TSV=""
BAM=""
OUTDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) VCF="$2"; shift 2 ;;
    --tsv) TSV="$2"; shift 2 ;;
    --bam) BAM="$2"; shift 2 ;;
    --output) OUTDIR="$2"; shift 2 ;;
    *) echo "Unknown option $1"; usage ;;
  esac
done

[[ -f "$VCF" ]] || { echo "VCF not found: $VCF"; exit 1; }
[[ -f "$TSV" ]] || { echo "TSV not found: $TSV"; exit 1; }
[[ -f "$BAM" ]] || { echo "BAM not found: $BAM"; exit 1; }
[[ -d "$OUTDIR" ]] || mkdir -p "$OUTDIR"

TMPDIR=$(mktemp -d)
echo "Using temporary directory: $TMPDIR"

ml samtools
ml bcftools

bcftools query -f '%CHROM\t%POS\t%ALT\n' "$VCF" |
while IFS=$'\t' read -r chrom pos alt; do
  echo -e "\n>>> Processing $chrom $pos $alt"

  safe_alt=$(echo "$alt" | sed 's/[^A-Za-z0-9._-]/_/g')
  safe_chrom=$(echo "$chrom" | sed 's/[^A-Za-z0-9._-]/_/g')
  outprefix="${safe_chrom}_${pos}_${safe_alt}"
  fullprefix="$TMPDIR/$outprefix"

  # Extract read names
  echo "Extracting read names..."
  grep -F "$chrom" "$TSV" | grep -F "$pos" | grep -F "$alt" | awk '{print $1}' > "${fullprefix}.txt"
  echo "Read names:"
  cat "${fullprefix}.txt"

  echo "Extracting read tags..."
  grep -F "$chrom" "$TSV" | grep -F "$pos" | grep -F "$alt" | awk '{print $1"\t"$3}' > "${fullprefix}.tags.txt"
  echo "Tags:"
  cat "${fullprefix}.tags.txt"

  # Extract reads from BAM
  echo "Extracting reads from BAM..."
  samtools view -@ 30 -N "${fullprefix}.txt" "$BAM" -b -o "${fullprefix}.bam"

  # Convert to SAM
  echo "Converting BAM to SAM..."
  samtools view -@ 30 -h "${fullprefix}.bam" > "${fullprefix}.sam"

  echo "Adding SM tags to SAM..."
  awk -v tagfile="${fullprefix}.tags.txt" '
    BEGIN {
      while ((getline < tagfile) > 0) {
        tags[$1] = $2
      }
    }
    {
      if ($1 ~ /^@/) {
        print $0
      } else {
        if ($1 in tags) {
          print $0 "\tSM:Z:" tags[$1]
        } else {
          print $0
        }
      }
    }
  ' "${fullprefix}.sam" > "${fullprefix}.tagged.sam"

  echo "Converting SAM back to BAM..."
  samtools view -@ 30 -bS "${fullprefix}.tagged.sam" > "${fullprefix}.tagged.bam"

  echo "Sorting and indexing..."
  samtools sort -@ 30 "${fullprefix}.tagged.bam" -o "${OUTDIR}/${outprefix}.sorted.bam"
  samtools index -@ 30 "${OUTDIR}/${outprefix}.sorted.bam"

  echo "✔ Written ${OUTDIR}/${outprefix}.sorted.bam"
done

rm -rf "$TMPDIR"
echo "🧹 Cleaned up temporary files"
