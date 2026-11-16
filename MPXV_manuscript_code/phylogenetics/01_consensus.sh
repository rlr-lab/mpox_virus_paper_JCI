#!/usr/bin/env bash
# ============================================================================
# Consensus genome generation from Nanopore reads (two-round reference mapping)
#
# Usage: bash 01_consensus.sh <MPXV_reference.fasta>
# Run in the directory containing the per-sample *.fastq.gz files.
# ============================================================================
set -euo pipefail

REF=${1:?Usage: bash 01_consensus.sh <MPXV_reference.fasta>}

for i in *.fastq.gz; do
  base=$(basename -s .fastq.gz "$i")

  # Round 1: map to reference, call consensus
  minimap2 -ax map-ont "$REF" "$i" > "${base}.sam"
  samtools view -bS -F 4 "${base}.sam" > "${base}.bam"
  samtools sort -o "${base}_sorted.bam" "${base}.bam"
  samtools index "${base}_sorted.bam"
  samtools consensus -m simple -A -q -d 1000 --threads 4 "${base}_sorted.bam" \
    -o "${base}_ConsensusR1.fasta"

  # Round 2: re-map to sample consensus, filter short alignments, re-call
  minimap2 -ax map-ont "${base}_ConsensusR1.fasta" "$i" > "${base}_ConMap.sam"
  samtools view -bS -F 4 "${base}_ConMap.sam" > "${base}_ConMap.bam"
  samtools view -e 'rlen>1000' -O BAM -o "${base}_ConMap_Filter.bam" "${base}_ConMap.bam"
  samtools sort -o "${base}_ConMap_sorted.bam" "${base}_ConMap_Filter.bam"
  samtools index "${base}_ConMap_sorted.bam"
  samtools consensus -A --threads 4 "${base}_ConMap_sorted.bam" \
    -o "${base}_ConMap_Consensus.fasta"
done

# Combine final consensus genomes for alignment
cat *_ConMap_Consensus.fasta > AllNanoporeMPXV_NM.fasta
