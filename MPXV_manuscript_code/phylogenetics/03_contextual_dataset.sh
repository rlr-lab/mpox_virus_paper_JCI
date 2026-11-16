#!/usr/bin/env bash
# ============================================================================
# Contextual dataset: QC of public US hMPXV genomes with Nextclade, joint
# alignment with the study genomes and maximum likelihood phylogeny.
#
# Usage: bash 03_contextual_dataset.sh
# Requires genbank_us.fasta (GenBank download), AllNanoporeMPXV_NM.fasta and
# NC_063383.1.fasta. Sequences passing QC (>= 99% genome coverage and overall
# QC score < 100 in nextclade.tsv) are saved as genbank_us_pass.fasta.
# ============================================================================
set -euo pipefail

# ---- Nextclade QC -------------------------------------------------------------
nextclade dataset get --name nextstrain/mpox/all-clades --tag 2024-11-19--14-18-53Z \
  --output-dir mpxv_dataset/

nextclade run --input-dataset mpxv_dataset/ \
  --output-tsv nextclade.tsv \
  genbank_us.fasta

# ---- Joint alignment and ML tree --------------------------------------------
cat AllNanoporeMPXV_NM.fasta genbank_us_pass.fasta > AllConsensusNM.fasta

mafft --auto --thread -1 --keeplength --addfragments AllConsensusNM.fasta \
  NC_063383.1.fasta > AllConsensusNM_Aligned.fasta

iqtree2 -s AllConsensusNM_Aligned.fasta -T AUTO --alrt 1000
