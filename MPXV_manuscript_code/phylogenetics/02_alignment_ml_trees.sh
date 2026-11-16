#!/usr/bin/env bash
# ============================================================================
# Alignment to the NC_063383.1 reference and maximum likelihood phylogenies
# of the study genomes.
#
# Usage: bash 02_alignment_ml_trees.sh
# Requires AllNanoporeMPXV_NM.fasta (from 01_consensus.sh) and NC_063383.1.fasta.
#
# Gblocks masking is done interactively in SeaView (Sites -> Gblocks) between
# the alignment and tree steps; see README.md.
# ============================================================================
set -euo pipefail

# ---- Alignment (keeping reference coordinates) ------------------------------
mafft --auto --thread -1 --keeplength --addfragments AllNanoporeMPXV_NM.fasta \
  NC_063383.1.fasta > AllNanoporeMPXV_NM_aligned.fasta

# ---- Maximum likelihood trees (ModelFinder/BIC, SH-aLRT 1000) ----------------
# Gblocks-trimmed alignment
iqtree2 -s AllNanoporeMPXV_NM_aligned_Gblocks_Clean_Aligned_GblockNoOut.fasta \
  -T AUTO --alrt 1000

# Untrimmed alignment
iqtree2 -s AllNanoporeMPXV_NM_aligned_NoOut.fasta -T AUTO --alrt 1000
