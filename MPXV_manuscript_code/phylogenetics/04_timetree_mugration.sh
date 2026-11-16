#!/usr/bin/env bash
# ============================================================================
# Time-scaled phylogeny and ancestral reconstruction of geographic region
# with TreeTime.
#
# Usage: bash 04_timetree_mugration.sh <dates.csv> <regions.csv>
# Requires AllConsensusNM_Aligned.fasta and its IQ-TREE .treefile
# (from 03_contextual_dataset.sh).
# ============================================================================
set -euo pipefail

DATES=${1:?Usage: bash 04_timetree_mugration.sh <dates.csv> <regions.csv>}
REGIONS=${2:?Usage: bash 04_timetree_mugration.sh <dates.csv> <regions.csv>}

# ---- Time tree: relaxed clock, skyline coalescent, 4-IQR clock filter -------
treetime --confidence --relax 1.0 0.5 \
  --aln AllConsensusNM_Aligned.fasta \
  --tree AllConsensusNM_Aligned.fasta.treefile \
  --dates "$DATES" \
  --coalescent skyline \
  --clock-filter 4 \
  --outdir timetree

# ---- Ancestral reconstruction of US Census Bureau region --------------------
treetime mugration --tree timetree/timetree.nexus \
  --states "$REGIONS" --attribute region \
  --outdir mugration
