# Phylogenetic analysis of hMPXV genomes

Commands used for genome assembly and phylogenetic analyses. Each numbered
script corresponds to the steps below; run them in order from the directory
holding the input files.

## Software

| Tool | Version |
|---|---|
| minimap2 | 2.26-r1175 |
| samtools | 1.18 (htslib 1.17) |
| MAFFT | 7.453 |
| SeaView (Gblocks) | 5.0.4 |
| IQ-TREE | 2.1.2 |
| Nextclade | 3.9.1 |
| TreeTime | 0.8.5 |

## 1. Consensus genome generation from Nanopore reads

`01_consensus.sh <MPXV_reference.fasta>`

1. Reads were mapped to the MPXV reference with minimap2, and a first-pass
   consensus was called with samtools.
2. Reads were then re-mapped to that sample-specific consensus, and
   alignments shorter than 1,000 bp were discarded.
3. A second, final consensus was called.

The final `*_ConMap_Consensus.fasta` genomes are concatenated into
`AllNanoporeMPXV_NM.fasta`.

## 2. Alignment

`02_alignment_ml_trees.sh`

Consensus genomes were added to the NC_063383.1 reference with MAFFT
`--addfragments`, keeping the reference coordinate system (`--keeplength`).

## 3. Masking of poorly aligned regions

Gblocks was run through the SeaView graphical interface (*Sites → Gblocks*).
This produced `AllNanoporeMPXV_NM_aligned_Gblocks_Clean_Aligned_GblockNoOut.fasta`.

## 4. Maximum likelihood phylogenies

`02_alignment_ml_trees.sh`

IQ-TREE trees were inferred from both the Gblocks-trimmed and the untrimmed
alignments. ModelFinder selected the substitution model by BIC, and branch
support was assessed with SH-aLRT (1,000 replicates). The selected models are
reported in the IQ-TREE `.log` files.

## 5. Contextual dataset

`03_contextual_dataset.sh`

- **Download:** 831 US hMPXV genomes collected from 1 January 2022 to
  18 April 2024 were downloaded from GenBank.
- **Lineage assignment and quality control:** Nextclade with the
  `nextstrain/mpox/all-clades` dataset (tag `2024-11-19--14-18-53Z`), keeping
  sequences with ≥ 99% genome coverage and an overall QC score < 100. Most
  excluded sequences had excess private mutations, insertions or deletions.
- **Tree:** the study and public sequences were aligned together against the
  reference, and a maximum likelihood tree was inferred as in step 4.

## 6. Time-scaled phylogeny

`04_timetree_mugration.sh <dates.csv> <regions.csv>`

TreeTime was used for root-to-tip regression, time tree estimation and
ancestral state reconstruction:

- relaxed clock (`--relax 1.0 0.5`);
- skyline coalescent prior;
- tips deviating more than 4 interquartile ranges from the clock regression
  were excluded (`--clock-filter 4`).

## 7. Ancestral reconstruction of geographic location

`04_timetree_mugration.sh`

US states were grouped into U.S. Census Bureau regions (census.gov), and the
region of each node was reconstructed with `treetime mugration`.

## Data availability

Whole-genome sequences generated in this study are available through NCBI
GenBank under accession numbers PZ534571–PZ534651.
