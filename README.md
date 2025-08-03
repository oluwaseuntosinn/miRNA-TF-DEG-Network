# miRNA-TF-DEG Network Analysis

This R script constructs a validated miRNA-target interaction network for a set of differentially expressed genes (DEGs), identifies highly interacting miRNAs, and explores shared TF-miRNA regulatory targets.

## Features

- Uses `multiMiR` to extract validated miRNA-target interactions
- Filters high-confidence interactions based on a scoring threshold
- Constructs and visualizes a directed igraph-based miRNA-target network
- Identifies transcription factor (TF) and miRNA crosstalk

## Requirements

- R packages: `multiMiR`, `igraph`, `dplyr`

## Running the Script

1. **Prepare your input files** and place them in the `data/` directory:
   - `miRNA_target_for_TF_miRNA.csv`
   - `tf_targets_for_TF_miRNA.csv`

   > ⚠️ These input files are not included in this repository. You must provide them yourself.

2. Run the script in R:
```R
source("miRNA_TF_DEG_Network.R")