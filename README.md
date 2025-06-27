# Taxonomic_-_BGC_Diversity
Microbial communities play a vital role in ecosystem function and stability, with natural products, which are specialized metabolites produced by organisms, being key to many of their interactions. These compounds are often encoded within biosynthetic gene clusters (BGCs) belonging to a gene cluster family (GCF), which have become a major focus in natural product discovery. Leveraging developments in metagenomics and genome mining, this project investigates the relationship between taxonomic diversity and biosynthetic diversity across different biomes using publicly available datasets from MGnify and the BGC Atlas.

This repository contains a collection of scripts for extracting, analyzing, and computing diversity metrics from processed genomic contig data. The pipeline handles downloading data from JSON files, counting GCF IDs, and computing alpha diversity statistics using both Bash and Python.

---

## Project Structure

```text
.
├── get_urls.py                  # Extracts processed contig download URLs from JSON assembly files
├── download_fasta.py           # Downloads contig FASTA files from extracted URLs
├── generate_pbs_richness.sh    # Creates PBS scripts to compute richness from .report files
├── count_gcf_ids.sh            # Counts distinct GCF IDs per sample
├── calculate_shannon.sh        # Computes Shannon-Wiener diversity index per sample
├── run_alpha_diversity.pbs     # PBS job to compute Shannon diversity from .bracken files
├── alpha_diversity.py          # Python script to compute alpha diversity (called by PBS job)
├── analysis.R                  # Initial R pipeline: genus-level taxonomy and preprocessing
├── rank_analysis.R             # R script extending analysis to family/species and comparisons
├── sample_taxa_analysis.R      # R script for comparative diversity between GCF and taxonomy
├── README.md                   # Project documentation (this file)
```

---

## Requirements

- **Python 3** with `pandas`, `requests`
- **GNU awk**
- **PBS/Torque scheduler** (or adapt scripts for SLURM, etc.)
- **R** with common tidyverse/data wrangling packages (e.g., `dplyr`, `ggplot2`)

---

### Python and Bash Components

| Script                     | Description                                                                                      |
| -------------------------- | ------------------------------------------------------------------------------------------------ |
| `get_urls.py`              | Parses assembly JSONs and retrieves download links for processed contigs                         |
| `download_fasta.py`        | Downloads `.fasta.gz` files from URLs                                                            |
| `generate_pbs_richness.sh` | Creates PBS scripts to compute family/genus/species richness from Kraken/Bracken `.report` files |
| `count_gcf_ids.sh`         | Counts distinct GCF IDs per sample, including unclassified `-1`                                  |
| `calculate_shannon.sh`     | Calculates Shannon-Wiener diversity index from GCF distributions                                 |
| `run_alpha_diversity.pbs`  | Batch submits diversity calculations using `.bracken` files                                      |
| `alpha_diversity.py`       | Computes Shannon index from `.bracken` files with error handling                                 |

---

### R Analysis Pipelines

| Script                   | Description                                                                                                                        |
| ------------------------ | ---------------------------------------------------------------------------------------------------------------------------------- |
| `analysis.R`             | Initial exploratory analysis using BGC Atlas metadata. Preprocesses input and examines genus-level taxonomy.                       |
| `rank_analysis.R`        | Extends taxonomy analysis to family and species ranks; compares diversity across taxonomic levels.                                 |
| `sample_taxa_analysis.R` | Performs comparative analysis between GCF-based diversity (from Bash scripts) and taxonomy-based diversity (from `.report` files). |

These R scripts are meant to be run interactively or in RStudio, and they assume that inputs such as `sample_gcf_counts.tsv`, `sample_gcf_shannon.tsv`, and richness `.tsv` outputs from PBS jobs are available.

---

## Notes

- Scripts assume column 2 = sample ID, column 14 = GCF ID in `regions_taxonomy.tsv`.
- PBS scripts are tuned for Torque environments — modify headers if using SLURM or others.
- Make sure the R scripts have access to the cleaned TSV outputs for downstream comparison.
- Diversity calculations in R use Shannon and richness indices consistent with Bash results.

---

## Contact

For questions, suggestions, or issues, feel free to open an issue or contact me.
