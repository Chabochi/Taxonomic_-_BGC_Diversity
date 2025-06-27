# Taxonomic_-_BGC_Diversity
Microbial communities play a vital role in ecosystem function and stability, with natural products, which are specialized metabolites produced by organisms, being key to many of their interactions. These compounds are often encoded within biosynthetic gene clusters (BGCs) belonging to a gene cluster family (GCF), which have become a major focus in natural product discovery. Leveraging developments in metagenomics and genome mining, this project investigates the relationship between taxonomic diversity and biosynthetic diversity across different biomes using publicly available datasets from MGnify and the BGC Atlas.

This repository contains a collection of scripts for extracting, analyzing, and computing diversity metrics from processed genomic contig data. The pipeline handles downloading data from JSON files, counting GCF IDs, and computing alpha diversity statistics using both Bash and Python.

---

## 📁 Project Structure

```text
.
├── get_urls.py                  # Extracts processed contig download URLs from JSON assembly files
├── download_fasta.py           # Downloads contig FASTA files from extracted URLs
├── generate_pbs_richness.sh    # Creates PBS scripts to compute richness from .report files
├── count_gcf_ids.sh            # Counts distinct GCF IDs per sample
├── calculate_shannon.sh        # Computes Shannon-Wiener diversity index per sample
├── run_alpha_diversity.pbs     # PBS job to compute Shannon diversity from .bracken files
├── alpha_diversity.py          # Python script to compute alpha diversity (called by PBS job)
├── README.md                   # Project documentation (this file)
```

---

## 🗏️ Requirements

- Python 3 (for `get_urls.py`, `download_fasta.py`, `alpha_diversity.py`)
- `pandas`, `requests` Python libraries
- `awk` (GNU Awk recommended)
- PBS/Torque job scheduler (for job scripts)
- Bracken `.report` and `.bracken` files for diversity estimation

---

## 🔧 Scripts Overview

### 1. `get_urls.py`

Extracts processed contig download links from a list of assembly JSON URLs and creates a `.tsv` mapping file.

**Usage**:

```bash
python get_urls.py assemblies.txt -n output.tsv
```

**Inputs**:

- `assemblies.txt`: Text file with one JSON URL per line

**Outputs**:

- TSV file with sample IDs and their corresponding download URLs

---

### 2. `download_fasta.py`

Downloads FASTA files using the mapping file generated above.

**Usage**:

```bash
python download_fasta.py "sample123\t<download_url>" -o output_folder
```

**Inputs**:

- A line from the TSV file: `<sample_id>\t<url>`

**Outputs**:

- Downloaded `.fasta.gz` file in specified directory

---

### 3. `generate_pbs_richness.sh`

Generates PBS job scripts that compute richness (family, genus, species) for each sample using `.report` files.

**Usage**:

```bash
./generate_pbs_richness.sh
```

**Inputs**:

- Subdirectories in `parent_dir` each containing a `.report` file

**Outputs**:

- PBS job scripts in `scripts/`, ready for submission

---

### 4. `count_gcf_ids.sh`

Counts distinct GCF IDs per sample, treating `-1` values as a special "unclassified" case.

**Usage**:

```bash
./count_gcf_ids.sh
```

**Inputs**:

- `regions_taxonomy.tsv`: Sample and GCF ID in columns 2 and 14

**Outputs**:

- `sample_gcf_counts.tsv`: Sample and GCF count

---

### 5. `calculate_shannon.sh`

Calculates Shannon-Wiener diversity index per sample based on GCF ID distributions.

**Usage**:

```bash
./calculate_shannon.sh
```

**Inputs**:

- `regions_taxonomy.tsv`: Sample and GCF ID in columns 2 and 14

**Outputs**:

- `sample_gcf_shannon.tsv`: Sample and Shannon index

---

### 6. `run_alpha_diversity.pbs`

PBS script to batch compute alpha diversity from `.bracken` files using `alpha_diversity.py`.

**Usage**: Submit to your cluster:

```bash
qsub run_alpha_diversity.pbs
```

**Outputs**:

- `alpha_diversity_results.csv`: Sample and Shannon diversity values
- `error_log.txt`: Samples skipped due to division errors

---

## Output Examples

### `output.tsv`

```tsv
SAMPLE_001    https://example.org/download/SAMPLE_001.fasta.gz
SAMPLE_002    https://example.org/download/SAMPLE_002.fasta.gz
```

### `sample_gcf_counts.tsv`

```tsv
SAMPLE_001    42
SAMPLE_002    37
```

### `sample_gcf_shannon.tsv`

```tsv
SAMPLE_001    3.154211
SAMPLE_002    2.837493
```

---

## Notes

- Scripts assume column 2 = sample ID, column 14 = GCF ID in `regions_taxonomy.tsv`.
- PBS scripts are tuned for Torque environments — modify headers if using SLURM or others.
- Make sure `alpha_diversity.py` is executable and in your PATH or working directory.

---

## Contact

For questions, suggestions, or issues, feel free to open an issue or contact me.
