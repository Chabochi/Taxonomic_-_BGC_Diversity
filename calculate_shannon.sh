#!/bin/bash

#
# Script to calculate Shannon-Wiener diversity index per sample based on GCF (Gene Cluster Family) assignments.
#
# This script reads a tab-separated taxonomy file (`regions_taxonomy.tsv`) and computes
# the Shannon-Wiener diversity index for each sample (assembly) based on the frequency
# of `bigslice_gcf_id` values (assumed to be in column 14).
#
# Behavior:
# - Counts occurrences of each `gcf_id` per sample.
# - Calculates relative abundance and applies the Shannon formula:
#     H = -∑ p_i * log(p_i)
# - Special case: if `gcf_id == -1`, the calculation assumes a probability of 1/total.
#
# Input:
#   - `regions_taxonomy.tsv`: TSV with sample IDs in column 2 and GCF IDs in column 14.
#
# Output:
#   - `sample_gcf_shannon.tsv`: TSV file with two columns: sample ID and Shannon index.
#
# Example Output Format:
#   SAMPLE_001    2.484921
#   SAMPLE_002    3.157298
#
# Requirements:
#   - `awk` (with support for `log` function, typically GNU awk)
#
# Usage:
#   ./calculate_shannon.sh


input_file="regions_taxonomy.tsv"
output_file="sample_gcf_shannon.tsv"

# Step 1: Count occurrences of bigslice_gcf_id per assembly
awk -F'\t' 'NR>1 {counts[$2][$14]++} 
    END {
        for (sample in counts) {
            total = 0
            for (gcf_id in counts[sample]) {
                total += counts[sample][gcf_id]
            }
            h = 0
            for (gcf_id in counts[sample]) {
                if (gcf_id == -1) {
                    p = 1 / total
                    h -= counts[sample][-1] * (p * log(p))
                } else {
                    p = counts[sample][gcf_id] / total
                    h -= p * log(p)
                }
            }
            printf "%s\t%.6f\n", sample, h
        }
    }' "$input_file" > "$output_file"

echo "Shannon-Wieder diversity calculated and saved to $output_file"
