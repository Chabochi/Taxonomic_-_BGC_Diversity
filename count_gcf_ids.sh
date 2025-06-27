#!/bin/bash
#
# Script to count distinct GCF (Gene Cluste Families) IDs per sample from a tab-separated taxonomy file.
#
# This script processes a `.tsv` file (`regions_taxonomy.tsv`) containing GCF identifiers
# and outputs the number of distinct GCF IDs associated with each sample.
#
# Logic:
#   - If column 14 is `-1`, the GCF ID is considered an "associated" (unclassified) GCF and counted directly.
#   - Otherwise, each unique non-`-1` GCF ID is counted once per sample.
#   - The final output is a TSV file with two columns: sample ID and distinct GCF count.
#
# Input:
#   - `regions_taxonomy.tsv`: TSV file where each row includes a sample ID in column 2 and a GCF ID in column 14.
#
# Output:
#   - `sample_gcf_counts.tsv`: TSV file with counts of distinct GCF IDs per sample.
#
# Example Output Format:
#   SAMPLE_001    37
#   SAMPLE_002    12
#
# Usage:
#   ./count_gcf_ids.sh
#
# Requirements:
#   - `awk`



input_file="regions_taxonomy.tsv"
output_file="sample_gcf_counts.tsv"

awk -F'\t' 'NR>1 {
    if ($14 == -1) {
        assoc_count[$2]++ # Count each occurrence of -1 as distinct
    } else {
        unique[$2][$14] = 1 # Track unique non -1 values
    }
} 
END {
    for (sample in unique) {
        distinct_count = assoc_count[sample]
        for (gcf_id in unique[sample]) {
            distinct_count++
        }

        printf "%s\t%d\n", sample, distinct_count
    }
}' "$input_file" > "$output_file"

echo "Distinct GCF IDs counted and saved to $output_file"
