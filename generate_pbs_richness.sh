#!/bin/bash

#
# Script to automatically generate PBS job scripts for computing taxonomic richness
# from `.report` files located within subdirectories of a parent folder.
#
# For each subdirectory containing a valid `.report` file (excluding *_bracken.report),
# the script:
#   - Extracts the base name of the report.
#   - Creates a corresponding PBS job script in the `scripts/` directory.
#   - The job script calculates family, genus, and species richness using `awk`.
#   - Outputs the result as a TSV file in the target `correct_richness/` folder.
#
# Input:
#   - Parent directory with subfolders, each potentially containing a `.report` file.
#
# Output:
#   - One PBS job script per report, saved in `scripts/`.
#
# Example Usage:
#   ./generate_pbs_richness.sh
#
# Notes:
#   - Update `parent_dir`, output paths, and `cd` command as needed to fit your directory structure.
#   - PBS directives assume a Torque/Portable Batch System environment.

# Directory containing the subfolders
parent_dir="/beegfs/work/tu_zxotj39/bgc/output"

# Create output directory for scripts if it doesn't exist
mkdir -p scripts

# Loop through each subfolder
for folder in "$parent_dir"/*; do
    # Check if it's a directory
    if [ -d "$folder" ]; then
        # Find the .report file in the current folder
        report_file=$(find "$folder" -maxdepth 1 -name "*.report" ! -name "*_bracken.report" | head -n 1)

        # Continue if no .report file is found
        if [ -z "$report_file" ]; then
            echo "No .report file found in $folder"
            continue
        fi

        # Extract the base name without extension
        base_name=$(basename "$report_file" .report)
        folder_name=$(basename "$folder")

        # Generate the PBS script
        cat <<EOF > scripts/${base_name}_script.sh
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash
#PBS -N ${base_name}
#PBS -j oe
#PBS -o ${base_name}_LOG

echo "Running Job:"

# Change to working directory
cd /beegfs/work/tu_zxotj39/final_analysis/richness

# Input report file
REPORT_FILE="$report_file"

# Output file path
OUTPUT_FILE="/beegfs/work/tu_zxotj39/final_analysis/richness/correct_richness/${folder_name}_richness.tsv"

# Check for required input
if [[ ! -f "\$REPORT_FILE" ]]; then
    echo "Error: Report file not found!"
    exit 1
fi

# Step 1: Calculate richness based on rank column
awk -F'\t' -v assembly="${base_name}" '
    NR > 1 {
        if (\$4 == "F") family[assembly]++
        if (\$4 == "G") genus[assembly]++
        if (\$4 == "S") species[assembly]++
    }
    END {
        print "assembly\tfamily_richness\tgenus_richness\tspecies_richness"
        f_count = family[assembly] ? family[assembly] : 0
        g_count = genus[assembly] ? genus[assembly] : 0
        s_count = species[assembly] ? species[assembly] : 0
        print assembly "\t" f_count "\t" g_count "\t" s_count
    }
' "\$REPORT_FILE" > "\$OUTPUT_FILE"

echo "Results saved to \$OUTPUT_FILE"
echo "Job finished"
EOF
    fi
done
