
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash
#PBS -N diversity
#PBS -j oe
#PBS -o DIVERSITY_LOG

#
# PBS Job Script to Compute Alpha Diversity from Bracken Files
#
# This script loops through all `.bracken` files in the specified directory,
# runs a Python script (`alpha_diversity.py`) to compute Shannon diversity,
# and logs results into a central CSV file.
#
# Features:
# - Outputs diversity values per sample to `alpha_diversity_results.csv`.
# - Handles ZeroDivisionError gracefully and logs those cases to `error_log.txt`.
# - Appends each result as: `sample_id,diversity_value`
#
# Input:
# - Folder of `.bracken` files: `/beegfs/work/tu_zxotj39/final_analysis/diversity/output/`
#
# Output:
# - `alpha_diversity_results.csv`: Sample IDs and corresponding diversity values.
# - `error_log.txt`: Log of samples that failed (e.g., division by zero).
#
# Requirements:
# - Python script `alpha_diversity.py` must be in the same or accessible directory.
#
# Usage:
# Submit this script as a PBS job to process alpha diversity on HPC clusters.

echo "Running Job:"

# Create or overwrite the CSV file and add header
output_file="/beegfs/work/tu_zxotj39/final_analysis/diversity/alpha_diversity_results.csv"
error_log="/beegfs/work/tu_zxotj39/final_analysis/diversity/error_log.txt"
#echo "filename,result" > "$output_file"
    
# Loop through all .bracken files in the current directory
for file in /beegfs/work/tu_zxotj39/final_analysis/diversity/output/*.bracken; do
        # Get the base name of the file (without the path)
        base_name=$(basename "$file" )
        sample_id=$(basename "$file" .bracken)        

        # Run the Python script and capture the result, including any errors
        result=$(python alpha_diversity.py -f "$file" -a Sh 2>&1)
        
        # Check if the result contains "ZeroDivisionError"
        if echo "$result" | grep -q "ZeroDivisionError"; then
            echo "$base_name,Error: Division by zero" >> "$error_log"
            echo "Skipped $base_name due to division by zero."
            continue
        fi
        
        # Extract only the number from the result
        number=$(echo "$result" | grep -oE '[0-9]+\.[0-9]+')
        
        # Append the base name and number to the CSV file
        echo "$sample_id,$number" >> "$output_file"
done

echo "Results saved to $output_file"
echo "Errors logged to $error_log"
