#!/bin/bash
set -euo pipefail

# =============================================================================
# ExpansionHunterDenovo Pipeline
# =============================================================================
# DBSCAN step 
# usage: Combine_RUN_DBSCAN.sh -m manifest_file
# =============================================================================

usage() {
    echo "Usage: $0 -m manifest_file -o output_dir -t threads"
    exit 1
}

while getopts "m:o:t:" opt; do
    case "$opt" in
        m) manifest_file=$OPTARG ;;
        o) output_dir=$OPTARG ;;
        t) threads=$OPTARG ;;
        *) usage ;;
    esac
done

if [ -z "$manifest_file" ] || [ -z "$output_dir" ]; then
    echo "Error: Manifest file and output directory are required"
    usage
fi

mkdir -p "$output_dir"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

sample_list="${output_dir}/sample_list.txt"
awk -F'\t' '{print $1}' "$manifest_file" > "$sample_list"

if [ ! -s "$sample_list" ]; then
    echo "Error: Derived sample list is empty"
    exit 1
fi

# Define the output file names with fixed names
combined_counts="${output_dir}/EHDn_DBSCAN.combinedCounts.json"
output_regions="${output_dir}/EHDn_DBSCAN.combinedCounts.bed"

# Paths to the helper scripts
combine_counts_script="${PIPELINE_DIR}/helper/Fazal2020Scripts/EHDn-v0.6.2_HelperScripts/combine_counts.py"
compare_anchored_irrs_script="${PIPELINE_DIR}/helper/Fazal2020Scripts/EHDn-v0.6.2_HelperScripts/compare_anchored_irrs.py"

[ -f "$combine_counts_script" ] || { echo "Missing $combine_counts_script"; exit 1; }
[ -f "$compare_anchored_irrs_script" ] || { echo "Missing $compare_anchored_irrs_script"; exit 1; }

echo "Running combine_counts.py..."
python3 "$combine_counts_script" --manifest "$manifest_file" --combinedCounts "$combined_counts"

echo "Running compare_anchored_irrs.py..."
python3 "$compare_anchored_irrs_script" \
    --manifest "$manifest_file" \
    --inputCounts "$combined_counts" \
    --outputRegions "$output_regions" \
    --minCount 2

r_script="${PIPELINE_DIR}/scripts/dbscan_only.R"
[ -f "$r_script" ] || { echo "Missing $r_script"; exit 1; }

echo "Running dbscan_only.R..."
Rscript "$r_script" "$output_dir" "$output_regions" "$sample_list" "$threads"

echo "Pipeline completed successfully"
echo "Results in: $output_dir"
