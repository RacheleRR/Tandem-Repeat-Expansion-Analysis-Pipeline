#!/bin/bash
set -euo pipefail


# =============================================================================
# ExpansionHunterDenovo Pipeline with GNU parallel support
# =============================================================================
# Usage:
# ./1_ehdn.sh -r /path/to/reference.fa -m /path/to/manifest.txt \
#             -o /path/to/output -e /path/to/ehdn \
#             [-p] [--max-parallel N]
# =============================================================================

# Default values
REF=""
MANIFEST_FILE=""
BASE_OUTPUT=""
EHDN_BIN=""
PARALLELIZE=false
MAX_PARALLEL=4

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -r) REF="$2"; shift 2 ;;
        -m) MANIFEST_FILE="$2"; shift 2 ;;
        -o) BASE_OUTPUT="$2"; shift 2 ;;
        -e) EHDN_BIN="$2"; shift 2 ;;
        -p) PARALLELIZE=true; shift ;;
        --max-parallel) MAX_PARALLEL="$2"; shift 2 ;;
        *) echo "Usage: $0 -r REFERENCE -m MANIFEST -o OUTPUT_DIR -e EHDN_BIN [-p] [--max-parallel N]"; exit 1 ;;
    esac
done

# Check mandatory arguments
if [[ -z "$REF" || -z "$MANIFEST_FILE" || -z "$BASE_OUTPUT" || -z "$EHDN_BIN" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 -r REFERENCE -m MANIFEST -o OUTPUT_DIR -e EHDN_BIN [-p] [--max-parallel N]"
    exit 1
fi

# Ensure output directory exists
mkdir -p "$BASE_OUTPUT"

# Initialize downstream manifest
MANIFEST_OUT="${BASE_OUTPUT}/ehdn_generated_manifest.tsv"
echo -e "SampleID\tStatus\tFilePath" > "$MANIFEST_OUT"

# Initialize pipeline log
LOG_FILE="${BASE_OUTPUT}/ehdn_pipeline.log"
echo "Starting EHDN pipeline at $(date)" > "$LOG_FILE"

# Function to process a single sample
process_sample() {
    local sample="$1"
    local status="$2"
    local bam_path="$3"

    if [[ -f "$bam_path" ]]; then
        OUTPUT_JSON="${BASE_OUTPUT}/${sample}.str_profile.json"

        echo "[$(date)] Running ExpansionHunterDenovo on $sample ($bam_path)" | tee -a "$LOG_FILE"
        "$EHDN_BIN" profile \
            --read "$bam_path" \
            --reference "$REF" \
            --min-anchor-mapq 50 \
            --max-irr-mapq 40 \
            --output-prefix "${BASE_OUTPUT}/${sample}" \
            >> "$LOG_FILE" 2>&1

        echo "[$(date)] Completed sample: $sample" | tee -a "$LOG_FILE"

        # Append to downstream manifest safely using flock
        (
            flock -x 200
            echo -e "${sample}\t${status}\t${OUTPUT_JSON}" >> "$MANIFEST_OUT"
        ) 200>"${BASE_OUTPUT}/manifest.lock"
    else
        echo "[$(date)] Warning: BAM not found for $sample at $bam_path" | tee -a "$LOG_FILE"
    fi
}

export -f process_sample
export EHDN_BIN REF BASE_OUTPUT LOG_FILE MANIFEST_OUT

# Read manifest lines into an array (skip header)
tail -n +2 "$MANIFEST_FILE" > "${BASE_OUTPUT}/samples.tmp"

if [[ "$PARALLELIZE" = true ]]; then
    echo "Running in parallel mode with max $MAX_PARALLEL jobs" | tee -a "$LOG_FILE"

    # Run samples in parallel using GNU parallel
    parallel -j "$MAX_PARALLEL" --colsep '\t' process_sample {1} {2} {3} :::: "${BASE_OUTPUT}/samples.tmp"
else
    # Serial execution
    while IFS=$'\t' read -r sample status bam_path; do
        process_sample "$sample" "$status" "$bam_path"
    done < "${BASE_OUTPUT}/samples.tmp"
fi

rm -f "${BASE_OUTPUT}/samples.tmp"

echo "EHDN pipeline finished at $(date)" | tee -a "$LOG_FILE"








