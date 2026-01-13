#!/bin/bash
set -euo pipefail

# =============================================================================
# ExpansionHunterDenovo Pipeline
# =============================================================================

# Usage: ./1_ehdn.sh -r /path/to/reference.fa -m /path/to/manifest.txt -o /path/to/output -e /path/to/ehdn
# =============================================================================

# Default values
REF=""
MANIFEST_FILE=""
BASE_OUTPUT=""
EHDN_BIN=""

# Parse command-line arguments
while getopts "n:r:m:o:e:" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        m) MANIFEST_FILE="$OPTARG" ;;
        o) BASE_OUTPUT="$OPTARG" ;;
        e) EHDN_BIN="$OPTARG" ;;
        *) echo "Usage: $0 -r REFERENCE -m MANIFEST -o OUTPUT_DIR -e EHDN_BIN"; exit 1 ;;
    esac
done

# Check mandatory arguments
if [[  -z "$REF" || -z "$MANIFEST_FILE" || -z "$BASE_OUTPUT" || -z "$EHDN_BIN" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 -r REFERENCE -m MANIFEST -o OUTPUT_DIR -e EHDN_BIN"
    exit 1
fi

# Ensure output directory exists
mkdir -p "$BASE_OUTPUT"

# Initialize downstream manifest
MANIFEST_OUT="${BASE_OUTPUT}/ehdn_generated_manifest.tsv"
echo -e "SampleID\tStatus\tFilePath" > "$MANIFEST_OUT"

# Initialize combined pipeline log
LOG_FILE="${BASE_OUTPUT}/ehdn_pipeline.log"
echo "Starting EHDN pipeline at $(date)" > "$LOG_FILE"

# Iterate through each sample in the manifest, skipping the header
tail -n +2 "$MANIFEST_FILE" | while IFS=$'\t' read -r sample status bam_path; do
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

        # Append to downstream manifest
        echo -e "${sample}\t${status}\t${OUTPUT_JSON}" >> "$MANIFEST_OUT"
    else
        echo "[$(date)] Warning: BAM not found for $sample at $bam_path" | tee -a "$LOG_FILE"
    fi
done

echo "EHDN pipeline finished at $(date)" | tee -a "$LOG_FILE"




