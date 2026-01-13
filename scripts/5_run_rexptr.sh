#!/usr/bin/env bash
#! 5_run_REXPTR_annotation.sh
set -euo pipefail

usage() {
    echo "Usage: $0 -i INPUT_TSV -o OUTPUT_DIR -l LIFTOVER -c CHAIN -r REXPTR_SCRIPT"
    exit 1
}

while getopts "i:o:l:c:r:" opt; do
    case "$opt" in
        i) INPUT_TSV="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        l) LIFTOVER_PATH="$OPTARG" ;;
        c) LIFTOVER_CHAIN="$OPTARG" ;;
        r) REXPTR_SCRIPT="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "${INPUT_TSV:-}" || -z "${OUTPUT_DIR:-}" || -z "${LIFTOVER_PATH:-}" || -z "${LIFTOVER_CHAIN:-}" || -z "${REXPTR_SCRIPT:-}" ]]; then
    echo "ERROR: Missing required arguments"
    usage
fi


# ============================================================================
# SETUP
# ============================================================================

echo "==========================================="
echo "RExPRT Annotation Pipeline - Step 5"
echo "==========================================="
echo ""
echo "Input file: $INPUT_TSV"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Check if input file exists
if [ ! -f "$INPUT_TSV" ]; then
    echo "ERROR: Input file not found: $INPUT_TSV"
    exit 1
fi

# Check if tools exist
if [ ! -f "$LIFTOVER_PATH" ]; then
    echo "ERROR: liftOver not found at: $LIFTOVER_PATH"
    exit 1
fi

if [ ! -f "$REXPTR_SCRIPT" ]; then
    echo "ERROR: RExPRT script not found at: $REXPTR_SCRIPT"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"


# Get base name for output files
BASE_NAME=$(basename "$INPUT_TSV" .tsv)

# Define intermediate file paths
BED_FILE="$OUTPUT_DIR/${BASE_NAME}.bed"
LIFTED_FILE="$OUTPUT_DIR/${BASE_NAME}_hg19.bed"
UNMAPPED_FILE="$OUTPUT_DIR/${BASE_NAME}_unmapped.bed"
REX_INPUT="$OUTPUT_DIR/${BASE_NAME}_rex_input.txt"
REX_OUTPUT="$OUTPUT_DIR/${BASE_NAME}_rex_input_TRsAnnotated_RExPRTscores.txt"

# ============================================================================
# STEP 1: PREPARE BED FILE
# ============================================================================

echo "=== Step 1: Preparing BED file ==="

Rscript -e "
library(dplyr)
library(readr)

# Read input
df <- read.delim('$INPUT_TSV', header = TRUE)
cat('Loaded', nrow(df), 'TRs\n')


# Create BED format
bed <- df %>%
  select(contig, start, end, motif) %>%
  mutate(
    contig = gsub('^chr', '', contig),
    contig = paste0('chr', contig),
    ID = paste(motif, contig, start, end, sep = '#')
  ) %>%
  select(contig, start, end, motif, ID)

# Write BED file
write.table(bed, '$BED_FILE', 
            quote = FALSE, sep = '\t', 
            col.names = FALSE, row.names = FALSE)

cat('Created BED file:', '$BED_FILE', '\n')
cat('Total regions:', nrow(bed), '\n')
"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create BED file"
    exit 1
fi

# ============================================================================
# STEP 2: LIFTOVER TO HG19
# ============================================================================

echo ""
echo "=== Step 2: Running liftOver (hg38 -> hg19) ==="

$LIFTOVER_PATH "$BED_FILE" "$LIFTOVER_CHAIN" "$LIFTED_FILE" "$UNMAPPED_FILE"

if [ $? -ne 0 ]; then
    echo "ERROR: liftOver failed"
    exit 1
fi

# Check results
LIFTED_COUNT=$(wc -l < "$LIFTED_FILE")
UNMAPPED_COUNT=$(wc -l < "$UNMAPPED_FILE")

echo "Successfully lifted: $LIFTED_COUNT regions"
echo "Unmapped: $UNMAPPED_COUNT regions"

if [ "$LIFTED_COUNT" -eq 0 ]; then
    echo "ERROR: No regions successfully lifted to hg19"
    exit 1
fi

# ============================================================================
# STEP 3: FORMAT FOR REXPTR
# ============================================================================

echo ""
echo "=== Step 3: Formatting for RExPRT ==="

Rscript -e "
library(dplyr)

# Read lifted BED file
lifted <- read.delim('$LIFTED_FILE', header = FALSE)
colnames(lifted) <- c('chr', 'start', 'end', 'motif', 'ID')

# Add sample IDs for RExPRT
lifted\$sampleID <- paste0('TR', seq_len(nrow(lifted)))

# Write RExPRT input
write.table(lifted, '$REX_INPUT',
            sep = '\t', quote = FALSE, row.names = FALSE)

cat('Created RExPRT input:', '$REX_INPUT', '\n')
cat('Total regions:', nrow(lifted), '\n')
"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to format RExPRT input"
    exit 1
fi

# ============================================================================
# STEP 4: RUN REXPTR
# ============================================================================

echo ""
echo "=== Step 4: Running RExPRT (this may take a while) ==="
echo "Started at: $(date)"

$REXPTR_SCRIPT "$REX_INPUT"

if [ $? -ne 0 ]; then
    echo "ERROR: RExPRT failed"
    exit 1
fi

echo "Completed at: $(date)"

# Check if output was created
if [ ! -f "$REX_OUTPUT" ]; then
    echo "ERROR: RExPRT output not found: $REX_OUTPUT"
    exit 1
fi

echo "RExPRT annotation complete!"
echo "Output file: $REX_OUTPUT"

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "==========================================="
echo "Step 5 Complete - RExPRT Raw Annotation"
echo "==========================================="
echo ""
echo "Files created:"
echo "  - BED file: $BED_FILE"
echo "  - Lifted (hg19): $LIFTED_FILE"
echo "  - Unmapped: $UNMAPPED_FILE"
echo "  - RExPRT input: $REX_INPUT"
echo "  - RExPRT output: $REX_OUTPUT"
echo ""
echo "Summary:"
echo "  - Input TRs: $(wc -l < "$INPUT_TSV" | tr -d ' ') (before filtering)"
echo "  - Lifted to hg19: $LIFTED_COUNT"
echo "  - Unmapped: $UNMAPPED_COUNT"
echo ""
echo "Next step: Run 5b_clean_REXPTR_results.R to merge with original data"
echo ""