#!/usr/bin/env bash
#05_annotate
set -euo pipefail

usage() {
    echo "Usage: $0 -i DBSCAN_TSV -o OUTPUT_DIR -b BUILDVER"
    exit 1
}

while getopts "i:o:b:" opt; do
    case "$opt" in
        i) DBSCAN_TSV="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        b) BUILDVER="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "${DBSCAN_TSV:-}" || -z "${OUTPUT_DIR:-}" || -z "${BUILDVER:-}" ]]; then
    echo "Error: missing required arguments"
    usage
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

REORDERED_TSV="${OUTPUT_DIR}/ehdn_DBSCAN_reorder.tsv"
ANNOTATED_TSV="${OUTPUT_DIR}/ehdn_DBSCAN_annotated.tsv"

mkdir -p "$OUTPUT_DIR"

r_script_prepare_annovar="${PIPELINE_DIR}/scripts/prepare_for_annovar.r"
sh_script_run_annovar="${PIPELINE_DIR}/program/ehdn/scripts/annotate_ehdn.sh"

Rscript "$r_script_prepare_annovar" \
    "$DBSCAN_TSV" \
    "$OUTPUT_DIR"

bash "$sh_script_run_annovar" \
    --ehdn-results "$REORDERED_TSV" \
    --ehdn-annotated-results "$ANNOTATED_TSV" \
    --annovar-annotate-variation helper/annotate_variation.pl \
    --annovar-humandb helper/humandb \
    --annovar-buildver "$BUILDVER"
