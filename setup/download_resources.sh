#!/usr/bin/env bash
set -euo pipefail

mkdir -p resources
cd resources

if [ ! -f "hg38.fa" ]; then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip hg38.fa.gz
fi

if [ ! -f "hg38ToHg19.over.chain.gz" ]; then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
fi