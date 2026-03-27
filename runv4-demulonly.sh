#!/bin/bash
# RNA-seq Analysis Pipeline v3.1 with UMI processing

# Parameter parsing
usage() {
    echo "Usage: $0 {h|m|r} [options]"
    echo "Options:"
    echo "  --r1 FILE         Read1 input (default: input/r1.fq.gz)"
    echo "  --r2 FILE         Read2 input (default: input/r2.fq.gz)"
    echo "  --sample FILE     Sample-barcode mapping (default: input/sample.csv)"
    echo "  --metadata FILE   Metadata file (default: input/metadata.csv)"
    exit 1
}

# Defaults
species=$1
shift
r1="input/r1.fq.gz"
r2="input/r2.fq.gz"
sample="input/sample.csv"
metadata="input/metadata.csv"

# Parse options
while [[ $# -gt 0 ]]; do
    case "$1" in
        --r1) r1="$2"; shift 2 ;;
        --r2) r2="$2"; shift 2 ;;
        --sample) sample="$2"; shift 2 ;;
        --metadata) metadata="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# Species index
case $species in
    h) index="genome/human/index" ;;
    m) index="genome/mouse/index" ;;
    r) index="genome/rat/index" ;;
    *) echo "Invalid species: $species"; exit 1 ;;
esac

case $species in
    h) gtf="genome/human/gtf/h.gtf.gz" ;;
    m) gtf="genome/mouse/gtf/m.gtf.gz" ;;
    r) gtf="genome/rat/gtf/r.gtf.gz" ;;
    *) echo "Invalid species: $species"; exit 1 ;;
esac

# 1. Demultiplex
echo "Demultiplexing..."
mkdir -p output/demuxed
idemux \
    --r1 "$r1" \
    --r2 "$r2" \
    --sample-sheet "$sample" \
    --out output/demuxed \
    --i1-start 11
