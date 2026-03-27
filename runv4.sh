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

# 1. Fastqc
mkdir -p output/QC
fastqc -o output/QC -t 8 $r1 $r2

# 2. Demultiplex
echo "Demultiplexing..."
mkdir -p output/demuxed
idemux \
    --r1 "$r1" \
    --r2 "$r2" \
    --sample-sheet "$sample" \
    --out output/demuxed \
    --i1-start 11

# 3. UMI Extraction
echo "Extracting UMIs..."
mkdir -p output/umi_extracted
for r2_file in output/demuxed/*_R2.fastq.gz; do
    base=$(basename "$r2_file" _R2.fastq.gz)
    r1_file="output/demuxed/${base}_R1.fastq.gz"
    umi_tools extract \
        -I $r2_file \
        -S "output/umi_extracted/${base}_R2.fastq.gz" \
        --read2-in $r1_file \
        --bc-pattern "NNNNNNNNNN" \
        --read2-out "output/umi_extracted/${base}_R1.fastq.gz" \
        --log "output/umi_extracted/${base}_umi.log"
done

# 4. Trimming
echo "Trimming..."
mkdir -p output/trimmed
for r1_file in output/umi_extracted/*_R1.fastq.gz; do
    base=$(basename "$r1_file" _R1.fastq.gz)
    cutadapt \
        -a CTGTCTCTTATACACATCT $r1_file|
    cutadapt \
        -a "A{15}" \
        --nextseq-trim=20 \
        -q 20 \
        -m 50 \
        - >  "output/trimmed/${base}_R1.fastq";
done

# 5. Alignment with STAR
echo "Aligning..."
mkdir -p output/aligned
STAR --genomeLoad LoadAndExit --genomeDir $index

for r1_file in output/trimmed/*_R1.fastq; do
    base=$(basename "$r1_file" _R1.fastq)
    outdir="output/aligned/$base"
    mkdir -p "$outdir"
    #mkdir -p STAR_temp  # Generate temp dir
    STAR \
        --runThreadN 40 \
        --genomeLoad LoadAndKeep \
        --genomeDir $index \
        --readFilesIn $r1_file \
        --outTmpDir /home/dingqi/3rnaseq/STARtemp/temp \
        --outFilterType BySJout \
        --outFilterMultimapNmax 200 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.6 \
        --outFilterScoreMinOverLread 0 \
        --outFilterMatchNminOverLread 0 \
        --outFilterMatchNmin 0 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI NM MD \
        --limitIObufferSize 200000000 200000000 \
        --limitOutSJcollapsed 5000000 \
        --limitBAMsortRAM 90000000000 \
        --outFileNamePrefix "$outdir/${base}_";
    rm -r /home/dingqi/3rnaseq/STARtemp/temp;
done

STAR --genomeLoad Remove --genomeDir "$index"

# 6. UMI Deduplication
echo "Deduplicating UMIs..."
for bam_file in output/aligned/*/*_Aligned.sortedByCoord.out.bam; do
    dir=$(dirname "$bam_file")
    base=$(basename "$bam_file" _Aligned.sortedByCoord.out.bam)
    samtools index "$bam_file"

    umi_tools dedup \
        --stdin "$bam_file" \
        --stdout "$dir/${base}_dedup.bam" \
        --method unique \
        --per-cell \

    samtools index "$dir/${base}_dedup.bam"
done

# 7. Quantification
echo "Counting reads..."
for bam_file in output/aligned/*/*_dedup.bam; do
    dir=$(dirname "$bam_file")
    base=$(basename "$bam_file" _dedup.bam)
    htseq-count \
        --format bam \
        "$bam_file" \
        "$gtf" \
        > "$dir/${base}_counts.tab"
done

# 8. Downstream analysis (uncomment when ready)
Rscript countsV3.R "$species" "$metadata"
Rscript 3RNA_analyzer.R "output/analysis/counts/counts.matrix.txt" "$metadata" "$species" "output/analysis"
Rscript upstream_analysis_summary.R "$metadata" "$sample"

#rm -r output/demuxed
#rm -r output/umi_extracted
#rm -r output/trimmed

echo "Pipeline completed!"