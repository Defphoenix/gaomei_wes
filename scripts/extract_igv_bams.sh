#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash extract_igv_bams.sh \
    --tumor-bam TUMOR.bam \
    --normal-bam NORMAL.bam \
    --loci-file loci.tsv \
    --output-dir igv_slices \
    [--flank 300] [--threads 2]

loci.tsv must contain: chrom, pos, ref, alt, gene
The first header line is optional. Coordinates are VCF-style 1-based positions.

The script only reads the input BAMs. It writes two small coordinate-sorted BAMs,
their BAI indexes, a BED file, and an IGV locus list in the output directory.
EOF
}

die() {
    printf '[ERROR] %s\n' "$*" >&2
    exit 1
}

tumor_bam=""
normal_bam=""
loci_file=""
output_dir=""
flank=300
threads=2

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tumor-bam) tumor_bam=${2:?}; shift 2 ;;
        --normal-bam) normal_bam=${2:?}; shift 2 ;;
        --loci-file) loci_file=${2:?}; shift 2 ;;
        --output-dir) output_dir=${2:?}; shift 2 ;;
        --flank) flank=${2:?}; shift 2 ;;
        --threads) threads=${2:?}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown option: $1" ;;
    esac
done

[[ -n "$tumor_bam" ]] || die "--tumor-bam is required"
[[ -n "$normal_bam" ]] || die "--normal-bam is required"
[[ -n "$loci_file" ]] || die "--loci-file is required"
[[ -n "$output_dir" ]] || die "--output-dir is required"
[[ -r "$tumor_bam" ]] || die "Tumor BAM is not readable: $tumor_bam"
[[ -r "$normal_bam" ]] || die "Normal BAM is not readable: $normal_bam"
[[ -r "$loci_file" ]] || die "Loci file is not readable: $loci_file"
[[ "$flank" =~ ^[0-9]+$ ]] || die "--flank must be a non-negative integer"
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"
command -v samtools >/dev/null 2>&1 || die "samtools is not available in PATH"

samtools_version=$(samtools --version 2>/dev/null | awk 'NR == 1 {print $2}' || true)
[[ "$samtools_version" =~ ^([0-9]+)\.([0-9]+) ]] \
    || die "Cannot detect samtools version. Activate the main WES environment (samtools >= 1.10)."
samtools_major=${BASH_REMATCH[1]}
samtools_minor=${BASH_REMATCH[2]}
if (( samtools_major < 1 || (samtools_major == 1 && samtools_minor < 10) )); then
    die "samtools >= 1.10 is required; detected $samtools_version. Activate the main WES environment."
fi
printf '[INFO] samtools version: %s\n' "$samtools_version"

has_bam_index() {
    local bam=$1
    [[ -r "${bam}.bai" || -r "${bam%.bam}.bai" || -r "${bam}.csi" || -r "${bam%.bam}.csi" ]]
}

has_bam_index "$tumor_bam" || die "Tumor BAM index is missing. Run: samtools index -@ $threads '$tumor_bam'"
has_bam_index "$normal_bam" || die "Normal BAM index is missing. Run: samtools index -@ $threads '$normal_bam'"
samtools quickcheck "$tumor_bam" || die "Tumor BAM failed samtools quickcheck"
samtools quickcheck "$normal_bam" || die "Normal BAM failed samtools quickcheck"

mkdir -p "$output_dir"
tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/wes-igv.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

samtools idxstats "$tumor_bam" | awk '$1 != "*" {print $1}' > "$tmp_dir/tumor.contigs"
samtools idxstats "$normal_bam" | awk '$1 != "*" {print $1}' > "$tmp_dir/normal.contigs"

resolve_contig() {
    local requested=$1
    local candidate
    for candidate in "$requested" "${requested#chr}" "chr${requested#chr}"; do
        if grep -Fqx "$candidate" "$tmp_dir/tumor.contigs" && grep -Fqx "$candidate" "$tmp_dir/normal.contigs"; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

regions_bed="$output_dir/igv_regions.bed"
locus_list="$output_dir/igv_loci.txt"
: > "$regions_bed"
: > "$locus_list"

while IFS=$'\t' read -r chrom pos ref alt gene _rest; do
    [[ -n "${chrom:-}" ]] || continue
    [[ "${chrom:0:1}" == "#" ]] && continue
    [[ "${pos:-}" =~ ^[0-9]+$ ]] || continue
    resolved=$(resolve_contig "$chrom") || die "Contig '$chrom' is not present in both BAM headers"
    start=$((pos - flank - 1))
    (( start < 0 )) && start=0
    end=$((pos + flank))
    label=${gene:-variant}_${resolved}_${pos}_${ref:-N}_${alt:-N}
    label=${label//[^A-Za-z0-9_.-]/_}
    printf '%s\t%s\t%s\t%s\n' "$resolved" "$start" "$end" "$label" >> "$regions_bed"
    printf '%s:%s\n' "$resolved" "$pos" >> "$locus_list"
done < "$loci_file"

[[ -s "$regions_bed" ]] || die "No valid loci were found in $loci_file"
LC_ALL=C sort -k1,1 -k2,2n -u "$regions_bed" -o "$regions_bed"

extract_one() {
    local role=$1
    local bam=$2
    local output_bam="$output_dir/${role}.igv.bam"
    printf '[INFO] Extracting %s reads -> %s\n' "$role" "$output_bam"
    samtools view -@ "$threads" -bh -M -L "$regions_bed" "$bam" \
        | samtools sort -@ "$threads" -o "$output_bam" -
    samtools index -@ "$threads" "$output_bam"
    samtools quickcheck "$output_bam" || die "Output BAM failed quickcheck: $output_bam"
}

extract_one tumor "$tumor_bam"
extract_one normal "$normal_bam"

printf '%s\n' \
    "IGV files are ready:" \
    "  $output_dir/tumor.igv.bam" \
    "  $output_dir/normal.igv.bam" \
    "  $regions_bed" \
    "  $locus_list" \
    "Load both BAM files in IGV, select the same reference genome, then paste a locus from igv_loci.txt into the search box." \
    > "$output_dir/README_IGV.txt"

cat "$output_dir/README_IGV.txt"
