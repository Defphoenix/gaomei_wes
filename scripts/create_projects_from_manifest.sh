#!/bin/bash
# Generate runnable WES projects from mixed_fastq_manifest.tsv.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

MANIFEST=""
OUT_BASE=""
COPY_MODE="link"
REFERENCE_DIR=""
REFERENCE_GENOME=""
INTERVAL_BED=""
ENV_ROOT=""
CONDA_BASE=""

usage() {
    cat <<EOF
Usage: bash scripts/create_projects_from_manifest.sh --manifest FILE --out-base DIR [options]

Required:
  --manifest FILE       mixed_fastq_manifest.tsv from make_mix_fastq.sh
  --out-base DIR        Output root for generated runnable projects.

Options:
  --copy-mode copy|link Default: ${COPY_MODE}
  --reference-dir DIR
  --reference-genome FA
  --interval-bed BED
  --env-root DIR
  --conda-base DIR
  -h, --help
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --manifest) MANIFEST="$2"; shift 2 ;;
        --out-base) OUT_BASE="$2"; shift 2 ;;
        --copy-mode) COPY_MODE="$2"; shift 2 ;;
        --reference-dir) REFERENCE_DIR="$2"; shift 2 ;;
        --reference-genome) REFERENCE_GENOME="$2"; shift 2 ;;
        --interval-bed) INTERVAL_BED="$2"; shift 2 ;;
        --env-root) ENV_ROOT="$2"; shift 2 ;;
        --conda-base) CONDA_BASE="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if [ -z "${MANIFEST}" ] || [ -z "${OUT_BASE}" ]; then
    usage >&2
    exit 1
fi

if [ ! -f "${MANIFEST}" ]; then
    echo "Manifest not found: ${MANIFEST}" >&2
    exit 1
fi

mkdir -p "${OUT_BASE}"

header="$(head -1 "${MANIFEST}")"
sample_col=$(awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="sample_id") print i}' <<< "${header}")
pair_col=$(awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="pair") print i}' <<< "${header}")
ratio_col=$(awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="ratio_a_percent") print i}' <<< "${header}")
r1_col=$(awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="fastq_r1") print i}' <<< "${header}")
r2_col=$(awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="fastq_r2") print i}' <<< "${header}")

if [ -z "${sample_col}" ] || [ -z "${r1_col}" ] || [ -z "${r2_col}" ]; then
    echo "Manifest must contain sample_id, fastq_r1 and fastq_r2 columns" >&2
    exit 1
fi

tail -n +2 "${MANIFEST}" | while IFS="$(printf '\t')" read -r line; do
    [ -n "${line}" ] || continue
    sample_id=$(awk -F'\t' -v c="${sample_col}" '{print $c}' <<< "${line}")
    pair=$(awk -F'\t' -v c="${pair_col:-0}" '{if(c>0) print $c; else print "NA"}' <<< "${line}")
    ratio=$(awk -F'\t' -v c="${ratio_col:-0}" '{if(c>0) print $c; else print "NA"}' <<< "${line}")
    fastq_r1=$(awk -F'\t' -v c="${r1_col}" '{print $c}' <<< "${line}")
    fastq_r2=$(awk -F'\t' -v c="${r2_col}" '{print $c}' <<< "${line}")

    [ -n "${sample_id}" ] || continue
    fastq_source="$(dirname "${fastq_r1}")"
    project_out="${OUT_BASE}/${sample_id}"

    cmd=(bash "${SCRIPT_DIR}/create_wes_project.sh"
        --fastq-source "${fastq_source}"
        --out-dir "${project_out}"
        --sample-id "${sample_id}"
        --copy-mode "${COPY_MODE}")

    [ -n "${REFERENCE_DIR}" ] && cmd+=(--reference-dir "${REFERENCE_DIR}")
    [ -n "${REFERENCE_GENOME}" ] && cmd+=(--reference-genome "${REFERENCE_GENOME}")
    [ -n "${INTERVAL_BED}" ] && cmd+=(--interval-bed "${INTERVAL_BED}")
    [ -n "${ENV_ROOT}" ] && cmd+=(--env-root "${ENV_ROOT}")
    [ -n "${CONDA_BASE}" ] && cmd+=(--conda-base "${CONDA_BASE}")

    echo "Generating project for ${sample_id} (${pair}, ${ratio}%)"
    "${cmd[@]}"
done

echo "All projects generated under: ${OUT_BASE}"
