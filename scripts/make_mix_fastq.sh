#!/bin/bash
# Generate mixed FASTQ benchmark datasets from two or more sample directories.

set -euo pipefail

BASE_DIR=""
OUT_BASE=""
RATIOS="5,10,15,20,30"
PAIRS="HG002:HG003,HG002:HG004"
SEED=100
SAMPLE_SPECS=""
TOTAL_FRACTION="1.0"

usage() {
    cat <<EOF
Usage: bash scripts/make_mix_fastq.sh --base-dir DIR --out-base DIR --sample NAME:DIR [--sample NAME:DIR ...] [options]

Required:
  --base-dir DIR       Benchmark root. Used only for relative sample dirs.
  --out-base DIR       Output root for mixed FASTQ.
  --sample NAME:DIR    Sample name and FASTQ directory. Can be repeated.

Options:
  --pairs A:B,C:D      Pairs to mix. Default: ${PAIRS}
  --ratios 5,10,15     Percent of sample A in each pair. Default: ${RATIOS}
  --total-fraction F   Extra downsample factor for the whole mix. Default: ${TOTAL_FRACTION}
                       Example: 0.1 means 1/10 total data; 0.01 means 1/100.
  --seed N             seqtk random seed. Default: ${SEED}
  -h, --help           Show this help.

Example:
  bash scripts/make_mix_fastq.sh \\
    --base-dir /PUBLIC/gomics/guofenghua/project/project_01_wespipline/hg002_benmark \\
    --out-base /PUBLIC/gomics/guofenghua/project/project_01_wespipline/hg002_benmark/moni \\
    --sample HG002:HG002_Sample_2A1 \\
    --sample HG003:HG003_Sample_3A1 \\
    --sample HG004:HG004_Sample_4A1
EOF
}

sample_dir() {
    local name="$1"
    local spec spec_name spec_dir
    OLD_IFS="${IFS}"
    IFS=','
    for spec in ${SAMPLE_SPECS}; do
        spec_name="${spec%%:*}"
        spec_dir="${spec#*:}"
        if [ "${spec_name}" = "${name}" ]; then
            if [ "${spec_dir#/}" != "${spec_dir}" ]; then
                echo "${spec_dir}"
            else
                echo "${BASE_DIR}/${spec_dir}"
            fi
            IFS="${OLD_IFS}"
            return 0
        fi
    done
    IFS="${OLD_IFS}"
    return 1
}

concat_sample_reads() {
    local dir="$1"
    local read_tag="$2"
    find "${dir}" -type f \( \
        -name "*_${read_tag}_*.fastq.gz" -o -name "*_${read_tag}.fastq.gz" -o \
        -name "*_${read_tag}_*.fq.gz" -o -name "*_${read_tag}.fq.gz" \
    \) | sort
}

stream_sample_reads() {
    local dir="$1"
    local read_tag="$2"
    local fq
    concat_sample_reads "${dir}" "${read_tag}" | while IFS= read -r fq; do
        [ -n "${fq}" ] || continue
        cat "${fq}"
    done
}

while [ $# -gt 0 ]; do
    case "$1" in
        --base-dir) BASE_DIR="$2"; shift 2 ;;
        --out-base) OUT_BASE="$2"; shift 2 ;;
        --sample)
            if [ -z "${SAMPLE_SPECS}" ]; then
                SAMPLE_SPECS="$2"
            else
                SAMPLE_SPECS="${SAMPLE_SPECS},$2"
            fi
            shift 2
            ;;
        --pairs) PAIRS="$2"; shift 2 ;;
        --ratios) RATIOS="$2"; shift 2 ;;
        --total-fraction) TOTAL_FRACTION="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if [ -z "${BASE_DIR}" ] || [ -z "${OUT_BASE}" ] || [ -z "${SAMPLE_SPECS}" ]; then
    usage >&2
    exit 1
fi

command -v seqtk >/dev/null 2>&1 || { echo "seqtk not found in PATH" >&2; exit 1; }
command -v gzip >/dev/null 2>&1 || { echo "gzip not found in PATH" >&2; exit 1; }
awk -v f="${TOTAL_FRACTION}" 'BEGIN{exit !(f > 0 && f <= 1)}' || {
    echo "--total-fraction must be > 0 and <= 1" >&2
    exit 1
}

mkdir -p "${OUT_BASE}"
MANIFEST="${OUT_BASE}/mixed_fastq_manifest.tsv"
echo -e "sample_id\tpair\tratio_a_percent\ttotal_fraction\tactual_fraction_a\tactual_fraction_b\tfastq_r1\tfastq_r2" > "${MANIFEST}"

OLD_IFS="${IFS}"
IFS=','
for pair in ${PAIRS}; do
    sample_a="${pair%%:*}"
    sample_b="${pair#*:}"
    path_a="$(sample_dir "${sample_a}")"
    path_b="$(sample_dir "${sample_b}")"

    if [ ! -d "${path_a}" ] || [ ! -d "${path_b}" ]; then
        echo "Missing sample directory for pair ${pair}: ${path_a} / ${path_b}" >&2
        exit 1
    fi

    for ratio in ${RATIOS}; do
        frac_a=$(awk -v r="${ratio}" -v t="${TOTAL_FRACTION}" 'BEGIN{printf "%.8f", (r/100)*t}')
        frac_b=$(awk -v r="${ratio}" -v t="${TOTAL_FRACTION}" 'BEGIN{printf "%.8f", (1-r/100)*t}')
        sample_id="mix_${sample_a}_${sample_b}_${ratio}pct"
        target_dir="${OUT_BASE}/${sample_id}"
        mkdir -p "${target_dir}"

        echo "========================================================="
        echo "Generating ${sample_id}: ${sample_a} ${ratio}% + ${sample_b} $((100 - ratio))%"
        echo "Total fraction: ${TOTAL_FRACTION}; actual sample fractions: ${sample_a}=${frac_a}, ${sample_b}=${frac_b}"

        r1_a_files="$(concat_sample_reads "${path_a}" "R1")"
        r1_b_files="$(concat_sample_reads "${path_b}" "R1")"
        r2_a_files="$(concat_sample_reads "${path_a}" "R2")"
        r2_b_files="$(concat_sample_reads "${path_b}" "R2")"

        if [ -z "${r1_a_files}" ] || [ -z "${r1_b_files}" ] || [ -z "${r2_a_files}" ] || [ -z "${r2_b_files}" ]; then
            echo "Missing R1/R2 FASTQ for ${sample_id}" >&2
            exit 1
        fi

        stream_sample_reads "${path_a}" "R1" | seqtk sample -s"${SEED}" - "${frac_a}" > "${target_dir}/tmp_A_R1.fq"
        stream_sample_reads "${path_b}" "R1" | seqtk sample -s"${SEED}" - "${frac_b}" > "${target_dir}/tmp_B_R1.fq"
        cat "${target_dir}/tmp_A_R1.fq" "${target_dir}/tmp_B_R1.fq" | gzip > "${target_dir}/mix_R1.fastq.gz"

        stream_sample_reads "${path_a}" "R2" | seqtk sample -s"${SEED}" - "${frac_a}" > "${target_dir}/tmp_A_R2.fq"
        stream_sample_reads "${path_b}" "R2" | seqtk sample -s"${SEED}" - "${frac_b}" > "${target_dir}/tmp_B_R2.fq"
        cat "${target_dir}/tmp_A_R2.fq" "${target_dir}/tmp_B_R2.fq" | gzip > "${target_dir}/mix_R2.fastq.gz"

        rm -f "${target_dir}/tmp_A_R1.fq" "${target_dir}/tmp_B_R1.fq" "${target_dir}/tmp_A_R2.fq" "${target_dir}/tmp_B_R2.fq"
        echo -e "${sample_id}\t${sample_a}_${sample_b}\t${ratio}\t${TOTAL_FRACTION}\t${frac_a}\t${frac_b}\t${target_dir}/mix_R1.fastq.gz\t${target_dir}/mix_R2.fastq.gz" >> "${MANIFEST}"
        echo "Done: ${target_dir}"
    done
done
IFS="${OLD_IFS}"

echo "All mixed FASTQ generated."
echo "Manifest: ${MANIFEST}"
