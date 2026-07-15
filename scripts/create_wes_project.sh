#!/bin/bash
# Generate a runnable WES project folder from one paired FASTQ source.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

FASTQ_SOURCE=""
MODE="single"
TUMOR_FASTQ_SOURCE=""
NORMAL_FASTQ_SOURCE=""
TUMOR_ID=""
NORMAL_ID=""
OUT_DIR=""
SAMPLE_ID=""
PROJECT_NAME="wes_project"
COPY_MODE="copy"
REFERENCE_DIR="/Users/mac/Documents/wes/reference_data"
REFERENCE_GENOME=""
INTERVAL_BED=""
CONDA_BASE="/Users/mac/anaconda3"
ENV_ROOT=""
MAIN_ENV_PREFIX_OVERRIDE=""
VEP_ENV_PREFIX_OVERRIDE=""
HLA_ENV_PREFIX_OVERRIDE=""
HLA_TYPING_ENV_PREFIX_OVERRIDE=""
CNV_ENV_PREFIX_OVERRIDE=""
SV_ENV_PREFIX_OVERRIDE=""
INCLUDE_TESTDATA=true

usage() {
    cat <<EOF
Usage:
  Single sample:
    bash scripts/create_wes_project.sh --mode single --fastq-source DIR --out-dir DIR [options]

  Tumor-normal pair:
    bash scripts/create_wes_project.sh --mode tumor-normal --tumor-fastq-source DIR --normal-fastq-source DIR --out-dir DIR [options]

Required:
  --fastq-source DIR      Directory containing one paired FASTQ set for single mode.
  --tumor-fastq-source DIR    Tumor FASTQ directory for tumor-normal mode.
  --normal-fastq-source DIR   Normal FASTQ directory for tumor-normal mode.
  --out-dir DIR           Target project directory to generate.

Options:
  --mode single|tumor-normal  Project mode. Default: ${MODE}
  --sample-id ID          Sample ID. Default: inferred from R1 filename.
  --tumor-id ID           Tumor sample ID for tumor-normal mode.
  --normal-id ID          Normal sample ID for tumor-normal mode.
  --project-name NAME     Project name. Default: ${PROJECT_NAME}
  --copy-mode copy|link   Copy FASTQ into target data/ or symlink them. Default: ${COPY_MODE}
  --reference-dir DIR     Reference root. Default: ${REFERENCE_DIR}
  --reference-genome FA   Reference FASTA. Default: reference-dir/hg38/Homo_sapiens_assembly38.fasta
  --interval-bed BED      WES target BED.
  --conda-base DIR        Conda base path. Default: ${CONDA_BASE}
  --env-root DIR          Env root containing big_wes_pipeline_env and wes_vep_env.
  --main-env-prefix DIR   Full path to big_wes_pipeline_env. Overrides --env-root for main tools.
  --vep-env-prefix DIR    Full path to wes_vep_env. Overrides --env-root for VEP.
  --hla-env-prefix DIR    Full path to wes_hla_env. Overrides --env-root for HLA.
  --hla-typing-env-prefix DIR Full path to wes_hla_typing_env. Overrides --env-root.
  --cnv-env-prefix DIR    Full path to wes_cnv_env. Overrides --env-root for CNVkit.
  --sv-env-prefix DIR     Full path to wes_sv_env. Overrides --env-root for Manta.
  --no-testdata           Do not copy bundled testdata/demo FASTQ.
  -h, --help              Show this help.

Examples:
  bash scripts/create_wes_project.sh \\
    --fastq-source /data/sample01 \\
    --out-dir /PUBLIC/project/sample01_run \\
    --sample-id sample01 \\
    --reference-dir /PUBLIC/reference/hg38_bundle \\
    --interval-bed /PUBLIC/reference/capture.bed
EOF
}

abs_path() {
    local path="$1"
    if [ -d "${path}" ]; then
        (cd "${path}" && pwd)
    else
        local dir base
        dir="$(dirname "${path}")"
        base="$(basename "${path}")"
        if [ ! -d "${dir}" ]; then
            echo "Parent directory not found for path: ${path}" >&2
            return 1
        fi
        echo "$(cd "${dir}" && pwd)/${base}"
    fi
}

require_dir() {
    local label="$1"
    local path="$2"
    if [ ! -d "${path}" ]; then
        echo "${label} directory not found: ${path}" >&2
        exit 1
    fi
}

require_file() {
    local label="$1"
    local path="$2"
    if [ ! -f "${path}" ]; then
        echo "${label} file not found: ${path}" >&2
        exit 1
    fi
}

find_r1() {
    local dir="$1"
    find "${dir}" -type f \( \
        -name '*_R1_*.fastq.gz' -o -name '*_R1.fastq.gz' -o -name '*_1.fastq.gz' -o \
        -name '*_R1_*.fq.gz' -o -name '*_R1.fq.gz' -o -name '*_1.fq.gz' \
    \) | sort | head -1
}

guess_r2() {
    local r1="$1"
    local candidate
    for candidate in \
        "${r1/_R1_/_R2_}" \
        "${r1/_R1./_R2.}" \
        "${r1/_R1/_R2}" \
        "${r1/_1.fastq.gz/_2.fastq.gz}" \
        "${r1/_1.fq.gz/_2.fq.gz}"; do
        if [ "${candidate}" != "${r1}" ] && [ -f "${candidate}" ]; then
            echo "${candidate}"
            return 0
        fi
    done
    return 1
}

infer_sample_id() {
    local r1_base="$1"
    r1_base="${r1_base%.fastq.gz}"
    r1_base="${r1_base%.fq.gz}"
    r1_base="${r1_base%%_R1_*}"
    r1_base="${r1_base%%_R1}"
    r1_base="${r1_base%%_1}"
    echo "${r1_base}"
}

copy_pipeline_code() {
    local target_pipeline="$1"
    mkdir -p "${target_pipeline}"

    if command -v rsync >/dev/null 2>&1; then
        local testdata_args=""
        if [ "${INCLUDE_TESTDATA}" != true ]; then
            testdata_args="--exclude testdata/"
        fi
        rsync -a \
            --exclude '.git/' \
            --exclude '.DS_Store' \
            --exclude 'results/' \
            --exclude 'logs/' \
            --exclude 'cache/' \
            --exclude '.runtime/' \
            --exclude 'tmp_*/' \
            ${testdata_args} \
            "${SOURCE_PROJECT_DIR}/" "${target_pipeline}/"
    else
        cp "${SOURCE_PROJECT_DIR}/run_pipeline.sh" "${target_pipeline}/"
        cp "${SOURCE_PROJECT_DIR}/config.sh" "${target_pipeline}/"
        cp "${SOURCE_PROJECT_DIR}/README.md" "${target_pipeline}/"
        [ ! -f "${SOURCE_PROJECT_DIR}/README_EN.md" ] || cp "${SOURCE_PROJECT_DIR}/README_EN.md" "${target_pipeline}/"
        cp "${SOURCE_PROJECT_DIR}"/*.yml "${target_pipeline}/"
        cp -R "${SOURCE_PROJECT_DIR}/scripts" "${target_pipeline}/"
        cp -R "${SOURCE_PROJECT_DIR}/docs" "${target_pipeline}/"
        if [ "${INCLUDE_TESTDATA}" = true ]; then
            cp -R "${SOURCE_PROJECT_DIR}/testdata" "${target_pipeline}/"
        fi
    fi
}

while [ $# -gt 0 ]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --fastq-source) FASTQ_SOURCE="$2"; shift 2 ;;
        --tumor-fastq-source) TUMOR_FASTQ_SOURCE="$2"; shift 2 ;;
        --normal-fastq-source) NORMAL_FASTQ_SOURCE="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --sample-id) SAMPLE_ID="$2"; shift 2 ;;
        --tumor-id) TUMOR_ID="$2"; shift 2 ;;
        --normal-id) NORMAL_ID="$2"; shift 2 ;;
        --project-name) PROJECT_NAME="$2"; shift 2 ;;
        --copy-mode) COPY_MODE="$2"; shift 2 ;;
        --reference-dir) REFERENCE_DIR="$2"; shift 2 ;;
        --reference-genome) REFERENCE_GENOME="$2"; shift 2 ;;
        --interval-bed) INTERVAL_BED="$2"; shift 2 ;;
        --conda-base) CONDA_BASE="$2"; shift 2 ;;
        --env-root) ENV_ROOT="$2"; shift 2 ;;
        --main-env-prefix) MAIN_ENV_PREFIX_OVERRIDE="$2"; shift 2 ;;
        --vep-env-prefix) VEP_ENV_PREFIX_OVERRIDE="$2"; shift 2 ;;
        --hla-env-prefix) HLA_ENV_PREFIX_OVERRIDE="$2"; shift 2 ;;
        --hla-typing-env-prefix) HLA_TYPING_ENV_PREFIX_OVERRIDE="$2"; shift 2 ;;
        --cnv-env-prefix) CNV_ENV_PREFIX_OVERRIDE="$2"; shift 2 ;;
        --sv-env-prefix) SV_ENV_PREFIX_OVERRIDE="$2"; shift 2 ;;
        --no-testdata) INCLUDE_TESTDATA=false; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if [ "${MODE}" = "tumor-normal" ] || [ "${MODE}" = "tumor_normal" ] || [ "${MODE}" = "paired" ]; then
    cmd=(bash "${SCRIPT_DIR}/create_tumor_normal_project.sh"
        --tumor-fastq-source "${TUMOR_FASTQ_SOURCE}"
        --normal-fastq-source "${NORMAL_FASTQ_SOURCE}"
        --out-dir "${OUT_DIR}"
        --project-name "${PROJECT_NAME}"
        --copy-mode "${COPY_MODE}"
        --reference-dir "${REFERENCE_DIR}"
        --conda-base "${CONDA_BASE}")
    [ -n "${TUMOR_ID}" ] && cmd+=(--tumor-id "${TUMOR_ID}")
    [ -n "${NORMAL_ID}" ] && cmd+=(--normal-id "${NORMAL_ID}")
    [ -n "${REFERENCE_GENOME}" ] && cmd+=(--reference-genome "${REFERENCE_GENOME}")
    [ -n "${INTERVAL_BED}" ] && cmd+=(--interval-bed "${INTERVAL_BED}")
    [ -n "${ENV_ROOT}" ] && cmd+=(--env-root "${ENV_ROOT}")
    [ -n "${MAIN_ENV_PREFIX_OVERRIDE}" ] && cmd+=(--main-env-prefix "${MAIN_ENV_PREFIX_OVERRIDE}")
    [ -n "${VEP_ENV_PREFIX_OVERRIDE}" ] && cmd+=(--vep-env-prefix "${VEP_ENV_PREFIX_OVERRIDE}")
    [ -n "${HLA_ENV_PREFIX_OVERRIDE}" ] && cmd+=(--hla-env-prefix "${HLA_ENV_PREFIX_OVERRIDE}")
    [ -n "${HLA_TYPING_ENV_PREFIX_OVERRIDE}" ] && cmd+=(--hla-typing-env-prefix "${HLA_TYPING_ENV_PREFIX_OVERRIDE}")
    [ -n "${CNV_ENV_PREFIX_OVERRIDE}" ] && cmd+=(--cnv-env-prefix "${CNV_ENV_PREFIX_OVERRIDE}")
    [ -n "${SV_ENV_PREFIX_OVERRIDE}" ] && cmd+=(--sv-env-prefix "${SV_ENV_PREFIX_OVERRIDE}")
    [ "${INCLUDE_TESTDATA}" = false ] && cmd+=(--no-testdata)
    exec "${cmd[@]}"
fi

if [ "${MODE}" != "single" ]; then
    echo "Unknown --mode: ${MODE}. Supported: single, tumor-normal" >&2
    exit 1
fi

if [ -z "${FASTQ_SOURCE}" ] || [ -z "${OUT_DIR}" ]; then
    usage >&2
    exit 1
fi

if [ "${COPY_MODE}" != "copy" ] && [ "${COPY_MODE}" != "link" ]; then
    echo "--copy-mode must be copy or link" >&2
    exit 1
fi

FASTQ_SOURCE="$(abs_path "${FASTQ_SOURCE}")"
OUT_DIR="$(mkdir -p "${OUT_DIR}" && abs_path "${OUT_DIR}")"
REFERENCE_DIR="$(abs_path "${REFERENCE_DIR}")"
require_dir "FASTQ source" "${FASTQ_SOURCE}"
require_dir "Reference root" "${REFERENCE_DIR}"

if [ -z "${REFERENCE_GENOME}" ]; then
    REFERENCE_GENOME="${REFERENCE_DIR}/hg38/Homo_sapiens_assembly38.fasta"
else
    REFERENCE_GENOME="$(abs_path "${REFERENCE_GENOME}")"
fi
require_file "Reference genome" "${REFERENCE_GENOME}"
if [ -n "${INTERVAL_BED}" ]; then
    INTERVAL_BED="$(abs_path "${INTERVAL_BED}")"
    require_file "Interval BED" "${INTERVAL_BED}"
fi

R1="$(find_r1 "${FASTQ_SOURCE}")"
if [ -z "${R1}" ]; then
    echo "No R1 FASTQ found in ${FASTQ_SOURCE}" >&2
    exit 1
fi

if ! R2="$(guess_r2 "${R1}")"; then
    echo "Cannot infer R2 for ${R1}" >&2
    exit 1
fi

if [ -z "${SAMPLE_ID}" ]; then
    SAMPLE_ID="$(infer_sample_id "$(basename "${R1}")")"
fi

PIPELINE_DIR="${OUT_DIR}/pipeline"
CONFIG_DIR="${OUT_DIR}/configs"
DATA_DIR="${OUT_DIR}/data/${SAMPLE_ID}"
RESULT_DIR="${OUT_DIR}/results/${SAMPLE_ID}"
mkdir -p "${CONFIG_DIR}" "${DATA_DIR}" "${RESULT_DIR}" "${OUT_DIR}/logs"

copy_pipeline_code "${PIPELINE_DIR}"

TARGET_R1="${DATA_DIR}/${SAMPLE_ID}_R1.fastq.gz"
TARGET_R2="${DATA_DIR}/${SAMPLE_ID}_R2.fastq.gz"
if [ "${COPY_MODE}" = "link" ]; then
    ln -sf "${R1}" "${TARGET_R1}"
    ln -sf "${R2}" "${TARGET_R2}"
else
    cp "${R1}" "${TARGET_R1}"
    cp "${R2}" "${TARGET_R2}"
fi

if [ -n "${ENV_ROOT}" ]; then
    MAIN_ENV_PREFIX="${ENV_ROOT}/big_wes_pipeline_env"
    VEP_ENV_PREFIX="${ENV_ROOT}/wes_vep_env"
    HLA_ENV_PREFIX="${ENV_ROOT}/wes_hla_env"
    HLA_TYPING_ENV_PREFIX="${ENV_ROOT}/wes_hla_typing_env"
    CNV_ENV_PREFIX="${ENV_ROOT}/wes_cnv_env"
    SV_ENV_PREFIX="${ENV_ROOT}/wes_sv_env"
else
    MAIN_ENV_PREFIX="${CONDA_BASE}/envs/big_wes_pipeline_env"
    VEP_ENV_PREFIX="${OUT_DIR}/.conda_envs/wes_vep_env"
    HLA_ENV_PREFIX="${OUT_DIR}/.conda_envs/wes_hla_env"
    HLA_TYPING_ENV_PREFIX="${OUT_DIR}/.conda_envs/wes_hla_typing_env"
    CNV_ENV_PREFIX="${OUT_DIR}/.conda_envs/wes_cnv_env"
    SV_ENV_PREFIX="${OUT_DIR}/.conda_envs/wes_sv_env"
fi
MAIN_ENV_PREFIX="${MAIN_ENV_PREFIX_OVERRIDE:-${MAIN_ENV_PREFIX}}"
VEP_ENV_PREFIX="${VEP_ENV_PREFIX_OVERRIDE:-${VEP_ENV_PREFIX}}"
HLA_ENV_PREFIX="${HLA_ENV_PREFIX_OVERRIDE:-${HLA_ENV_PREFIX}}"
HLA_TYPING_ENV_PREFIX="${HLA_TYPING_ENV_PREFIX_OVERRIDE:-${HLA_TYPING_ENV_PREFIX}}"
CNV_ENV_PREFIX="${CNV_ENV_PREFIX_OVERRIDE:-${CNV_ENV_PREFIX}}"
SV_ENV_PREFIX="${SV_ENV_PREFIX_OVERRIDE:-${SV_ENV_PREFIX}}"

CONFIG_FILE="${CONFIG_DIR}/${SAMPLE_ID}.config.sh"
cat > "${CONFIG_FILE}" <<EOF
#!/bin/bash
# Auto-generated by scripts/create_wes_project.sh

PIPELINE_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")/../pipeline" && pwd)"
source "\${PIPELINE_DIR}/config.sh"

PROJECT_NAME="${PROJECT_NAME}"
PROJECT_DIR="\${PIPELINE_DIR}"
SAMPLE_ID="${SAMPLE_ID}"
SAMPLE_TYPE="germline"

CONDA_BASE="${CONDA_BASE}"
MAIN_ENV_PREFIX="${MAIN_ENV_PREFIX}"
VEP_ENV_PREFIX="${VEP_ENV_PREFIX}"
HLA_ENV_PREFIX="${HLA_ENV_PREFIX}"
HLA_TYPING_ENV_PREFIX="${HLA_TYPING_ENV_PREFIX}"
CNV_ENV_PREFIX="${CNV_ENV_PREFIX}"
SV_ENV_PREFIX="${SV_ENV_PREFIX}"
PIPELINE_EXTRA_PATHS="\${MAIN_ENV_PREFIX}/bin:\${VEP_ENV_PREFIX}/bin:\${HLA_ENV_PREFIX}/bin:\${HLA_TYPING_ENV_PREFIX}/bin:\${CNV_ENV_PREFIX}/bin:\${SV_ENV_PREFIX}/bin"
export PATH="\${PIPELINE_EXTRA_PATHS}:\${PATH}"
export VEP_ENV="\${VEP_ENV_PREFIX}"
PIPELINE_JAVA_HOME="\${MAIN_ENV_PREFIX}"
export JAVA_HOME="\${PIPELINE_JAVA_HOME}"

RAW_DATA_DIR="${DATA_DIR}"
FASTQ_R1="${TARGET_R1}"
FASTQ_R2="${TARGET_R2}"

REFERENCE_GENOME="${REFERENCE_GENOME}"
REFERENCE_DICT="${REFERENCE_GENOME%.fasta}.dict"
REFERENCE_BWA_INDEX="\${REFERENCE_GENOME}"
REFERENCE_GENOME_VERSION="hg38"

DBSNP_VCF="${REFERENCE_DIR}/dbsnp_146.hg38.vcf.gz"
DBSNP_VCF_INDEX="\${DBSNP_VCF}.tbi"
MILLS_VCF="${REFERENCE_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
THOUSAND_G_VCF="${REFERENCE_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
VEP_CACHE_DIR="${REFERENCE_DIR}/vep_cache"
VEP_FASTA="\${REFERENCE_GENOME}"
NEOANTIGEN_PROTEIN_FASTA="${REFERENCE_DIR}/protein/protein.fa"
RUN_HLA_TYPING="auto"
HLA_LA_GRAPH_DIR="${REFERENCE_DIR}/hla/PRG_MHC_GRCh38_withIMGT"

INTERVAL_FILE="${INTERVAL_BED}"
CNVKIT_TARGET_BED="${INTERVAL_BED}"
TMB_EFFECTIVE_CODING_BED="${INTERVAL_BED}"

RESULT_DIR="${RESULT_DIR}"
DIR_FASTQC="\${RESULT_DIR}/fastqc"
DIR_TRIMMED="\${RESULT_DIR}/trimmed"
DIR_ALIGNED="\${RESULT_DIR}/aligned"
DIR_POSTQC="\${RESULT_DIR}/post_align_qc"
DIR_BQSR="\${RESULT_DIR}/bqsr"
DIR_VARIANTS="\${RESULT_DIR}/variants"
DIR_ANNOTATION="\${RESULT_DIR}/annotation"
DIR_CNV="\${RESULT_DIR}/cnv"
DIR_SV="\${RESULT_DIR}/sv"
DIR_MSI="\${RESULT_DIR}/msi"
DIR_COVERAGE="\${RESULT_DIR}/coverage"
DIR_TMB="\${RESULT_DIR}/tmb"
DIR_NEOANTIGEN="\${RESULT_DIR}/neoantigen"
DIR_HLA_TYPING="\${RESULT_DIR}/hla_typing"
DIR_SUMMARY="\${RESULT_DIR}/summary"
DIR_MULTIQC="\${RESULT_DIR}/multiqc"
DIR_LOGS="${OUT_DIR}/logs"

RUN_BQSR=false
RUN_SNPEFF=false
RUN_VEP=true
RUN_NEOANTIGEN=true
RUN_HLA_BINDING=false
RUN_CNV=false
RUN_SV=false
RUN_MSI=false
RUN_COVERAGE=true
RUN_TMB=false

SKIP_BQSR=true
SKIP_SNPEFF=true
SKIP_VEP=false
SKIP_NEOANTIGEN=true
SKIP_CNV=true
SKIP_SV=true
SKIP_MSI=true
SKIP_MULTIQC=true
MIN_DISK_GB=1
EOF

RUNNER="${OUT_DIR}/run_${SAMPLE_ID}.sh"
cat > "${RUNNER}" <<EOF
#!/bin/bash
set -euo pipefail
cd "${PIPELINE_DIR}"
bash run_pipeline.sh --config "${CONFIG_FILE}"
EOF
chmod +x "${RUNNER}"

PROJECT_RUNNER="${OUT_DIR}/run_pipeline.sh"
cat > "${PROJECT_RUNNER}" <<EOF
#!/bin/bash
set -euo pipefail

PROJECT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${PIPELINE_DIR}"
CONFIG_FILE="${CONFIG_FILE}"
SAMPLE_RUNNER="${RUNNER}"

usage() {
    cat <<USAGE
Usage: bash run_pipeline.sh [all|check|status|list|step N|from N]

Default:
  bash run_pipeline.sh              Run the generated single-sample workflow.

Debug:
  bash run_pipeline.sh check
  bash run_pipeline.sh status
  bash run_pipeline.sh step 4
  bash run_pipeline.sh from 4
USAGE
}

run_internal() {
    (cd "\${PIPELINE_DIR}" && bash run_pipeline.sh --config "\${CONFIG_FILE}" "\$@")
}

case "\${1:-all}" in
    all|run)
        bash "\${SAMPLE_RUNNER}"
        ;;
    check|status|list|dry-run)
        run_internal "\$1" "\${2:-}"
        ;;
    step|from)
        if [ -z "\${2:-}" ]; then
            usage >&2
            exit 1
        fi
        run_internal "\$1" "\$2"
        ;;
    --help|-h|help)
        usage
        ;;
    *)
        echo "Unknown mode: \$1" >&2
        usage >&2
        exit 1
        ;;
esac
EOF
chmod +x "${PROJECT_RUNNER}"

cat > "${OUT_DIR}/README_RUN_${SAMPLE_ID}.md" <<EOF
# ${SAMPLE_ID} WES run

Generated project:

- Pipeline: ${PIPELINE_DIR}
- Config: ${CONFIG_FILE}
- Project runner: ${PROJECT_RUNNER}
- FASTQ R1: ${TARGET_R1}
- FASTQ R2: ${TARGET_R2}
- Results: ${RESULT_DIR}

Run:

\`\`\`bash
cd "${OUT_DIR}"
bash run_pipeline.sh
\`\`\`

Single-step debug:

\`\`\`bash
cd "${OUT_DIR}"
bash run_pipeline.sh status
bash run_pipeline.sh step 7c
\`\`\`
EOF

echo "Project generated:"
echo "  out_dir: ${OUT_DIR}"
echo "  sample: ${SAMPLE_ID}"
echo "  config: ${CONFIG_FILE}"
echo "  runner: ${RUNNER}"
echo "  project_runner: ${PROJECT_RUNNER}"
