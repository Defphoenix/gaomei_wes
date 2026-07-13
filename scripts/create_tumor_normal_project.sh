#!/bin/bash
# Generate a runnable tumor-normal WES project.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

TUMOR_FASTQ_SOURCE=""
NORMAL_FASTQ_SOURCE=""
OUT_DIR=""
TUMOR_ID=""
NORMAL_ID=""
PROJECT_NAME="wes_tumor_normal_project"
COPY_MODE="copy"
REFERENCE_DIR="/Users/mac/Documents/wes/reference_data"
REFERENCE_GENOME=""
INTERVAL_BED=""
CONDA_BASE="/Users/mac/anaconda3"
ENV_ROOT=""
MAIN_ENV_PREFIX_OVERRIDE=""
VEP_ENV_PREFIX_OVERRIDE=""
HLA_ENV_PREFIX_OVERRIDE=""
INCLUDE_TESTDATA=true

usage() {
    cat <<EOF
Usage: bash scripts/create_tumor_normal_project.sh --tumor-fastq-source DIR --normal-fastq-source DIR --out-dir DIR [options]

Required:
  --tumor-fastq-source DIR    Tumor paired FASTQ directory.
  --normal-fastq-source DIR   Matched normal paired FASTQ directory.
  --out-dir DIR               Target project directory.

Options:
  --tumor-id ID               Tumor sample ID. Default: inferred from tumor R1.
  --normal-id ID              Normal sample ID. Default: inferred from normal R1.
  --project-name NAME         Project name. Default: ${PROJECT_NAME}
  --copy-mode copy|link       Copy FASTQ or symlink. Default: ${COPY_MODE}
  --reference-dir DIR         Reference root. Default: ${REFERENCE_DIR}
  --reference-genome FA       Reference FASTA. Default: reference-dir/hg38/Homo_sapiens_assembly38.fasta
  --interval-bed BED          WES target BED.
  --conda-base DIR            Conda base path. Default: ${CONDA_BASE}
  --env-root DIR              Env root containing big_wes_pipeline_env and wes_vep_env.
  --main-env-prefix DIR       Full path to big_wes_pipeline_env. Overrides --env-root for main tools.
  --vep-env-prefix DIR        Full path to wes_vep_env. Overrides --env-root for VEP.
  --hla-env-prefix DIR        Full path to wes_hla_env. Overrides --env-root for HLA.
  --no-testdata               Do not copy bundled testdata/demo FASTQ.
  -h, --help                  Show this help.
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
        --tumor-fastq-source) TUMOR_FASTQ_SOURCE="$2"; shift 2 ;;
        --normal-fastq-source) NORMAL_FASTQ_SOURCE="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
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
        --no-testdata) INCLUDE_TESTDATA=false; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if [ -z "${TUMOR_FASTQ_SOURCE}" ] || [ -z "${NORMAL_FASTQ_SOURCE}" ] || [ -z "${OUT_DIR}" ]; then
    usage >&2
    exit 1
fi

if [ "${COPY_MODE}" != "copy" ] && [ "${COPY_MODE}" != "link" ]; then
    echo "--copy-mode must be copy or link" >&2
    exit 1
fi

TUMOR_FASTQ_SOURCE="$(abs_path "${TUMOR_FASTQ_SOURCE}")"
NORMAL_FASTQ_SOURCE="$(abs_path "${NORMAL_FASTQ_SOURCE}")"
OUT_DIR="$(mkdir -p "${OUT_DIR}" && abs_path "${OUT_DIR}")"
REFERENCE_DIR="$(abs_path "${REFERENCE_DIR}")"
require_dir "Tumor FASTQ source" "${TUMOR_FASTQ_SOURCE}"
require_dir "Normal FASTQ source" "${NORMAL_FASTQ_SOURCE}"
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

TUMOR_R1="$(find_r1 "${TUMOR_FASTQ_SOURCE}")"
NORMAL_R1="$(find_r1 "${NORMAL_FASTQ_SOURCE}")"
[ -n "${TUMOR_R1}" ] || { echo "No tumor R1 FASTQ found: ${TUMOR_FASTQ_SOURCE}" >&2; exit 1; }
[ -n "${NORMAL_R1}" ] || { echo "No normal R1 FASTQ found: ${NORMAL_FASTQ_SOURCE}" >&2; exit 1; }
TUMOR_R2="$(guess_r2 "${TUMOR_R1}")"
NORMAL_R2="$(guess_r2 "${NORMAL_R1}")"

[ -n "${TUMOR_ID}" ] || TUMOR_ID="$(infer_sample_id "$(basename "${TUMOR_R1}")")"
[ -n "${NORMAL_ID}" ] || NORMAL_ID="$(infer_sample_id "$(basename "${NORMAL_R1}")")"

PAIR_ID="${TUMOR_ID}_vs_${NORMAL_ID}"
PIPELINE_DIR="${OUT_DIR}/pipeline"
CONFIG_DIR="${OUT_DIR}/configs"
TUMOR_DATA_DIR="${OUT_DIR}/data/${TUMOR_ID}"
NORMAL_DATA_DIR="${OUT_DIR}/data/${NORMAL_ID}"
TUMOR_RESULT_DIR="${OUT_DIR}/results/${TUMOR_ID}"
NORMAL_RESULT_DIR="${OUT_DIR}/results/${NORMAL_ID}"
SOMATIC_RESULT_DIR="${OUT_DIR}/results/${PAIR_ID}"
mkdir -p "${CONFIG_DIR}" "${TUMOR_DATA_DIR}" "${NORMAL_DATA_DIR}" "${TUMOR_RESULT_DIR}" "${NORMAL_RESULT_DIR}" "${SOMATIC_RESULT_DIR}" "${OUT_DIR}/logs"

copy_pipeline_code "${PIPELINE_DIR}"

TUMOR_TARGET_R1="${TUMOR_DATA_DIR}/${TUMOR_ID}_R1.fastq.gz"
TUMOR_TARGET_R2="${TUMOR_DATA_DIR}/${TUMOR_ID}_R2.fastq.gz"
NORMAL_TARGET_R1="${NORMAL_DATA_DIR}/${NORMAL_ID}_R1.fastq.gz"
NORMAL_TARGET_R2="${NORMAL_DATA_DIR}/${NORMAL_ID}_R2.fastq.gz"

if [ "${COPY_MODE}" = "link" ]; then
    ln -sf "${TUMOR_R1}" "${TUMOR_TARGET_R1}"
    ln -sf "${TUMOR_R2}" "${TUMOR_TARGET_R2}"
    ln -sf "${NORMAL_R1}" "${NORMAL_TARGET_R1}"
    ln -sf "${NORMAL_R2}" "${NORMAL_TARGET_R2}"
else
    cp "${TUMOR_R1}" "${TUMOR_TARGET_R1}"
    cp "${TUMOR_R2}" "${TUMOR_TARGET_R2}"
    cp "${NORMAL_R1}" "${NORMAL_TARGET_R1}"
    cp "${NORMAL_R2}" "${NORMAL_TARGET_R2}"
fi

if [ -n "${ENV_ROOT}" ]; then
    MAIN_ENV_PREFIX="${ENV_ROOT}/big_wes_pipeline_env"
    VEP_ENV_PREFIX="${ENV_ROOT}/wes_vep_env"
    HLA_ENV_PREFIX="${ENV_ROOT}/wes_hla_env"
else
    MAIN_ENV_PREFIX="${CONDA_BASE}/envs/big_wes_pipeline_env"
    VEP_ENV_PREFIX="${OUT_DIR}/.conda_envs/wes_vep_env"
    HLA_ENV_PREFIX="${OUT_DIR}/.conda_envs/wes_hla_env"
fi
MAIN_ENV_PREFIX="${MAIN_ENV_PREFIX_OVERRIDE:-${MAIN_ENV_PREFIX}}"
VEP_ENV_PREFIX="${VEP_ENV_PREFIX_OVERRIDE:-${VEP_ENV_PREFIX}}"
HLA_ENV_PREFIX="${HLA_ENV_PREFIX_OVERRIDE:-${HLA_ENV_PREFIX}}"

write_common_config() {
    local sample_id="$1"
    local sample_type="$2"
    local r1="$3"
    local r2="$4"
    local result_dir="$5"
    local config_file="$6"

    cat > "${config_file}" <<EOF
#!/bin/bash
# Auto-generated tumor-normal config.

PIPELINE_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")/../pipeline" && pwd)"
source "\${PIPELINE_DIR}/config.sh"

PROJECT_NAME="${PROJECT_NAME}"
PROJECT_DIR="\${PIPELINE_DIR}"
ANALYSIS_MODE="tumor_normal"
SAMPLE_ID="${sample_id}"
SAMPLE_TYPE="${sample_type}"

CONDA_BASE="${CONDA_BASE}"
MAIN_ENV_PREFIX="${MAIN_ENV_PREFIX}"
VEP_ENV_PREFIX="${VEP_ENV_PREFIX}"
HLA_ENV_PREFIX="${HLA_ENV_PREFIX}"
PIPELINE_EXTRA_PATHS="\${MAIN_ENV_PREFIX}/bin:\${VEP_ENV_PREFIX}/bin:\${HLA_ENV_PREFIX}/bin"
export PATH="\${PIPELINE_EXTRA_PATHS}:\${PATH}"
export VEP_ENV="\${VEP_ENV_PREFIX}"
PIPELINE_JAVA_HOME="\${MAIN_ENV_PREFIX}"
export JAVA_HOME="\${PIPELINE_JAVA_HOME}"

RAW_DATA_DIR="$(dirname "${r1}")"
FASTQ_R1="${r1}"
FASTQ_R2="${r2}"

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
INTERVAL_FILE="${INTERVAL_BED}"
CNVKIT_TARGET_BED="${INTERVAL_BED}"

RESULT_DIR="${result_dir}"
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
DIR_SUMMARY="\${RESULT_DIR}/summary"
DIR_MULTIQC="\${RESULT_DIR}/multiqc"
DIR_LOGS="${OUT_DIR}/logs"

RUN_BQSR=false
RUN_SNPEFF=false
RUN_VEP=false
RUN_NEOANTIGEN=false
RUN_CNV=false
RUN_SV=false
RUN_MSI=false
RUN_COVERAGE=true
RUN_TMB=false

SKIP_BQSR=true
SKIP_VARIANT_CALLING=true
SKIP_VARIANT_FILTER=true
SKIP_SNPEFF=true
SKIP_VEP=true
SKIP_NEOANTIGEN=true
SKIP_CNV=true
SKIP_SV=true
SKIP_MSI=true
SKIP_MULTIQC=true
MIN_DISK_GB=1
EOF
}

TUMOR_CONFIG="${CONFIG_DIR}/${TUMOR_ID}.align.config.sh"
NORMAL_CONFIG="${CONFIG_DIR}/${NORMAL_ID}.align.config.sh"
SOMATIC_CONFIG="${CONFIG_DIR}/${PAIR_ID}.somatic.config.sh"
write_common_config "${TUMOR_ID}" "tumor" "${TUMOR_TARGET_R1}" "${TUMOR_TARGET_R2}" "${TUMOR_RESULT_DIR}" "${TUMOR_CONFIG}"
write_common_config "${NORMAL_ID}" "normal" "${NORMAL_TARGET_R1}" "${NORMAL_TARGET_R2}" "${NORMAL_RESULT_DIR}" "${NORMAL_CONFIG}"

cat > "${SOMATIC_CONFIG}" <<EOF
#!/bin/bash
# Auto-generated tumor-normal somatic config.

PIPELINE_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")/../pipeline" && pwd)"
source "\${PIPELINE_DIR}/config.sh"

PROJECT_NAME="${PROJECT_NAME}"
PROJECT_DIR="\${PIPELINE_DIR}"
ANALYSIS_MODE="tumor_normal"
SAMPLE_ID="${PAIR_ID}"
SAMPLE_TYPE="tumor"
TUMOR_SAMPLE_ID="${TUMOR_ID}"
NORMAL_SAMPLE_ID="${NORMAL_ID}"

CONDA_BASE="${CONDA_BASE}"
MAIN_ENV_PREFIX="${MAIN_ENV_PREFIX}"
VEP_ENV_PREFIX="${VEP_ENV_PREFIX}"
HLA_ENV_PREFIX="${HLA_ENV_PREFIX}"
PIPELINE_EXTRA_PATHS="\${MAIN_ENV_PREFIX}/bin:\${VEP_ENV_PREFIX}/bin:\${HLA_ENV_PREFIX}/bin"
export PATH="\${PIPELINE_EXTRA_PATHS}:\${PATH}"
export VEP_ENV="\${VEP_ENV_PREFIX}"
PIPELINE_JAVA_HOME="\${MAIN_ENV_PREFIX}"
export JAVA_HOME="\${PIPELINE_JAVA_HOME}"

TUMOR_BAM="${TUMOR_RESULT_DIR}/aligned/${TUMOR_ID}.dedup.bam"
NORMAL_BAM="${NORMAL_RESULT_DIR}/aligned/${NORMAL_ID}.dedup.bam"
FASTQ_R1="${TUMOR_TARGET_R1}"
FASTQ_R2="${TUMOR_TARGET_R2}"

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
INTERVAL_FILE="${INTERVAL_BED}"
CNVKIT_TARGET_BED="${INTERVAL_BED}"

RESULT_DIR="${SOMATIC_RESULT_DIR}"
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
DIR_SUMMARY="\${RESULT_DIR}/summary"
DIR_MULTIQC="\${RESULT_DIR}/multiqc"
DIR_LOGS="${OUT_DIR}/logs"

CALLER_MODE="mutect2"
RUN_BQSR=false
RUN_SNPEFF=false
RUN_VEP=true
RUN_NEOANTIGEN=true
RUN_HLA_BINDING=false
RUN_CNV=true
RUN_SV=false
RUN_MSI=true
RUN_COVERAGE=true
RUN_TMB=true

SKIP_FASTQC=true
SKIP_TRIM=true
SKIP_ALIGN=true
SKIP_SORT=true
SKIP_MARKDUP=true
SKIP_POSTQC=true
SKIP_BQSR=true
SKIP_VARIANT_CALLING=false
SKIP_VARIANT_FILTER=false
SKIP_SNPEFF=true
SKIP_VEP=false
SKIP_NEOANTIGEN=false
SKIP_CNV=false
SKIP_SV=true
SKIP_MSI=false
SKIP_COVERAGE=false
SKIP_MULTIQC=true
MIN_DISK_GB=1
EOF

RUNNER="${OUT_DIR}/run_${PAIR_ID}.sh"
cat > "${RUNNER}" <<EOF
#!/bin/bash
set -euo pipefail
cd "${PIPELINE_DIR}"

usage() {
    cat <<USAGE
Usage: bash \$(basename "\$0") [all|all-core|normal|normal5|tumor|tumor5|somatic|somatic-core|from-normal5|from-tumor|from-tumor5|from-somatic|from-somatic-core]

  all            Run normal alignment, tumor alignment, then full somatic steps.
  all-core       Run normal alignment, tumor alignment, then core somatic steps.
  normal         Run normal check + steps 1-5 only.
  normal5        Re-run normal step 5 only.
  tumor          Run tumor check + steps 1-5 only.
  tumor5         Re-run tumor step 5 only.
  somatic        Run full somatic steps only: Mutect2/filter/VEP/neoantigen/CNV/MSI/coverage/TMB/summary.
  somatic-core   Run core somatic steps only: Mutect2/filter/coverage.
  from-normal5   Re-run normal step 5, then tumor alignment and full somatic steps.
  from-tumor     Run tumor alignment and full somatic steps.
  from-tumor5    Re-run tumor step 5, then full somatic steps.
  from-somatic   Same as somatic.
  from-somatic-core Same as somatic-core.
USAGE
}

run_normal_all() {
    echo "[1/3] Align normal: ${NORMAL_ID}"
    bash run_pipeline.sh --config "${NORMAL_CONFIG}" check
    bash run_pipeline.sh --config "${NORMAL_CONFIG}" step 1
    bash run_pipeline.sh --config "${NORMAL_CONFIG}" step 2
    bash run_pipeline.sh --config "${NORMAL_CONFIG}" step 3
    bash run_pipeline.sh --config "${NORMAL_CONFIG}" step 4
    bash run_pipeline.sh --config "${NORMAL_CONFIG}" step 5
}

run_normal_step5() {
    echo "[resume] Normal step 5: ${NORMAL_ID}"
    bash run_pipeline.sh --config "${NORMAL_CONFIG}" step 5
}

run_tumor_all() {
    echo "[2/3] Align tumor: ${TUMOR_ID}"
    bash run_pipeline.sh --config "${TUMOR_CONFIG}" check
    bash run_pipeline.sh --config "${TUMOR_CONFIG}" step 1
    bash run_pipeline.sh --config "${TUMOR_CONFIG}" step 2
    bash run_pipeline.sh --config "${TUMOR_CONFIG}" step 3
    bash run_pipeline.sh --config "${TUMOR_CONFIG}" step 4
    bash run_pipeline.sh --config "${TUMOR_CONFIG}" step 5
}

run_tumor_step5() {
    echo "[resume] Tumor step 5: ${TUMOR_ID}"
    bash run_pipeline.sh --config "${TUMOR_CONFIG}" step 5
}

run_somatic() {
    echo "[3/3] Full somatic analysis: ${PAIR_ID}"
    bash run_pipeline.sh --config "${SOMATIC_CONFIG}" step 6
    bash run_pipeline.sh --config "${SOMATIC_CONFIG}" step 7
    bash run_pipeline.sh --config "${SOMATIC_CONFIG}" step 7c
    bash run_pipeline.sh --config "${SOMATIC_CONFIG}" step 7d
    bash run_pipeline.sh --config "${SOMATIC_CONFIG}" step 8
    bash run_pipeline.sh --config "${SOMATIC_CONFIG}" step 9
    bash run_pipeline.sh --config "${SOMATIC_CONFIG}" step 11
    bash run_pipeline.sh --config "${SOMATIC_CONFIG}" step 12
    bash run_pipeline.sh --config "${SOMATIC_CONFIG}" step 13
}

run_somatic_core() {
    echo "[3/3] Core somatic analysis: ${PAIR_ID}"
    bash run_pipeline.sh --config "${SOMATIC_CONFIG}" step 6
    bash run_pipeline.sh --config "${SOMATIC_CONFIG}" step 7
    bash run_pipeline.sh --config "${SOMATIC_CONFIG}" step 11
}

case "\${1:-all}" in
    all)
        run_normal_all
        run_tumor_all
        run_somatic
        ;;
    all-core)
        run_normal_all
        run_tumor_all
        run_somatic_core
        ;;
    normal)
        run_normal_all
        ;;
    normal5)
        run_normal_step5
        ;;
    tumor)
        run_tumor_all
        ;;
    tumor5)
        run_tumor_step5
        ;;
    somatic|from-somatic)
        run_somatic
        ;;
    somatic-core|from-somatic-core)
        run_somatic_core
        ;;
    from-normal5)
        run_normal_step5
        run_tumor_all
        run_somatic
        ;;
    from-tumor)
        run_tumor_all
        run_somatic
        ;;
    from-tumor5)
        run_tumor_step5
        run_somatic
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
chmod +x "${RUNNER}"

PROJECT_RUNNER="${OUT_DIR}/run_pipeline.sh"
cat > "${PROJECT_RUNNER}" <<EOF
#!/bin/bash
set -euo pipefail

PROJECT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="${PIPELINE_DIR}"
NORMAL_CONFIG="${NORMAL_CONFIG}"
TUMOR_CONFIG="${TUMOR_CONFIG}"
SOMATIC_CONFIG="${SOMATIC_CONFIG}"
PAIR_RUNNER="${RUNNER}"

usage() {
    cat <<USAGE
Usage: bash run_pipeline.sh [all-core|all|full|normal|tumor|somatic-core|somatic|check|status|list|step ROLE N|from ROLE N]

Default:
  bash run_pipeline.sh              Run full: normal + tumor + Mutect2/filter/VEP/neoantigen/CNV/MSI/coverage/TMB/summary.

Project modes:
  all-core                          Run normal, tumor, then core somatic analysis.
  all | full                        Run normal, tumor, then full somatic analysis.
  normal                            Run normal FASTQ to BAM.
  tumor                             Run tumor FASTQ to BAM.
  somatic-core                      Run Mutect2/filter/coverage only.
  somatic                           Run full somatic analysis only.

Debug modes:
  check                             Check normal, tumor and somatic configs.
  status [ROLE]                     Show status for normal/tumor/somatic.
  list                              List internal pipeline steps.
  step normal 4                     Run one normal step.
  step tumor 5                      Run one tumor step.
  step somatic 6                    Run one somatic step.
  from normal 4                     Continue normal from step 4.
  from tumor 4                      Continue tumor from step 4.
  from somatic 6                    Continue somatic from step 6.
USAGE
}

config_for_role() {
    case "\$1" in
        normal) echo "\${NORMAL_CONFIG}" ;;
        tumor) echo "\${TUMOR_CONFIG}" ;;
        somatic|pair) echo "\${SOMATIC_CONFIG}" ;;
        *) echo "Unknown role: \$1" >&2; return 1 ;;
    esac
}

run_internal() {
    local config="\$1"
    shift
    (cd "\${PIPELINE_DIR}" && bash run_pipeline.sh --config "\${config}" "\$@")
}

case "\${1:-all}" in
    all-core|core)
        bash "\${PAIR_RUNNER}" all-core
        ;;
    all|full)
        bash "\${PAIR_RUNNER}" all
        ;;
    normal|tumor|somatic|somatic-core)
        bash "\${PAIR_RUNNER}" "\$1"
        ;;
    check)
        run_internal "\${NORMAL_CONFIG}" check
        run_internal "\${TUMOR_CONFIG}" check
        run_internal "\${SOMATIC_CONFIG}" check
        ;;
    status)
        role="\${2:-somatic}"
        config="\$(config_for_role "\${role}")"
        run_internal "\${config}" status
        ;;
    list)
        run_internal "\${SOMATIC_CONFIG}" list
        ;;
    step|from)
        action="\$1"
        role="\${2:-}"
        step_id="\${3:-}"
        if [ -z "\${role}" ] || [ -z "\${step_id}" ]; then
            usage >&2
            exit 1
        fi
        config="\$(config_for_role "\${role}")"
        run_internal "\${config}" "\${action}" "\${step_id}"
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

cat > "${OUT_DIR}/README_RUN_${PAIR_ID}.md" <<EOF
# ${PAIR_ID} tumor-normal WES run

- Tumor config: ${TUMOR_CONFIG}
- Normal config: ${NORMAL_CONFIG}
- Somatic config: ${SOMATIC_CONFIG}
- Runner: ${RUNNER}
- Project runner: ${PROJECT_RUNNER}

Run:

\`\`\`bash
cd "${OUT_DIR}"
bash run_pipeline.sh
\`\`\`

Core test run:

\`\`\`bash
cd "${OUT_DIR}"
bash run_pipeline.sh all-core
\`\`\`

Full somatic run:

\`\`\`bash
cd "${OUT_DIR}"
bash run_pipeline.sh full
\`\`\`

Single-step debug:

\`\`\`bash
cd "${OUT_DIR}"
bash run_pipeline.sh step normal 4
bash run_pipeline.sh step tumor 5
bash run_pipeline.sh step somatic 6
\`\`\`

Somatic outputs:

\`\`\`text
${SOMATIC_RESULT_DIR}/variants/${PAIR_ID}.mutect2.filtered.vcf.gz
${SOMATIC_RESULT_DIR}/variants/${PAIR_ID}.mutect2.pass.vcf.gz
${SOMATIC_RESULT_DIR}/annotation/${PAIR_ID}.vep.vcf.gz
\`\`\`
EOF

echo "Tumor-normal project generated:"
echo "  out_dir: ${OUT_DIR}"
echo "  pair: ${PAIR_ID}"
echo "  tumor_config: ${TUMOR_CONFIG}"
echo "  normal_config: ${NORMAL_CONFIG}"
echo "  somatic_config: ${SOMATIC_CONFIG}"
echo "  runner: ${RUNNER}"
echo "  project_runner: ${PROJECT_RUNNER}"
