#!/bin/bash
# Create conda/mamba environments used by the WES pipeline.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

ENV_ROOT="${ENV_ROOT:-${PROJECT_DIR}/.conda_envs}"
MAMBA_BIN="${MAMBA_BIN:-}"
CREATE_HLA="${CREATE_HLA:-false}"
CREATE_HLA_TYPING="${CREATE_HLA_TYPING:-false}"
CREATE_CNV="${CREATE_CNV:-false}"
CREATE_SV="${CREATE_SV:-false}"
FETCH_MHCFLURRY_MODELS="${FETCH_MHCFLURRY_MODELS:-false}"
VERIFY_INSTALL="${VERIFY_INSTALL:-true}"
WRITE_MANIFESTS="${WRITE_MANIFESTS:-true}"
UPDATE_EXISTING="${UPDATE_EXISTING:-false}"
SOLVER_ARGS="${SOLVER_ARGS:---override-channels --channel-priority flexible}"
CLEAN_INCOMPLETE="${CLEAN_INCOMPLETE:-false}"

usage() {
    cat <<EOF
Usage: bash scripts/create_conda_envs.sh [options]

Options:
  --env-root DIR       Environment prefix root. Default: ${ENV_ROOT}
  --mamba-bin CMD      mamba/micromamba/conda command. Default: auto-detect mamba > micromamba > conda
  --solver-args ARGS   Extra solver args. Default: ${SOLVER_ARGS}
  --clean-incomplete   Remove incomplete env dirs before recreating.
  --update-existing    Update existing prefixes from current YML files and prune removed packages.
  --with-hla           Also create mhcflurry HLA environment.
  --with-hla-typing    Also create the Linux HLA*LA typing environment.
  --with-cnv           Also create the isolated CNVkit environment.
  --with-sv            Also create the isolated Manta environment (Linux only).
  --all-optional       Create binding, HLA typing, CNVkit and Manta environments.
  --fetch-hla-models   After --with-hla, fetch MHCflurry class-I presentation models.
  --no-verify          Skip post-install tool checks.
  --no-manifests       Do not write resolved package/version manifests.
  -h, --help           Show this help.

Examples:
  bash scripts/create_conda_envs.sh --env-root /PUBLIC/envs/wes --with-hla --with-cnv
  bash scripts/create_conda_envs.sh --mamba-bin micromamba --all-optional --fetch-hla-models
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --env-root) ENV_ROOT="$2"; shift 2 ;;
        --mamba-bin) MAMBA_BIN="$2"; shift 2 ;;
        --solver-args) SOLVER_ARGS="$2"; shift 2 ;;
        --clean-incomplete) CLEAN_INCOMPLETE=true; shift ;;
        --update-existing) UPDATE_EXISTING=true; shift ;;
        --with-hla) CREATE_HLA=true; shift ;;
        --with-hla-typing) CREATE_HLA_TYPING=true; shift ;;
        --with-cnv) CREATE_CNV=true; shift ;;
        --with-sv) CREATE_SV=true; shift ;;
        --all-optional) CREATE_HLA=true; CREATE_HLA_TYPING=true; CREATE_CNV=true; CREATE_SV=true; shift ;;
        --fetch-hla-models) FETCH_MHCFLURRY_MODELS=true; CREATE_HLA=true; shift ;;
        --no-verify) VERIFY_INSTALL=false; shift ;;
        --no-manifests) WRITE_MANIFESTS=false; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

detect_mamba_bin() {
    if [ -n "${MAMBA_BIN}" ]; then
        command -v "${MAMBA_BIN}" >/dev/null 2>&1 || {
            echo "[ERROR] command not found: ${MAMBA_BIN}" >&2
            return 1
        }
        "${MAMBA_BIN}" --version >/dev/null 2>&1 || {
            echo "[ERROR] solver exists but cannot run: ${MAMBA_BIN}" >&2
            echo "        Repair it or choose another solver with --mamba-bin." >&2
            return 1
        }
        return 0
    fi

    if command -v mamba >/dev/null 2>&1 && mamba --version >/dev/null 2>&1; then
        MAMBA_BIN="mamba"
    elif command -v micromamba >/dev/null 2>&1 && micromamba --version >/dev/null 2>&1; then
        MAMBA_BIN="micromamba"
    elif command -v conda >/dev/null 2>&1 && conda --version >/dev/null 2>&1; then
        MAMBA_BIN="conda"
    else
        echo "[ERROR] mamba/micromamba/conda not found in PATH" >&2
        return 1
    fi
}

run_create() {
    local name="$1"
    local yml="$2"
    local prefix="${ENV_ROOT}/${name}"

    if [ -f "${prefix}/conda-meta/history" ]; then
        if [ "${UPDATE_EXISTING}" = true ]; then
            echo "[UPDATE] ${name}: ${prefix}"
            if [ "${MAMBA_BIN##*/}" = "micromamba" ]; then
                "${MAMBA_BIN}" update -y ${SOLVER_ARGS} -p "${prefix}" -f "${PROJECT_DIR}/${yml}" --prune
            elif [ "${MAMBA_BIN##*/}" = "conda" ]; then
                CONDA_OVERRIDE_CHANNELS=true CONDA_CHANNELS=conda-forge,bioconda CONDA_CHANNEL_PRIORITY=flexible \
                    "${MAMBA_BIN}" env update -p "${prefix}" -f "${PROJECT_DIR}/${yml}" --prune
            else
                "${MAMBA_BIN}" env update -y ${SOLVER_ARGS} -p "${prefix}" -f "${PROJECT_DIR}/${yml}" --prune
            fi
        else
            echo "[SKIP] ${name} exists: ${prefix} (use --update-existing to apply YML changes)"
        fi
        return 0
    fi

    if [ -d "${prefix}" ]; then
        if [ "${CLEAN_INCOMPLETE}" = true ]; then
            echo "[CLEAN] removing incomplete env dir: ${prefix}"
            rm -rf "${prefix}"
        else
            echo "[ERROR] incomplete env dir exists: ${prefix}" >&2
            echo "        Retry with --clean-incomplete, or remove it manually:" >&2
            echo "        rm -rf '${prefix}'" >&2
            return 1
        fi
    fi

    echo "[CREATE] ${name}: ${prefix}"
    if [ "${MAMBA_BIN##*/}" = "micromamba" ]; then
        "${MAMBA_BIN}" create -y ${SOLVER_ARGS} -p "${prefix}" -f "${PROJECT_DIR}/${yml}"
    elif [ "${MAMBA_BIN##*/}" = "conda" ]; then
        CONDA_OVERRIDE_CHANNELS=true CONDA_CHANNELS=conda-forge,bioconda CONDA_CHANNEL_PRIORITY=flexible \
            "${MAMBA_BIN}" env create -p "${prefix}" -f "${PROJECT_DIR}/${yml}"
    else
        "${MAMBA_BIN}" env create -y ${SOLVER_ARGS} -p "${prefix}" -f "${PROJECT_DIR}/${yml}"
    fi
}

run_in_env() {
    local prefix="$1"
    shift
    if [ "${MAMBA_BIN##*/}" = "micromamba" ]; then
        "${MAMBA_BIN}" run -p "${prefix}" "$@"
    elif [ "${MAMBA_BIN##*/}" = "conda" ]; then
        "${MAMBA_BIN}" run -p "${prefix}" "$@"
    else
        "${MAMBA_BIN}" run -p "${prefix}" "$@"
    fi
}

verify_tool() {
    local prefix="$1"
    local label="$2"
    shift 2
    printf '[CHECK] %-18s' "${label}"
    if run_in_env "${prefix}" "$@" >/tmp/wes_env_check.$$ 2>&1; then
        echo "OK"
    else
        echo "FAILED"
        sed 's/^/        /' /tmp/wes_env_check.$$ >&2 || true
        rm -f /tmp/wes_env_check.$$
        return 1
    fi
    rm -f /tmp/wes_env_check.$$
}

write_manifest() {
    local name="$1"
    local yml="$2"
    local prefix="${ENV_ROOT}/${name}"
    local manifest_dir="${ENV_ROOT}/manifests"

    [ "${WRITE_MANIFESTS}" = true ] || return 0
    mkdir -p "${manifest_dir}"
    cp "${PROJECT_DIR}/${yml}" "${manifest_dir}/${name}.requested.yml"
    if ! "${MAMBA_BIN}" list -p "${prefix}" --explicit > "${manifest_dir}/${name}.explicit.txt"; then
        echo "[WARN] solver does not support explicit export; continuing without lock file" >&2
        rm -f "${manifest_dir}/${name}.explicit.txt"
    fi
    if ! "${MAMBA_BIN}" list -p "${prefix}" > "${manifest_dir}/${name}.packages.txt"; then
        echo "[WARN] unable to export package list for ${name}" >&2
        rm -f "${manifest_dir}/${name}.packages.txt"
        return 0
    fi
    echo "[MANIFEST] ${manifest_dir}/${name}.packages.txt"
}

detect_mamba_bin

if [ "${CREATE_HLA_TYPING}" = true ] && [ "$(uname -s)" != "Linux" ]; then
    echo "[ERROR] HLA*LA Bioconda package is Linux-only; create wes_hla_typing_env on the analysis server." >&2
    exit 1
fi

echo "[INFO] project: ${PROJECT_DIR}"
echo "[INFO] env root: ${ENV_ROOT}"
echo "[INFO] solver: ${MAMBA_BIN}"
echo "[INFO] solver args: ${SOLVER_ARGS}"

mkdir -p "${ENV_ROOT}"

run_create "big_wes_pipeline_env" "wes_env_version01.yml"
run_create "wes_vep_env" "wes_vep_env.yml"
write_manifest "big_wes_pipeline_env" "wes_env_version01.yml"
write_manifest "wes_vep_env" "wes_vep_env.yml"

if [ "${CREATE_HLA}" = true ]; then
    run_create "wes_hla_env" "wes_hla_env.yml"
    write_manifest "wes_hla_env" "wes_hla_env.yml"
    if [ "${FETCH_MHCFLURRY_MODELS}" = true ]; then
        echo "[FETCH] MHCflurry models_class1_presentation"
        run_in_env "${ENV_ROOT}/wes_hla_env" mhcflurry-downloads fetch models_class1_presentation
    fi
fi

if [ "${CREATE_HLA_TYPING}" = true ]; then
    run_create "wes_hla_typing_env" "wes_hla_typing_env.yml"
    write_manifest "wes_hla_typing_env" "wes_hla_typing_env.yml"
fi

if [ "${CREATE_CNV}" = true ]; then
    run_create "wes_cnv_env" "wes_cnv_env.yml"
    write_manifest "wes_cnv_env" "wes_cnv_env.yml"
fi

if [ "${CREATE_SV}" = true ]; then
    run_create "wes_sv_env" "wes_sv_env.yml"
    write_manifest "wes_sv_env" "wes_sv_env.yml"
fi

ENV_SH="${ENV_ROOT}/env.sh"
cat > "${ENV_SH}" <<EOF
# Source this file before running or testing the WES pipeline.
# Usage:
#   source "${ENV_SH}"

export MAIN_ENV_PREFIX="${ENV_ROOT}/big_wes_pipeline_env"
export VEP_ENV_PREFIX="${ENV_ROOT}/wes_vep_env"
export HLA_ENV_PREFIX="${ENV_ROOT}/wes_hla_env"
export HLA_TYPING_ENV_PREFIX="${ENV_ROOT}/wes_hla_typing_env"
export CNV_ENV_PREFIX="${ENV_ROOT}/wes_cnv_env"
export SV_ENV_PREFIX="${ENV_ROOT}/wes_sv_env"
export JAVA_HOME="\${MAIN_ENV_PREFIX}"
export VEP_ENV="\${VEP_ENV_PREFIX}"
export PATH="\${MAIN_ENV_PREFIX}/bin:\${VEP_ENV_PREFIX}/bin:\${HLA_ENV_PREFIX}/bin:\${HLA_TYPING_ENV_PREFIX}/bin:\${CNV_ENV_PREFIX}/bin:\${SV_ENV_PREFIX}/bin:\${PATH}"
EOF

if [ "${VERIFY_INSTALL}" = true ]; then
    echo
    echo "Verifying core tools..."
    verify_tool "${ENV_ROOT}/big_wes_pipeline_env" "python" python --version
    verify_tool "${ENV_ROOT}/big_wes_pipeline_env" "java" java -version
    verify_tool "${ENV_ROOT}/big_wes_pipeline_env" "gatk" gatk --version
    verify_tool "${ENV_ROOT}/big_wes_pipeline_env" "bwa" bash -c 'bwa 2>&1 | grep -qi version'
    verify_tool "${ENV_ROOT}/big_wes_pipeline_env" "samtools" samtools --version
    verify_tool "${ENV_ROOT}/big_wes_pipeline_env" "bcftools" bcftools --version
    verify_tool "${ENV_ROOT}/big_wes_pipeline_env" "fastp" fastp --version
    verify_tool "${ENV_ROOT}/big_wes_pipeline_env" "fastqc" fastqc --version
    verify_tool "${ENV_ROOT}/big_wes_pipeline_env" "bedtools" bedtools --version
    verify_tool "${ENV_ROOT}/big_wes_pipeline_env" "picard" picard -h
    verify_tool "${ENV_ROOT}/big_wes_pipeline_env" "mosdepth" mosdepth --version
    verify_tool "${ENV_ROOT}/big_wes_pipeline_env" "msisensor-pro" bash -c 'msisensor-pro 2>&1 | grep -qi msisensor'
    verify_tool "${ENV_ROOT}/big_wes_pipeline_env" "multiqc" multiqc --version
    verify_tool "${ENV_ROOT}/big_wes_pipeline_env" "snpeff" snpEff -version
    verify_tool "${ENV_ROOT}/wes_vep_env" "vep" bash -c 'vep --help >/dev/null'
    if [ "${CREATE_HLA}" = true ]; then
        verify_tool "${ENV_ROOT}/wes_hla_env" "mhcflurry" mhcflurry-predict --help
    fi
    if [ "${CREATE_HLA_TYPING}" = true ]; then
        verify_tool "${ENV_ROOT}/wes_hla_typing_env" "HLA-LA" bash -c 'command -v HLA-LA.pl >/dev/null'
    fi
    if [ "${CREATE_CNV}" = true ]; then
        verify_tool "${ENV_ROOT}/wes_cnv_env" "cnvkit" cnvkit.py version
    fi
    if [ "${CREATE_SV}" = true ]; then
        verify_tool "${ENV_ROOT}/wes_sv_env" "manta" configManta.py --help
    fi
fi

cat <<EOF

Environment creation finished.

Use these paths in config.sh:
  MAIN_ENV_PREFIX="${ENV_ROOT}/big_wes_pipeline_env"
  VEP_ENV_PREFIX="${ENV_ROOT}/wes_vep_env"
  HLA_ENV_PREFIX="${ENV_ROOT}/wes_hla_env"
  HLA_TYPING_ENV_PREFIX="${ENV_ROOT}/wes_hla_typing_env"
  CNV_ENV_PREFIX="${ENV_ROOT}/wes_cnv_env"
  SV_ENV_PREFIX="${ENV_ROOT}/wes_sv_env"

Because these environments are prefix-based, activate them by full path:
  mamba activate "${ENV_ROOT}/big_wes_pipeline_env"
  mamba activate "${ENV_ROOT}/wes_vep_env"

Or source the helper before running/testing pipeline tools:
  source "${ENV_SH}"
EOF
