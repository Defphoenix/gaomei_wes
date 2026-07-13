#!/bin/bash
# Create conda/mamba environments used by the WES pipeline.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

ENV_ROOT="${ENV_ROOT:-${PROJECT_DIR}/.conda_envs}"
MAMBA_BIN="${MAMBA_BIN:-mamba}"
CREATE_HLA="${CREATE_HLA:-false}"
SOLVER_ARGS="${SOLVER_ARGS:---override-channels --channel-priority flexible}"
CLEAN_INCOMPLETE="${CLEAN_INCOMPLETE:-false}"

usage() {
    cat <<EOF
Usage: bash scripts/create_conda_envs.sh [options]

Options:
  --env-root DIR       Environment prefix root. Default: ${ENV_ROOT}
  --mamba-bin CMD      mamba/micromamba/conda command. Default: ${MAMBA_BIN}
  --solver-args ARGS   Extra solver args. Default: ${SOLVER_ARGS}
  --clean-incomplete   Remove incomplete env dirs before recreating.
  --with-hla           Also create mhcflurry HLA environment.
  -h, --help           Show this help.

Examples:
  bash scripts/create_conda_envs.sh --env-root /PUBLIC/envs/wes
  bash scripts/create_conda_envs.sh --mamba-bin micromamba --with-hla
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --env-root) ENV_ROOT="$2"; shift 2 ;;
        --mamba-bin) MAMBA_BIN="$2"; shift 2 ;;
        --solver-args) SOLVER_ARGS="$2"; shift 2 ;;
        --clean-incomplete) CLEAN_INCOMPLETE=true; shift ;;
        --with-hla) CREATE_HLA=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

run_create() {
    local name="$1"
    local yml="$2"
    local prefix="${ENV_ROOT}/${name}"

    if [ -f "${prefix}/conda-meta/history" ]; then
        echo "[SKIP] ${name} exists: ${prefix}"
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
    if [ "${MAMBA_BIN}" = "micromamba" ]; then
        "${MAMBA_BIN}" create -y ${SOLVER_ARGS} -p "${prefix}" -f "${PROJECT_DIR}/${yml}"
    else
        "${MAMBA_BIN}" env create ${SOLVER_ARGS} -p "${prefix}" -f "${PROJECT_DIR}/${yml}"
    fi
}

mkdir -p "${ENV_ROOT}"

run_create "big_wes_pipeline_env" "wes_env_version01.yml"
run_create "wes_vep_env" "wes_vep_env.yml"

if [ "${CREATE_HLA}" = true ]; then
    run_create "wes_hla_env" "wes_hla_env.yml"
fi

ENV_SH="${ENV_ROOT}/env.sh"
cat > "${ENV_SH}" <<EOF
# Source this file before running or testing the WES pipeline.
# Usage:
#   source "${ENV_SH}"

export MAIN_ENV_PREFIX="${ENV_ROOT}/big_wes_pipeline_env"
export VEP_ENV_PREFIX="${ENV_ROOT}/wes_vep_env"
export HLA_ENV_PREFIX="${ENV_ROOT}/wes_hla_env"
export JAVA_HOME="\${MAIN_ENV_PREFIX}"
export VEP_ENV="\${VEP_ENV_PREFIX}"
export PATH="\${MAIN_ENV_PREFIX}/bin:\${VEP_ENV_PREFIX}/bin:\${HLA_ENV_PREFIX}/bin:\${PATH}"
EOF

cat <<EOF

Environment creation finished.

Use these paths in config.sh:
  MAIN_ENV_PREFIX="${ENV_ROOT}/big_wes_pipeline_env"
  VEP_ENV_PREFIX="${ENV_ROOT}/wes_vep_env"
  HLA_ENV_PREFIX="${ENV_ROOT}/wes_hla_env"

Because these environments are prefix-based, activate them by full path:
  mamba activate "${ENV_ROOT}/big_wes_pipeline_env"
  mamba activate "${ENV_ROOT}/wes_vep_env"

Or source the helper before running/testing pipeline tools:
  source "${ENV_SH}"
EOF
