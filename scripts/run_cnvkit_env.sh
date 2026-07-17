#!/bin/bash
# Run CNVkit from its dedicated conda/mamba environment.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEFAULT_ENV_ROOT="${ENV_ROOT:-${PROJECT_DIR}/.conda_envs}"
CNV_ENV="${CNV_ENV:-${CNV_ENV_PREFIX:-${DEFAULT_ENV_ROOT}/wes_cnv_env}}"
CNVKIT_BIN="${CNV_ENV}/bin/cnvkit.py"

if [ ! -x "${CNVKIT_BIN}" ]; then
    echo "CNVkit executable not found: ${CNVKIT_BIN}" >&2
    echo "Set CNV_ENV_PREFIX or CNV_ENV to the wes_cnv_env prefix." >&2
    exit 127
fi

export PATH="${CNV_ENV}/bin:${PATH}"
export R_HOME="${CNV_ENV}/lib/R"

exec "${CNVKIT_BIN}" "$@"
