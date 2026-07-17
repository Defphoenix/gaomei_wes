#!/bin/bash
# Run SnpEff with the Java runtime from its dedicated environment.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEFAULT_ENV_ROOT="${ENV_ROOT:-${PROJECT_DIR}/.conda_envs}"
SNPEFF_ENV="${SNPEFF_ENV:-${SNPEFF_ENV_PREFIX:-${DEFAULT_ENV_ROOT}/wes_snpeff_env}}"
SNPEFF_BIN="${SNPEFF_ENV}/bin/snpEff"

if [ ! -x "${SNPEFF_BIN}" ] && [ -n "${MAIN_ENV_PREFIX:-}" ]; then
    SNPEFF_ENV="$(dirname "${MAIN_ENV_PREFIX}")/wes_snpeff_env"
    SNPEFF_BIN="${SNPEFF_ENV}/bin/snpEff"
fi

if [ ! -x "${SNPEFF_BIN}" ]; then
    echo "SnpEff executable not found: ${SNPEFF_BIN}" >&2
    echo "Set SNPEFF_ENV_PREFIX or SNPEFF_ENV to the wes_snpeff_env prefix." >&2
    exit 127
fi

export JAVA_HOME="${SNPEFF_ENV}"
export PATH="${SNPEFF_ENV}/bin:${PATH}"

if [ "${1:-}" = "--version" ]; then
    shift
    set -- -version "$@"
fi

exec "${SNPEFF_BIN}" "$@"
