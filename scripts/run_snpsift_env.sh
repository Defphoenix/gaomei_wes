#!/bin/bash
# Run SnpSift with the Java runtime from the dedicated SnpEff environment.

set -euo pipefail

SNPEFF_ENV="${SNPEFF_ENV:-${SNPEFF_ENV_PREFIX:-/Users/mac/Documents/wes/.conda_envs/wes_snpeff_env}}"
SNPSIFT_BIN="${SNPEFF_ENV}/bin/SnpSift"

if [ ! -x "${SNPSIFT_BIN}" ] && [ -n "${MAIN_ENV_PREFIX:-}" ]; then
    SNPEFF_ENV="$(dirname "${MAIN_ENV_PREFIX}")/wes_snpeff_env"
    SNPSIFT_BIN="${SNPEFF_ENV}/bin/SnpSift"
fi

if [ ! -x "${SNPSIFT_BIN}" ]; then
    echo "SnpSift executable not found: ${SNPSIFT_BIN}" >&2
    echo "Set SNPEFF_ENV_PREFIX or SNPEFF_ENV to the wes_snpeff_env prefix." >&2
    exit 127
fi

export JAVA_HOME="${SNPEFF_ENV}"
export PATH="${SNPEFF_ENV}/bin:${PATH}"

exec "${SNPSIFT_BIN}" "$@"
