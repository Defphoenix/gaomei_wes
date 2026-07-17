#!/bin/bash
# Run Ensembl VEP from its dedicated conda/mamba environment.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEFAULT_ENV_ROOT="${ENV_ROOT:-${PROJECT_DIR}/.conda_envs}"
VEP_ENV="${VEP_ENV:-${VEP_ENV_PREFIX:-${DEFAULT_ENV_ROOT}/wes_vep_env}}"
VEP_BIN="${VEP_ENV}/bin/vep"

if [ ! -x "${VEP_BIN}" ]; then
    echo "VEP executable not found: ${VEP_BIN}" >&2
    echo "Set VEP_ENV_PREFIX or VEP_ENV to the wes_vep_env prefix." >&2
    exit 127
fi

unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
unset PERL_MB_OPT
unset PERL_MM_OPT

export LC_ALL=C
export LANG=C
export LC_CTYPE=C
export PATH="${VEP_ENV}/bin:${PATH}"

perl5lib=""
for perl_dir in \
    "${VEP_ENV}/lib/perl5/site_perl" \
    "${VEP_ENV}"/lib/perl5/*/site_perl \
    "${VEP_ENV}/lib/perl5/vendor_perl" \
    "${VEP_ENV}"/lib/perl5/*/vendor_perl \
    "${VEP_ENV}/lib/perl5/core_perl" \
    "${VEP_ENV}"/lib/perl5/*/core_perl; do
    [ -d "${perl_dir}" ] || continue
    if [ -z "${perl5lib}" ]; then
        perl5lib="${perl_dir}"
    else
        perl5lib="${perl5lib}:${perl_dir}"
    fi
done
[ -z "${perl5lib}" ] || export PERL5LIB="${perl5lib}"

exec "${VEP_BIN}" "$@"
