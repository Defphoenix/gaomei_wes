#!/bin/bash
# Run Ensembl VEP from the dedicated local conda/micromamba environment.

set -euo pipefail

VEP_ENV="${VEP_ENV:-/Users/mac/Documents/wes/.conda_envs/wes_vep_env}"

unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
unset PERL_MB_OPT
unset PERL_MM_OPT

export LC_ALL=C
export LANG=C
export LC_CTYPE=C
export PATH="${VEP_ENV}/bin:${PATH}"
export PERL5LIB="${VEP_ENV}/lib/perl5/site_perl:${VEP_ENV}/lib/perl5/5.32/site_perl:${VEP_ENV}/lib/perl5/vendor_perl:${VEP_ENV}/lib/perl5/5.32/vendor_perl:${VEP_ENV}/lib/perl5/core_perl:${VEP_ENV}/lib/perl5/5.32/core_perl"

exec "${VEP_ENV}/bin/vep" "$@"
