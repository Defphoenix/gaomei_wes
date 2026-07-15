#!/bin/bash

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${PIPELINE_DIR}/config.sh"

PROJECT_DIR="${PIPELINE_DIR}"
SAMPLE_ID="mock_pair"
SAMPLE_TYPE="tumor"
TUMOR_SAMPLE_ID="tumor01"
NORMAL_SAMPLE_ID="normal01"
CALLER_MODE="mutect2"
RUN_TMB=true
SKIP_TMB=false
TMB_VEP_VCF="${PIPELINE_DIR}/testdata/tmb/mock_vep_somatic.vcf"
TMB_EFFECTIVE_CODING_BED="${PIPELINE_DIR}/testdata/tmb/mock_effective_coding.bed"
TMB_DENOMINATOR_VALIDATED=true
DIR_TMB="${TMPDIR:-/tmp}/gaomei_wes_tmb_test"
DIR_LOGS="${DIR_TMB}/logs"
MIN_DISK_GB=0
