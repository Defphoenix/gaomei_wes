#!/bin/bash

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${PIPELINE_DIR}/config.sh"

PROJECT_DIR="${PIPELINE_DIR}"
SAMPLE_ID="mock_tumor_vs_normal"
SAMPLE_TYPE="tumor"
CALLER_MODE="mutect2"
TEST_ROOT="${CNV_TEST_ROOT:-${TMPDIR:-/tmp}/gaomei_wes_cnv_matched_test}"
RESULT_DIR="${TEST_ROOT}"
DIR_CNV="${TEST_ROOT}/cnv"
DIR_ALIGNED="${TEST_ROOT}/aligned"
DIR_BQSR="${TEST_ROOT}/bqsr"
DIR_LOGS="${TEST_ROOT}/logs"
mkdir -p "${DIR_CNV}" "${DIR_ALIGNED}" "${DIR_BQSR}" "${DIR_LOGS}"

TUMOR_BAM="${DIR_BQSR}/tumor.bqsr.bam"
NORMAL_BAM="${DIR_BQSR}/normal.bqsr.bam"
printf 'mock tumor bam\n' > "${TUMOR_BAM}"
printf 'mock normal bam\n' > "${NORMAL_BAM}"
printf 'mock tumor index\n' > "${TUMOR_BAM}.bai"
printf 'mock normal index\n' > "${NORMAL_BAM}.bai"

REFERENCE_GENOME="${PIPELINE_DIR}/testdata/cnv/mock_reference.fa"
INTERVAL_FILE="${PIPELINE_DIR}/testdata/cnv/mock_targets.bed"
CNVKIT_TARGET_BED="${INTERVAL_FILE}"
CNVKIT_REFERENCE=""
CNV_METHOD="cnvkit"
CNV_REQUIRE_REFERENCE=true
CNVKIT_PROCESSES=1
RUN_CNV=true
SKIP_CNV=false
TOOL_CNVKIT="${PIPELINE_DIR}/testdata/cnv/mock_cnvkit.sh"
TOOL_SAMTOOLS="${PIPELINE_DIR}/testdata/cnv/mock_samtools.sh"
