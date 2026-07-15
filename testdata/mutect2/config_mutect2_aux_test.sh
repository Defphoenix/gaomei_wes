#!/bin/bash

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${PIPELINE_DIR}/config.sh"

PROJECT_DIR="${PIPELINE_DIR}"
SAMPLE_ID="mock_pair"
SAMPLE_TYPE="tumor"
TUMOR_SAMPLE_ID="tumor01"
NORMAL_SAMPLE_ID="normal01"
CALLER_MODE="mutect2"
TEST_ROOT="${TMPDIR:-/tmp}/gaomei_wes_mutect2_aux_test"
DIR_VARIANTS="${TEST_ROOT}/variants"
DIR_ALIGNED="${TEST_ROOT}/aligned"
DIR_BQSR="${TEST_ROOT}/bqsr"
DIR_LOGS="${TEST_ROOT}/logs"
mkdir -p "${DIR_VARIANTS}" "${DIR_ALIGNED}" "${DIR_LOGS}"
printf 'mock raw vcf\n' > "${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.raw.vcf.gz"
printf 'mock f1r2\n' > "${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.f1r2.tar.gz"
printf 'mock tumor bam\n' > "${DIR_ALIGNED}/tumor01.dedup.bam"
printf 'mock normal bam\n' > "${DIR_ALIGNED}/normal01.dedup.bam"

TUMOR_BAM="${DIR_ALIGNED}/tumor01.dedup.bam"
NORMAL_BAM="${DIR_ALIGNED}/normal01.dedup.bam"
REFERENCE_GENOME="${PIPELINE_DIR}/testdata/mutect2/mock_reference.fa"
INTERVAL_FILE="${PIPELINE_DIR}/testdata/mutect2/mock_targets.bed"
MUTECT2_COMMON_VARIANTS_VCF="${PIPELINE_DIR}/testdata/mutect2/mock_common.vcf.gz"
RUN_MUTECT2_ORIENTATION_MODEL=true
RUN_MUTECT2_CONTAMINATION=true
MUTECT2_REQUIRE_AUXILIARY=true
SKIP_VARIANT_FILTER=false
TOOL_GATK="bash ${PIPELINE_DIR}/testdata/mutect2/mock_gatk.sh"
TOOL_BCFTOOLS="bash ${PIPELINE_DIR}/testdata/mutect2/mock_bcftools.sh"
GATK_JAVA_MEM="1g"

