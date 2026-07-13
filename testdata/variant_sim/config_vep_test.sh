#!/bin/bash
# Minimal config for testing VEP 115 with local GRCh38 cache.

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${PROJECT_DIR}/testdata/variant_sim/config_variant_sim.sh"

RUN_VEP=true
SKIP_VEP=false
VEP_INPUT_VCF="${PROJECT_DIR}/testdata/variant_sim/results/variants/sim100.pass.vcf.gz"
VEP_CACHE_DIR="/Users/mac/Documents/wes/reference_data/vep_cache"
VEP_CACHE_VERSION="115"
VEP_FASTA="/Users/mac/Documents/wes/reference_data/hg38/Homo_sapiens_assembly38.fasta"
TOOL_VEP="${PROJECT_DIR}/scripts/run_vep_env.sh"
GATK_THREADS=1
