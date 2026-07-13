#!/bin/bash
# Minimal config for CNV/MSI smoke tests using the simulated 100-read BAM.

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${PROJECT_DIR}/testdata/variant_sim/config_variant_sim.sh"

RUN_CNV=true
SKIP_CNV=false
CNV_METHOD="depth"
CNVKIT_TARGET_BED="${SIM_DIR}/${SAMPLE_ID}_target.bed"
TOOL_MOSDEPTH="mosdepth"

RUN_MSI=true
SKIP_MSI=false
TOOL_MSISENSOR2="msisensor-pro"
MSISENSOR_MODE="smoke"
MSI_SMOKE_TEST=true
