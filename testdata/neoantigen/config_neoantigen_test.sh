#!/bin/bash
# Minimal config for testing only the neoantigen step.

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${PROJECT_DIR}/config.sh"

SAMPLE_ID="mock"
RESULT_DIR="${PROJECT_DIR}/testdata/neoantigen/out_step"
DIR_ANNOTATION="${RESULT_DIR}/annotation"
DIR_NEOANTIGEN="${RESULT_DIR}/neoantigen"
DIR_LOGS="${RESULT_DIR}/logs"

RUN_NEOANTIGEN=true
SKIP_NEOANTIGEN=false
NEOANTIGEN_VEP_VCF="${PROJECT_DIR}/testdata/neoantigen/mock_vep.vcf"
NEOANTIGEN_PROTEIN_FASTA="${PROJECT_DIR}/testdata/neoantigen/mock_proteins.fa"
NEOANTIGEN_PEPTIDE_LENGTHS="8,9"
NEOANTIGEN_PEPTIDE_FLANK=12
RUN_HLA_BINDING=true
HLA_ALLELES="HLA-A*02:01,HLA-B*07:02"
HLA_BINDING_TOOL="simple"
HLA_BINDING_PREDICTION_THRESHOLD_NM=500
