#!/bin/bash
# Minimal config for a 100-read local variant-calling smoke test.

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${PROJECT_DIR}/config.sh"

SAMPLE_ID="sim100"
SAMPLE_TYPE="germline"
CALLER_MODE="haplotypecaller"

SIM_DIR="${PROJECT_DIR}/testdata/variant_sim"
RESULT_DIR="${SIM_DIR}/results"

FASTQ_R1="${SIM_DIR}/${SAMPLE_ID}_R1.fastq.gz"
FASTQ_R2="${SIM_DIR}/${SAMPLE_ID}_R2.fastq.gz"

REFERENCE_GENOME="/Users/mac/Documents/wes/reference_data/hg38/Homo_sapiens_assembly38.fasta"
REFERENCE_BWA_INDEX="${REFERENCE_GENOME}"
REFERENCE_DICT="/Users/mac/Documents/wes/reference_data/hg38/Homo_sapiens_assembly38.dict"
REFERENCE_GENOME_VERSION="hg38"
INTERVAL_FILE="${SIM_DIR}/${SAMPLE_ID}_target.bed"

TOOL_FASTQC="fastqc"
TOOL_FASTP="fastp"
TOOL_BWA="bwa"
TOOL_SAMTOOLS="samtools"
TOOL_GATK="gatk"
TOOL_BCFTOOLS="bcftools"
TOOL_BEDTOOLS="bedtools"
TOOL_PICARD="picard"

BWA_THREADS=2
SAMTOOLS_THREADS=2
GATK_THREADS=1
GATK_JAVA_MEM="4g"
USE_FASTP=false
BWA_MEM_PARAMS="-M"

DIR_FASTQC="${RESULT_DIR}/fastqc"
DIR_TRIMMED="${RESULT_DIR}/trimmed"
DIR_ALIGNED="${RESULT_DIR}/aligned"
DIR_POSTQC="${RESULT_DIR}/post_align_qc"
DIR_BQSR="${RESULT_DIR}/bqsr"
DIR_VARIANTS="${RESULT_DIR}/variants"
DIR_ANNOTATION="${RESULT_DIR}/annotation"
DIR_CNV="${RESULT_DIR}/cnv"
DIR_SV="${RESULT_DIR}/sv"
DIR_MSI="${RESULT_DIR}/msi"
DIR_COVERAGE="${RESULT_DIR}/coverage"
DIR_TMB="${RESULT_DIR}/tmb"
DIR_NEOANTIGEN="${RESULT_DIR}/neoantigen"
DIR_SUMMARY="${RESULT_DIR}/summary"
DIR_MULTIQC="${RESULT_DIR}/multiqc"
DIR_LOGS="${SIM_DIR}/logs"

SKIP_FASTQC=true
SKIP_TRIM=true
SKIP_ALIGN=false
SKIP_SORT=false
SKIP_MARKDUP=true
SKIP_POSTQC=true
SKIP_BQSR=true
SKIP_VARIANT_CALLING=false
SKIP_VARIANT_FILTER=false
SKIP_SNPEFF=true
SKIP_VEP=true
SKIP_NEOANTIGEN=true
SKIP_CNV=true
SKIP_SV=true
SKIP_MSI=true
SKIP_COVERAGE=true
SKIP_TMB=true
SKIP_SUMMARY=true
SKIP_MULTIQC=true

RUN_BQSR=false
RUN_SNPEFF=false
RUN_VEP=false
RUN_NEOANTIGEN=false
RUN_CNV=false
RUN_SV=false
RUN_MSI=false
RUN_COVERAGE=false
RUN_TMB=false
