#!/bin/bash
#===============================================================================
# config_test.sh - 测试用配置文件 (chr21小片段)
#
# 说明: 用于小规模测试的配置文件
#       使用chr21:9000000-12000000区域 (3Mb) 进行测试
#       下载完数据后，复制此文件覆盖config.sh即可测试
#
# 用法:
#   1. 先运行: bash prepare_data.sh all
#   2. 复制:   cp config_test.sh config.sh
#   3. 运行:   bash run_pipeline.sh check
#   4. 运行:   bash run_pipeline.sh
#===============================================================================

#---------------------------------------
# 项目与样本 (测试用)
#---------------------------------------
PROJECT_NAME="mutation_pipeline_test"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLE_ID="test_sample"
SAMPLE_TYPE="germline"                      # 测试用胚系样本

NORMAL_SAMPLE_ID=""
NORMAL_BAM=""

#---------------------------------------
# 原始数据 (测试数据)
#---------------------------------------
RAW_DATA_DIR="${PROJECT_DIR}/testdata"
FASTQ_R1="${RAW_DATA_DIR}/test_R1.fastq.gz"
FASTQ_R2="${RAW_DATA_DIR}/test_R2.fastq.gz"

#---------------------------------------
# 参考基因组 (本地下载路径)
#---------------------------------------
REFERENCE_GENOME="${PROJECT_DIR}/reference/hg38.fa"
REFERENCE_DICT="${PROJECT_DIR}/reference/hg38.dict"
REFERENCE_BWA_INDEX="${PROJECT_DIR}/reference/hg38.fa"
REFERENCE_GENOME_VERSION="hg38"

#---------------------------------------
# 数据库 (本地下载路径)
#---------------------------------------
DBSNP_VCF="${PROJECT_DIR}/database/dbsnp_146.hg38.vcf.gz"
DBSNP_VCF_INDEX="${PROJECT_DIR}/database/dbsnp_146.hg38.vcf.gz.tbi"
HAPMAP_VCF="${PROJECT_DIR}/database/hapmap_3.3.hg38.vcf.gz"
OMNI_VCF="${PROJECT_DIR}/database/1000G_omni2.5.hg38.vcf.gz"
THOUSAND_G_VCF="${PROJECT_DIR}/database/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
MILLS_VCF="${PROJECT_DIR}/database/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN_INDELS_VCF="${PROJECT_DIR}/database/Homo_sapiens_assembly38.known_indels.vcf.gz"

#---------------------------------------
# 区间文件 (chr21测试区域)
#---------------------------------------
INTERVAL_FILE="${PROJECT_DIR}/testdata/test_target.bed"
INTERVAL_LIST=""
EXOME_KIT="test_chr21"

#---------------------------------------
# SnpEff
#---------------------------------------
SNPEFF_ENV_PREFIX="${SNPEFF_ENV_PREFIX:-${PROJECT_DIR}/.conda_envs/wes_snpeff_env}"
export SNPEFF_ENV="${SNPEFF_ENV:-${SNPEFF_ENV_PREFIX}}"
SNPEFF_DB="GRCh38.105"
SNPEFF_DATA_DIR="${PROJECT_DIR}/snpeff_db"
SNPEFF_CONFIG="${PROJECT_DIR}/snpeff_db/snpEff.config"

#---------------------------------------
# VEP
#---------------------------------------
VEP_CACHE_DIR="${PROJECT_DIR}/vep_cache"
VEP_SPECIES="homo_sapiens"
VEP_ASSEMBLY="GRCh38"
VEP_CACHE_VERSION="110"

#---------------------------------------
# CNVkit
#---------------------------------------
CNVKIT_REFERENCE=""
CNVKIT_TARGET_BED="${PROJECT_DIR}/testdata/test_target.bed"
CNVKIT_ANTITARGET_BED=""

#---------------------------------------
# MSIsensor2
#---------------------------------------
MSISENSOR2_LIST=""

#---------------------------------------
# Manta
#---------------------------------------
MANTA_CONFIG=""

#---------------------------------------
# 工具路径
#---------------------------------------
TOOL_FASTQC="fastqc"
TOOL_TRIMMOMATIC="trimmomatic"
TOOL_FASTP="fastp"
TOOL_MULTIQC="multiqc"
TOOL_BWA="bwa"
TOOL_SAMTOOLS="samtools"
TOOL_PICARD="picard"
TOOL_BEDTOOLS="bedtools"
TOOL_GATK="gatk"
TOOL_BCFTOOLS="bcftools"
TOOL_SNPEFF="${PROJECT_DIR}/scripts/run_snpeff_env.sh"
TOOL_VEP="vep"
TOOL_SNPSIFT="${PROJECT_DIR}/scripts/run_snpsift_env.sh"
TOOL_CNVKIT="cnvkit"
TOOL_MANTA="/path/to/manta/bin/configManta.py"
TOOL_MSISENSOR2="msisensor2"
TOOL_MOSDEPTH="mosdepth"
TOOL_QUALIMAP="qualimap"

#---------------------------------------
# 分析参数 (测试用较小参数加速)
#---------------------------------------
TRIM_QUALITY=20
TRIM_MIN_LENGTH=30
TRIM_ADAPTER="AGATCGGAAGAGC"
USE_FASTP=true

BWA_THREADS=4
BWA_MEM_PARAMS="-M"
SAMTOOLS_THREADS=2

GATK_JAVA_MEM="8g"
GATK_THREADS=2
CALLER_MODE="haplotypecaller"

PANEL_OF_NORMALS=""
PON_INDEX=""
MUTECT2_LEARNING_DATA=""

RUN_BQSR=false                            # 测试时可关闭BQSR加速

HC_SNP_QUAL_SCORE="QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
HC_INDEL_QUAL_SCORE="QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"

RUN_SNPEFF=true
SNPEFF_EXTRA_PARAMS=""
RUN_VEP=true
VEP_PLUGINS=""
VEP_EXTRA_PARAMS=""

RUN_CNV=false                             # 测试时可关闭CNV
RUN_SV=false                              # 测试时可关闭SV
RUN_MSI=false                             # 测试时可关闭MSI
RUN_COVERAGE=true
MOSDEPTH_THREADS=2
COVERAGE_TARGETS="10,30,50,100"
RUN_TMB=true
TMB_CODING_SIZE=3                         # 测试区域约3Mb
TMB_MIN_QUAL=20
TMB_MIN_AF=0.05
TMB_INCLUDE_TYPES="missense,nonsense,frameshift,inframe,splice"

#---------------------------------------
# 输出目录
#---------------------------------------
RESULT_DIR="${PROJECT_DIR}/results"
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
DIR_SUMMARY="${RESULT_DIR}/summary"
DIR_MULTIQC="${RESULT_DIR}/multiqc"
DIR_LOGS="${PROJECT_DIR}/logs"

#---------------------------------------
# 运行控制 (测试模式: 关闭不需要的模块)
#---------------------------------------
SKIP_FASTQC=false
SKIP_TRIM=false
SKIP_ALIGN=false
SKIP_SORT=false
SKIP_MARKDUP=false
SKIP_POSTQC=false
SKIP_BQSR=true                            # 测试跳过BQSR
SKIP_VARIANT_CALLING=false
SKIP_VARIANT_FILTER=false
SKIP_SNPEFF=true                          # 测试先跳过注释
SKIP_VEP=true                             # 测试先跳过注释
SKIP_CNV=true                             # 测试跳过CNV
SKIP_SV=true                              # 测试跳过SV
SKIP_MSI=true                             # 测试跳过MSI
SKIP_COVERAGE=false
SKIP_TMB=false
SKIP_SUMMARY=false
SKIP_MULTIQC=true                         # 测试跳过MultiQC

CLEAN_INTERMEDIATE=false
