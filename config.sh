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
ANALYSIS_MODE="${ANALYSIS_MODE:-single}"     # single / tumor_normal / multi

NORMAL_SAMPLE_ID=""
NORMAL_BAM=""
TUMOR_SAMPLE_ID="${TUMOR_SAMPLE_ID:-${SAMPLE_ID}}"
TUMOR_BAM="${TUMOR_BAM:-}"

#---------------------------------------
# 原始数据 (测试数据)
#---------------------------------------
RAW_DATA_DIR="/PUBLIC/gomics/guofenghua/ns_workflow/test-datasets-modules/data/genomics/homo_sapiens/illumina/fastq"
FASTQ_R1="${RAW_DATA_DIR}/test_1.fastq.gz"
FASTQ_R2="${RAW_DATA_DIR}/test_2.fastq.gz"

#---------------------------------------
# 参考基因组 (本地下载路径)
#---------------------------------------
REFERENCE_GENOME="${REFERENCE_GENOME:-/Users/mac/Documents/wes/reference_data/hg38/Homo_sapiens_assembly38.fasta}"
REFERENCE_DICT="${REFERENCE_DICT:-/Users/mac/Documents/wes/reference_data/hg38/Homo_sapiens_assembly38.dict}"
REFERENCE_BWA_INDEX="${REFERENCE_BWA_INDEX:-${REFERENCE_GENOME}}"
REFERENCE_GENOME_VERSION="hg38"

#---------------------------------------
# 数据库 (本地下载路径)
#---------------------------------------
DBSNP_VCF="${DBSNP_VCF:-/Users/mac/Documents/wes/reference_data/dbsnp_146.hg38.vcf.gz}"
DBSNP_VCF_INDEX="${DBSNP_VCF_INDEX:-/Users/mac/Documents/wes/reference_data/dbsnp_146.hg38.vcf.gz.tbi}"
MILLS_VCF="${MILLS_VCF:-/Users/mac/Documents/wes/reference_data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz}"
THOUSAND_G_VCF="${THOUSAND_G_VCF:-/Users/mac/Documents/wes/reference_data/1000G_phase1.snps.high_confidence.hg38.vcf.gz}"

#---------------------------------------
# 区间文件 (chr21测试区域)
#---------------------------------------
INTERVAL_FILE="/PUBLIC/gomics/guofenghua/ns_workflow/test-datasets-modules/data/genomics/homo_sapiens/genome/chr21/sequence/multi_intervals.bed"
CNVKIT_TARGET_BED="/PUBLIC/gomics/guofenghua/ns_workflow/test-datasets-modules/data/genomics/homo_sapiens/genome/chr21/sequence/multi_intervals.bed"

#---------------------------------------
# SnpEff
#---------------------------------------
SNPEFF_DB="GRCh38.105"
SNPEFF_DATA_DIR="${PROJECT_DIR}/snpeff_db"
SNPEFF_CONFIG="${PROJECT_DIR}/snpeff_db/snpEff.config"

#---------------------------------------
# VEP
#---------------------------------------
VEP_INPUT_VCF="${VEP_INPUT_VCF:-}"         # 留空时按 CALLER_MODE 自动选择过滤VCF
VEP_CACHE_DIR="${VEP_CACHE_DIR:-/Users/mac/Documents/wes/reference_data/vep_cache}"
VEP_SPECIES="homo_sapiens"
VEP_ASSEMBLY="GRCh38"
VEP_CACHE_VERSION="115"
VEP_FASTA="${VEP_FASTA:-${REFERENCE_GENOME}}"

#---------------------------------------
# CNVkit
#---------------------------------------
CNVKIT_REFERENCE=""
CNVKIT_TARGET_BED="${PROJECT_DIR}/testdata/test_target.bed"
CNVKIT_ANTITARGET_BED=""
CNV_METHOD="auto"                         # auto/cnvkit/depth
CNV_NEUTRAL_LOW=0.75
CNV_NEUTRAL_HIGH=1.25
CNVKIT_METHOD="hybrid"
CNVKIT_SEGMENT_PARAMS=""
CNVKIT_CALL_PARAMS=""

#---------------------------------------
# MSIsensor-pro
#---------------------------------------
MSISENSOR2_LIST=""
MSISENSOR_BASELINE=""
MSISENSOR_MODE="auto"                     # auto/paired/tumor_only/smoke
MSISENSOR_THREADS=2
MSI_MIN_LENGTH=10
MSI_MIN_HOMOPOLYMER=1
MSI_MIN_COVERAGE=15
MSI_SMOKE_TEST=true                       # 无normal/baseline时生成可读summary，不冒充正式MSI结果

#---------------------------------------
# Manta
#---------------------------------------
MANTA_CONFIG=""
MANTA_THREADS=2

#---------------------------------------
# 运行环境与工具路径
# 说明:
#   1. 优先在这里配置服务器上的 conda 环境和软件路径。
#   2. 如果工具在 PATH 中，保留工具名即可；如果服务器路径不同，改成绝对路径。
#   3. ${VAR:-default} 写法允许在命令行或外层调度系统中临时覆盖。
#---------------------------------------

CONDA_BASE="${CONDA_BASE:-/Users/mac/anaconda3}"
MAIN_ENV_PREFIX="${MAIN_ENV_PREFIX:-${CONDA_BASE}/envs/big_wes_pipeline_env}"
VEP_ENV_PREFIX="${VEP_ENV_PREFIX:-/Users/mac/Documents/wes/.conda_envs/wes_vep_env}"
HLA_ENV_PREFIX="${HLA_ENV_PREFIX:-/Users/mac/Documents/wes/.conda_envs/wes_hla_env}"

PIPELINE_EXTRA_PATHS="${PIPELINE_EXTRA_PATHS:-${MAIN_ENV_PREFIX}/bin:${VEP_ENV_PREFIX}/bin:${HLA_ENV_PREFIX}/bin}"
export PATH="${PIPELINE_EXTRA_PATHS}:${PATH}"
export VEP_ENV="${VEP_ENV:-${VEP_ENV_PREFIX}}"
PIPELINE_JAVA_HOME="${PIPELINE_JAVA_HOME:-${MAIN_ENV_PREFIX}}"
export JAVA_HOME="${PIPELINE_JAVA_HOME}"
PIPELINE_LOCALE="${PIPELINE_LOCALE:-C}"
export LC_ALL="${PIPELINE_LOCALE}"
export LANG="${PIPELINE_LOCALE}"
export LC_CTYPE="${PIPELINE_LOCALE}"
MIN_DISK_GB="${MIN_DISK_GB:-100}"

TOOL_PYTHON="${TOOL_PYTHON:-python3}"
TOOL_FASTQC="${TOOL_FASTQC:-fastqc}"
TOOL_TRIMMOMATIC="${TOOL_TRIMMOMATIC:-trimmomatic}"
TOOL_FASTP="${TOOL_FASTP:-fastp}"
TOOL_MULTIQC="${TOOL_MULTIQC:-multiqc}"
TOOL_BWA="${TOOL_BWA:-bwa}"
TOOL_SAMTOOLS="${TOOL_SAMTOOLS:-samtools}"
TOOL_PICARD="${TOOL_PICARD:-picard}"
TOOL_BEDTOOLS="${TOOL_BEDTOOLS:-bedtools}"
TOOL_GATK="${TOOL_GATK:-gatk}"
TOOL_BCFTOOLS="${TOOL_BCFTOOLS:-bcftools}"
TOOL_BGZIP="${TOOL_BGZIP:-bgzip}"
TOOL_TABIX="${TOOL_TABIX:-tabix}"
TOOL_SNPEFF="${TOOL_SNPEFF:-snpEff}"
TOOL_VEP="${TOOL_VEP:-${PROJECT_DIR}/scripts/run_vep_env.sh}"
TOOL_NETMHCPAN="${TOOL_NETMHCPAN:-netMHCpan}"
TOOL_MHCFLURRY="${TOOL_MHCFLURRY:-mhcflurry-predict}"
TOOL_SNPSIFT="${TOOL_SNPSIFT:-SnpSift}"
TOOL_CNVKIT="${TOOL_CNVKIT:-cnvkit.py}"
TOOL_MANTA="${TOOL_MANTA:-configManta.py}"
TOOL_MSISENSOR2="${TOOL_MSISENSOR2:-msisensor-pro}"
TOOL_MOSDEPTH="${TOOL_MOSDEPTH:-mosdepth}"
TOOL_QUALIMAP="${TOOL_QUALIMAP:-qualimap}"

CLINVAR_VCF="${CLINVAR_VCF:-}"
COSMIC_VCF="${COSMIC_VCF:-}"

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
GERMLINE_RESOURCE_VCF="${GERMLINE_RESOURCE_VCF:-}"
GERMLINE_RESOURCE_INDEX="${GERMLINE_RESOURCE_INDEX:-}"
MUTECT2_CONTAMINATION_TABLE="${MUTECT2_CONTAMINATION_TABLE:-}"
MUTECT2_EXTRA_PARAMS="${MUTECT2_EXTRA_PARAMS:-}"
FILTER_MUTECT_EXTRA_PARAMS="${FILTER_MUTECT_EXTRA_PARAMS:-}"

RUN_BQSR=false                            # 测试时可关闭BQSR加速

HC_SNP_QUAL_SCORE="QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
HC_INDEL_QUAL_SCORE="QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"

RUN_SNPEFF=true
SNPEFF_EXTRA_PARAMS=""
RUN_VEP=true
VEP_PLUGINS=""
VEP_EXTRA_PARAMS=""

#---------------------------------------
# 新抗原候选肽
#---------------------------------------
RUN_NEOANTIGEN=true
NEOANTIGEN_VEP_VCF=""                     # 留空时默认使用 ${DIR_ANNOTATION}/${SAMPLE_ID}.vep.vcf.gz
NEOANTIGEN_PROTEIN_FASTA="${PROJECT_DIR}/reference/protein.fa"
NEOANTIGEN_PEPTIDE_LENGTHS="8,9,10,11"    # 历史变量名；现在表示突变位点上下游flank aa列表。8 => 8+1+8=17aa窗口
NEOANTIGEN_PEPTIDE_FLANK=30               # 保留给扩展/兜底窗口
HLA_ALLELES=""                            # 例如: HLA-A*02:01,HLA-B*07:02
RUN_HLA_BINDING=false                     # 需要本地安装 netMHCpan 并配置 HLA_ALLELES
HLA_BINDING_TOOL="auto"                   # auto/netmhcpan/mhcflurry/simple
HLA_BINDING_PREDICTION_THRESHOLD_NM=500

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
DIR_NEOANTIGEN="${RESULT_DIR}/neoantigen"
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
SKIP_NEOANTIGEN=true                      # 默认跳过；VEP完成后可单步运行 step 7d
SKIP_CNV=true                             # 测试跳过CNV
SKIP_SV=true                              # 测试跳过SV
SKIP_MSI=true                             # 测试跳过MSI
SKIP_COVERAGE=false
SKIP_TMB=false
SKIP_SUMMARY=false
SKIP_MULTIQC=true                         # 测试跳过MultiQC

CLEAN_INTERMEDIATE=false
