#!/bin/bash
#===============================================================================
# utils.sh - 公共工具函数 (完整版)
#
# 说明: 提供日志记录、环境检查、步骤执行等通用功能
#       被 run_pipeline.sh 和各子脚本 source 引用
#       包含所有分析模块的工具检查
#===============================================================================

#---------------------------------------
# 配置默认值兜底
#---------------------------------------
apply_runtime_defaults() {
    CONDA_BASE="${CONDA_BASE:-/Users/mac/anaconda3}"
    MAIN_ENV_PREFIX="${MAIN_ENV_PREFIX:-${CONDA_BASE}/envs/big_wes_pipeline_env}"
    VEP_ENV_PREFIX="${VEP_ENV_PREFIX:-/Users/mac/Documents/wes/.conda_envs/wes_vep_env}"
    SNPEFF_ENV_PREFIX="${SNPEFF_ENV_PREFIX:-/Users/mac/Documents/wes/.conda_envs/wes_snpeff_env}"
    HLA_ENV_PREFIX="${HLA_ENV_PREFIX:-/Users/mac/Documents/wes/.conda_envs/wes_hla_env}"
    HLA_TYPING_ENV_PREFIX="${HLA_TYPING_ENV_PREFIX:-/Users/mac/Documents/wes/.conda_envs/wes_hla_typing_env}"
    CNV_ENV_PREFIX="${CNV_ENV_PREFIX:-/Users/mac/Documents/wes/.conda_envs/wes_cnv_env}"
    SV_ENV_PREFIX="${SV_ENV_PREFIX:-/Users/mac/Documents/wes/.conda_envs/wes_sv_env}"
    PIPELINE_EXTRA_PATHS="${PIPELINE_EXTRA_PATHS:-${MAIN_ENV_PREFIX}/bin:${VEP_ENV_PREFIX}/bin:${HLA_ENV_PREFIX}/bin:${HLA_TYPING_ENV_PREFIX}/bin:${CNV_ENV_PREFIX}/bin:${SV_ENV_PREFIX}/bin}"
    export PATH="${PIPELINE_EXTRA_PATHS}:${PATH}"
    export VEP_ENV="${VEP_ENV:-${VEP_ENV_PREFIX}}"
    export SNPEFF_ENV="${SNPEFF_ENV:-${SNPEFF_ENV_PREFIX}}"
    PIPELINE_JAVA_HOME="${PIPELINE_JAVA_HOME:-${MAIN_ENV_PREFIX}}"
    export JAVA_HOME="${PIPELINE_JAVA_HOME}"
    PIPELINE_LOCALE="${PIPELINE_LOCALE:-C}"
    export LC_ALL="${PIPELINE_LOCALE}"
    export LANG="${PIPELINE_LOCALE}"
    export LC_CTYPE="${PIPELINE_LOCALE}"

    TOOL_PYTHON="${TOOL_PYTHON:-python3}"
    TOOL_BGZIP="${TOOL_BGZIP:-bgzip}"
    TOOL_TABIX="${TOOL_TABIX:-tabix}"
    TOOL_SNPEFF="${TOOL_SNPEFF:-${PROJECT_DIR}/scripts/run_snpeff_env.sh}"
    TOOL_SNPSIFT="${TOOL_SNPSIFT:-${PROJECT_DIR}/scripts/run_snpsift_env.sh}"
    TOOL_VEP="${TOOL_VEP:-${PROJECT_DIR}/scripts/run_vep_env.sh}"
    TOOL_MANTA="${TOOL_MANTA:-configManta.py}"
}

apply_runtime_defaults

#---------------------------------------
# 颜色定义
#---------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m'

#---------------------------------------
# 日志函数
#---------------------------------------
log_info() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $*"
    echo -e "${GREEN}${msg}${NC}"
    echo "$msg" >> "${DIR_LOGS:-./logs}/pipeline.log"
}

log_warn() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $*"
    echo -e "${YELLOW}${msg}${NC}"
    echo "$msg" >> "${DIR_LOGS:-./logs}/pipeline.log"
}

log_error() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*"
    echo -e "${RED}${msg}${NC}" >&2
    echo "$msg" >> "${DIR_LOGS:-./logs}/pipeline.log"
}

log_step() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] [STEP] ========== $* =========="
    echo -e "${CYAN}${msg}${NC}"
    echo "$msg" >> "${DIR_LOGS:-./logs}/pipeline.log"
}

log_module() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] [MODULE] >>> $* <<<"
    echo -e "${MAGENTA}${msg}${NC}"
    echo "$msg" >> "${DIR_LOGS:-./logs}/pipeline.log"
}

#---------------------------------------
# 运行命令并记录日志
#---------------------------------------
run_cmd() {
    local desc="$1"
    shift
    local cmd="$*"

    log_info "执行: ${desc}"
    log_info "命令: ${cmd}"

    local start_time=$(date +%s)
    bash -o pipefail -c "${cmd}" 2>&1 | tee -a "${DIR_LOGS:-./logs}/pipeline.log"
    local exit_code=${PIPESTATUS[0]}
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    local minutes=$(( elapsed / 60 ))
    local seconds=$(( elapsed % 60 ))

    if [ ${exit_code} -eq 0 ]; then
        log_info "完成: ${desc} (耗时: ${minutes}分${seconds}秒)"
    else
        log_error "失败: ${desc} (退出码: ${exit_code}, 耗时: ${minutes}分${seconds}秒)"
        return ${exit_code}
    fi
}

#---------------------------------------
# 检查文件是否存在
#---------------------------------------
check_file() {
    local desc="$1"
    local filepath="$2"

    if [ ! -f "${filepath}" ]; then
        log_error "${desc} 不存在: ${filepath}"
        return 1
    fi
    log_info "检查通过: ${desc} -> ${filepath}"
    return 0
}

#---------------------------------------
# 检查工具是否可用
#---------------------------------------
check_tool() {
    local name="$1"
    local cmd="$2"

    if ! command -v "${cmd}" &>/dev/null; then
        log_error "工具未找到: ${name} (${cmd})，请检查安装和PATH配置"
        return 1
    fi
    local version
    version=$(${cmd} --version 2>&1 | head -1)
    log_info "工具检查通过: ${name} -> ${version}"
    return 0
}

#---------------------------------------
# 检查可选工具 (不存在只警告)
#---------------------------------------
check_tool_optional() {
    local name="$1"
    local cmd="$2"

    if ! command -v "${cmd}" &>/dev/null; then
        log_warn "可选工具未找到: ${name} (${cmd})，相关功能将被跳过"
        return 1
    fi
    local version
    version=$(${cmd} --version 2>&1 | head -1)
    log_info "可选工具通过: ${name} -> ${version}"
    return 0
}

#---------------------------------------
# 打印当前运行环境配置
#---------------------------------------
print_runtime_config() {
    log_step "运行环境配置"
    log_info "CONFIG_FILE: ${CONFIG_FILE:-未设置}"
    log_info "MAIN_ENV_PREFIX: ${MAIN_ENV_PREFIX:-未设置}"
    log_info "VEP_ENV_PREFIX: ${VEP_ENV_PREFIX:-未设置}"
    log_info "SNPEFF_ENV_PREFIX: ${SNPEFF_ENV_PREFIX:-未设置}"
    log_info "HLA_ENV_PREFIX: ${HLA_ENV_PREFIX:-未设置}"
    log_info "HLA_TYPING_ENV_PREFIX: ${HLA_TYPING_ENV_PREFIX:-未设置}"
    log_info "CNV_ENV_PREFIX: ${CNV_ENV_PREFIX:-未设置}"
    log_info "SV_ENV_PREFIX: ${SV_ENV_PREFIX:-未设置}"
    log_info "PIPELINE_EXTRA_PATHS: ${PIPELINE_EXTRA_PATHS:-未设置}"
    log_info "JAVA_HOME: ${JAVA_HOME:-未设置}"
    if command -v java >/dev/null 2>&1; then
        log_info "JAVA_VERSION: $(java -version 2>&1 | head -1)"
    fi
    log_info "TOOL_PYTHON: ${TOOL_PYTHON:-python3}"
    log_info "TOOL_VEP: ${TOOL_VEP:-vep}"
}

#---------------------------------------
# 检查所有必需工具 (核心流程)
#---------------------------------------
check_core_tools() {
    local has_error=0

    log_step "检查核心工具"

    check_tool "Python3" "${TOOL_PYTHON:-python3}" || has_error=1
    check_tool "FastQC" "${TOOL_FASTQC}" || has_error=1
    check_tool "BWA" "${TOOL_BWA}" || has_error=1
    check_tool "Samtools" "${TOOL_SAMTOOLS}" || has_error=1
    check_tool "GATK" "${TOOL_GATK}" || has_error=1
    check_tool "BEDtools" "${TOOL_BEDTOOLS}" || has_error=1
    check_tool "BCFtools" "${TOOL_BCFTOOLS}" || has_error=1
    check_tool_optional "bgzip" "${TOOL_BGZIP:-bgzip}" || true
    check_tool_optional "tabix" "${TOOL_TABIX:-tabix}" || true

    # 修剪工具
    if [ "${USE_FASTP}" = true ]; then
        check_tool "fastp" "${TOOL_FASTP}" || has_error=1
    else
        check_tool "Trimmomatic" "${TOOL_TRIMMOMATIC}" || has_error=1
    fi

    # Picard
    if ! check_tool_optional "Picard" "${TOOL_PICARD}"; then
        log_warn "Picard 未找到，将使用 samtools markdup 替代"
    fi

    if [ ${has_error} -ne 0 ]; then
        log_error "核心工具缺失，请安装后重试"
        return 1
    fi

    log_info "核心工具检查通过!"
    return 0
}

#---------------------------------------
# 检查注释与分析工具 (可选)
#---------------------------------------
check_analysis_tools() {
    log_step "检查分析工具 (可选模块)"

    # SnpEff
    if [ "${RUN_SNPEFF}" = true ] && [ "${SKIP_SNPEFF}" = false ]; then
        check_tool_optional "SnpEff" "${TOOL_SNPEFF}" || log_warn "SnpEff不可用，注释步骤将跳过"
    fi

    # VEP
    if [ "${RUN_VEP}" = true ] && [ "${SKIP_VEP}" = false ]; then
        check_tool_optional "VEP" "${TOOL_VEP}" || log_warn "VEP不可用，注释步骤将跳过"
    fi

    # CNV
    if [ "${RUN_CNV}" = true ] && [ "${SKIP_CNV}" = false ]; then
        if ! check_tool_optional "CNVkit" "${TOOL_CNVKIT}"; then
            check_tool_optional "mosdepth" "${TOOL_MOSDEPTH}" || log_warn "CNVkit和mosdepth均不可用，CNV分析将失败"
        fi
    fi

    # Manta
    if [ "${RUN_SV}" = true ] && [ "${SKIP_SV}" = false ]; then
        check_tool_optional "Manta" "python3" || log_warn "Manta需要Python3环境"
        if [ -f "${TOOL_MANTA}" ]; then
            log_info "Manta配置: ${TOOL_MANTA}"
        else
            log_warn "Manta脚本未找到: ${TOOL_MANTA}，SV分析将跳过"
        fi
    fi

    # MSIsensor-pro
    if [ "${RUN_MSI}" = true ] && [ "${SKIP_MSI}" = false ]; then
        check_tool_optional "MSIsensor-pro" "${TOOL_MSISENSOR2}" || log_warn "MSIsensor-pro不可用，MSI检测将跳过"
    fi

    # mosdepth
    if [ "${RUN_COVERAGE}" = true ] && [ "${SKIP_COVERAGE}" = false ]; then
        check_tool_optional "mosdepth" "${TOOL_MOSDEPTH}" || log_warn "mosdepth不可用，覆盖度分析将跳过"
    fi

    # MultiQC
    if [ "${SKIP_MULTIQC}" = false ]; then
        check_tool_optional "MultiQC" "${TOOL_MULTIQC}" || log_warn "MultiQC不可用，汇总QC报告将跳过"
    fi

    # SnpSift
    check_tool_optional "SnpSift" "${TOOL_SNPSIFT}" || true

    # Qualimap
    check_tool_optional "Qualimap" "${TOOL_QUALIMAP}" || true

    # HLA binding predictors
    if [ "${RUN_HLA_BINDING:-false}" != false ]; then
        if command -v "${TOOL_NETMHCPAN:-netMHCpan}" &>/dev/null; then
            check_tool_optional "netMHCpan" "${TOOL_NETMHCPAN:-netMHCpan}" || true
        elif command -v "${TOOL_MHCFLURRY:-mhcflurry-predict}" &>/dev/null; then
            check_tool_optional "MHCflurry" "${TOOL_MHCFLURRY:-mhcflurry-predict}" || true
        else
            log_warn "netMHCpan和mhcflurry-predict均不可用；auto模式会停止，simple仅可显式用于连通性测试"
        fi
    fi

    if [ "${RUN_HLA_TYPING:-false}" != false ]; then
        check_tool_optional "HLA*LA" "${TOOL_HLA_LA:-HLA-LA.pl}" || \
            log_warn "HLA*LA不可用；RUN_HLA_TYPING=auto时将跳过，required模式会终止"
    fi

    return 0
}

#---------------------------------------
# 检查所有工具 (核心 + 分析)
#---------------------------------------
check_all_tools() {
    check_core_tools || return 1
    check_analysis_tools
    return 0
}

#---------------------------------------
# 检查参考文件
#---------------------------------------
check_reference() {
    local has_error=0

    log_step "检查参考文件"

    check_file "参考基因组" "${REFERENCE_GENOME}" || has_error=1
    check_file "参考基因组索引" "${REFERENCE_GENOME}.fai" || has_error=1
    check_file "参考基因组dict" "${REFERENCE_DICT}" || has_error=1

    if [ ${has_error} -ne 0 ]; then
        log_error "参考文件缺失，请检查 config.sh 中的路径配置"
        return 1
    fi

    # 可选参考文件检查
    if [ -n "${DBSNP_VCF}" ] && [ "${DBSNP_VCF}" != "/path/to/dbsnp_146.hg38.vcf.gz" ]; then
        check_file "dbSNP" "${DBSNP_VCF}" || log_warn "dbSNP未配置，BQSR和注释可能受影响"
    fi

    log_info "参考文件检查通过!"
    return 0
}

#---------------------------------------
# 检查输入数据
#---------------------------------------
check_input() {
    local has_error=0

    log_step "检查输入数据"

    check_file "FASTQ R1" "${FASTQ_R1}" || has_error=1
    check_file "FASTQ R2" "${FASTQ_R2}" || has_error=1

    if [ ${has_error} -ne 0 ]; then
        log_error "输入数据缺失，请检查 config.sh 中的FASTQ路径"
        return 1
    fi

    log_info "输入数据检查通过!"
    return 0
}

#---------------------------------------
# 创建输出目录
#---------------------------------------
create_output_dirs() {
    log_step "创建输出目录"

    local dirs=(
        "${DIR_FASTQC}"
        "${DIR_TRIMMED}"
        "${DIR_ALIGNED}"
        "${DIR_POSTQC}"
        "${DIR_BQSR}"
        "${DIR_VARIANTS}"
        "${DIR_ANNOTATION}"
        "${DIR_CNV}"
        "${DIR_SV}"
        "${DIR_MSI}"
        "${DIR_COVERAGE}"
        "${DIR_TMB}"
        "${DIR_NEOANTIGEN:-${RESULT_DIR}/neoantigen}"
        "${DIR_SUMMARY}"
        "${DIR_MULTIQC}"
        "${DIR_LOGS}"
    )

    for dir in "${dirs[@]}"; do
        mkdir -p "${dir}"
        log_info "目录: ${dir}"
    done

    log_info "输出目录创建完成"
}

#---------------------------------------
# 打印流程横幅
#---------------------------------------
print_banner() {
    echo -e "${BLUE}"
    echo "╔══════════════════════════════════════════════════════════════════╗"
    echo "║       突变分析流程 (Mutation Analysis Pipeline) - 完整版       ║"
    echo "║                                                                  ║"
    echo "║  FASTQ → QC → Trim → Align → BQSR → Call → Filter → Annotate  ║"
    echo "║  CNV → SV → MSI → Coverage → TMB → MultiQC → Report           ║"
    echo "╚══════════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
    echo "  样本ID:      ${SAMPLE_ID}"
    echo "  样本类型:    ${SAMPLE_TYPE}"
    echo "  检测模式:    ${CALLER_MODE}"
    echo "  参考基因组:  ${REFERENCE_GENOME}"
    echo "  项目目录:    ${PROJECT_DIR}"
    echo "  开始时间:    $(date '+%Y-%m-%d %H:%M:%S')"
    echo "──────────────────────────────────────────────────────────────────"
    echo ""
    echo "  分析模块:"
    echo "    [质控] FastQC + fastp/Trimmomatic"
    echo "    [比对] BWA-MEM + Samtools + Picard"
    echo "    [校准] GATK BQSR"
    echo "    [变异] GATK ${CALLER_MODE}"
    echo "    [注释] SnpEff + VEP"
    echo "    [免疫] 新抗原候选肽 + HLA结合预测"
    echo "    [CNV]  CNVkit / mosdepth depth-ratio"
    echo "    [SV]   Manta"
    echo "    [MSI]  MSIsensor-pro"
    echo "    [覆盖] mosdepth + bedtools"
    echo "    [TMB]  自定义计算"
    echo "    [汇总] MultiQC"
    echo ""
}

#---------------------------------------
# 打印流程结束信息
#---------------------------------------
print_finish() {
    local start_ts="$1"
    local end_ts=$(date +%s)
    local elapsed=$(( end_ts - start_ts ))
    local hours=$(( elapsed / 3600 ))
    local minutes=$(( (elapsed % 3600) / 60 ))

    echo ""
    echo -e "${GREEN}╔══════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║                    流程运行完成!                                ║${NC}"
    echo -e "${GREEN}╚══════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    echo "  总耗时:   ${hours}小时${minutes}分钟"
    echo "  结果目录: ${RESULT_DIR}"
    echo "  日志文件: ${DIR_LOGS}/pipeline.log"
    echo "  QC报告:   ${DIR_MULTIQC}/multiqc_report.html"
    echo "  完成时间: $(date '+%Y-%m-%d %H:%M:%S')"
    echo ""
    echo "  主要结果文件:"
    echo "    变异VCF:    ${DIR_ANNOTATION}/${SAMPLE_ID}*.annotated.vcf.gz"
    echo "    CNV结果:    ${DIR_CNV}/${SAMPLE_ID}*.cnr / ${DIR_CNV}/${SAMPLE_ID}*.cns"
    echo "    SV结果:     ${DIR_SV}/results/variants/*.vcf.gz"
    echo "    MSI结果:    ${DIR_MSI}/${SAMPLE_ID}_msi_result.txt"
    echo "    覆盖度:     ${DIR_COVERAGE}/${SAMPLE_ID}*.mosdepth.summary.txt"
    echo "    TMB:        ${DIR_TMB}/${SAMPLE_ID}_tmb_result.txt"
    echo "    新抗原:     ${DIR_NEOANTIGEN:-${RESULT_DIR}/neoantigen}/${SAMPLE_ID}_neoantigen_peptides.fa"
    echo ""
}

#---------------------------------------
# 检查磁盘空间 (单位: GB)
#---------------------------------------
check_disk_space() {
    local path="$1"
    local required_gb="${2:-10}"

    local available_kb
    available_kb=$(df -k "${path}" | tail -1 | awk '{print $4}')
    local available_gb=$(( available_kb / 1024 / 1024 ))

    if [ ${available_gb} -lt ${required_gb} ]; then
        log_warn "磁盘空间不足! 需要约 ${required_gb}GB, 可用 ${available_gb}GB (${path})"
        return 1
    fi

    log_info "磁盘空间充足: 可用 ${available_gb}GB (${path})"
    return 0
}

#---------------------------------------
# 获取最终BAM路径 (自动判断是否经过BQSR)
#---------------------------------------
get_final_bam() {
    if [ -n "${TUMOR_BAM:-}" ] && [ -f "${TUMOR_BAM}" ]; then
        echo "${TUMOR_BAM}"
        return 0
    fi

    local bqsr_bam="${DIR_BQSR}/${SAMPLE_ID}.bqsr.bam"
    local dedup_bam="${DIR_ALIGNED}/${SAMPLE_ID}.dedup.bam"
    local sorted_bam="${DIR_ALIGNED}/${SAMPLE_ID}.sorted.bam"

    if [ -f "${bqsr_bam}" ]; then
        echo "${bqsr_bam}"
    elif [ -f "${dedup_bam}" ]; then
        echo "${dedup_bam}"
    elif [ -f "${sorted_bam}" ]; then
        echo "${sorted_bam}"
    else
        log_error "未找到可用的BAM文件!"
        return 1
    fi
}
