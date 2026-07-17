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
    local runtime_project_dir="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
    local runtime_env_root="${ENV_ROOT:-${runtime_project_dir}/.conda_envs}"
    CONDA_BASE="${CONDA_BASE:-${HOME}/miniforge3}"
    MAIN_ENV_PREFIX="${MAIN_ENV_PREFIX:-${runtime_env_root}/big_wes_pipeline_env}"
    VEP_ENV_PREFIX="${VEP_ENV_PREFIX:-${runtime_env_root}/wes_vep_env}"
    SNPEFF_ENV_PREFIX="${SNPEFF_ENV_PREFIX:-${runtime_env_root}/wes_snpeff_env}"
    HLA_ENV_PREFIX="${HLA_ENV_PREFIX:-${runtime_env_root}/wes_hla_env}"
    HLA_TYPING_ENV_PREFIX="${HLA_TYPING_ENV_PREFIX:-${runtime_env_root}/wes_hla_typing_env}"
    CNV_ENV_PREFIX="${CNV_ENV_PREFIX:-${runtime_env_root}/wes_cnv_env}"
    SV_ENV_PREFIX="${SV_ENV_PREFIX:-${runtime_env_root}/wes_sv_env}"
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
    TOOL_CNVKIT="${TOOL_CNVKIT:-${PROJECT_DIR}/scripts/run_cnvkit_env.sh}"
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
    log_info "SAMPLE: ${SAMPLE_ID:-未设置} (${SAMPLE_TYPE:-未设置}, ${ANALYSIS_MODE:-未设置})"
    log_info "REFERENCE_DIR: ${REFERENCE_DIR:-未设置}"
    log_info "REFERENCE_GENOME: ${REFERENCE_GENOME:-未设置}"
    log_info "INTERVAL_FILE: ${INTERVAL_FILE:-未设置}"
    log_info "BQSR: RUN=${RUN_BQSR:-false}, SKIP=${SKIP_BQSR:-true}"
    log_info "HLA_TYPING: RUN=${RUN_HLA_TYPING:-false}, SKIP=${SKIP_HLA_TYPING:-false}"
    log_info "VEP: RUN=${RUN_VEP:-false}, SKIP=${SKIP_VEP:-true}"
    log_info "NEOANTIGEN: RUN=${RUN_NEOANTIGEN:-false}, SKIP=${SKIP_NEOANTIGEN:-true}"
    log_info "CNV/MSI/TMB: ${RUN_CNV:-false}/${RUN_MSI:-false}/${RUN_TMB:-false}"
}

#---------------------------------------
# 检查配置开关和关键资源是否自洽
#---------------------------------------
check_config_consistency() {
    local has_error=0
    log_step "检查配置一致性"

    if [ "${SKIP_ALIGN:-false}" = false ] && \
       { [ "${RUN_BQSR:-false}" != true ] || [ "${SKIP_BQSR:-true}" = true ]; }; then
        log_warn "当前样本会执行比对，但BQSR已关闭 (RUN_BQSR=${RUN_BQSR:-false}, SKIP_BQSR=${SKIP_BQSR:-true})"
    fi

    if [ "${RUN_BQSR:-false}" = true ] && [ "${SKIP_BQSR:-true}" = false ]; then
        if [ -z "${DBSNP_VCF:-}" ] || [ ! -f "${DBSNP_VCF}" ]; then
            log_error "BQSR已启用，但DBSNP_VCF不存在: ${DBSNP_VCF:-未配置}"
            has_error=1
        elif [ ! -f "${DBSNP_VCF_INDEX:-${DBSNP_VCF}.tbi}" ] && [ ! -f "${DBSNP_VCF}.idx" ]; then
            log_error "BQSR dbSNP缺少索引: ${DBSNP_VCF_INDEX:-${DBSNP_VCF}.tbi} 或 ${DBSNP_VCF}.idx"
            has_error=1
        fi
        local known_indels="${KNOWN_INDELS_VCF:-${MILLS_VCF:-}}"
        if [ -z "${known_indels}" ] || [ ! -f "${known_indels}" ]; then
            log_warn "BQSR known-indels未配置；建议设置 KNOWN_INDELS_VCF/MILLS_VCF"
        fi
    fi

    if [ "${RUN_VEP:-false}" = true ] && [ "${SKIP_VEP:-true}" = false ]; then
        if [ -z "${VEP_CACHE_DIR:-}" ] || [ ! -d "${VEP_CACHE_DIR}" ]; then
            log_error "VEP已启用，但VEP_CACHE_DIR不存在: ${VEP_CACHE_DIR:-未配置}"
            has_error=1
        fi
    fi

    if [ "${RUN_NEOANTIGEN:-false}" = true ] && [ "${SKIP_NEOANTIGEN:-true}" = false ]; then
        if [ -z "${NEOANTIGEN_PROTEIN_FASTA:-}" ] || [ ! -f "${NEOANTIGEN_PROTEIN_FASTA}" ]; then
            log_error "新抗原已启用，但蛋白FASTA不存在: ${NEOANTIGEN_PROTEIN_FASTA:-未配置}"
            has_error=1
        fi
    fi

    if [ "${RUN_HLA_TYPING:-false}" != false ] && [ "${SKIP_HLA_TYPING:-false}" = false ]; then
        if [ -z "${HLA_LA_GRAPH_DIR:-}" ] || [ ! -f "${HLA_LA_GRAPH_DIR}/serializedGRAPH" ]; then
            if [ "${HLA_TYPING_REQUIRED:-false}" = true ]; then
                log_error "HLA分型为必需，但graph未完成prepare: ${HLA_LA_GRAPH_DIR:-未配置}"
                has_error=1
            else
                log_warn "HLA*LA graph未完成prepare，auto模式将跳过: ${HLA_LA_GRAPH_DIR:-未配置}"
            fi
        fi
    fi

    if [ "${RUN_MSI:-false}" = true ] && [ "${SKIP_MSI:-true}" = false ] && \
       { [ -z "${MSISENSOR2_LIST:-}" ] || [ ! -f "${MSISENSOR2_LIST}" ]; } && \
       [ "${MSISENSOR_SCAN_REFERENCE:-false}" != true ]; then
        log_warn "MSI已启用但没有MSISENSOR2_LIST；当前只会生成smoke/NOT_RUN结果"
    fi

    if [ "${CALLER_MODE:-}" = mutect2 ] && [ "${SKIP_VARIANT_CALLING:-false}" = false ]; then
        [ -n "${GERMLINE_RESOURCE_VCF:-}" ] && [ -f "${GERMLINE_RESOURCE_VCF}" ] || \
            log_warn "Mutect2未配置AF-only germline resource"
        [ -n "${PANEL_OF_NORMALS:-}" ] && [ -f "${PANEL_OF_NORMALS}" ] || \
            log_warn "Mutect2未配置Panel of Normals"
        if [ "${RUN_MUTECT2_CONTAMINATION:-true}" = true ] && \
           { [ -z "${MUTECT2_COMMON_VARIANTS_VCF:-}" ] || [ ! -f "${MUTECT2_COMMON_VARIANTS_VCF}" ]; }; then
            if [ "${MUTECT2_REQUIRE_AUXILIARY:-false}" = true ]; then
                log_error "Mutect2污染估计为必需，但common-variants VCF缺失"
                has_error=1
            else
                log_warn "Mutect2 common-variants VCF缺失，污染估计将跳过"
            fi
        fi
        if [ -n "${MUTECT2_MAX_READS_PER_ALIGNMENT_START:-}" ] && \
           [ "${MUTECT2_MAX_READS_PER_ALIGNMENT_START}" -lt 20 ] 2>/dev/null; then
            log_warn "Mutect2 max-reads-per-alignment-start=${MUTECT2_MAX_READS_PER_ALIGNMENT_START}偏低，可能损失低VAF突变；正式分析建议使用50"
        fi
    fi

    if [ "${RUN_CNV:-false}" = true ] && [ "${SKIP_CNV:-true}" = false ]; then
        local has_cnv_baseline=false
        if [ -n "${CNVKIT_REFERENCE:-}" ] && [ -s "${CNVKIT_REFERENCE}" ]; then
            has_cnv_baseline=true
        elif [ -n "${NORMAL_BAM:-}" ] && bam_is_complete "${NORMAL_BAM}"; then
            has_cnv_baseline=true
        fi
        if [ "${has_cnv_baseline}" != true ]; then
            if [ "${CNV_REQUIRE_REFERENCE:-false}" = true ]; then
                log_error "CNV要求matched/pooled normal reference，但当前未找到"
                has_error=1
            else
                log_warn "CNV缺少matched/pooled normal reference；auto模式只会输出depth_qc，不会输出正式CNV calls"
            fi
        fi
    fi

    if [ "${RUN_TMB:-false}" = true ] && [ "${SKIP_TMB:-true}" = false ] && \
       [ "${TMB_DENOMINATOR_VALIDATED:-false}" != true ]; then
        log_warn "TMB分母BED尚未经方法学验证，结果应标记为研发用"
    fi

    if [ ${has_error} -ne 0 ]; then
        return 1
    fi
    log_info "配置一致性检查通过"
    return 0
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
    echo "║       WES分析流程 (WES Analysis Pipeline) v${PIPELINE_VERSION:-1.0.0}             ║"
    echo "║                                                                  ║"
    echo "║  FASTQ → QC → Trim → Align → BQSR → Call → Filter → Annotate  ║"
    echo "║  CNV → SV → MSI → Coverage → TMB → MultiQC → Report           ║"
    echo "╚══════════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
    echo "  样本ID:      ${SAMPLE_ID}"
    echo "  流程版本:    ${PIPELINE_VERSION:-1.0.0}"
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
    echo "    [CNV]  CNVkit matched/reference CNV；mosdepth depth QC"
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
bam_is_complete() {
    local bam="$1"
    [ -s "${bam}" ] || return 1
    "${TOOL_SAMTOOLS:-samtools}" quickcheck "${bam}" >/dev/null 2>&1
}

get_final_bam() {
    if [ -n "${TUMOR_BAM:-}" ] && bam_is_complete "${TUMOR_BAM}"; then
        echo "${TUMOR_BAM}"
        return 0
    fi

    local bqsr_bam="${DIR_BQSR}/${SAMPLE_ID}.bqsr.bam"
    local dedup_bam="${DIR_ALIGNED}/${SAMPLE_ID}.dedup.bam"
    local sorted_bam="${DIR_ALIGNED}/${SAMPLE_ID}.sorted.bam"

    if bam_is_complete "${bqsr_bam}"; then
        echo "${bqsr_bam}"
    elif bam_is_complete "${dedup_bam}"; then
        echo "${dedup_bam}"
    elif bam_is_complete "${sorted_bam}"; then
        echo "${sorted_bam}"
    else
        log_error "未找到完整且可读取的BAM文件!"
        return 1
    fi
}
