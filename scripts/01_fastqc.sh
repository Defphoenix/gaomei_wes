#!/bin/bash
#===============================================================================
# 01_fastqc.sh - 原始FASTQ数据质控
#
# 说明: 使用FastQC对原始测序数据进行质量评估
#       生成HTML报告，用于评估数据质量
#===============================================================================

set -euo pipefail

# 加载配置和工具函数
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/../config.sh}"
source "${CONFIG_FILE}"
source "${SCRIPT_DIR}/utils.sh"

#---------------------------------------
# 主函数
#---------------------------------------
main() {
    log_step "步骤1: FastQC 原始数据质控"

    # 检查是否跳过
    if [ "${SKIP_FASTQC}" = true ]; then
        log_warn "跳过 FastQC (SKIP_FASTQC=true)"
        return 0
    fi

    # 检查输入文件
    check_file "FASTQ R1" "${FASTQ_R1}" || return 1
    check_file "FASTQ R2" "${FASTQ_R2}" || return 1

    # 确保输出目录存在
    mkdir -p "${DIR_FASTQC}"

    # 检查是否已有结果 (避免重复运行)
    local r1_html="${DIR_FASTQC}/$(basename ${FASTQ_R1%.fastq.gz})_fastqc.html"
    local r2_html="${DIR_FASTQC}/$(basename ${FASTQ_R2%.fastq.gz})_fastqc.html"

    if [ -f "${r1_html}" ] && [ -f "${r2_html}" ]; then
        log_warn "FastQC结果已存在，跳过运行"
        log_info "R1报告: ${r1_html}"
        log_info "R2报告: ${r2_html}"
        return 0
    fi

    # 运行FastQC - R1
    run_cmd "FastQC R1 质控" \
        "${TOOL_FASTQC}" \
        --threads 4 \
        --outdir "${DIR_FASTQC}" \
        --extract \
        "${FASTQ_R1}"

    # 运行FastQC - R2
    run_cmd "FastQC R2 质控" \
        "${TOOL_FASTQC}" \
        --threads 4 \
        --outdir "${DIR_FASTQC}" \
        --extract \
        "${FASTQ_R2}"

    log_info "FastQC质控完成! 结果目录: ${DIR_FASTQC}"
    log_info "请查看HTML报告评估数据质量"
}

# 执行主函数
main "$@"
