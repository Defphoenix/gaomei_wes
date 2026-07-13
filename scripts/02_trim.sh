#!/bin/bash
#===============================================================================
# 02_trim.sh - 数据修剪 (Adapter去除 & 低质量过滤)
#
# 说明: 支持两种修剪工具:
#       1. fastp (推荐) - 快速一体化质控和修剪
#       2. Trimmomatic  - 经典修剪工具
#       通过 config.sh 中 USE_FASTP 参数切换
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/../config.sh}"
source "${CONFIG_FILE}"
source "${SCRIPT_DIR}/utils.sh"

#---------------------------------------
# fastp 修剪
#---------------------------------------
run_fastp() {
    local output_r1="${DIR_TRIMMED}/${SAMPLE_ID}_trimmed_R1.fastq.gz"
    local output_r2="${DIR_TRIMMED}/${SAMPLE_ID}_trimmed_R2.fastq.gz"
    local html_report="${DIR_TRIMMED}/${SAMPLE_ID}_fastp.html"
    local json_report="${DIR_TRIMMED}/${SAMPLE_ID}_fastp.json"

    run_cmd "fastp 数据修剪和质控" \
        "${TOOL_FASTP}" \
        --in1 ${FASTQ_R1} \
        --in2 ${FASTQ_R2} \
        --out1 ${output_r1} \
        --out2 ${output_r2} \
        --html ${html_report} \
        --json ${json_report} \
        --thread ${BWA_THREADS} \
        --qualified_quality_phred ${TRIM_QUALITY} \
        --length_required ${TRIM_MIN_LENGTH} \
        --detect_adapter_for_pe \
        --adapter_sequence ${TRIM_ADAPTER} \
        --adapter_sequence_r2 ${TRIM_ADAPTER} \
        --correction \

    log_info "fastp修剪完成!"
    log_info "质控报告: ${html_report}"
}

#---------------------------------------
# Trimmomatic 修剪
#---------------------------------------
run_trimmomatic() {
    local output_r1="${DIR_TRIMMED}/${SAMPLE_ID}_trimmed_R1.fastq.gz"
    local output_r1_unpaired="${DIR_TRIMMED}/${SAMPLE_ID}_trimmed_R1_unpaired.fastq.gz"
    local output_r2="${DIR_TRIMMED}/${SAMPLE_ID}_trimmed_R2.fastq.gz"
    local output_r2_unpaired="${DIR_TRIMMED}/${SAMPLE_ID}_trimmed_R2_unpaired.fastq.gz"

    run_cmd "Trimmomatic 数据修剪" \
        "${TOOL_TRIMMOMATIC}" PE \
        -threads ${BWA_THREADS} \
        -phred33 \
        ${FASTQ_R1} ${FASTQ_R2} \
        ${output_r1} ${output_r1_unpaired} \
        ${output_r2} ${output_r2_unpaired} \
        ILLUMINACLIP:${TRIM_ADAPTER}:2:30:10 \
        LEADING:${TRIM_QUALITY} \
        TRAILING:${TRIM_QUALITY} \
        SLIDINGWINDOW:4:${TRIM_QUALITY} \
        MINLEN:${TRIM_MIN_LENGTH}

    # 清理未配对文件 (通常不需要)
    if [ -f "${output_r1_unpaired}" ] && [ ! -s "${output_r1_unpaired}" ]; then
        rm -f "${output_r1_unpaired}"
    fi
    if [ -f "${output_r2_unpaired}" ] && [ ! -s "${output_r2_unpaired}" ]; then
        rm -f "${output_r2_unpaired}"
    fi

    log_info "Trimmomatic修剪完成!"
}

#---------------------------------------
# 主函数
#---------------------------------------
main() {
    log_step "步骤2: 数据修剪 (Adapter去除 & 质控过滤)"

    # 检查是否跳过
    if [ "${SKIP_TRIM}" = true ]; then
        log_warn "跳过数据修剪 (SKIP_TRIM=true)"
        return 0
    fi

    # 检查输入文件
    check_file "FASTQ R1" "${FASTQ_R1}" || return 1
    check_file "FASTQ R2" "${FASTQ_R2}" || return 1

    # 确保输出目录存在
    mkdir -p "${DIR_TRIMMED}"

    # 根据配置选择修剪工具
    if [ "${USE_FASTP}" = true ]; then
        log_info "使用 fastp 进行数据修剪"
        run_fastp
    else
        log_info "使用 Trimmomatic 进行数据修剪"
        run_trimmomatic
    fi

    log_info "数据修剪完成! 结果目录: ${DIR_TRIMMED}"
}

main "$@"
