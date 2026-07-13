#!/bin/bash
#===============================================================================
# 03_align.sh - 序列比对 (BWA-MEM)
#
# 说明: 使用BWA-MEM将修剪后的FASTQ比对到参考基因组
#       输出SAM格式的比对结果
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/../config.sh}"
source "${CONFIG_FILE}"
source "${SCRIPT_DIR}/utils.sh"

#---------------------------------------
# 主函数
#---------------------------------------
main() {
    log_step "步骤3: BWA-MEM 序列比对"

    # 检查是否跳过
    if [ "${SKIP_ALIGN}" = true ]; then
        log_warn "跳过序列比对 (SKIP_ALIGN=true)"
        return 0
    fi

    # 确定输入FASTQ (修剪后或原始)
    local input_r1 input_r2
    if [ "${SKIP_TRIM}" = false ]; then
        # 使用修剪后的数据
        if [ "${USE_FASTP}" = true ]; then
            input_r1="${DIR_TRIMMED}/${SAMPLE_ID}_trimmed_R1.fastq.gz"
            input_r2="${DIR_TRIMMED}/${SAMPLE_ID}_trimmed_R2.fastq.gz"
        else
            input_r1="${DIR_TRIMMED}/${SAMPLE_ID}_trimmed_R1.fastq.gz"
            input_r2="${DIR_TRIMMED}/${SAMPLE_ID}_trimmed_R2.fastq.gz"
        fi
    else
        # 使用原始数据
        input_r1="${FASTQ_R1}"
        input_r2="${FASTQ_R2}"
    fi

    # 检查输入文件
    check_file "输入FASTQ R1" "${input_r1}" || return 1
    check_file "输入FASTQ R2" "${input_r2}" || return 1
    check_file "参考基因组" "${REFERENCE_BWA_INDEX}" || return 1

    # 确保输出目录存在
    mkdir -p "${DIR_ALIGNED}"

    # 输出SAM文件路径
    local output_sam="${DIR_ALIGNED}/${SAMPLE_ID}.sam"

    # BWA-MEM 比对
    # -M: 将较短split hits标记为secondary (Picard兼容)
    # -R: 添加read group信息
    # -t: 线程数
    # -Y: 使用soft clipping for supplementary alignments
    local read_group="@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tLB:lib1\tPL:ILLUMINA"

    log_info "执行: BWA-MEM 比对 (原生模式，防止引号丢失)"
    ${TOOL_BWA} mem \
        ${BWA_MEM_PARAMS} \
        -t ${BWA_THREADS} \
        -R "${read_group}" \
        -Y \
        "${REFERENCE_BWA_INDEX}" \
        "${input_r1}" \
        "${input_r2}" \
        > "${output_sam}"

    # 检查SAM文件是否生成成功
    if [ ! -s "${output_sam}" ]; then
        log_error "BWA比对失败，SAM文件为空"
        return 1
    fi

    # 统计比对信息
    local total_reads
    total_reads=$(grep -c "^@" "${output_sam}" || echo "N/A")
    log_info "SAM文件生成: ${output_sam}"

    log_info "BWA-MEM比对完成!"
}

main "$@"
