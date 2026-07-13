#!/bin/bash
#===============================================================================
# 04_sort_index.sh - BAM排序和索引
#
# 说明: 将SAM转换为排序后的BAM文件，并建立索引
#       排序后的BAM是后续分析的基础输入
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
    log_step "步骤4: BAM排序和索引"

    # 检查是否跳过
    if [ "${SKIP_SORT}" = true ]; then
        log_warn "跳过BAM排序 (SKIP_SORT=true)"
        return 0
    fi

    # 输入SAM文件
    local input_sam="${DIR_ALIGNED}/${SAMPLE_ID}.sam"
    check_file "SAM文件" "${input_sam}" || return 1

    # 输出文件路径
    local sorted_bam="${DIR_ALIGNED}/${SAMPLE_ID}.sorted.bam"
    local sorted_bam_index="${sorted_bam}.bai"

    # 检查是否已有排序结果
    if [ -f "${sorted_bam}" ] && [ -f "${sorted_bam_index}" ]; then
        log_warn "排序BAM和索引已存在，跳过"
        return 0
    fi

    # 步骤4a: SAM -> 排序BAM
    # samtools sort: 按坐标排序
    # -@: 线程数
    # -m: 每线程最大内存
    run_cmd "SAM转换为排序BAM" \
        "${TOOL_SAMTOOLS} sort" \
        -@ ${SAMTOOLS_THREADS} \
        -m 2G \
        -o "${sorted_bam}" \
        "${input_sam}"

    # 步骤4b: 建立BAM索引
    # 索引文件用于快速随机访问BAM的特定区域
    run_cmd "建立BAM索引" \
        "${TOOL_SAMTOOLS} index" \
        -@ ${SAMTOOLS_THREADS} \
        "${sorted_bam}"

    # 步骤4c: 统计比对信息
    run_cmd "统计比对信息" \
        "${TOOL_SAMTOOLS} flagstat" \
        -@ ${SAMTOOLS_THREADS} \
        "${sorted_bam}" \
        "> ${DIR_ALIGNED}/${SAMPLE_ID}.flagstat.txt"

    # 统计目标区域覆盖度 (如果有BED文件)
    if [ -f "${INTERVAL_FILE}" ]; then
        run_cmd "统计目标区域覆盖度" \
            "${TOOL_SAMTOOLS} bedcov" \
            "${INTERVAL_FILE}" \
            "${sorted_bam}" \
            "> ${DIR_ALIGNED}/${SAMPLE_ID}.bedcov.txt" || \
            log_warn "bedcov统计失败，可能BED文件格式不兼容"
    fi

    # 删除原始SAM文件以节省空间
    if [ -f "${sorted_bam}" ] && [ -f "${sorted_bam_index}" ]; then
        rm -f "${input_sam}"
        log_info "已删除原始SAM文件以节省空间"
    fi

    log_info "BAM排序和索引完成!"
    log_info "排序BAM: ${sorted_bam}"
    log_info "索引文件: ${sorted_bam_index}"
}

main "$@"
