#!/bin/bash
#===============================================================================
# 05_mark_duplicates.sh - 标记/去除PCR重复
#
# 说明: PCR扩增过程中产生的重复reads需要被标记
#       支持 Picard MarkDuplicates 或 samtools markdup
#       标记后的BAM用于后续变异检测
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
    log_step "步骤5: 标记PCR重复 (Mark Duplicates)"

    # 检查是否跳过
    if [ "${SKIP_MARKDUP}" = true ]; then
        log_warn "跳过标记重复 (SKIP_MARKDUP=true)"
        return 0
    fi

    # 输入BAM
    local input_bam="${DIR_ALIGNED}/${SAMPLE_ID}.sorted.bam"
    check_file "排序BAM" "${input_bam}" || return 1

    # 输出文件
    local dedup_bam="${DIR_ALIGNED}/${SAMPLE_ID}.dedup.bam"
    local metrics_file="${DIR_ALIGNED}/${SAMPLE_ID}.dup_metrics.txt"

    # 检查是否已有结果
    if [ -f "${dedup_bam}" ] && [ -f "${dedup_bam}.bai" ]; then
        log_warn "去重BAM已存在，跳过"
        return 0
    fi

    # 根据工具可用性选择标记方法
    if command -v "${TOOL_PICARD}" &>/dev/null; then
        #---------------------------------------
        # 方法1: Picard MarkDuplicates (推荐)
        #---------------------------------------
        log_info "使用 Picard MarkDuplicates"

        run_cmd "Picard 标记PCR重复" \
            "${TOOL_PICARD} MarkDuplicates" \
            INPUT=${input_bam} \
            OUTPUT=${dedup_bam} \
            METRICS_FILE=${metrics_file} \
            VALIDATION_STRINGENCY=LENIENT \
            REMOVE_DUPLICATES=false \
            CREATE_INDEX=true \
            MAX_RECORDS_IN_RAM=500000

    else
        #---------------------------------------
        # 方法2: samtools markdup (备用)
        #---------------------------------------
        log_info "Picard不可用，使用 samtools markdup"

        # samtools markdup 需要name-sorted BAM
        local namesorted_bam="${DIR_ALIGNED}/${SAMPLE_ID}.namesorted.bam"
        local fixmate_bam="${DIR_ALIGNED}/${SAMPLE_ID}.fixmate.bam"
        local fixmate_sorted_bam="${DIR_ALIGNED}/${SAMPLE_ID}.fixmate.sorted.bam"

        # 先转为name-sorted
        run_cmd "转换为name-sorted BAM" \
            "${TOOL_SAMTOOLS} sort" \
            -n \
            -@ ${SAMTOOLS_THREADS} \
            -o "${namesorted_bam}" \
            "${input_bam}"

        run_cmd "samtools fixmate 添加mate tags" \
            "${TOOL_SAMTOOLS} fixmate" \
            -m \
            "${namesorted_bam}" \
            "${fixmate_bam}"

        run_cmd "fixmate BAM坐标排序" \
            "${TOOL_SAMTOOLS} sort" \
            -@ ${SAMTOOLS_THREADS} \
            -o "${fixmate_sorted_bam}" \
            "${fixmate_bam}"

        local markdup_tag_params=()
        local markdup_help
        markdup_help="$("${TOOL_SAMTOOLS}" markdup --help 2>&1 || true)"
        if printf '%s\n' "${markdup_help}" | grep -q -- '--write-tags'; then
            markdup_tag_params+=(--write-tags)
        elif printf '%s\n' "${markdup_help}" | grep -q -- '--duplicate-count'; then
            markdup_tag_params+=(--duplicate-count)
        else
            log_warn "当前 samtools markdup 不支持 --write-tags/--duplicate-count，仅标记duplicate flag"
        fi

        # 标记重复
        run_cmd "samtools markdup 标记重复" \
            "${TOOL_SAMTOOLS} markdup" \
            --threads ${SAMTOOLS_THREADS} \
            "${markdup_tag_params[@]}" \
            "${fixmate_sorted_bam}" \
            "${dedup_bam}"

        run_cmd "去重BAM建立索引" \
            "${TOOL_SAMTOOLS} index" \
            -@ ${SAMTOOLS_THREADS} \
            "${dedup_bam}"

        # 生成简单的重复统计
        echo "PCR duplicate statistics (samtools markdup)" > "${metrics_file}"
        echo "See BAM tags for duplicate information" >> "${metrics_file}"

        # 清理临时文件
        rm -f "${namesorted_bam}" "${fixmate_bam}" "${fixmate_sorted_bam}"
    fi

    # 验证输出
    if [ ! -f "${dedup_bam}" ]; then
        log_error "去重BAM文件未生成!"
        return 1
    fi

    # 统计去重后的比对信息
    run_cmd "去重后比对统计" \
        "${TOOL_SAMTOOLS} flagstat" \
        "${dedup_bam}" \
        "> ${DIR_ALIGNED}/${SAMPLE_ID}.dedup.flagstat.txt"

    log_info "标记PCR重复完成!"
    log_info "去重BAM: ${dedup_bam}"
    log_info "重复统计: ${metrics_file}"
}

main "$@"
