#!/bin/bash
#===============================================================================
# 05b_post_align_qc.sh - 比对后质量控制
#
# 说明: 使用Picard和Samtools收集比对后的质量指标
#       包括: 插入片段大小、文库复杂度、目标区域覆盖指标等
#       这些指标用于评估测序和捕获实验的质量
#
# 输入: 去重后的BAM文件
# 输出: 各类QC指标文件 (txt/pdf)
# 依赖: Picard, Samtools, Qualimap(可选)
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
    log_step "步骤5b: 比对后质量控制 (Post-Alignment QC)"

    if [ "${SKIP_POSTQC}" = true ]; then
        log_warn "跳过比对后QC (SKIP_POSTQC=true)"
        return 0
    fi

    # 确定输入BAM
    local input_bam
    if [ -f "${DIR_ALIGNED}/${SAMPLE_ID}.dedup.bam" ]; then
        input_bam="${DIR_ALIGNED}/${SAMPLE_ID}.dedup.bam"
    elif [ -f "${DIR_ALIGNED}/${SAMPLE_ID}.sorted.bam" ]; then
        input_bam="${DIR_ALIGNED}/${SAMPLE_ID}.sorted.bam"
    else
        log_error "未找到输入BAM文件!"
        return 1
    fi

    check_file "输入BAM" "${input_bam}" || return 1
    mkdir -p "${DIR_POSTQC}"

    #---------------------------------------
    # 1. Picard CollectInsertSizeMetrics - 插入片段大小分布
    #---------------------------------------
    if command -v "${TOOL_PICARD}" &>/dev/null; then
        run_cmd "收集插入片段大小分布" \
            "${TOOL_PICARD} CollectInsertSizeMetrics" \
            INPUT=${input_bam} \
            OUTPUT=${DIR_POSTQC}/${SAMPLE_ID}.insert_size_metrics \
            HISTOGRAM_FILE=${DIR_POSTQC}/${SAMPLE_ID}.insert_size_histogram.pdf \
            VALIDATION_STRINGENCY=LENIENT \
            MINIMUM_PERCENTAGE=0.05

        #---------------------------------------
        # 2. Picard CollectMultipleMetrics - 多种QC指标
        #---------------------------------------
        run_cmd "收集多种比对指标" \
            "${TOOL_PICARD} CollectMultipleMetrics" \
            INPUT=${input_bam} \
            OUTPUT=${DIR_POSTQC}/${SAMPLE_ID}.multiple_metrics \
            VALIDATION_STRINGENCY=LENIENT \
            PROGRAM=null \
            PROGRAM=CollectAlignmentSummaryMetrics \
            PROGRAM=CollectBaseDistributionByCycle \
            PROGRAM=CollectGcBiasMetrics \
            PROGRAM=MeanQualityByCycle \
            PROGRAM=QualityScoreDistribution \
            REFERENCE_SEQUENCE=${REFERENCE_GENOME}

        #---------------------------------------
        # 3. Picard EstimateLibraryComplexity - 文库复杂度估计
        #---------------------------------------
        run_cmd "估计文库复杂度" \
            "${TOOL_PICARD} EstimateLibraryComplexity" \
            INPUT=${input_bam} \
            OUTPUT=${DIR_POSTQC}/${SAMPLE_ID}.library_complexity.txt \
            VALIDATION_STRINGENCY=LENIENT
    fi

    #---------------------------------------
    # 4. Samtools stats - 详细比对统计
    #---------------------------------------
    run_cmd "Samtools详细比对统计" \
        "${TOOL_SAMTOOLS} stats" \
        -@ ${SAMTOOLS_THREADS} \
        "${input_bam}" \
        "> ${DIR_POSTQC}/${SAMPLE_ID}.samtools_stats.txt"

    #---------------------------------------
    # 5. 外显子组覆盖指标 (如果是外显子组数据)
    #---------------------------------------
    if command -v "${TOOL_PICARD}" &>/dev/null; then
        if [ -n "${INTERVAL_FILE}" ] && [ -f "${INTERVAL_FILE}" ]; then
            run_cmd "收集外显子组覆盖指标 (HsMetrics)" \
                "${TOOL_PICARD} CollectHsMetrics" \
                INPUT=${input_bam} \
                OUTPUT=${DIR_POSTQC}/${SAMPLE_ID}.hs_metrics.txt \
                VALIDATION_STRINGENCY=LENIENT \
                REFERENCE_SEQUENCE=${REFERENCE_GENOME} \
                BAI_INTERVAL_FILE=${INTERVAL_FILE} \
                TARGET_INTERVALS=${INTERVAL_FILE}
        fi

        #---------------------------------------
        # 6. CollectWgsMetrics - 全基因组覆盖指标
        #---------------------------------------
        run_cmd "收集WGS覆盖指标" \
            "${TOOL_PICARD} CollectWgsMetrics" \
            INPUT=${input_bam} \
            OUTPUT=${DIR_POSTQC}/${SAMPLE_ID}.wgs_metrics.txt \
            VALIDATION_STRINGENCY=LENIENT \
            REFERENCE_SEQUENCE=${REFERENCE_GENOME} \
            MINIMUM_BASE_QUALITY=0 \
            MINIMUM_MAPPING_QUALITY=0
    fi

    #---------------------------------------
    # 7. Qualimap bamqc (可选)
    #---------------------------------------
    if command -v "${TOOL_QUALIMAP}" &>/dev/null; then
        local qualimap_extra=""
        if [ -n "${INTERVAL_FILE}" ] && [ -f "${INTERVAL_FILE}" ]; then
            qualimap_extra="-gff ${INTERVAL_FILE}"
        fi

        run_cmd "Qualimap BAM质控" \
            "${TOOL_QUALIMAP} bamqc" \
            -bam ${input_bam} \
            ${qualimap_extra} \
            -outdir ${DIR_POSTQC}/qualimap_${SAMPLE_ID} \
            -nt ${SAMTOOLS_THREADS} \
            -outformat HTML
    fi

    log_info "比对后QC完成! 结果目录: ${DIR_POSTQC}"
}

main "$@"
