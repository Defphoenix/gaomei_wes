#!/bin/bash
#===============================================================================
# 05c_bqsr.sh - 碱基质量分数重校准 (Base Quality Score Recalibration)
#
# 说明: GATK BQSR 校正测序仪系统误差导致的碱基质量偏差
#       步骤:
#         1. BaseRecalibrator - 构建重校准模型
#         2. ApplyBQSR - 应用校正到BAM
#       注意: 需要dbSNP和已知indels数据库
#
# 输入: 去重后的BAM
# 输出: 重校准后的BAM + 重校准报告
# 依赖: GATK4
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
    log_step "步骤5c: 碱基质量重校准 (BQSR)"

    if [ "${SKIP_BQSR}" = true ] || [ "${RUN_BQSR}" = false ]; then
        log_warn "跳过BQSR (SKIP_BQSR=${SKIP_BQSR}, RUN_BQSR=${RUN_BQSR})"
        return 0
    fi

    # 输入BAM (去重后)
    local input_bam="${DIR_ALIGNED}/${SAMPLE_ID}.dedup.bam"
    if [ ! -f "${input_bam}" ]; then
        input_bam="${DIR_ALIGNED}/${SAMPLE_ID}.sorted.bam"
    fi
    check_file "输入BAM" "${input_bam}" || return 1

    # 检查已知位点数据库
    if [ ! -f "${DBSNP_VCF}" ]; then
        log_error "dbSNP未配置或不存在: ${DBSNP_VCF}"
        log_error "BQSR需要dbSNP数据库，请从GATK Resource Bundle下载"
        return 1
    fi

    mkdir -p "${DIR_BQSR}"

    # 输出文件
    local recal_table="${DIR_BQSR}/${SAMPLE_ID}.recal_data.table"
    local bqsr_bam="${DIR_BQSR}/${SAMPLE_ID}.bqsr.bam"
    local bqsr_report_before="${DIR_BQSR}/${SAMPLE_ID}.recal_report_before.pdf"
    local bqsr_report_after="${DIR_BQSR}/${SAMPLE_ID}.recal_report_after.pdf"

    # 构建已知位点参数
    local known_sites_params="-known-sites ${DBSNP_VCF}"
    if [ -n "${KNOWN_INDELS_VCF}" ] && [ -f "${KNOWN_INDELS_VCF}" ]; then
        known_sites_params="${known_sites_params} -known-sites ${KNOWN_INDELS_VCF}"
    fi

    #---------------------------------------
    # 步骤1: BaseRecalibrator - 构建重校准模型
    #---------------------------------------
    # 分析BAM中的碱基质量与已知位点的差异
    # 生成重校准表 (recalibration table)
    run_cmd "BaseRecalibrator 构建重校准模型" \
        "${TOOL_GATK} --java-options \"-Xmx${GATK_JAVA_MEM}\" BaseRecalibrator" \
        -R "${REFERENCE_GENOME}" \
        -I "${input_bam}" \
        ${known_sites_params} \
        -O "${recal_table}"

    #---------------------------------------
    # 步骤2: GatherBQSRReports (可选, 用于分片分析)
    #---------------------------------------
    # 单样本通常不需要

    #---------------------------------------
    # 步骤3: AnalyzeCovariates - 生成重校准报告 (Before)
    #---------------------------------------
    run_cmd "生成重校准前报告" \
        "${TOOL_GATK} --java-options \"-Xmx${GATK_JAVA_MEM}\" AnalyzeCovariates" \
        -before "${recal_table}" \
        -plots "${bqsr_report_before}" || \
        log_warn "AnalyzeCovariates生成报告失败 (非致命)"

    #---------------------------------------
    # 步骤4: ApplyBQSR - 应用重校准
    #---------------------------------------
    # 根据重校准表修正BAM中的碱基质量分数
    run_cmd "ApplyBQSR 应用碱基质量重校准" \
        "${TOOL_GATK} --java-options \"-Xmx${GATK_JAVA_MEM}\" ApplyBQSR" \
        -R "${REFERENCE_GENOME}" \
        -I "${input_bam}" \
        --bqsr-recal-file "${recal_table}" \
        -O "${bqsr_bam}"

    #---------------------------------------
    # 步骤5: 对BQSR后的BAM建立索引
    #---------------------------------------
    run_cmd "BQSR BAM建立索引" \
        "${TOOL_SAMTOOLS} index" \
        -@ ${SAMTOOLS_THREADS} \
        "${bqsr_bam}"

    #---------------------------------------
    # 步骤6: 生成重校准后报告 (用于对比)
    #---------------------------------------
    run_cmd "构建重校准后模型" \
        "${TOOL_GATK} --java-options \"-Xmx${GATK_JAVA_MEM}\" BaseRecalibrator" \
        -R "${REFERENCE_GENOME}" \
        -I "${bqsr_bam}" \
        ${known_sites_params} \
        -O "${DIR_BQSR}/${SAMPLE_ID}.recal_data_after.table" || \
        log_warn "重校准后模型构建失败 (非致命)"

    if [ -f "${DIR_BQSR}/${SAMPLE_ID}.recal_data_after.table" ]; then
        run_cmd "生成重校准前后对比报告" \
            "${TOOL_GATK} --java-options \"-Xmx${GATK_JAVA_MEM}\" AnalyzeCovariates" \
            -before "${recal_table}" \
            -after "${DIR_BQSR}/${SAMPLE_ID}.recal_data_after.table" \
            -plots "${bqsr_report_after}" || \
            log_warn "对比报告生成失败 (非致命)"
    fi

    # 验证输出
    if [ ! -f "${bqsr_bam}" ]; then
        log_error "BQSR BAM文件未生成!"
        return 1
    fi

    log_info "BQSR完成!"
    log_info "重校准BAM: ${bqsr_bam}"
    log_info "重校准表: ${recal_table}"
    log_info "校准报告: ${bqsr_report_before}"
}

main "$@"
