#!/bin/bash
#===============================================================================
# 11_coverage.sh - 覆盖度分析 (mosdepth + bedtools)
#
# 说明: 使用mosdepth进行快速覆盖度分析，统计目标区域深度及覆盖比例
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
    log_step "步骤11: 覆盖度分析 (mosdepth + bedtools)"

    if [ "${SKIP_COVERAGE:-false}" = true ] || [ "${RUN_COVERAGE:-true}" = false ]; then
        log_warn "跳过覆盖度分析"
        return 0
    fi

    # 确定输入BAM；tumor-normal汇总配置没有自身的aligned目录，优先使用肿瘤BAM。
    local input_bam
    if [ -n "${TUMOR_BAM:-}" ] && [ -f "${TUMOR_BAM}" ]; then
        input_bam="${TUMOR_BAM}"
        log_info "使用配置指定的Tumor BAM进行覆盖度分析: ${input_bam}"
    else
        input_bam=$(get_final_bam) || return 1
    fi
    check_file "输入BAM" "${input_bam}" || return 1

    mkdir -p "${DIR_COVERAGE}"

    # 输出文件前缀
    local prefix="${DIR_COVERAGE}/${SAMPLE_ID}"

    #---------------------------------------
    # 步骤1: mosdepth 全局覆盖度
    #---------------------------------------
    # 构建阈值参数 (例如: -T 10,30,50,100)
    local threshold_params=""
    if [ -n "${COVERAGE_TARGETS:-}" ]; then
        threshold_params="-T ${COVERAGE_TARGETS// /,}"
    fi

    log_info "执行: mosdepth 覆盖度分析"
    
    # 调整参数顺序: 选项在前，<prefix> 和 <bam> 在最后，确保 mosdepth 正确解析
    local mosdepth_cmd=(${TOOL_MOSDEPTH} --threads ${MOSDEPTH_THREADS} --mapq 20 --no-per-base)
    
    if [ -n "${INTERVAL_FILE:-}" ] && [ -f "${INTERVAL_FILE}" ]; then
        mosdepth_cmd+=(--by "${INTERVAL_FILE}")
        log_info "使用目标区域BED: ${INTERVAL_FILE}"
    fi
    
    if [ -n "${threshold_params}" ]; then
        mosdepth_cmd+=(${threshold_params})
    fi
    
    # 执行命令
    "${mosdepth_cmd[@]}" "${prefix}" "${input_bam}"

    #---------------------------------------
    # 步骤2: 解析mosdepth结果
    #---------------------------------------
    local summary_file="${prefix}.mosdepth.summary.txt"
    if [ -f "${summary_file}" ]; then
        log_info "mosdepth汇总统计:"
        cat "${summary_file}" | sed 's/^/  /'
    fi

    #---------------------------------------
    # 步骤3: bedtools 目标区域覆盖统计
    #---------------------------------------
    if [ -n "${INTERVAL_FILE:-}" ] && [ -f "${INTERVAL_FILE}" ]; then
        log_info "执行: bedtools 目标区域覆盖统计"
        ${TOOL_BEDTOOLS} coverage -a "${INTERVAL_FILE}" -b "${input_bam}" -sorted > "${DIR_COVERAGE}/${SAMPLE_ID}.bedtools_coverage.txt"

        log_info "执行: bedtools 目标区域read计数"
        ${TOOL_BEDTOOLS} multicov -bams "${input_bam}" -bed "${INTERVAL_FILE}" > "${DIR_COVERAGE}/${SAMPLE_ID}.multicov.txt" || \
            log_warn "bedtools multicov失败 (非致命)"
    fi

    #---------------------------------------
    # 步骤4: 生成覆盖度汇总报告
    #---------------------------------------
    {
        echo "覆盖度分析报告 - ${SAMPLE_ID}"
        echo "================================"
        echo "分析时间: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "目标区域: ${INTERVAL_FILE:-全基因组}"
        echo ""

        if [ -f "${summary_file}" ]; then
            echo "【mosdepth覆盖度统计】"
            cat "${summary_file}"; echo ""
        fi

        local bedtools_cov="${DIR_COVERAGE}/${SAMPLE_ID}.bedtools_coverage.txt"
        if [ -f "${bedtools_cov}" ]; then
            echo "【bedtools目标区域覆盖】"
            awk 'BEGIN{OFS="\t"} {bases+=$4; len+=$3} END{if(len>0) printf "目标区域平均深度: %.1fx\n",bases/len}' "${bedtools_cov}"
        fi
    } > "${DIR_COVERAGE}/${SAMPLE_ID}_coverage_report.txt"

    log_info "覆盖度分析完成! 报告: ${DIR_COVERAGE}/${SAMPLE_ID}_coverage_report.txt"
}

main "$@"
