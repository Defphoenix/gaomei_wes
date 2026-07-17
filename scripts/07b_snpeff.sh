#!/bin/bash
#===============================================================================
# 07b_snpeff.sh - SnpEff 变异功能注释
#
# 说明: 使用SnpEff对变异位点进行功能注释
#       注释内容包括: 基因名、转录本、氨基酸变化、影响程度等
#       同时使用SnpSift进行数据库注释 (dbSNP, ClinVar, COSMIC等)
#
# 输入: 过滤后的VCF文件
# 输出: 注释后的VCF + 基因级别统计 (csv/html)
# 依赖: SnpEff, SnpSift (可选)
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
    log_step "步骤7b: SnpEff 变异功能注释"

    if [ "${SKIP_SNPEFF}" = true ] || [ "${RUN_SNPEFF}" = false ]; then
        log_warn "跳过SnpEff注释 (SKIP_SNPEFF=${SKIP_SNPEFF})"
        return 0
    fi

    # 确定输入VCF
    local input_vcf=""
    if [ "${CALLER_MODE}" = "haplotypecaller" ] || [ "${CALLER_MODE}" = "hc" ]; then
        input_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.pass.vcf.gz"
    elif [ "${CALLER_MODE}" = "mutect2" ] || [ "${CALLER_MODE}" = "mt2" ]; then
        input_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.pass.vcf.gz"
    fi

    if [ -z "${input_vcf}" ] || [ ! -f "${input_vcf}" ]; then
        log_error "未找到过滤后的VCF文件!"
        return 1
    fi

    check_file "输入VCF" "${input_vcf}" || return 1
    mkdir -p "${DIR_ANNOTATION}"

    # 输出文件
    local snpeff_vcf="${DIR_ANNOTATION}/${SAMPLE_ID}.snpeff.vcf.gz"
    local snpeff_report="${DIR_ANNOTATION}/${SAMPLE_ID}.snpeff_summary.html"
    local snpeff_genes="${DIR_ANNOTATION}/${SAMPLE_ID}.snpeff_genes.csv"

    #---------------------------------------
    # 步骤1: SnpEff 功能注释
    #---------------------------------------
    # -v: 详细输出
    # -stats: 生成统计报告
    # -gene: 只输出基因级别统计
    # -canon: 只使用canonical转录本
    run_cmd "SnpEff 功能注释 (${SNPEFF_DB})" \
        "${TOOL_SNPEFF}" ann \
        -v \
        -dataDir ${SNPEFF_DATA_DIR} \
        -config ${SNPEFF_CONFIG} \
        -stats "${snpeff_report}" \
        -gene "${snpeff_genes}" \
        -canon \
        ${SNPEFF_EXTRA_PARAMS:-} \
        ${SNPEFF_DB} \
        "${input_vcf}" \
        | "${TOOL_BGZIP:-bgzip}" > "${snpeff_vcf}" || \
    run_cmd "SnpEff 功能注释 (备用方式)" \
        "${TOOL_SNPEFF}" ann \
        -v \
        -dataDir ${SNPEFF_DATA_DIR} \
        ${SNPEFF_EXTRA_PARAMS:-} \
        ${SNPEFF_DB} \
        "${input_vcf}" \
        > "${DIR_ANNOTATION}/${SAMPLE_ID}.snpeff.vcf"

    # 建立索引
    if [ -f "${snpeff_vcf}" ]; then
        run_cmd "SnpEff VCF建立索引" \
            "${TOOL_BCFTOOLS} index" -t "${snpeff_vcf}"
    fi

    #---------------------------------------
    # 步骤2: SnpSift 数据库注释 (可选)
    #---------------------------------------
    if command -v "${TOOL_SNPSIFT}" &>/dev/null; then
        local current_vcf="${snpeff_vcf}"
        if [ ! -f "${current_vcf}" ]; then
            current_vcf="${DIR_ANNOTATION}/${SAMPLE_ID}.snpeff.vcf"
        fi

        # dbSNP注释
        if [ -f "${DBSNP_VCF}" ]; then
            run_cmd "SnpSift dbSNP注释" \
                "${TOOL_SNPSIFT} dbnsfp" \
                -v \
                -db "${DBSNP_VCF}" \
                "${current_vcf}" \
                > "${DIR_ANNOTATION}/${SAMPLE_ID}.snpsift_dbsnp.vcf" || \
                log_warn "SnpSift dbSNP注释失败 (非致命)"
        fi

        # ClinVar注释 (如果有数据库)
        local clinvar_vcf="${CLINVAR_VCF:-}"
        if [ -f "${clinvar_vcf}" ]; then
            run_cmd "SnpSift ClinVar注释" \
                "${TOOL_SNPSIFT} annotate" \
                -v \
                "${clinvar_vcf}" \
                "${current_vcf}" \
                > "${DIR_ANNOTATION}/${SAMPLE_ID}.snpsift_clinvar.vcf" || \
                log_warn "SnpSift ClinVar注释失败 (非致命)"
        fi

        # COSMIC注释 (如果有数据库, 肿瘤分析重要)
        local cosmic_vcf="${COSMIC_VCF:-}"
        if [ -f "${cosmic_vcf}" ]; then
            run_cmd "SnpSift COSMIC注释" \
                "${TOOL_SNPSIFT} annotate" \
                -v \
                "${cosmic_vcf}" \
                "${current_vcf}" \
                > "${DIR_ANNOTATION}/${SAMPLE_ID}.snpsift_cosmic.vcf" || \
                log_warn "SnpSift COSMIC注释失败 (非致命)"
        fi
    fi

    #---------------------------------------
    # 步骤3: 提取注释统计
    #---------------------------------------
    if [ -f "${snpeff_vcf}" ]; then
        # 统计变异类型分布
        run_cmd "统计SnpEff变异类型" \
            "${TOOL_BCFTOOLS} query" \
            -f '%INFO/ANN\n' "${snpeff_vcf}" 2>/dev/null | \
            tr ',' '\n' | \
            awk -F'|' '{print $2}' | \
            sort | uniq -c | sort -rn | \
            head -20 > "${DIR_ANNOTATION}/${SAMPLE_ID}.snpeff_variant_types.txt" || \
            log_warn "SnpEff变异类型统计失败"
    fi

    log_info "SnpEff注释完成!"
    log_info "注释VCF: ${snpeff_vcf}"
    log_info "注释报告: ${snpeff_report}"
}

main "$@"
