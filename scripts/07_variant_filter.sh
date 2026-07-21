#!/bin/bash
#===============================================================================
# 07_variant_filter.sh - 变异过滤 (完全修复版)
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/../config.sh}"
source "${CONFIG_FILE}"
source "${SCRIPT_DIR}/utils.sh"

#---------------------------------------
# HaplotypeCaller 变异硬过滤
#---------------------------------------
filter_haplotypecaller() {
    local input_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.raw.vcf.gz"
    local snp_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.snps.vcf.gz"
    local indel_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.indels.vcf.gz"
    local filtered_snp="${DIR_VARIANTS}/${SAMPLE_ID}.snps.filtered.vcf.gz"
    local filtered_indel="${DIR_VARIANTS}/${SAMPLE_ID}.indels.filtered.vcf.gz"
    local final_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.filtered.vcf.gz"

    check_file "原始VCF" "${input_vcf}" || return 1

    # 1. 分离 SNP 和 InDel (修正了参数名为 -select-type)
    log_info "执行: 分离SNP位点"
    ${TOOL_GATK} --java-options "-Xmx${GATK_JAVA_MEM}" SelectVariants \
        -R "${REFERENCE_GENOME}" \
        -V "${input_vcf}" \
        -select-type SNP \
        -O "${snp_vcf}"

    log_info "执行: 分离InDel位点"
    ${TOOL_GATK} --java-options "-Xmx${GATK_JAVA_MEM}" SelectVariants \
        -R "${REFERENCE_GENOME}" \
        -V "${input_vcf}" \
        -select-type INDEL \
        -O "${indel_vcf}"

    # 2. 执行硬过滤 (使用原生的 GATK 调用方式，避开 Shell 管道符号解析错误)
    log_info "执行: SNP硬过滤"
    ${TOOL_GATK} --java-options "-Xmx${GATK_JAVA_MEM}" VariantFiltration \
        -R "${REFERENCE_GENOME}" \
        -V "${snp_vcf}" \
        -O "${filtered_snp}" \
        --filter-expression "${HC_SNP_QUAL_SCORE}" \
        --filter-name "SNP_FILTER"

    log_info "执行: InDel硬过滤"
    ${TOOL_GATK} --java-options "-Xmx${GATK_JAVA_MEM}" VariantFiltration \
        -R "${REFERENCE_GENOME}" \
        -V "${indel_vcf}" \
        -O "${filtered_indel}" \
        --filter-expression "${HC_INDEL_QUAL_SCORE}" \
        --filter-name "INDEL_FILTER"

    # 3. 合并结果
    log_info "执行: 合并SNP和InDel结果"
    ${TOOL_GATK} --java-options "-Xmx${GATK_JAVA_MEM}" MergeVcfs \
        -I "${filtered_snp}" \
        -I "${filtered_indel}" \
        -O "${final_vcf}"

    # 4. 提取 PASS 位点
    local pass_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.pass.vcf.gz"
    log_info "执行: 提取PASS位点"
    ${TOOL_BCFTOOLS} view -f PASS -Oz -o "${pass_vcf}" "${final_vcf}"

    # 5. 建立索引并清理
    ${TOOL_BCFTOOLS} index -t "${pass_vcf}"
    rm -f "${snp_vcf}" "${indel_vcf}" "${filtered_snp}" "${filtered_indel}" "${final_vcf}"
    
    log_info "变异过滤完成! 最终结果: ${pass_vcf}"
}

#---------------------------------------
# Mutect2 体细胞变异过滤
#---------------------------------------
build_mutect2_orientation_model() {
    local f1r2_tar="$1"
    local orientation_model="$2"

    if [ "${RUN_MUTECT2_ORIENTATION_MODEL:-true}" != true ]; then
        log_warn "未启用Mutect2方向偏倚模型"
        return 0
    fi
    if [ ! -s "${f1r2_tar}" ]; then
        log_warn "未找到F1R2计数文件，无法建立方向偏倚模型: ${f1r2_tar}"
        [ "${MUTECT2_REQUIRE_AUXILIARY:-false}" = true ] && return 1
        return 0
    fi

    if ! run_cmd "Mutect2学习read-orientation偏倚模型" \
        "${TOOL_GATK} --java-options \"-Xmx${GATK_JAVA_MEM}\" LearnReadOrientationModel" \
        -I "${f1r2_tar}" \
        -O "${orientation_model}"; then
        rm -f "${orientation_model}"
        [ "${MUTECT2_REQUIRE_AUXILIARY:-false}" = true ] && return 1
        log_warn "方向偏倚模型生成失败，继续执行基础FilterMutectCalls"
        return 0
    fi
    check_file "方向偏倚模型" "${orientation_model}" || return 1
}

resolve_mutect2_common_variants() {
    local candidate
    for candidate in \
        "${MUTECT2_COMMON_VARIANTS_VCF:-}" \
        "${REFERENCE_DIR:-}/mutect2/common_biallelic.hg38.vcf.gz" \
        "${REFERENCE_DIR:-}/mutect2/small_exac_common_3.hg38.vcf.gz" \
        "${REFERENCE_DIR:-}/small_exac_common_3.hg38.vcf.gz"; do
        if [ -n "${candidate}" ] && [ -f "${candidate}" ]; then
            echo "${candidate}"
            return 0
        fi
    done
    return 1
}

build_mutect2_contamination_model() {
    local tumor_bam="$1"
    local contamination_table="$2"
    local segmentation_table="$3"
    local common_vcf=""
    local tumor_pileups="${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.tumor.pileups.table"
    local normal_pileups="${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.normal.pileups.table"
    local interval_params=""

    if [ "${RUN_MUTECT2_CONTAMINATION:-true}" != true ]; then
        log_warn "未启用Mutect2污染估计"
        return 0
    fi
    if ! common_vcf=$(resolve_mutect2_common_variants); then
        log_warn "未找到Mutect2 common biallelic variants资源，跳过污染估计"
        log_warn "请配置 MUTECT2_COMMON_VARIANTS_VCF，例如small_exac_common_3.hg38.vcf.gz"
        [ "${MUTECT2_REQUIRE_AUXILIARY:-false}" = true ] && return 1
        return 0
    fi
    if [ ! -f "${common_vcf}.tbi" ] && [ ! -f "${common_vcf}.idx" ]; then
        log_warn "Mutect2 common variants缺少索引: ${common_vcf}.tbi 或 ${common_vcf}.idx"
        [ "${MUTECT2_REQUIRE_AUXILIARY:-false}" = true ] && return 1
        return 0
    fi

    if [ -n "${INTERVAL_LIST:-}" ] && [ -f "${INTERVAL_LIST}" ]; then
        interval_params="-L ${INTERVAL_LIST} -L ${common_vcf} --interval-set-rule INTERSECTION"
    elif [ -n "${INTERVAL_FILE:-}" ] && [ -f "${INTERVAL_FILE}" ]; then
        interval_params="-L ${INTERVAL_FILE} -L ${common_vcf} --interval-set-rule INTERSECTION"
    else
        interval_params="-L ${common_vcf}"
    fi

    log_info "使用污染估计位点集: ${common_vcf}"
    if ! run_cmd "计算Tumor common-variant pileup" \
        "${TOOL_GATK} --java-options \"-Xmx${GATK_JAVA_MEM}\" GetPileupSummaries" \
        -R "${REFERENCE_GENOME}" \
        -I "${tumor_bam}" \
        -V "${common_vcf}" \
        -O "${tumor_pileups}" \
        ${interval_params}; then
        rm -f "${tumor_pileups}" "${contamination_table}" "${segmentation_table}"
        [ "${MUTECT2_REQUIRE_AUXILIARY:-false}" = true ] && return 1
        log_warn "Tumor pileup生成失败，跳过污染估计"
        return 0
    fi

    local matched_param=""
    if [ -n "${NORMAL_BAM:-}" ] && [ -f "${NORMAL_BAM}" ]; then
        if ! run_cmd "计算Normal common-variant pileup" \
            "${TOOL_GATK} --java-options \"-Xmx${GATK_JAVA_MEM}\" GetPileupSummaries" \
            -R "${REFERENCE_GENOME}" \
            -I "${NORMAL_BAM}" \
            -V "${common_vcf}" \
            -O "${normal_pileups}" \
            ${interval_params}; then
            rm -f "${normal_pileups}" "${contamination_table}" "${segmentation_table}"
            [ "${MUTECT2_REQUIRE_AUXILIARY:-false}" = true ] && return 1
            log_warn "Normal pileup生成失败，跳过污染估计"
            return 0
        fi
        matched_param="--matched ${normal_pileups}"
    fi

    if ! run_cmd "Mutect2计算样本污染比例" \
        "${TOOL_GATK} --java-options \"-Xmx${GATK_JAVA_MEM}\" CalculateContamination" \
        -I "${tumor_pileups}" \
        ${matched_param} \
        -O "${contamination_table}" \
        --tumor-segmentation "${segmentation_table}"; then
        rm -f "${contamination_table}" "${segmentation_table}"
        [ "${MUTECT2_REQUIRE_AUXILIARY:-false}" = true ] && return 1
        log_warn "污染比例计算失败，继续执行不带污染模型的FilterMutectCalls"
        return 0
    fi

    check_file "污染估计表" "${contamination_table}" || return 1
    check_file "污染分段表" "${segmentation_table}" || return 1
}

filter_mutect2() {
    local input_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.raw.vcf.gz"
    local f1r2_tar="${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.f1r2.tar.gz"
    local orientation_model="${MUTECT2_ORIENTATION_MODEL:-${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.read-orientation-model.tar.gz}"
    local contamination_table="${MUTECT2_CONTAMINATION_TABLE:-${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.contamination.table}"
    local segmentation_table="${MUTECT2_SEGMENTATION_TABLE:-${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.segments.table}"
    local filtered_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.filtered.vcf.gz"
    local pass_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.pass.vcf.gz"
    local tumor_bam=""

    check_file "Mutect2原始VCF" "${input_vcf}" || return 1

    if [ -n "${TUMOR_BAM:-}" ] && [ -f "${TUMOR_BAM}" ]; then
        tumor_bam="${TUMOR_BAM}"
    else
        tumor_bam=$(get_final_bam) || return 1
    fi
    check_file "Tumor BAM" "${tumor_bam}" || return 1

    build_mutect2_orientation_model "${f1r2_tar}" "${orientation_model}" || return 1
    build_mutect2_contamination_model "${tumor_bam}" "${contamination_table}" "${segmentation_table}" || return 1

    local contamination_param=""
    if [ -f "${contamination_table}" ]; then
        contamination_param="--contamination-table ${contamination_table}"
        if [ -f "${segmentation_table}" ]; then
            contamination_param="${contamination_param} --tumor-segmentation ${segmentation_table}"
        fi
        log_info "使用污染估计表: ${contamination_table}"
    fi

    local orientation_param=""
    if [ -f "${orientation_model}" ]; then
        orientation_param="--ob-priors ${orientation_model}"
        log_info "使用方向偏倚模型: ${orientation_model}"
    fi

    local interval_param=""
    local interval_padding_param=""
    if [ -n "${INTERVAL_LIST:-}" ] && [ -f "${INTERVAL_LIST:-}" ]; then
        interval_param="-L ${INTERVAL_LIST}"
    elif [ -n "${INTERVAL_FILE:-}" ] && [ -f "${INTERVAL_FILE:-}" ]; then
        interval_param="-L ${INTERVAL_FILE}"
    fi
    if [ -n "${interval_param}" ]; then
        if [[ ! "${MUTECT2_INTERVAL_PADDING:-100}" =~ ^[0-9]+$ ]]; then
            log_error "MUTECT2_INTERVAL_PADDING必须为非负整数: ${MUTECT2_INTERVAL_PADDING:-}"
            return 1
        fi
        interval_padding_param="--interval-padding ${MUTECT2_INTERVAL_PADDING:-100}"
    fi

    log_info "执行: Mutect2 FilterMutectCalls"
    ${TOOL_GATK} --java-options "-Xmx${GATK_JAVA_MEM}" FilterMutectCalls \
        -R "${REFERENCE_GENOME}" \
        -V "${input_vcf}" \
        -O "${filtered_vcf}" \
        ${interval_param} \
        ${interval_padding_param} \
        ${contamination_param} \
        ${orientation_param} \
        ${FILTER_MUTECT_EXTRA_PARAMS:-}

    check_file "Mutect2过滤VCF" "${filtered_vcf}" || return 1

    log_info "执行: 提取Mutect2 PASS位点"
    ${TOOL_BCFTOOLS} view -f PASS -Oz -o "${pass_vcf}" "${filtered_vcf}"
    ${TOOL_BCFTOOLS} index -t "${pass_vcf}"

    log_info "Mutect2过滤完成! 过滤VCF: ${filtered_vcf}"
    log_info "Mutect2 PASS结果: ${pass_vcf}"
}

#---------------------------------------
# 主函数
#---------------------------------------
main() {
    log_step "步骤7: 变异过滤 (Variant Filtering)"
    if [ "${SKIP_VARIANT_FILTER:-false}" = true ]; then return 0; fi
    mkdir -p "${DIR_VARIANTS}"

    case "${CALLER_MODE}" in
        haplotypecaller|hc)
            filter_haplotypecaller
            ;;
        mutect2|mt2)
            filter_mutect2
            ;;
        *)
            log_error "未知的CALLER_MODE: ${CALLER_MODE}"
            log_error "支持: haplotypecaller / mutect2"
            return 1
            ;;
    esac
}

main "$@"
