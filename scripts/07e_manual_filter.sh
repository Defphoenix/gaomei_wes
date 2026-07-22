#!/bin/bash
# Apply configurable evidence and population-frequency thresholds after VEP.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/../config.sh}"
source "${CONFIG_FILE}"
source "${SCRIPT_DIR}/utils.sh"

MANUAL_FILTER_TMP_NORMALIZED_VCF=""
MANUAL_FILTER_TMP_PLAIN_VCF=""

cleanup_manual_filter_tmp() {
    if [ -n "${MANUAL_FILTER_TMP_NORMALIZED_VCF}" ]; then
        rm -f -- \
            "${MANUAL_FILTER_TMP_NORMALIZED_VCF}" \
            "${MANUAL_FILTER_TMP_NORMALIZED_VCF}.tbi" \
            "${MANUAL_FILTER_TMP_NORMALIZED_VCF}.csi"
    fi
    if [ -n "${MANUAL_FILTER_TMP_PLAIN_VCF}" ]; then
        rm -f -- "${MANUAL_FILTER_TMP_PLAIN_VCF}"
    fi
}

main() {
    log_step "步骤7e: VEP后人工硬过滤"

    if [ "${SKIP_MANUAL_FILTER:-false}" = true ] || [ "${RUN_MANUAL_FILTER:-true}" = false ]; then
        log_warn "跳过VEP后人工硬过滤"
        return 0
    fi
    if [ "${CALLER_MODE:-}" != mutect2 ] && [ "${CALLER_MODE:-}" != mt2 ]; then
        log_warn "人工体细胞过滤只适用于Mutect2，当前CALLER_MODE=${CALLER_MODE:-未配置}"
        return 0
    fi
    if [ -z "${MANUAL_FILTER_VEP_VCF:-}" ] && \
       { [ "${SKIP_VEP:-false}" = true ] || [ "${RUN_VEP:-true}" = false ]; }; then
        log_warn "VEP未启用且未指定MANUAL_FILTER_VEP_VCF，跳过人工硬过滤"
        return 0
    fi

    local input_vcf="${MANUAL_FILTER_VEP_VCF:-${DIR_ANNOTATION}/${SAMPLE_ID}.vep.vcf.gz}"
    local normalized_vcf="${DIR_ANNOTATION}/${SAMPLE_ID}.vep.manual_filter.normalized.tmp.vcf.gz"
    local plain_vcf="${DIR_ANNOTATION}/${SAMPLE_ID}.vep.manual_filtered.tmp.vcf"
    local output_vcf="${DIR_ANNOTATION}/${SAMPLE_ID}.vep.manual_filtered.vcf.gz"
    local audit_tsv="${DIR_ANNOTATION}/${SAMPLE_ID}.vep.manual_filter_audit.tsv"
    local summary_json="${DIR_ANNOTATION}/${SAMPLE_ID}.vep.manual_filter_summary.json"

    MANUAL_FILTER_TMP_NORMALIZED_VCF="${normalized_vcf}"
    MANUAL_FILTER_TMP_PLAIN_VCF="${plain_vcf}"

    check_file "VEP注释VCF" "${input_vcf}" || return 1
    check_tool "BCFtools" "${TOOL_BCFTOOLS}" || return 1
    check_tool "Python3" "${TOOL_PYTHON:-python3}" || return 1
    mkdir -p "${DIR_ANNOTATION}"

    trap cleanup_manual_filter_tmp EXIT

    run_cmd "VEP VCF拆分为双等位记录" \
        "${TOOL_BCFTOOLS} norm" \
        -m -any \
        -f "${REFERENCE_GENOME}" \
        -Oz \
        -o "${normalized_vcf}" \
        "${input_vcf}"

    run_cmd "按人工阈值筛选VEP注释变异" \
        "${TOOL_PYTHON:-python3}" "${SCRIPT_DIR}/manual_filter_vep.py" \
        --input-vcf "${normalized_vcf}" \
        --output-vcf "${plain_vcf}" \
        --audit-tsv "${audit_tsv}" \
        --summary-json "${summary_json}" \
        --tumor-sample "${TUMOR_SAMPLE_ID:-${SAMPLE_ID}}" \
        --normal-sample "${NORMAL_SAMPLE_ID:-}" \
        --min-tlod "${MANUAL_FILTER_MIN_TLOD:-6.3}" \
        --min-tumor-dp "${MANUAL_FILTER_MIN_TUMOR_DP:-20}" \
        --min-tumor-alt-reads "${MANUAL_FILTER_MIN_TUMOR_ALT_READS:-5}" \
        --min-tumor-af "${MANUAL_FILTER_MIN_TUMOR_AF:-0.02}" \
        --min-normal-dp "${MANUAL_FILTER_MIN_NORMAL_DP:-20}" \
        --max-normal-alt-reads "${MANUAL_FILTER_MAX_NORMAL_ALT_READS:-2}" \
        --max-normal-af "${MANUAL_FILTER_MAX_NORMAL_AF:-0.02}" \
        --max-population-af "${MANUAL_FILTER_MAX_POPULATION_AF:-0.001}" \
        --population-fields "${MANUAL_FILTER_POPULATION_AF_FIELDS:-MAX_AF,gnomADe_AF,gnomADg_AF,AF}"

    run_cmd "压缩人工过滤VCF" \
        "${TOOL_BCFTOOLS} view" -Oz -o "${output_vcf}" "${plain_vcf}"
    run_cmd "人工过滤VCF建立索引" \
        "${TOOL_BCFTOOLS} index" -f -t "${output_vcf}"

    log_info "VEP后人工硬过滤完成"
    log_info "过滤VCF: ${output_vcf}"
    log_info "逐位点审计表: ${audit_tsv}"
    log_info "参数与计数: ${summary_json}"
}

main "$@"
