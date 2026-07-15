#!/bin/bash
# Strict TMB calculation from VEP CSQ, paired sample evidence and an effective coding denominator.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/../config.sh}"
source "${CONFIG_FILE}"
source "${SCRIPT_DIR}/utils.sh"

main() {
    log_step "步骤12: 基于VEP注释的严格TMB计算"
    if [ "${SKIP_TMB:-false}" = true ] || [ "${RUN_TMB:-false}" = false ]; then
        log_warn "跳过TMB计算"
        return 0
    fi

    local input_vcf="${TMB_VEP_VCF:-${DIR_ANNOTATION}/${SAMPLE_ID}.vep.vcf.gz}"
    local coding_bed="${TMB_EFFECTIVE_CODING_BED:-}"
    local accepted="${DIR_TMB}/${SAMPLE_ID}_tmb_accepted_variants.tsv"
    local rejected="${DIR_TMB}/${SAMPLE_ID}_tmb_rejected_variants.tsv"
    local summary_json="${DIR_TMB}/${SAMPLE_ID}_tmb_summary.json"
    local summary_tsv="${DIR_TMB}/${SAMPLE_ID}_tmb_summary.tsv"
    local report="${DIR_TMB}/${SAMPLE_ID}_tmb_result.txt"
    local bed_param=""

    mkdir -p "${DIR_TMB}"
    check_file "VEP注释VCF" "${input_vcf}" || return 1
    check_tool "Python3" "${TOOL_PYTHON:-python3}" || return 1
    if [ -n "${coding_bed}" ]; then
        check_file "TMB有效编码区BED" "${coding_bed}" || return 1
        bed_param="--effective-coding-bed ${coding_bed}"
        if [ "${TMB_DENOMINATOR_VALIDATED:-false}" != true ]; then
            log_warn "当前BED尚未标记为经过panel验证的有效编码区；TMB仅作为研发结果"
        fi
    else
        log_warn "未配置TMB_EFFECTIVE_CODING_BED，使用配置的固定分母 ${TMB_CODING_SIZE:-0} Mb"
    fi

    run_cmd "严格筛选VEP体细胞突变并计算TMB" \
        "${TOOL_PYTHON:-python3}" "${SCRIPT_DIR}/tmb_from_vep.py" \
        --vcf "${input_vcf}" \
        --accepted "${accepted}" \
        --rejected "${rejected}" \
        --summary-json "${summary_json}" \
        --summary-tsv "${summary_tsv}" \
        --tumor-sample "${TUMOR_SAMPLE_ID:-${SAMPLE_ID}}" \
        --normal-sample "${NORMAL_SAMPLE_ID:-}" \
        ${bed_param} \
        --denominator-mb "${TMB_CODING_SIZE:-0}" \
        --denominator-validated "${TMB_DENOMINATOR_VALIDATED:-false}" \
        --min-qual "${TMB_MIN_QUAL:-0}" \
        --min-tlod "${TMB_MIN_TLOD:-6.3}" \
        --min-tumor-dp "${TMB_MIN_TUMOR_DP:-20}" \
        --min-tumor-alt-reads "${TMB_MIN_TUMOR_ALT_READS:-5}" \
        --min-tumor-af "${TMB_MIN_AF:-0.05}" \
        --min-normal-dp "${TMB_MIN_NORMAL_DP:-10}" \
        --max-normal-alt-reads "${TMB_MAX_NORMAL_ALT_READS:-2}" \
        --max-normal-af "${TMB_MAX_NORMAL_AF:-0.02}" \
        --max-population-af "${TMB_MAX_POPULATION_AF:-0.001}" \
        --consequences "${TMB_VEP_CONSEQUENCES}" \
        --population-fields "${TMB_POPULATION_AF_FIELDS:-MAX_AF,gnomADe_AF,gnomADg_AF,AF}"

    "${TOOL_PYTHON:-python3}" - "${summary_json}" > "${report}" <<'PY'
import json
import sys
data = json.load(open(sys.argv[1], encoding="utf-8"))
tmb = float(data["tmb_mutations_per_mb"])
level = "TMB-H" if tmb >= 10 else ("TMB-中等" if tmb >= 5 else "TMB-L")
print("TMB（肿瘤突变负荷）严格计算结果")
print("=" * 42)
print(f"输入VCF: {data['input_vcf']}")
print(f"肿瘤样本: {data['tumor_sample']}")
print(f"配对正常样本: {data['normal_sample'] or '无'}")
print(f"纳入突变数: {data['accepted_variants']}")
print(f"排除突变数: {data['rejected_variants']}")
print(f"有效编码区: {data['denominator_mb']:.6f} Mb ({data['denominator_source']})")
print(f"分母验证状态: {'已验证' if data['denominator_validated'] else '未验证（研发用途）'}")
print(f"TMB: {tmb:.4f} mutations/Mb")
print(f"分层: {level}")
print("说明: 分层阈值仅作流程展示，临床使用必须按癌种、panel和验证方案校准。")
PY

    log_info "TMB结果: ${report}"
    log_info "纳入明细: ${accepted}"
    log_info "排除及原因: ${rejected}"
    cat "${report}"
}

main "$@"
