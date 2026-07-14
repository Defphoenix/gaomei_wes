#!/bin/bash
#===============================================================================
# 09_msi.sh - MSIsensor-pro 微卫星不稳定性检测
#
# 支持:
#   1. tumor-normal paired: msisensor-pro msi
#   2. tumor-only/single sample: msisensor-pro pro
#   3. 无可用位点列表时: smoke summary，不冒充正式 MSI 结果
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/../config.sh}"
source "${CONFIG_FILE}"
source "${SCRIPT_DIR}/utils.sh"

resolve_msi_list() {
    local out_list="${DIR_MSI}/msisensor_reference.list"

    if [ -n "${MSISENSOR2_LIST:-}" ] && [ -f "${MSISENSOR2_LIST}" ]; then
        echo "${MSISENSOR2_LIST}"
        return 0
    fi

    if [ "${MSISENSOR_SCAN_REFERENCE:-false}" = true ]; then
        echo "${out_list}"
        return 0
    fi

    echo ""
    return 0
}

write_smoke_summary() {
    local input_bam="$1"
    local reason="$2"
    local summary="${DIR_MSI}/${SAMPLE_ID}_msi_summary.txt"
    local result="${DIR_MSI}/${SAMPLE_ID}_msi_result.txt"
    local call="${DIR_MSI}/${SAMPLE_ID}_msi_call.tsv"

    {
        echo -e "sample\tmode\tstatus\treason\tbam"
        echo -e "${SAMPLE_ID}\tsmoke\tNOT_RUN\t${reason}\t${input_bam}"
    } > "${result}"

    {
        echo -e "sample\tmode\tmsi_status\tmsi_score_percent\tunstable_sites\ttotal_sites\tlow_threshold_percent\thigh_threshold_percent\treason\tresult_files"
        echo -e "${SAMPLE_ID}\tsmoke\tNOT_DETERMINED\t\t\t\t${MSI_STATUS_LOW_THRESHOLD:-3}\t${MSI_STATUS_HIGH_THRESHOLD:-20}\t${reason}\t${result}"
    } > "${call}"

    {
        echo "MSI检测摘要 - ${SAMPLE_ID}"
        echo "================================"
        echo "模式: smoke"
        echo "状态: NOT_RUN"
        echo "原因: ${reason}"
        echo ""
        echo "说明: 这是流程连通性结果，不是正式MSI判定。"
        echo "正式运行需要提供 MSISENSOR2_LIST，或设置 MSISENSOR_SCAN_REFERENCE=true 先扫描参考基因组。"
        echo "标准化判定表: ${call}"
    } > "${summary}"
}

summarize_msi() {
    local output_prefix="$1"
    local mode="$2"
    local summary="${DIR_MSI}/${SAMPLE_ID}_msi_summary.txt"
    local result_file=""

    if [ -f "${output_prefix}" ]; then
        result_file="${output_prefix}"
    elif [ -f "${output_prefix}.tsv" ]; then
        result_file="${output_prefix}.tsv"
    elif [ -f "${output_prefix}_dis" ]; then
        result_file="${output_prefix}_dis"
    else
        result_file="${output_prefix}"
    fi

    {
        echo "MSI检测摘要 - ${SAMPLE_ID}"
        echo "================================"
        echo "检测时间: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "模式: ${mode}"
        echo "结果前缀: ${output_prefix}"
        echo ""
        echo "结果文件:"
        ls -1 "${output_prefix}"* 2>/dev/null | sed 's/^/  /' || echo "  未找到结果文件"
        echo ""
        if [ -f "${result_file}" ]; then
            echo "结果预览:"
            sed -n '1,20p' "${result_file}" | sed 's/^/  /'
        fi
    } > "${summary}"
}

write_msi_call() {
    local output_prefix="$1"
    local mode="$2"
    local call="${DIR_MSI}/${SAMPLE_ID}_msi_call.tsv"
    local call_summary="${DIR_MSI}/${SAMPLE_ID}_msi_call_summary.txt"

    run_cmd "解析MSI结果并生成标准化判定表" \
        "${TOOL_PYTHON:-python3}" "${SCRIPT_DIR}/msi_call_from_msisensor.py" \
        --prefix "${output_prefix}" \
        --sample "${SAMPLE_ID}" \
        --mode "${mode}" \
        --out "${call}" \
        --summary "${call_summary}" \
        --low-threshold "${MSI_STATUS_LOW_THRESHOLD:-3}" \
        --high-threshold "${MSI_STATUS_HIGH_THRESHOLD:-20}"
}

main() {
    log_step "步骤9: MSIsensor-pro 微卫星不稳定性检测"

    if [ "${SKIP_MSI:-true}" = true ] || [ "${RUN_MSI:-false}" = false ]; then
        log_warn "跳过MSI检测 (SKIP_MSI=${SKIP_MSI:-true}, RUN_MSI=${RUN_MSI:-false})"
        return 0
    fi

    local input_bam
    input_bam=$(get_final_bam) || return 1
    check_file "输入BAM" "${input_bam}" || return 1
    check_file "输入BAM索引" "${input_bam}.bai" || return 1

    mkdir -p "${DIR_MSI}"

    local mode="${MSISENSOR_MODE:-auto}"
    if [ "${mode}" = "smoke" ]; then
        write_smoke_summary "${input_bam}" "forced_smoke_mode"
        log_info "MSI smoke summary已生成: ${DIR_MSI}/${SAMPLE_ID}_msi_summary.txt"
        return 0
    fi

    check_tool "msisensor-pro" "${TOOL_MSISENSOR2:-msisensor-pro}" || return 1

    local msi_list
    msi_list=$(resolve_msi_list)
    if [ "${MSISENSOR_SCAN_REFERENCE:-false}" = true ] && [ -n "${msi_list}" ] && [ ! -f "${msi_list}" ]; then
        run_cmd "MSIsensor-pro 扫描微卫星位点" \
            "${TOOL_MSISENSOR2:-msisensor-pro}" scan \
            -d "${REFERENCE_GENOME}" \
            -o "${msi_list}" \
            -l "${MSI_MIN_LENGTH:-10}"
    fi

    if [ -z "${msi_list}" ] || [ ! -f "${msi_list}" ]; then
        if [ "${MSI_SMOKE_TEST:-true}" = true ]; then
            log_warn "未提供MSI位点列表，生成smoke summary"
            write_smoke_summary "${input_bam}" "missing_msisensor_list"
            return 0
        fi
        log_error "未提供MSI位点列表: 请配置 MSISENSOR2_LIST 或设置 MSISENSOR_SCAN_REFERENCE=true"
        return 1
    fi

    local output_prefix="${DIR_MSI}/${SAMPLE_ID}_msi"

    if [ "${mode}" = "auto" ]; then
        if [ "${SAMPLE_TYPE:-}" = "tumor" ] && [ -n "${NORMAL_BAM:-}" ] && [ -f "${NORMAL_BAM}" ]; then
            mode="paired"
        else
            mode="tumor_only"
        fi
    fi

    case "${mode}" in
        paired)
            check_file "Normal BAM" "${NORMAL_BAM}" || return 1
            check_file "Normal BAM索引" "${NORMAL_BAM}.bai" || return 1
            run_cmd "MSIsensor-pro paired MSI" \
                "${TOOL_MSISENSOR2:-msisensor-pro}" msi \
                -d "${msi_list}" \
                -n "${NORMAL_BAM}" \
                -t "${input_bam}" \
                -o "${output_prefix}" \
                -e "${INTERVAL_FILE:-}" \
                -c "${MSI_MIN_COVERAGE:-15}" \
                -b "${MSISENSOR_THREADS:-2}"
            ;;
        tumor_only)
            run_cmd "MSIsensor-pro tumor-only MSI" \
                "${TOOL_MSISENSOR2:-msisensor-pro}" pro \
                -d "${msi_list}" \
                -t "${input_bam}" \
                -o "${output_prefix}" \
                -e "${INTERVAL_FILE:-}" \
                -i "${MSISENSOR_TUMOR_ONLY_THRESHOLD:-0.1}" \
                -c "${MSI_MIN_COVERAGE:-15}" \
                -b "${MSISENSOR_THREADS:-2}"
            ;;
        *)
            log_error "未知 MSISENSOR_MODE: ${mode}"
            return 1
            ;;
    esac

    summarize_msi "${output_prefix}" "${mode}"
    write_msi_call "${output_prefix}" "${mode}"

    log_info "MSI检测完成!"
    log_info "MSI摘要: ${DIR_MSI}/${SAMPLE_ID}_msi_summary.txt"
    log_info "MSI标准化判定表: ${DIR_MSI}/${SAMPLE_ID}_msi_call.tsv"
}

main "$@"
