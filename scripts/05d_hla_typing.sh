#!/bin/bash
# HLA*LA typing from an indexed WES BAM. Full calls and binding-compatible alleles are separate outputs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/../config.sh}"
source "${CONFIG_FILE}"
source "${SCRIPT_DIR}/utils.sh"

main() {
    log_step "步骤5d: HLA*LA高分辨率HLA分型"

    if [ "${SKIP_HLA_TYPING:-false}" = true ] || [ "${RUN_HLA_TYPING:-auto}" = false ]; then
        log_warn "跳过HLA分型"
        return 0
    fi

    local graph_dir="${HLA_LA_GRAPH_DIR:-}"
    local graph_name="${HLA_LA_GRAPH_NAME:-PRG_MHC_GRCh38_withIMGT}"
    local hla_bin="${TOOL_HLA_LA:-HLA-LA.pl}"
    local required="${HLA_TYPING_REQUIRED:-false}"
    if [ -z "${graph_dir}" ] || [ ! -d "${graph_dir}" ]; then
        log_warn "HLA*LA graph目录不存在: ${graph_dir:-未配置}"
        log_warn "请准备PRG_MHC_GRCh38_withIMGT图数据并设置 HLA_LA_GRAPH_DIR"
        [ "${required}" = true ] && return 1
        return 0
    fi
    if ! command -v "${hla_bin}" >/dev/null 2>&1; then
        log_warn "HLA*LA未安装或不在PATH中: ${hla_bin}"
        [ "${required}" = true ] && return 1
        return 0
    fi

    local hla_prefix="${HLA_TYPING_ENV_PREFIX:-}"
    if [ -z "${hla_prefix}" ]; then
        hla_prefix=$(cd "$(dirname "$(command -v "${hla_bin}")")/.." && pwd)
    fi
    local hla_install_root="${HLA_LA_INSTALL_ROOT:-${hla_prefix}/opt/hla-la}"
    local installed_graph="${hla_install_root}/graphs/${graph_name}"
    mkdir -p "${hla_install_root}/graphs" || {
        log_error "无法创建HLA*LA graph目录: ${hla_install_root}/graphs"
        return 1
    }
    if [ ! -e "${installed_graph}" ]; then
        ln -s "${graph_dir}" "${installed_graph}"
        log_info "已将参考graph链接到HLA*LA安装目录: ${installed_graph}"
    fi
    if [ ! -f "${installed_graph}/serializedGRAPH" ]; then
        log_warn "HLA*LA graph尚未prepare，缺少: ${installed_graph}/serializedGRAPH"
        log_warn "请先运行 scripts/prepare_hlala_graph.sh"
        [ "${required}" = true ] && return 1
        return 0
    fi

    local input_bam="${HLA_TYPING_BAM:-}"
    if [ -z "${input_bam}" ]; then
        input_bam=$(get_final_bam) || return 1
    fi
    check_file "HLA分型BAM" "${input_bam}" || return 1
    if [ ! -f "${input_bam}.bai" ]; then
        "${TOOL_SAMTOOLS}" index "${input_bam}"
    fi

    local working_dir="${DIR_HLA_TYPING}/working"
    local raw_result="${working_dir}/${SAMPLE_ID}/hla/R1_bestguess_G.txt"
    local normalized="${DIR_HLA_TYPING}/${SAMPLE_ID}_hla_typing.tsv"
    local binding_alleles="${DIR_HLA_TYPING}/${SAMPLE_ID}_hla_binding_alleles.txt"
    mkdir -p "${DIR_HLA_TYPING}" "${working_dir}"

    run_cmd "HLA*LA WES HLA分型" \
        "${hla_bin}" \
        --BAM "${input_bam}" \
        --graph "${graph_name}" \
        --sampleID "${SAMPLE_ID}" \
        --maxThreads "${HLA_TYPING_THREADS:-${GATK_THREADS:-4}}" \
        --workingDir "${working_dir}"
    check_file "HLA*LA G-group结果" "${raw_result}" || return 1

    run_cmd "标准化HLA*LA结果" \
        "${TOOL_PYTHON:-python3}" "${SCRIPT_DIR}/parse_hlala_results.py" \
        --input "${raw_result}" \
        --output "${normalized}" \
        --alleles-output "${binding_alleles}" \
        --sample "${SAMPLE_ID}"

    log_info "HLA完整/G-group分型: ${normalized}"
    log_info "HLA-I binding兼容等位基因: ${binding_alleles}"
}

main "$@"
