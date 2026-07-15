#!/bin/bash
#===============================================================================
# 07d_neoantigen.sh - 从 VEP 注释 VCF 生成新抗原候选肽，并可选运行 HLA 结合预测
#
# 输入: VEP 注释后的 VCF + 蛋白 FASTA
# 输出: peptide FASTA + manifest TSV + 可选 netMHCpan 结果
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/../config.sh}"
source "${CONFIG_FILE}"
source "${SCRIPT_DIR}/utils.sh"

main() {
    log_step "步骤7d: 新抗原候选肽生成与HLA结合预测"

    if [ "${SKIP_NEOANTIGEN:-true}" = true ] || [ "${RUN_NEOANTIGEN:-false}" = false ]; then
        log_warn "跳过新抗原分析 (SKIP_NEOANTIGEN=${SKIP_NEOANTIGEN:-true}, RUN_NEOANTIGEN=${RUN_NEOANTIGEN:-false})"
        return 0
    fi

    local input_vcf="${NEOANTIGEN_VEP_VCF:-}"
    if [ -z "${input_vcf}" ]; then
        input_vcf="${DIR_ANNOTATION}/${SAMPLE_ID}.vep.vcf.gz"
    fi

    local protein_fasta="${NEOANTIGEN_PROTEIN_FASTA:-${PROJECT_DIR}/reference/protein.fa}"
    local output_fasta="${DIR_NEOANTIGEN}/${SAMPLE_ID}_neoantigen_peptides.fa"
    local manifest="${DIR_NEOANTIGEN}/${SAMPLE_ID}_neoantigen_manifest.tsv"
    local protein_table="${DIR_NEOANTIGEN}/${SAMPLE_ID}_variant_proteins.tsv"
    local peptide_table="${DIR_NEOANTIGEN}/${SAMPLE_ID}_neoantigen_peptides.tsv"
    local fasta_dir="${DIR_NEOANTIGEN}/fasta_by_mer"
    local binding_out="${DIR_NEOANTIGEN}/${SAMPLE_ID}_hla_binding.tsv"

    mkdir -p "${DIR_NEOANTIGEN}" "${DIR_LOGS}"

    check_file "VEP注释VCF" "${input_vcf}" || return 1
    check_file "蛋白FASTA" "${protein_fasta}" || return 1
    check_tool "Python3" "${TOOL_PYTHON:-python3}" || return 1

    run_cmd "生成新抗原候选肽FASTA" \
        "${TOOL_PYTHON:-python3}" "${SCRIPT_DIR}/neoantigen_from_vep.py" \
        --vep-vcf "${input_vcf}" \
        --protein-fasta "${protein_fasta}" \
        --output-fasta "${output_fasta}" \
        --manifest "${manifest}" \
        --protein-table "${protein_table}" \
        --peptide-table "${peptide_table}" \
        --fasta-dir "${fasta_dir}" \
        ${NEOANTIGEN_ANNOVAR_TXT:+--annovar-txt "${NEOANTIGEN_ANNOVAR_TXT}"} \
        --sample "${SAMPLE_ID}" \
        --lengths "${NEOANTIGEN_PEPTIDE_LENGTHS:-8,9,10,11,12,13,14,15}" \
        --flank "${NEOANTIGEN_PEPTIDE_FLANK:-30}"

    if [ ! -s "${output_fasta}" ]; then
        log_warn "没有生成候选肽；请查看 manifest 中的 skipped/no_novel_peptide 原因: ${manifest}"
        return 0
    fi

    if [ "${RUN_HLA_BINDING:-false}" = false ]; then
        log_info "HLA结合预测未启用，候选肽FASTA已生成: ${output_fasta}"
        log_info "候选肽按mer拆分目录: ${fasta_dir}"
        log_info "突变蛋白明细表: ${protein_table}"
        log_info "候选肽明细表: ${peptide_table}"
        log_info "如需预测，配置 HLA_ALLELES 并设置 RUN_HLA_BINDING=true"
        return 0
    fi

    local resolved_hla_alleles="${HLA_ALLELES:-}"
    local typing_alleles_file="${HLA_TYPING_ALLELES_FILE:-${DIR_HLA_TYPING:-${RESULT_DIR}/hla_typing}/${SAMPLE_ID}_hla_binding_alleles.txt}"
    if [ -z "${resolved_hla_alleles}" ] && [ -s "${typing_alleles_file}" ]; then
        resolved_hla_alleles=$(tr -d '[:space:]' < "${typing_alleles_file}")
        log_info "使用自动HLA分型的binding等位基因: ${typing_alleles_file}"
    fi

    if [ -z "${resolved_hla_alleles}" ]; then
        if [ "${RUN_HLA_BINDING:-false}" = true ]; then
            log_error "RUN_HLA_BINDING=true 但未配置HLA_ALLELES且无自动分型结果"
            return 1
        fi
        log_warn "没有可用HLA等位基因，自动跳过HLA结合预测"
        return 0
    fi

    if [ "${RUN_HLA_BINDING:-false}" = auto ] && \
       ! command -v "${TOOL_NETMHCPAN:-netMHCpan}" >/dev/null 2>&1 && \
       ! command -v "${TOOL_MHCFLURRY:-mhcflurry-predict}" >/dev/null 2>&1; then
        log_warn "未安装netMHCpan或MHCflurry，自动跳过HLA结合预测"
        return 0
    fi

    local allele_file="${DIR_NEOANTIGEN}/${SAMPLE_ID}_hla_alleles.txt"
    tr ',' '\n' <<< "${resolved_hla_alleles}" > "${allele_file}"

    run_cmd "HLA结合预测 (${HLA_BINDING_TOOL:-auto})" \
        "${TOOL_PYTHON:-python3}" "${SCRIPT_DIR}/hla_binding_predict.py" \
        --peptides "${output_fasta}" \
        --alleles "${resolved_hla_alleles}" \
        --output "${binding_out}" \
        --tool "${HLA_BINDING_TOOL:-auto}" \
        --netmhcpan-bin "${TOOL_NETMHCPAN:-netMHCpan}" \
        --mhcflurry-bin "${TOOL_MHCFLURRY:-mhcflurry-predict}" \
        --threshold-nm "${HLA_BINDING_PREDICTION_THRESHOLD_NM:-500}"

    log_info "新抗原分析完成!"
    log_info "候选肽FASTA: ${output_fasta}"
    log_info "候选肽按mer拆分目录: ${fasta_dir}"
    log_info "候选肽manifest: ${manifest}"
    log_info "突变蛋白明细表: ${protein_table}"
    log_info "候选肽明细表: ${peptide_table}"
    log_info "HLA结合预测: ${binding_out}"
}

main "$@"
