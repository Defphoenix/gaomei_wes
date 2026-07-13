#!/bin/bash
#===============================================================================
# 07c_vep.sh - Ensembl VEP 变异注释
#
# 说明: 使用Ensembl VEP (Variant Effect Predictor)进行变异注释
#       VEP提供更丰富的功能注释，包括:
#       - 基因/转录本/蛋白质影响
#       - SIFT/PolyPhen预测
#       - CADD/REVEL评分
#       - 已有数据库交叉引用 (ClinVar, COSMIC, gnomAD等)
#
# 输入: 过滤后的VCF文件
# 输出: 注释VCF + TSV表格 + HTML报告
# 依赖: Ensembl VEP
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/../config.sh}"
source "${CONFIG_FILE}"
source "${SCRIPT_DIR}/utils.sh"

resolve_vep_cache_dir() {
    local species="${VEP_SPECIES:-homo_sapiens}"
    local version="${VEP_CACHE_VERSION:-}"
    local assembly="${VEP_ASSEMBLY:-GRCh38}"
    local candidate
    local candidates=()

    [ -n "${VEP_CACHE_DIR:-}" ] && candidates+=("${VEP_CACHE_DIR}")
    [ -n "${REFERENCE_DIR:-}" ] && candidates+=("${REFERENCE_DIR}/vep_cache" "${REFERENCE_DIR}/.vep" "${REFERENCE_DIR}")
    [ -n "${HOME:-}" ] && candidates+=("${HOME}/.vep")

    for candidate in "${candidates[@]}"; do
        [ -n "${candidate}" ] || continue
        if [ -d "${candidate}/${species}" ]; then
            echo "${candidate}"
            return 0
        fi
        if [ "$(basename "${candidate}")" = "${species}" ] && [ -d "${candidate}" ]; then
            dirname "${candidate}"
            return 0
        fi
    done

    log_error "未找到VEP离线缓存目录。需要存在类似: <cache_root>/${species}/${version}_${assembly}"
    log_error "当前检查过: ${candidates[*]:-未配置}"
    log_error "请在config中设置 VEP_CACHE_DIR 为包含 ${species}/ 的父目录，例如: /PUBLIC/.../reference_data/vep_cache"
    return 1
}

#---------------------------------------
# 主函数
#---------------------------------------
main() {
    log_step "步骤7c: Ensembl VEP 变异注释"

    if [ "${SKIP_VEP}" = true ] || [ "${RUN_VEP}" = false ]; then
        log_warn "跳过VEP注释 (SKIP_VEP=${SKIP_VEP})"
        return 0
    fi

    # 确定输入VCF
    local input_vcf="${VEP_INPUT_VCF:-}"
    if [ -n "${input_vcf}" ]; then
        log_info "使用配置指定的VEP输入VCF: ${input_vcf}"
    elif [ "${CALLER_MODE}" = "haplotypecaller" ] || [ "${CALLER_MODE}" = "hc" ]; then
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
    local vep_vcf="${DIR_ANNOTATION}/${SAMPLE_ID}.vep.vcf.gz"
    local vep_tsv="${DIR_ANNOTATION}/${SAMPLE_ID}.vep.tsv"
    local vep_html="${DIR_ANNOTATION}/${SAMPLE_ID}.vep_summary.html"
    local vep_warnings="${DIR_ANNOTATION}/${SAMPLE_ID}.vep_warnings.txt"

    local first_variant=""
    first_variant=$("${TOOL_BCFTOOLS}" view -H "${input_vcf}" 2>/dev/null | awk 'NR==1{print; exit}' || true)
    if [ -z "${first_variant}" ]; then
        log_warn "输入VCF没有变异记录，跳过VEP实质注释并生成空结果: ${input_vcf}"
        run_cmd "生成空VEP VCF" \
            "${TOOL_BCFTOOLS} view" \
            -Oz \
            -o "${vep_vcf}" \
            "${input_vcf}"
        run_cmd "空VEP VCF建立索引" \
            "${TOOL_BCFTOOLS} index" -t "${vep_vcf}"
        {
            echo -e "Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence"
        } > "${vep_tsv}"
        {
            echo "<html><body><h1>VEP skipped</h1><p>No variant records in input VCF.</p></body></html>"
        } > "${vep_html}"
        log_info "VEP空结果生成完成: ${vep_vcf}"
        return 0
    fi

    # 构建VEP参数
    local resolved_cache_dir
    resolved_cache_dir=$(resolve_vep_cache_dir) || return 1
    local vep_cache_params="--dir_cache ${resolved_cache_dir}"
    log_info "使用VEP cache目录: ${resolved_cache_dir}"

    local vep_fasta_params=""
    if [ -n "${VEP_FASTA:-}" ] && [ -f "${VEP_FASTA}" ]; then
        vep_fasta_params="--fasta ${VEP_FASTA}"
    fi

    # 插件参数
    local plugin_params=""
    if [ -n "${VEP_PLUGINS}" ]; then
        # 根据配置的插件添加参数
        IFS=',' read -ra PLUGINS <<< "${VEP_PLUGINS}"
        for plugin in "${PLUGINS[@]}"; do
            case ${plugin} in
                dbNSFP)
                    plugin_params="${plugin_params} --plugin dbNSFP"
                    ;;
                REVEL)
                    plugin_params="${plugin_params} --plugin REVEL"
                    ;;
                CADD)
                    plugin_params="${plugin_params} --plugin CADD"
                    ;;
            esac
        done
    fi

    #---------------------------------------
    # 步骤1: VEP 核心注释
    #---------------------------------------
    # --everything: 启用所有注释 (SIFT, PolyPhen, CADD等)
    # --vcf: 输出VCF格式
    # --tsv: 同时输出TSV表格
    # --html: 生成HTML统计报告
    # --fork: 并行线程数
    # --cache: 使用离线缓存
    # --assembly: 基因组版本
    run_cmd "VEP 变异注释 (${VEP_ASSEMBLY})" \
        "${TOOL_VEP}" \
        --input_file "${input_vcf}" \
        --format vcf \
        --output_file "${vep_vcf}" \
        --vcf \
        --compress_output bgzip \
        --stats_file "${vep_html}" \
        --force_overwrite \
        --species "${VEP_SPECIES}" \
        --assembly "${VEP_ASSEMBLY}" \
        --cache \
        --offline \
        ${vep_cache_params} \
        --cache_version "${VEP_CACHE_VERSION}" \
        ${vep_fasta_params} \
        --fork ${GATK_THREADS} \
        --everything \
        --no_progress \
        --check_existing \
        --clinvar \
        --pubmed \
        --per_gene \
        --symbol \
        --canonical \
        --biotype \
        --protein \
        --uniprot \
        --hgvs \
        --hgvsg \
        --variant_class \
        ${VEP_EXTRA_PARAMS} \
        ${plugin_params} || \
    {
        log_warn "VEP完整注释失败，尝试简化模式"
        run_cmd "VEP 简化注释" \
            "${TOOL_VEP}" \
            --input_file "${input_vcf}" \
            --format vcf \
            --output_file "${vep_vcf}" \
            --vcf \
            --compress_output bgzip \
            --stats_file "${vep_html}" \
            --force_overwrite \
            --species "${VEP_SPECIES}" \
            --assembly "${VEP_ASSEMBLY}" \
            --cache \
            --offline \
            ${vep_cache_params} \
            --cache_version "${VEP_CACHE_VERSION}" \
            ${vep_fasta_params} \
            --fork ${GATK_THREADS} \
            --symbol \
            --canonical \
            --biotype
    }

    # 建立索引
    if [ -f "${vep_vcf}" ]; then
        run_cmd "VEP VCF建立索引" \
            "${TOOL_BCFTOOLS} index" -t "${vep_vcf}"
    fi

    #---------------------------------------
    # 步骤2: 生成TSV表格 (便于查看)
    #---------------------------------------
    if [ -f "${vep_vcf}" ]; then
        run_cmd "VEP生成TSV表格" \
            "${TOOL_VEP}" \
            --input_file "${input_vcf}" \
            --format vcf \
            --output_file "${vep_tsv}" \
            --force_overwrite \
            --species "${VEP_SPECIES}" \
            --assembly "${VEP_ASSEMBLY}" \
            --cache \
            --offline \
            ${vep_cache_params} \
            --cache_version "${VEP_CACHE_VERSION}" \
            ${vep_fasta_params} \
            --everything \
            --no_progress || \
            log_warn "VEP TSV生成失败 (非致命)"
    fi

    #---------------------------------------
    # 步骤3: 合并SnpEff和VEP注释 (可选)
    #---------------------------------------
    local snpeff_vcf="${DIR_ANNOTATION}/${SAMPLE_ID}.snpeff.vcf.gz"
    if [ -f "${snpeff_vcf}" ] && [ -f "${vep_vcf}" ]; then
        log_info "SnpEff和VEP注释均已完成"
        log_info "SnpEff结果: ${snpeff_vcf}"
        log_info "VEP结果: ${vep_vcf}"
        log_info "可使用 SnpSift 或 bcftools annotate 合并注释"
    fi

    log_info "VEP注释完成!"
    log_info "注释VCF: ${vep_vcf}"
    log_info "TSV表格: ${vep_tsv}"
    log_info "HTML报告: ${vep_html}"
}

main "$@"
