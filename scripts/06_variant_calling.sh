#!/bin/bash
#===============================================================================
# 06_variant_calling.sh - 变异检测
#
# 说明: 使用GATK进行变异检测
#       支持两种模式:
#       1. HaplotypeCaller - 胚系突变检测 (germline)
#       2. Mutect2         - 体细胞突变检测 (somatic, 需tumor-normal配对)
#       通过 config.sh 中 CALLER_MODE 参数切换
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/../config.sh}"
source "${CONFIG_FILE}"
source "${SCRIPT_DIR}/utils.sh"

#---------------------------------------
# GATK HaplotypeCaller 胚系变异检测
#---------------------------------------
run_haplotypecaller() {
    local input_bam="$1"
    local output_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.raw.vcf.gz"
    local output_idx="${output_vcf}.tbi"

    # 构建区间参数 (使用 :- 语法防止变量未定义报错)
    local interval_param=""
    if [ -n "${INTERVAL_LIST:-}" ] && [ -f "${INTERVAL_LIST:-}" ]; then
        interval_param="-L ${INTERVAL_LIST}"
        log_info "使用interval_list限制区域: ${INTERVAL_LIST}"
    elif [ -n "${INTERVAL_FILE:-}" ] && [ -f "${INTERVAL_FILE:-}" ]; then
        interval_param="-L ${INTERVAL_FILE}"
        log_info "使用BED文件限制区域: ${INTERVAL_FILE}"
    fi

    # GATK HaplotypeCaller
    # -R: 参考基因组
    # -I: 输入BAM
    # -O: 输出VCF
    # -L: 限定分析区域 (外显子组分析推荐)
    # --native-pair-hmm-threads: 线程数
    run_cmd "GATK HaplotypeCaller 胚系变异检测" \
        "${TOOL_GATK} --java-options \"-Xmx${GATK_JAVA_MEM}\" HaplotypeCaller" \
        -R "${REFERENCE_GENOME}" \
        -I "${input_bam}" \
        -O "${output_vcf}" \
        ${interval_param} \
        --native-pair-hmm-threads ${GATK_THREADS} \
        --dont-use-soft-clipped-bases \
        -stand-call-conf 30 \
        2>&1

    # 检查输出
    if [ ! -f "${output_vcf}" ]; then
        log_error "HaplotypeCaller 运行失败，VCF未生成"
        return 1
    fi

    log_info "HaplotypeCaller完成! 原始VCF: ${output_vcf}"
}

#---------------------------------------
# GATK Mutect2 体细胞变异检测
#---------------------------------------
run_mutect2() {
    local input_bam="$1"
    local output_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.raw.vcf.gz"
    local f1r2_tar="${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.f1r2.tar.gz"

    # 构建区间参数 (使用 :- 语法防止变量未定义报错)
    local interval_param=""
    if [ -n "${INTERVAL_LIST:-}" ] && [ -f "${INTERVAL_LIST:-}" ]; then
        interval_param="-L ${INTERVAL_LIST}"
    elif [ -n "${INTERVAL_FILE:-}" ] && [ -f "${INTERVAL_FILE:-}" ]; then
        interval_param="-L ${INTERVAL_FILE}"
    fi

    # Panel of Normals (可选，提高特异性)
    local pon_param=""
    if [ -n "${PANEL_OF_NORMALS:-}" ] && [ -f "${PANEL_OF_NORMALS:-}" ]; then
        pon_param="--panel-of-normals ${PANEL_OF_NORMALS}"
        log_info "使用Panel of Normals: ${PANEL_OF_NORMALS}"
    fi

    local germline_param=""
    if [ -n "${GERMLINE_RESOURCE_VCF:-}" ] && [ -f "${GERMLINE_RESOURCE_VCF:-}" ]; then
        germline_param="--germline-resource ${GERMLINE_RESOURCE_VCF}"
        log_info "使用germline resource: ${GERMLINE_RESOURCE_VCF}"
    fi

    local normal_param=""
    if [ -n "${NORMAL_BAM:-}" ] && bam_is_complete "${NORMAL_BAM}"; then
        check_file "配对Normal BAM" "${NORMAL_BAM}" || return 1
        if [ ! -f "${NORMAL_BAM}.bai" ]; then
            log_info "未检测到Normal BAM索引，正在使用 samtools index 补齐..."
            "${TOOL_SAMTOOLS}" index "${NORMAL_BAM}"
        fi
        check_file "配对Normal BAM索引" "${NORMAL_BAM}.bai" || return 1
        if [ -z "${NORMAL_SAMPLE_ID:-}" ]; then
            log_error "配置了 NORMAL_BAM 但 NORMAL_SAMPLE_ID 为空；Mutect2配对模式需要normal样本名"
            return 1
        fi
        normal_param="-I ${NORMAL_BAM} -normal ${NORMAL_SAMPLE_ID}"
        log_info "Mutect2配对模式: tumor=${SAMPLE_ID}, normal=${NORMAL_SAMPLE_ID}"
    else
        log_warn "未配置配对Normal BAM，Mutect2将以tumor-only模式运行"
    fi

    local max_reads_param=""
    if [ -n "${MUTECT2_MAX_READS_PER_ALIGNMENT_START:-}" ]; then
        if [[ ! "${MUTECT2_MAX_READS_PER_ALIGNMENT_START}" =~ ^[0-9]+$ ]]; then
            log_error "MUTECT2_MAX_READS_PER_ALIGNMENT_START必须为非负整数: ${MUTECT2_MAX_READS_PER_ALIGNMENT_START}"
            return 1
        fi
        max_reads_param="--max-reads-per-alignment-start ${MUTECT2_MAX_READS_PER_ALIGNMENT_START}"
        log_info "Mutect2每个alignment-start最大reads: ${MUTECT2_MAX_READS_PER_ALIGNMENT_START}"
    fi

    # GATK Mutect2
    # 体细胞突变检测，适用于肿瘤样本
    # 如果有配对正常样本，可用 -I normal 指定
    run_cmd "GATK Mutect2 体细胞变异检测" \
        "${TOOL_GATK} --java-options \"-Xmx${GATK_JAVA_MEM}\" Mutect2" \
        -R "${REFERENCE_GENOME}" \
        -I "${input_bam}" \
        ${normal_param} \
        -O "${output_vcf}" \
        --f1r2-tar-gz "${f1r2_tar}" \
        ${interval_param} \
        ${pon_param} \
        ${germline_param} \
        ${max_reads_param} \
        --minimum-mapping-quality 30 \
        --native-pair-hmm-threads ${GATK_THREADS} \
        ${MUTECT2_EXTRA_PARAMS:-} \
        2>&1

    # 检查输出
    if [ ! -f "${output_vcf}" ]; then
        log_error "Mutect2 运行失败，VCF未生成"
        return 1
    fi

    log_info "Mutect2完成! 原始VCF: ${output_vcf}"
    log_info "Mutect2 F1R2: ${f1r2_tar}"
}

#---------------------------------------
# 主函数
#---------------------------------------
main() {
    log_step "步骤6: 变异检测 (Variant Calling)"

    # 检查是否跳过
    if [ "${SKIP_VARIANT_CALLING}" = true ]; then
        log_warn "跳过变异检测 (SKIP_VARIANT_CALLING=true)"
        return 0
    fi

    # 确定输入BAM (优先使用去重后的BAM)
    local input_bam
    if [ -n "${TUMOR_BAM:-}" ] && bam_is_complete "${TUMOR_BAM}"; then
        input_bam="${TUMOR_BAM}"
        log_info "使用配置指定的Tumor BAM: ${input_bam}"
    elif bam_is_complete "${DIR_ALIGNED}/${SAMPLE_ID}.dedup.bam"; then
        input_bam="${DIR_ALIGNED}/${SAMPLE_ID}.dedup.bam"
        log_info "使用去重后BAM: ${input_bam}"
    elif bam_is_complete "${DIR_ALIGNED}/${SAMPLE_ID}.sorted.bam"; then
        input_bam="${DIR_ALIGNED}/${SAMPLE_ID}.sorted.bam"
        log_info "使用排序BAM (未去重): ${input_bam}"
    else
        log_error "未找到完整且通过samtools quickcheck的输入BAM文件!"
        return 1
    fi

    # ──────────────────────────────────────────────────────────
    # ✨ 【自动修复1】在检查索引之前，主动用 samtools 补齐标准索引
    # ──────────────────────────────────────────────────────────
    if [ ! -f "${input_bam}.bai" ]; then
        log_info "未检测到标准的 .bam.bai 索引，正在使用 samtools index 补齐..."
        ${TOOL_SAMTOOLS} index "${input_bam}"
    fi

    check_file "输入BAM" "${input_bam}" || return 1
    check_file "输入BAM索引" "${input_bam}.bai" || return 1
    check_file "参考基因组" "${REFERENCE_GENOME}" || return 1

    # 确保输出目录存在
    mkdir -p "${DIR_VARIANTS}"

    # 根据配置选择变异检测方法
    case "${CALLER_MODE}" in
        haplotypecaller|hc)
            log_info "变异检测模式: HaplotypeCaller (胚系突变)"
            run_haplotypecaller "${input_bam}"
            ;;
        mutect2|mt2)
            log_info "变异检测模式: Mutect2 (体细胞突变)"
            run_mutect2 "${input_bam}"
            ;;
        *)
            log_error "未知的CALLER_MODE: ${CALLER_MODE}"
            log_error "支持: haplotypecaller / mutect2"
            return 1
            ;;
    esac

    log_info "变异检测完成!"
}

main "$@"
