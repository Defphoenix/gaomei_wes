#!/bin/bash
#===============================================================================
# run_pipeline.sh - WES/WGS 突变分析流程主控脚本
#
# 支持全流程运行、断点续跑、单步运行、环境检查和关键结果状态查看。
# 兼容 macOS 自带 Bash 3.x，不依赖 Bash 4 关联数组。
#===============================================================================

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${PROJECT_DIR}/config.sh"

show_short_usage() {
    echo "用法: bash run_pipeline.sh [--config FILE] [check|list|status|dry-run|step N|through N|from N]"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --help|-h|help)
            ACTION="help"
            shift
            break
            ;;
        --*)
            echo "未知参数: $1" >&2
            show_short_usage
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

if [ ! -f "${CONFIG_FILE}" ]; then
    echo "配置文件不存在: ${CONFIG_FILE}" >&2
    exit 1
fi

source "${CONFIG_FILE}"
source "${PROJECT_DIR}/scripts/utils.sh"

PIPELINE_START_TIME=$(date +%s)
STEP_ORDER="1 2 3 4 5 5b 5c 5d 6 7 7b 7c 7e 7d 8 9 10 11 12 13"

step_script() {
    case "$1" in
        1) echo "01_fastqc.sh" ;;
        2) echo "02_trim.sh" ;;
        3) echo "03_align.sh" ;;
        4) echo "04_sort_index.sh" ;;
        5) echo "05_mark_duplicates.sh" ;;
        5b) echo "05b_post_align_qc.sh" ;;
        5c) echo "05c_bqsr.sh" ;;
        5d) echo "05d_hla_typing.sh" ;;
        6) echo "06_variant_calling.sh" ;;
        7) echo "07_variant_filter.sh" ;;
        7b) echo "07b_snpeff.sh" ;;
        7c) echo "07c_vep.sh" ;;
        7e) echo "07e_manual_filter.sh" ;;
        7d) echo "07d_neoantigen.sh" ;;
        8) echo "08_cnv.sh" ;;
        9) echo "09_msi.sh" ;;
        10) echo "10_sv.sh" ;;
        11) echo "11_coverage.sh" ;;
        12) echo "12_tmb.sh" ;;
        13) echo "13_final_summary.sh" ;;
        *) return 1 ;;
    esac
}

step_name() {
    case "$1" in
        1) echo "FastQC 原始数据质控" ;;
        2) echo "数据修剪 (fastp/Trimmomatic)" ;;
        3) echo "BWA-MEM 序列比对" ;;
        4) echo "BAM排序和索引" ;;
        5) echo "标记PCR重复" ;;
        5b) echo "比对后QC统计 (Picard/Samtools)" ;;
        5c) echo "碱基质量重校准 (BQSR)" ;;
        5d) echo "HLA*LA高分辨率HLA分型" ;;
        6) echo "变异检测 (GATK HC/Mutect2)" ;;
        7) echo "变异过滤" ;;
        7b) echo "SnpEff 功能注释" ;;
        7c) echo "VEP 功能注释" ;;
        7e) echo "VEP后人工体细胞硬过滤" ;;
        7d) echo "新抗原候选肽生成与HLA结合预测" ;;
        8) echo "CNV拷贝数分析 (CNVkit) / 深度QC (mosdepth)" ;;
        9) echo "MSI微卫星不稳定性检测 (MSIsensor-pro)" ;;
        10) echo "结构变异检测 (Manta)" ;;
        11) echo "覆盖度分析 (mosdepth/bedtools)" ;;
        12) echo "TMB肿瘤突变负荷计算" ;;
        13) echo "MultiQC汇总报告" ;;
        *) return 1 ;;
    esac
}

step_count() {
    set -- ${STEP_ORDER}
    echo "$#"
}

usage() {
    local total
    total=$(step_count)
    show_short_usage
    echo ""
    echo "命令:"
    echo "  check             检查运行环境和输入"
    echo "  list              列出流程步骤"
    echo "  status [N]        查看关键输出，或查看单步输出"
    echo "  step <N>          只运行指定步骤，例如 7d"
    echo "  through <N>       从步骤1运行到指定步骤，例如预处理 through 5d"
    echo "  from <N>          从指定步骤继续运行"
    echo "  dry-run           只检查，不执行分析"
    echo "  无命令            运行完整流程"
    echo ""
    echo "步骤列表 (共${total}步):"
    list_steps
}

list_steps() {
    local i=1
    echo "────────────────────────────────────────────────────────"
    for step_id in ${STEP_ORDER}; do
        printf "  %2d. [%-3s] %s\n" "${i}" "${step_id}" "$(step_name "${step_id}")"
        i=$((i + 1))
    done
    echo "────────────────────────────────────────────────────────"
}

do_check() {
    print_banner
    mkdir -p "${DIR_LOGS}"

    local has_error=0
    log_step "环境检查"

    print_runtime_config
    check_config_consistency || has_error=1
    check_core_tools || has_error=1
    check_analysis_tools || true
    check_reference || has_error=1
    check_input || has_error=1
    check_disk_space "${PROJECT_DIR}" "${MIN_DISK_GB:-100}" || has_error=1

    if [ -n "${INTERVAL_FILE:-}" ] && [ "${INTERVAL_FILE}" != "/path/to/target_regions.bed" ]; then
        check_file "区间BED文件" "${INTERVAL_FILE}" || has_error=1
    else
        log_warn "区间文件未配置"
    fi

    echo ""
    if [ ${has_error} -eq 0 ]; then
        log_info "环境检查通过!"
        return 0
    fi

    log_error "环境检查存在问题，请修复后重试"
    return 1
}

run_step() {
    local step_id="$1"
    local script_name
    if ! script_name=$(step_script "${step_id}"); then
        log_error "未知步骤: ${step_id}"
        list_steps
        return 1
    fi

    local script_path="${PROJECT_DIR}/scripts/${script_name}"
    if [ ! -f "${script_path}" ]; then
        log_error "脚本不存在: ${script_path}"
        return 1
    fi

    mkdir -p "${DIR_LOGS}"
    log_step "运行步骤 ${step_id}: $(step_name "${step_id}")"
    CONFIG_FILE="${CONFIG_FILE}" bash "${script_path}"
}

run_full_pipeline() {
    print_banner
    do_check || {
        log_error "环境检查未通过，流程终止"
        return 1
    }

    local total i step_id
    total=$(step_count)
    i=1
    log_step "开始运行完整流程 (${total}步)"

    for step_id in ${STEP_ORDER}; do
        echo ""
        echo -e "${BLUE}──────────────────────────────────────────────────────────────${NC}"
        echo -e "${BLUE}  步骤 ${i}/${total} [${step_id}]: $(step_name "${step_id}")${NC}"
        echo -e "${BLUE}──────────────────────────────────────────────────────────────${NC}"

        if ! run_step "${step_id}"; then
            log_error "步骤 [${step_id}] ($(step_name "${step_id}")) 运行失败!"
            log_error "修复后可使用 'bash run_pipeline.sh from ${step_id}' 继续"
            return 1
        fi
        i=$((i + 1))
    done

    print_finish "${PIPELINE_START_TIME}"
}

run_from_step() {
    local start_id="$1"
    local found=false
    local step_id

    print_banner
    do_check || {
        log_error "环境检查未通过"
        return 1
    }

    for step_id in ${STEP_ORDER}; do
        if [ "${step_id}" = "${start_id}" ]; then
            found=true
        fi
        if [ "${found}" = true ]; then
            if ! run_step "${step_id}"; then
                log_error "步骤 [${step_id}] ($(step_name "${step_id}")) 运行失败!"
                return 1
            fi
        fi
    done

    if [ "${found}" = false ]; then
        log_error "未知步骤: ${start_id}"
        list_steps
        return 1
    fi

    print_finish "${PIPELINE_START_TIME}"
}

run_through_step() {
    local end_id="$1"
    local found=false
    local step_id

    print_banner
    do_check || {
        log_error "环境检查未通过"
        return 1
    }

    for step_id in ${STEP_ORDER}; do
        if ! run_step "${step_id}"; then
            log_error "步骤 [${step_id}] ($(step_name "${step_id}")) 运行失败!"
            return 1
        fi
        if [ "${step_id}" = "${end_id}" ]; then
            found=true
            break
        fi
    done

    if [ "${found}" = false ]; then
        log_error "未知结束步骤: ${end_id}"
        list_steps
        return 1
    fi

    print_finish "${PIPELINE_START_TIME}"
}

show_file_status() {
    local label="$1"
    local filepath="$2"
    if [ -f "${filepath}" ]; then
        printf "  [OK]   %-18s %s\n" "${label}" "${filepath}"
    else
        printf "  [MISS] %-18s %s\n" "${label}" "${filepath}"
    fi
}

show_status() {
    local step_id="${1:-all}"
    echo "样本: ${SAMPLE_ID}"
    echo "结果目录: ${RESULT_DIR}"
    echo ""

    case "${step_id}" in
        5d)
            show_file_status "HLA typing" "${DIR_HLA_TYPING}/${SAMPLE_ID}_hla_typing.tsv"
            show_file_status "Binding alleles" "${DIR_HLA_TYPING}/${SAMPLE_ID}_hla_binding_alleles.txt"
            ;;
        6)
            show_file_status "HC raw VCF" "${DIR_VARIANTS}/${SAMPLE_ID}.raw.vcf.gz"
            show_file_status "Mutect2 raw" "${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.raw.vcf.gz"
            show_file_status "Mutect2 F1R2" "${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.f1r2.tar.gz"
            ;;
        7)
            show_file_status "HC PASS VCF" "${DIR_VARIANTS}/${SAMPLE_ID}.pass.vcf.gz"
            show_file_status "Mutect2 filtered" "${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.filtered.vcf.gz"
            show_file_status "Mutect2 PASS" "${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.pass.vcf.gz"
            show_file_status "Orientation model" "${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.read-orientation-model.tar.gz"
            show_file_status "Contamination" "${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.contamination.table"
            show_file_status "Tumor segments" "${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.segments.table"
            ;;
        7c)
            show_file_status "VEP VCF" "${DIR_ANNOTATION}/${SAMPLE_ID}.vep.vcf.gz"
            show_file_status "VEP TSV" "${DIR_ANNOTATION}/${SAMPLE_ID}.vep.tsv"
            ;;
        7e)
            show_file_status "人工过滤VEP VCF" "${DIR_ANNOTATION}/${SAMPLE_ID}.vep.manual_filtered.vcf.gz"
            show_file_status "人工过滤审计表" "${DIR_ANNOTATION}/${SAMPLE_ID}.vep.manual_filter_audit.tsv"
            show_file_status "人工过滤摘要" "${DIR_ANNOTATION}/${SAMPLE_ID}.vep.manual_filter_summary.json"
            ;;
        8)
            show_file_status "CNV CNR" "${DIR_CNV}/${SAMPLE_ID}.cnr"
            show_file_status "CNV call" "${DIR_CNV}/${SAMPLE_ID}.call.cns"
            show_file_status "Depth QC" "${DIR_CNV}/${SAMPLE_ID}.depth_qc.tsv"
            show_file_status "CNV summary" "${DIR_CNV}/${SAMPLE_ID}_cnv_summary.txt"
            ;;
        9)
            show_file_status "MSI result" "${DIR_MSI}/${SAMPLE_ID}_msi_result.txt"
            show_file_status "MSI summary" "${DIR_MSI}/${SAMPLE_ID}_msi_summary.txt"
            ;;
        12)
            show_file_status "TMB report" "${DIR_TMB}/${SAMPLE_ID}_tmb_result.txt"
            show_file_status "TMB accepted" "${DIR_TMB}/${SAMPLE_ID}_tmb_accepted_variants.tsv"
            show_file_status "TMB rejected" "${DIR_TMB}/${SAMPLE_ID}_tmb_rejected_variants.tsv"
            show_file_status "TMB JSON" "${DIR_TMB}/${SAMPLE_ID}_tmb_summary.json"
            ;;
        7d)
            show_file_status "neo FASTA" "${DIR_NEOANTIGEN}/${SAMPLE_ID}_neoantigen_peptides.fa"
            show_file_status "neo manifest" "${DIR_NEOANTIGEN}/${SAMPLE_ID}_neoantigen_manifest.tsv"
            show_file_status "HLA binding" "${DIR_NEOANTIGEN}/${SAMPLE_ID}_hla_binding.tsv"
            ;;
        all)
            show_file_status "dedup BAM" "${DIR_ALIGNED}/${SAMPLE_ID}.dedup.bam"
            show_file_status "HLA typing" "${DIR_HLA_TYPING}/${SAMPLE_ID}_hla_typing.tsv"
            show_file_status "PASS VCF" "${DIR_VARIANTS}/${SAMPLE_ID}.pass.vcf.gz"
            show_file_status "Mutect2 PASS" "${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.pass.vcf.gz"
            show_file_status "Contamination" "${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.contamination.table"
            show_file_status "VEP VCF" "${DIR_ANNOTATION}/${SAMPLE_ID}.vep.vcf.gz"
            show_file_status "人工过滤VEP VCF" "${DIR_ANNOTATION}/${SAMPLE_ID}.vep.manual_filtered.vcf.gz"
            show_file_status "neo FASTA" "${DIR_NEOANTIGEN}/${SAMPLE_ID}_neoantigen_peptides.fa"
            show_file_status "CNV call" "${DIR_CNV}/${SAMPLE_ID}.call.cns"
            show_file_status "CNV depth QC" "${DIR_CNV}/${SAMPLE_ID}.depth_qc.tsv"
            show_file_status "MSI summary" "${DIR_MSI}/${SAMPLE_ID}_msi_summary.txt"
            show_file_status "TMB" "${DIR_TMB}/${SAMPLE_ID}_tmb_result.txt"
            show_file_status "summary" "${DIR_SUMMARY}/${SAMPLE_ID}_final_report.txt"
            ;;
        *)
            log_warn "目前 status 内置 all、5d、6、7、7c、7e、7d、8、9、12；步骤 ${step_id} 可直接检查对应结果目录"
            ;;
    esac
}

main() {
    mkdir -p "${DIR_LOGS}"
    local action="${1:-full}"

    case "${action}" in
        check) do_check ;;
        step)
            if [ -z "${2:-}" ]; then
                log_error "请指定步骤编号"
                usage
                exit 1
            fi
            run_step "$2"
            ;;
        from)
            if [ -z "${2:-}" ]; then
                log_error "请指定起始步骤"
                usage
                exit 1
            fi
            run_from_step "$2"
            ;;
        through)
            if [ -z "${2:-}" ]; then
                log_error "请指定结束步骤"
                usage
                exit 1
            fi
            run_through_step "$2"
            ;;
        list) list_steps ;;
        status) show_status "${2:-all}" ;;
        dry-run)
            log_info "模拟运行模式"
            do_check
            log_info "模拟完成"
            ;;
        help|--help|-h) usage ;;
        full) run_full_pipeline ;;
        *)
            log_error "未知操作: ${action}"
            usage
            exit 1
            ;;
    esac
}

main "$@"
