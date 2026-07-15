#!/bin/bash
#===============================================================================
# 10_sv.sh - Manta 结构变异检测
#
# 说明: 使用Manta检测结构变异 (SV)
#       包括: 缺失、插入、倒位、易位、重复等
#       支持单样本和tumor-normal配对模式
#
# 输入: 去重后的BAM
# 输出: 候选SV VCF + 体细胞SV VCF
# 依赖: Manta
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
    log_step "步骤10: Manta 结构变异检测"

    if [ "${SKIP_SV}" = true ] || [ "${RUN_SV}" = false ]; then
        log_warn "跳过SV检测 (SKIP_SV=${SKIP_SV})"
        return 0
    fi

    # 确定输入BAM
    local input_bam
    input_bam=$(get_final_bam) || return 1
    check_file "输入BAM" "${input_bam}" || return 1

    mkdir -p "${DIR_SV}"

    # Manta工作目录
    local manta_workdir="${DIR_SV}/manta_work"

    #---------------------------------------
    # 步骤1: 配置Manta分析流程
    #---------------------------------------
    local config_cmd=""
    if [ -f "${TOOL_MANTA}" ]; then
        config_cmd="${TOOL_MANTA}"
    elif command -v "${TOOL_MANTA:-configManta.py}" &>/dev/null; then
        config_cmd="${TOOL_MANTA:-configManta.py}"
    else
        log_error "Manta未找到! 请检查 TOOL_MANTA 配置"
        return 1
    fi

    # 构建配置参数
    local config_params="--bam=${input_bam} --referenceFasta=${REFERENCE_GENOME} --runDir=${manta_workdir}"

    # 额外配置
    if [ -n "${MANTA_CONFIG}" ] && [ -f "${MANTA_CONFIG}" ]; then
        config_params="${config_params} --config=${MANTA_CONFIG}"
    fi

    # 胚系模式 vs 体细胞模式
    if [ "${SAMPLE_TYPE}" = "tumor" ] && [ -n "${NORMAL_BAM}" ] && [ -f "${NORMAL_BAM}" ]; then
        # 体细胞模式
        config_params="${config_params} --tumorBam=${input_bam} --normalBam=${NORMAL_BAM}"
        log_info "Manta体细胞模式 (tumor-normal配对)"
    else
        # 胚系模式
        config_params="--bam=${input_bam} --referenceFasta=${REFERENCE_GENOME} --runDir=${manta_workdir}"
        log_info "Manta胚系模式 (单样本)"
    fi

    #---------------------------------------
    # 步骤2: 配置Manta
    #---------------------------------------
    run_cmd "Manta 配置分析流程" \
        ${config_cmd} ${config_params}

    #---------------------------------------
    # 步骤3: 运行Manta
    #---------------------------------------
    run_cmd "Manta 运行结构变异检测" \
        "${manta_workdir}/runWorkflow.py" \
        -m local \
        -j ${MANTA_THREADS:-2}

    #---------------------------------------
    # 步骤4: 整理结果
    #---------------------------------------
    local manta_results="${manta_workdir}/results/variants"

    # 胚系候选SV
    if [ -f "${manta_results}/candidateSV.vcf.gz" ]; then
        cp "${manta_results}/candidateSV.vcf.gz" "${DIR_SV}/${SAMPLE_ID}.candidateSV.vcf.gz"
        cp "${manta_results}/candidateSV.vcf.gz.tbi" "${DIR_SV}/${SAMPLE_ID}.candidateSV.vcf.gz.tbi" 2>/dev/null || true
        log_info "候选SV: ${DIR_SV}/${SAMPLE_ID}.candidateSV.vcf.gz"
    fi

    # 胚系diploid SV
    if [ -f "${manta_results}/diploidSV.vcf.gz" ]; then
        cp "${manta_results}/diploidSV.vcf.gz" "${DIR_SV}/${SAMPLE_ID}.diploidSV.vcf.gz"
        cp "${manta_results}/diploidSV.vcf.gz.tbi" "${DIR_SV}/${SAMPLE_ID}.diploidSV.vcf.gz.tbi" 2>/dev/null || true
        log_info "Diploid SV: ${DIR_SV}/${SAMPLE_ID}.diploidSV.vcf.gz"
    fi

    # 体细胞SV (如果有)
    if [ -f "${manta_results}/somaticSV.vcf.gz" ]; then
        cp "${manta_results}/somaticSV.vcf.gz" "${DIR_SV}/${SAMPLE_ID}.somaticSV.vcf.gz"
        cp "${manta_results}/somaticSV.vcf.gz.tbi" "${DIR_SV}/${SAMPLE_ID}.somaticSV.vcf.gz.tbi" 2>/dev/null || true
        log_info "体细胞SV: ${DIR_SV}/${SAMPLE_ID}.somaticSV.vcf.gz"
    fi

    #---------------------------------------
    # 步骤5: SV统计
    #---------------------------------------
    {
        echo "结构变异检测结果 - ${SAMPLE_ID}"
        echo "================================"
        echo "检测工具: Manta"
        echo "样本类型: ${SAMPLE_TYPE}"
        echo ""

        for vcf_file in "${DIR_SV}"/*.vcf.gz; do
            if [ -f "${vcf_file}" ]; then
                local name=$(basename "${vcf_file}")
                local count=$(${TOOL_BCFTOOLS} view -H "${vcf_file}" 2>/dev/null | wc -l || echo "N/A")
                echo "文件: ${name}"
                echo "  SV数量: ${count}"

                # 按SV类型统计
                echo "  SV类型分布:"
                ${TOOL_BCFTOOLS} view -H "${vcf_file}" 2>/dev/null | \
                    grep -oP 'SVTYPE=\w+' | sort | uniq -c | sort -rn | \
                    while read cnt type; do
                        echo "    ${type}: ${cnt}"
                    done || echo "    无法统计"
                echo ""
            fi
        done
    } > "${DIR_SV}/${SAMPLE_ID}_sv_summary.txt"

    log_info "Manta SV检测完成!"
    log_info "结果目录: ${DIR_SV}"
    cat "${DIR_SV}/${SAMPLE_ID}_sv_summary.txt"
}

main "$@"
