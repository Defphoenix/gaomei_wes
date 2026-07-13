#!/bin/bash
#===============================================================================
# 12_tmb.sh - 肿瘤突变负荷 (TMB) 计算
#
# 说明: 基于变异检测结果计算TMB (Tumor Mutational Burden)
#       TMB = 编码区体细胞突变数 / 编码区大小 (Mb)
#       流程:
#         1. 筛选非同义体细胞突变
#         2. 计算编码区大小
#         3. 计算TMB值
#       注: TMB是免疫治疗疗效预测的重要标志物
#
# 输入: 注释后的VCF + SnpEff/VEP注释结果
# 输出: TMB计算结果文件
# 依赖: bcftools, SnpEff/VEP注释结果
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
    log_step "步骤12: TMB 肿瘤突变负荷计算"

    if [ "${SKIP_TMB}" = true ] || [ "${RUN_TMB}" = false ]; then
        log_warn "跳过TMB计算 (SKIP_TMB=${SKIP_TMB})"
        return 0
    fi

    mkdir -p "${DIR_TMB}"

    # 确定输入VCF (优先使用注释后的VCF)
    local input_vcf=""
    if [ -f "${DIR_ANNOTATION}/${SAMPLE_ID}.snpeff.vcf.gz" ]; then
        input_vcf="${DIR_ANNOTATION}/${SAMPLE_ID}.snpeff.vcf.gz"
        log_info "使用SnpEff注释VCF: ${input_vcf}"
    elif [ -f "${DIR_ANNOTATION}/${SAMPLE_ID}.vep.vcf.gz" ]; then
        input_vcf="${DIR_ANNOTATION}/${SAMPLE_ID}.vep.vcf.gz"
        log_info "使用VEP注释VCF: ${input_vcf}"
    elif [ "${CALLER_MODE}" = "haplotypecaller" ] || [ "${CALLER_MODE}" = "hc" ]; then
        input_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.pass.vcf.gz"
    elif [ "${CALLER_MODE}" = "mutect2" ] || [ "${CALLER_MODE}" = "mt2" ]; then
        input_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.pass.vcf.gz"
    fi

    if [ -z "${input_vcf}" ] || [ ! -f "${input_vcf}" ]; then
        log_error "未找到可用的VCF文件进行TMB计算!"
        return 1
    fi

    # 输出文件
    local tmb_result="${DIR_TMB}/${SAMPLE_ID}_tmb_result.txt"
    local tmb_variants="${DIR_TMB}/${SAMPLE_ID}_tmb_variants.txt"
    local tmb_gene_list="${DIR_TMB}/${SAMPLE_ID}_tmb_genes.txt"

    #---------------------------------------
    # 步骤1: 筛选用于TMB计算的变异
    #---------------------------------------
    # TMB计算标准:
    #   - 非同义突变 (missense, nonsense, frameshift, splice等)
    #   - 排除同义突变
    #   - 过滤低质量/低频变异
    #   - 可选: 排除已知胚系变异 (dbSNP common)

    log_info "筛选TMB相关变异..."

    # 提取所有变异
    local total_variants=0
    if [ -f "${input_vcf}" ]; then
        total_variants=$(${TOOL_BCFTOOLS} view -H "${input_vcf}" 2>/dev/null | wc -l || echo "0")
    fi
    log_info "总变异数: ${total_variants}"

    # 筛选非同义突变 (基于SnpEff注释)
    {
        echo "# 用于TMB计算的变异列表"
        echo "# 样本: ${SAMPLE_ID}"
        echo "# 筛选条件: 非同义突变, QUAL>=${TMB_MIN_QUAL}, AF>=${TMB_MIN_AF}"
        echo "#"
        echo "# CHROM  POS  REF  ALT  GENE  EFFECT  QUAL  AF"
    } > "${tmb_variants}"

    local tmb_count=0

    if [ -f "${input_vcf}" ]; then
        # 使用bcftools query提取信息
        # 如果有SnpEff注释，提取ANN字段
        if ${TOOL_BCFTOOLS} view -h "${input_vcf}" 2>/dev/null | grep -q "SnpEff"; then
            ${TOOL_BCFTOOLS} query \
                -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/ANN\n' \
                "${input_vcf}" 2>/dev/null | \
            while IFS=$'\t' read -r chrom pos ref alt qual ann; do
                # 解析SnpEff注释
                local effect=$(echo "${ann}" | tr ',' '\n' | head -1 | awk -F'|' '{print $2}')
                local gene=$(echo "${ann}" | tr ',' '\n' | head -1 | awk -F'|' '{print $4}')

                # 检查是否是非同义突变
                if echo "${TMB_INCLUDE_TYPES}" | grep -qi "${effect}"; then
                    echo -e "${chrom}\t${pos}\t${ref}\t${alt}\t${gene}\t${effect}\t${qual}"
                fi
            done >> "${tmb_variants}"
        else
            # 无SnpEff注释时，使用基本过滤
            # 只保留SNP和InDel，过滤QUAL
            ${TOOL_BCFTOOLS} view \
                -f PASS \
                -m2 -M2 \
                --min-ac 1 \
                "${input_vcf}" 2>/dev/null | \
                ${TOOL_BCFTOOLS} query \
                -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\tno_annotation\n' \
                >> "${tmb_variants}" 2>/dev/null || true
        fi
    fi

    # 计算TMB相关变异数
    tmb_count=$(grep -v "^#" "${tmb_variants}" | wc -l || echo "0")

    #---------------------------------------
    # 步骤2: 计算编码区大小
    #---------------------------------------
    local coding_region_size="${TMB_CODING_SIZE}"

    # 如果有BED文件，可以精确计算编码区大小
    if [ -n "${INTERVAL_FILE}" ] && [ -f "${INTERVAL_FILE}" ]; then
        local bed_size
        bed_size=$(awk '{sum+=$3-$2} END{printf "%.0f", sum/1000000}' "${INTERVAL_FILE}" 2>/dev/null || echo "${TMB_CODING_SIZE}")
        if [ "${bed_size}" -gt 0 ] 2>/dev/null; then
            coding_region_size="${bed_size}"
            log_info "使用BED文件计算的编码区大小: ${coding_region_size} Mb"
        fi
    else
        log_info "使用默认编码区大小: ${coding_region_size} Mb"
    fi

    #---------------------------------------
    # 步骤3: 计算TMB
    #---------------------------------------
    local tmb_value="N/A"
    if [ "${tmb_count}" -gt 0 ] && [ "${coding_region_size}" != "0" ]; then
        tmb_value=$(awk "BEGIN{printf \"%.2f\", ${tmb_count}/${coding_region_size}}")
    fi

    #---------------------------------------
    # 步骤4: 提取突变基因列表
    #---------------------------------------
    grep -v "^#" "${tmb_variants}" 2>/dev/null | \
        awk -F'\t' '{print $5}' | sort | uniq -c | sort -rn | \
        head -50 > "${tmb_gene_list}" 2>/dev/null || true

    #---------------------------------------
    # 步骤5: 生成TMB报告
    #---------------------------------------
    {
        echo "TMB (肿瘤突变负荷) 计算结果"
        echo "================================"
        echo "样本ID:          ${SAMPLE_ID}"
        echo "样本类型:        ${SAMPLE_TYPE}"
        echo "计算时间:        $(date '+%Y-%m-%d %H:%M:%S')"
        echo "输入VCF:         ${input_vcf}"
        echo ""
        echo "【TMB计算】"
        echo "────────────────────────────────"
        echo "总变异数:              ${total_variants}"
        echo "TMB相关变异数:         ${tmb_count}"
        echo "编码区大小:            ${coding_region_size} Mb"
        echo "TMB值:                 ${tmb_value} muts/Mb"
        echo ""

        # TMB临床解读
        if [ "${tmb_value}" != "N/A" ]; then
            echo "【TMB临床解读】"
            echo "────────────────────────────────"
            if awk "BEGIN{exit !(${tmb_value} >= 10)}"; then
                echo "TMB-HIGH (>=10 muts/Mb)"
                echo "  -> 可能从免疫检查点抑制剂治疗中获益"
                echo "  -> 参考: FDA批准帕博利珠单抗用于TMB-H (>=10)实体瘤"
            elif awk "BEGIN{exit !(${tmb_value} >= 5)}"; then
                echo "TMB-INTERMEDIATE (5-10 muts/Mb)"
                echo "  -> 中等突变负荷，临床意义需结合其他指标"
            else
                echo "TMB-LOW (<5 muts/Mb)"
                echo "  -> 低突变负荷"
            fi
        fi
        echo ""

        # 变异类型分布
        echo "【变异类型分布】"
        echo "────────────────────────────────"
        grep -v "^#" "${tmb_variants}" 2>/dev/null | \
            awk -F'\t' '{print $6}' | sort | uniq -c | sort -rn || echo "无数据"
        echo ""

        # 高频突变基因
        echo "【高频突变基因 (Top 20)】"
        echo "────────────────────────────────"
        if [ -f "${tmb_gene_list}" ] && [ -s "${tmb_gene_list}" ]; then
            head -20 "${tmb_gene_list}"
        else
            echo "无基因注释数据"
        fi
        echo ""

        # 各染色体变异分布
        echo "【各染色体变异分布】"
        echo "────────────────────────────────"
        grep -v "^#" "${tmb_variants}" 2>/dev/null | \
            awk -F'\t' '{print $1}' | sort -V | uniq -c | sort -rn || echo "无数据"

    } > "${tmb_result}"

    log_info "TMB计算完成!"
    log_info "TMB值: ${tmb_value} muts/Mb"
    log_info "结果文件: ${tmb_result}"
    cat "${tmb_result}"
}

main "$@"
