#!/bin/bash
#===============================================================================
# 13_final_summary.sh - MultiQC 汇总报告 + 最终统计
#
# 说明: 使用MultiQC汇总所有QC报告
#       整合所有分析模块的结果，生成最终的综合报告
#       包括: 质控、比对、变异、新抗原、CNV、SV、MSI、覆盖度、TMB
#
# 输入: 各步骤的输出文件和QC指标
# 输出: MultiQC HTML报告 + 综合分析报告
# 依赖: MultiQC
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
    log_step "步骤13: MultiQC汇总报告 + 最终统计"

    mkdir -p "${DIR_SUMMARY}"
    mkdir -p "${DIR_MULTIQC}"

    #---------------------------------------
    # 步骤1: MultiQC 汇总所有QC报告
    #---------------------------------------
    if [ "${SKIP_MULTIQC}" = false ]; then
        if command -v "${TOOL_MULTIQC}" &>/dev/null; then
            log_info "运行MultiQC汇总..."

            # MultiQC扫描所有结果目录
            run_cmd "MultiQC 汇总QC报告" \
                "${TOOL_MULTIQC}" \
                "${DIR_FASTQC}" \
                "${DIR_ALIGNED}" \
                "${DIR_POSTQC}" \
                "${DIR_VARIANTS}" \
                "${DIR_COVERAGE}" \
                "${DIR_BQSR}" \
                "${DIR_ANNOTATION}" \
                -o "${DIR_MULTIQC}" \
                -n "${SAMPLE_ID}_multiqc" \
                --title "突变分析报告 - ${SAMPLE_ID}" \
                --force \
                -m fastqc \
                -m samtools \
                -m picard \
                -m qualimap \
                -m mosdepth \
                -m bcftools \
                2>/dev/null || \
            {
                log_warn "MultiQC部分模块运行失败，尝试简化模式"
                run_cmd "MultiQC 简化模式" \
                    "${TOOL_MULTIQC}" \
                    "${DIR_FASTQC}" \
                    "${DIR_ALIGNED}" \
                    "${DIR_POSTQC}" \
                    -o "${DIR_MULTIQC}" \
                    -n "${SAMPLE_ID}_multiqc" \
                    --force || \
                    log_warn "MultiQC运行失败"
            }

            log_info "MultiQC报告: ${DIR_MULTIQC}/${SAMPLE_ID}_multiqc_report.html"
        else
            log_warn "MultiQC不可用，跳过汇总QC报告"
        fi
    fi

    #---------------------------------------
    # 步骤2: 生成综合分析报告
    #---------------------------------------
    local final_report="${DIR_SUMMARY}/${SAMPLE_ID}_final_report.txt"

    {
        echo "=========================================================================="
        echo "           突变分析综合报告 - ${SAMPLE_ID}"
        echo "           生成时间: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "=========================================================================="
        echo ""
        echo "  样本信息:"
        echo "    样本ID:    ${SAMPLE_ID}"
        echo "    样本类型:  ${SAMPLE_TYPE}"
        echo "    检测模式:  ${CALLER_MODE}"
        echo "    参考基因组: ${REFERENCE_GENOME_VERSION}"
        echo ""

        #=======================================
        # 模块1: 测序质控
        #=======================================
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "  [模块1] 测序质控"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo ""

        # FastQC结果
        echo "  FastQC:"
        local fq_r1="${DIR_FASTQC}/$(basename ${FASTQ_R1%.fastq.gz})_fastqc.html"
        local fq_r2="${DIR_FASTQC}/$(basename ${FASTQ_R2%.fastq.gz})_fastqc.html"
        [ -f "${fq_r1}" ] && echo "    R1报告: ${fq_r1}" || echo "    R1报告: 未生成"
        [ -f "${fq_r2}" ] && echo "    R2报告: ${fq_r2}" || echo "    R2报告: 未生成"

        # fastp结果
        if [ "${USE_FASTP}" = true ]; then
            local fastp_json="${DIR_TRIMMED}/${SAMPLE_ID}_fastp.json"
            [ -f "${fastp_json}" ] && echo "  fastp报告: ${fastp_json}"
        fi
        echo ""

        #=======================================
        # 模块2: 比对统计
        #=======================================
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "  [模块2] 比对统计"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo ""

        local flagstat="${DIR_ALIGNED}/${SAMPLE_ID}.flagstat.txt"
        if [ -f "${flagstat}" ]; then
            echo "  比对统计:"
            cat "${flagstat}" | sed 's/^/    /'
        fi

        # 插入片段大小
        local insert_metrics="${DIR_POSTQC}/${SAMPLE_ID}.insert_size_metrics"
        if [ -f "${insert_metrics}" ]; then
            echo ""
            echo "  插入片段大小:"
            grep -v "^#" "${insert_metrics}" | head -1 | awk '{print "    中位数: "$5" bp"}' 2>/dev/null || true
        fi

        # 重复率
        local dup_metrics="${DIR_ALIGNED}/${SAMPLE_ID}.dup_metrics.txt"
        if [ -f "${dup_metrics}" ]; then
            local dup_pct
            dup_pct=$(grep -v "^#" "${dup_metrics}" 2>/dev/null | awk '{print $8}' | tail -1)
            [ -n "${dup_pct}" ] && echo "  PCR重复率: ${dup_pct}"
        fi
        echo ""

        #=======================================
        # 模块3: BQSR
        #=======================================
        if [ "${RUN_BQSR}" = true ] && [ "${SKIP_BQSR}" = false ]; then
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "  [模块3] 碱基质量重校准 (BQSR)"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo ""
            local bqsr_bam="${DIR_BQSR}/${SAMPLE_ID}.bqsr.bam"
            [ -f "${bqsr_bam}" ] && echo "  BQSR BAM: ${bqsr_bam}" || echo "  BQSR: 未完成"
            echo ""
        fi

        #=======================================
        # 模块4: 变异检测与注释
        #=======================================
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "  [模块4] 变异检测与注释"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo ""

        # 变异统计
        local final_vcf=""
        if [ "${CALLER_MODE}" = "haplotypecaller" ] || [ "${CALLER_MODE}" = "hc" ]; then
            final_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.pass.vcf.gz"
        elif [ "${CALLER_MODE}" = "mutect2" ] || [ "${CALLER_MODE}" = "mt2" ]; then
            final_vcf="${DIR_VARIANTS}/${SAMPLE_ID}.mutect2.pass.vcf.gz"
        fi

        if [ -f "${final_vcf}" ]; then
            local total_var snp_count indel_count
            total_var=$(${TOOL_BCFTOOLS} view -H "${final_vcf}" 2>/dev/null | wc -l || echo "N/A")
            snp_count=$(${TOOL_BCFTOOLS} view -H -v snps "${final_vcf}" 2>/dev/null | wc -l || echo "N/A")
            indel_count=$(${TOOL_BCFTOOLS} view -H -v indels "${final_vcf}" 2>/dev/null | wc -l || echo "N/A")

            echo "  变异检测 (${CALLER_MODE}):"
            echo "    总变异数: ${total_var}"
            echo "    SNP:      ${snp_count}"
            echo "    InDel:    ${indel_count}"
        fi

        # 注释结果
        local snpeff_vcf="${DIR_ANNOTATION}/${SAMPLE_ID}.snpeff.vcf.gz"
        local vep_vcf="${DIR_ANNOTATION}/${SAMPLE_ID}.vep.vcf.gz"
        [ -f "${snpeff_vcf}" ] && echo "  SnpEff注释: ${snpeff_vcf}"
        [ -f "${vep_vcf}" ] && echo "  VEP注释: ${vep_vcf}"
        echo ""

        #=======================================
        # 模块4b: 新抗原
        #=======================================
        if [ "${RUN_NEOANTIGEN:-false}" = true ] && [ "${SKIP_NEOANTIGEN:-true}" = false ]; then
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "  [模块4b] 新抗原候选肽"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo ""
            local neo_fasta="${DIR_NEOANTIGEN}/${SAMPLE_ID}_neoantigen_peptides.fa"
            local neo_manifest="${DIR_NEOANTIGEN}/${SAMPLE_ID}_neoantigen_manifest.tsv"
            local neo_proteins="${DIR_NEOANTIGEN}/${SAMPLE_ID}_variant_proteins.tsv"
            local neo_peptides="${DIR_NEOANTIGEN}/${SAMPLE_ID}_neoantigen_peptides.tsv"
            local neo_by_mer="${DIR_NEOANTIGEN}/fasta_by_mer"
            local hla_binding="${DIR_NEOANTIGEN}/${SAMPLE_ID}_hla_binding.tsv"
            if [ -f "${neo_fasta}" ]; then
                local neo_count
                neo_count=$(grep -c "^>" "${neo_fasta}" 2>/dev/null || echo "0")
                echo "  候选肽数量: ${neo_count}"
                echo "  候选肽FASTA: ${neo_fasta}"
                [ -f "${neo_manifest}" ] && echo "  候选肽manifest: ${neo_manifest}"
                [ -f "${neo_proteins}" ] && echo "  突变蛋白明细表: ${neo_proteins}"
                [ -f "${neo_peptides}" ] && echo "  候选肽明细表: ${neo_peptides}"
                [ -d "${neo_by_mer}" ] && echo "  按mer拆分FASTA目录: ${neo_by_mer}"
                [ -f "${hla_binding}" ] && echo "  HLA结合预测: ${hla_binding}"
            else
                echo "  新抗原分析: 未完成"
            fi
            echo ""
        fi

        #=======================================
        # 模块5: CNV分析
        #=======================================
        if [ "${RUN_CNV}" = true ] && [ "${SKIP_CNV}" = false ]; then
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "  [模块5] 拷贝数变异 (CNV)"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo ""
            local cnv_summary="${DIR_CNV}/${SAMPLE_ID}_cnv_summary.txt"
            if [ -f "${cnv_summary}" ]; then
                cat "${cnv_summary}" | sed 's/^/  /'
            else
                echo "  CNV分析: 未完成"
            fi
            echo ""
        fi

        #=======================================
        # 模块6: 结构变异
        #=======================================
        if [ "${RUN_SV}" = true ] && [ "${SKIP_SV}" = false ]; then
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "  [模块6] 结构变异 (SV)"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo ""
            local sv_summary="${DIR_SV}/${SAMPLE_ID}_sv_summary.txt"
            if [ -f "${sv_summary}" ]; then
                cat "${sv_summary}" | sed 's/^/  /'
            else
                echo "  SV检测: 未完成"
            fi
            echo ""
        fi

        #=======================================
        # 模块7: MSI检测
        #=======================================
        if [ "${RUN_MSI}" = true ] && [ "${SKIP_MSI}" = false ]; then
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "  [模块7] 微卫星不稳定性 (MSI)"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo ""
            local msi_call="${DIR_MSI}/${SAMPLE_ID}_msi_call.tsv"
            local msi_call_summary="${DIR_MSI}/${SAMPLE_ID}_msi_call_summary.txt"
            local msi_summary="${DIR_MSI}/${SAMPLE_ID}_msi_summary.txt"
            if [ -f "${msi_call_summary}" ]; then
                cat "${msi_call_summary}" | sed 's/^/  /'
                [ -f "${msi_call}" ] && echo "  标准化判定表: ${msi_call}"
            elif [ -f "${msi_summary}" ]; then
                cat "${msi_summary}" | sed 's/^/  /'
            else
                echo "  MSI检测: 未完成"
            fi
            echo ""
        fi

        #=======================================
        # 模块8: 覆盖度
        #=======================================
        if [ "${RUN_COVERAGE}" = true ] && [ "${SKIP_COVERAGE}" = false ]; then
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "  [模块8] 覆盖度分析"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo ""
            local cov_report="${DIR_COVERAGE}/${SAMPLE_ID}_coverage_report.txt"
            if [ -f "${cov_report}" ]; then
                cat "${cov_report}" | sed 's/^/  /'
            else
                echo "  覆盖度分析: 未完成"
            fi
            echo ""
        fi

        #=======================================
        # 模块9: TMB
        #=======================================
        if [ "${RUN_TMB}" = true ] && [ "${SKIP_TMB}" = false ]; then
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "  [模块9] 肿瘤突变负荷 (TMB)"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo ""
            local tmb_result="${DIR_TMB}/${SAMPLE_ID}_tmb_result.txt"
            if [ -f "${tmb_result}" ]; then
                # 只提取关键信息
                grep -E "TMB值|TMB-HIGH|TMB-LOW|TMB-INTERMEDIATE|TMB相关变异数|编码区大小" "${tmb_result}" | \
                    sed 's/^/  /'
            else
                echo "  TMB计算: 未完成"
            fi
            echo ""
        fi

        #=======================================
        # 输出文件清单
        #=======================================
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "  输出文件清单"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo ""
        echo "  结果根目录: ${RESULT_DIR}"
        echo ""

        for dir_name in fastqc trimmed aligned bqsr variants annotation neoantigen cnv sv msi coverage tmb summary multiqc; do
            local dir_path="${RESULT_DIR}/${dir_name}"
            if [ -d "${dir_path}" ]; then
                echo "  [${dir_name}]"
                ls -lh "${dir_path}"/*.* 2>/dev/null | awk '{print "    "$NF" ("$5")"}' | head -10 || echo "    (空)"
                echo ""
            fi
        done

        echo "=========================================================================="
        echo "  报告生成完毕"
        echo "  MultiQC报告: ${DIR_MULTIQC}/${SAMPLE_ID}_multiqc_report.html"
        echo "  综合报告:    ${final_report}"
        echo "=========================================================================="

    } > "${final_report}"

    log_info "综合报告已生成: ${final_report}"
    echo ""
    cat "${final_report}"
}

main "$@"
