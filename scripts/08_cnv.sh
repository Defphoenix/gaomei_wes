#!/bin/bash
#===============================================================================
# 08_cnv.sh - CNV copy-number analysis
#
# CNVkit is used only when a pooled reference or matched-normal BAM is
# available. Mosdepth fallback is explicitly reported as depth QC and never
# emits CNVkit-like calls.
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/../config.sh}"
source "${CONFIG_FILE}"
source "${SCRIPT_DIR}/utils.sh"

cnvkit_available() {
    [ -x "${TOOL_CNVKIT:-}" ] || command -v "${TOOL_CNVKIT:-cnvkit.py}" >/dev/null 2>&1
}

cnvkit_baseline_available() {
    if [ -n "${CNVKIT_REFERENCE:-}" ] && [ -s "${CNVKIT_REFERENCE}" ]; then
        return 0
    fi
    [ -n "${NORMAL_BAM:-}" ] && bam_is_complete "${NORMAL_BAM}"
}

choose_cnv_method() {
    local requested="${CNV_METHOD:-auto}"
    case "${requested}" in
        cnvkit) echo "cnvkit" ;;
        depth|depth_qc) echo "depth_qc" ;;
        auto)
            if cnvkit_available && cnvkit_baseline_available; then
                echo "cnvkit"
            else
                echo "depth_qc"
            fi
            ;;
        *) echo "invalid" ;;
    esac
}

write_depth_qc_summary() {
    local depth_qc="$1"
    local target_bed="$2"
    local summary="${DIR_CNV}/${SAMPLE_ID}_cnv_summary.txt"
    {
        echo "CNV/深度质控摘要 - ${SAMPLE_ID}"
        echo "================================"
        echo "方法: mosdepth depth_qc"
        echo "目标区域: ${target_bed}"
        echo "深度质控表: ${depth_qc}"
        echo "状态: EXPLORATORY_DEPTH_QC_ONLY"
        echo ""
        echo "警告: 本结果仅表示相对样本均值的覆盖深度异常。"
        echo "未使用matched/pooled normal、GC/可比对性校正或分段模型，不能作为正式CNV调用。"
        echo "不会生成伪装成CNVkit结果的CNR/CNS/VCF文件。"
    } > "${summary}"
}

run_depth_qc() {
    local input_bam="$1"
    local target_bed="$2"

    check_tool "mosdepth" "${TOOL_MOSDEPTH:-mosdepth}" || return 1
    log_warn "CNVkit或normal reference不可用；仅生成mosdepth探索性深度QC，不进行CNV calling"

    local prefix="${DIR_CNV}/${SAMPLE_ID}.cnv_depth"
    local regions_gz="${prefix}.regions.bed.gz"
    local depth_qc="${DIR_CNV}/${SAMPLE_ID}.depth_qc.tsv"

    run_cmd "mosdepth CNV目标区域深度QC" \
        "${TOOL_MOSDEPTH:-mosdepth}" \
        -t "${MOSDEPTH_THREADS:-2}" \
        -n \
        --by "${target_bed}" \
        "${prefix}" \
        "${input_bam}"

    check_file "mosdepth区域深度" "${regions_gz}" || return 1

    gzip -cd "${regions_gz}" | awk '
        BEGIN{OFS="\t"; print "chromosome","start","end","region","depth","ratio_to_weighted_mean","log2_to_weighted_mean","qc_flag"}
        {
            rows[NR]=$0;
            len=$3-$2;
            if (len > 0) {weighted_sum += $4*len; total_length += len}
            n++
        }
        END{
            mean=(total_length > 0 ? weighted_sum/total_length : 0);
            if (mean <= 0) mean=1;
            for (i=1; i<=n; i++) {
                split(rows[i], a, "\t");
                ratio=a[4]/mean;
                log2=(ratio>0 ? log(ratio)/log(2) : -10);
                region=(a[5] != "" ? a[5] : "bin_"i);
                flag=(ratio < 0.75 ? "LOW_DEPTH" : (ratio > 1.25 ? "HIGH_DEPTH" : "EXPECTED_DEPTH"));
                print a[1],a[2],a[3],region,a[4],ratio,log2,flag;
            }
        }' > "${depth_qc}"

    write_depth_qc_summary "${depth_qc}" "${target_bed}"
}

write_cnvkit_summary() {
    local call_file="$1"
    local target_bed="$2"
    local reference_cnn="$3"
    local baseline_type="$4"
    local summary="${DIR_CNV}/${SAMPLE_ID}_cnv_summary.txt"

    {
        echo "CNV分析摘要 - ${SAMPLE_ID}"
        echo "================================"
        echo "方法: CNVkit 0.9.x"
        echo "基线类型: ${baseline_type}"
        echo "目标区域: ${target_bed}"
        echo "CNVkit reference: ${reference_cnn}"
        echo "调用结果: ${call_file}"
        echo "状态: CNVKIT_CALL_COMPLETE"
        echo ""
        echo "片段统计:"
        awk -F '\t' '
            NR==1 {
                for (i=1; i<=NF; i++) if ($i=="cn") cn_col=i;
                next
            }
            cn_col>0 && $cn_col!="" {
                total++;
                if ($cn_col < 2) loss++;
                else if ($cn_col > 2) gain++;
                else neutral++;
            }
            END {
                print "  总片段: " total+0;
                print "  loss: " loss+0;
                print "  neutral: " neutral+0;
                print "  gain: " gain+0;
            }' "${call_file}"
        echo ""
        echo "说明: matched-normal结果优于单样本深度比例，但正式报告仍需同平台normal队列、纯度/倍性评估和方法学验证。"
    } > "${summary}"
}

run_cnvkit_cnv() {
    local tumor_bam="$1"
    local target_bed="$2"
    local cnvkit="${TOOL_CNVKIT:-${SCRIPT_DIR}/run_cnvkit_env.sh}"

    if ! "${cnvkit}" version >/dev/null 2>&1; then
        log_error "CNVkit不可用: ${cnvkit}"
        return 1
    fi
    log_info "CNVkit检查通过: $("${cnvkit}" version 2>&1 | head -1)"

    local normalized_targets="${DIR_CNV}/${SAMPLE_ID}.targets.bed"
    local antitarget_bed="${CNVKIT_ANTITARGET_BED:-${DIR_CNV}/${SAMPLE_ID}.antitarget.bed}"
    local tumor_target_cnn="${DIR_CNV}/${SAMPLE_ID}.tumor.targetcoverage.cnn"
    local tumor_antitarget_cnn="${DIR_CNV}/${SAMPLE_ID}.tumor.antitargetcoverage.cnn"
    local normal_target_cnn="${DIR_CNV}/${SAMPLE_ID}.normal.targetcoverage.cnn"
    local normal_antitarget_cnn="${DIR_CNV}/${SAMPLE_ID}.normal.antitargetcoverage.cnn"
    local reference_cnn="${CNVKIT_REFERENCE:-${DIR_CNV}/${SAMPLE_ID}.matched_normal_reference.cnn}"
    local cnr_file="${DIR_CNV}/${SAMPLE_ID}.cnr"
    local cns_file="${DIR_CNV}/${SAMPLE_ID}.cns"
    local call_file="${DIR_CNV}/${SAMPLE_ID}.call.cns"
    local baseline_type="pooled_or_prebuilt_reference"

    run_cmd "CNVkit规范化target BED" \
        "${cnvkit}" target "${target_bed}" --split -o "${normalized_targets}"

    if [ -z "${CNVKIT_ANTITARGET_BED:-}" ] || [ ! -s "${CNVKIT_ANTITARGET_BED:-}" ]; then
        run_cmd "CNVkit生成antitarget BED" \
            "${cnvkit}" antitarget "${normalized_targets}" -o "${antitarget_bed}"
    fi
    check_file "CNVkit target BED" "${normalized_targets}" || return 1
    check_file "CNVkit antitarget BED" "${antitarget_bed}" || return 1

    run_cmd "CNVkit tumor target coverage" \
        "${cnvkit}" coverage "${tumor_bam}" "${normalized_targets}" \
        -p "${CNVKIT_PROCESSES:-2}" -o "${tumor_target_cnn}"
    run_cmd "CNVkit tumor antitarget coverage" \
        "${cnvkit}" coverage "${tumor_bam}" "${antitarget_bed}" \
        -p "${CNVKIT_PROCESSES:-2}" -o "${tumor_antitarget_cnn}"

    if [ -z "${CNVKIT_REFERENCE:-}" ] || [ ! -s "${CNVKIT_REFERENCE:-}" ]; then
        if [ -z "${NORMAL_BAM:-}" ] || ! bam_is_complete "${NORMAL_BAM}"; then
            log_error "CNVkit没有预构建reference，且matched normal BAM不可用"
            return 1
        fi
        baseline_type="single_matched_normal"
        check_file "CNV matched normal BAM" "${NORMAL_BAM}" || return 1
        run_cmd "CNVkit normal target coverage" \
            "${cnvkit}" coverage "${NORMAL_BAM}" "${normalized_targets}" \
            -p "${CNVKIT_PROCESSES:-2}" -o "${normal_target_cnn}"
        run_cmd "CNVkit normal antitarget coverage" \
            "${cnvkit}" coverage "${NORMAL_BAM}" "${antitarget_bed}" \
            -p "${CNVKIT_PROCESSES:-2}" -o "${normal_antitarget_cnn}"
        run_cmd "CNVkit构建matched-normal reference" \
            "${cnvkit}" reference "${normal_target_cnn}" "${normal_antitarget_cnn}" \
            --fasta "${REFERENCE_GENOME}" -o "${reference_cnn}"
    fi
    check_file "CNVkit reference" "${reference_cnn}" || return 1

    run_cmd "CNVkit校正肿瘤覆盖度" \
        "${cnvkit}" fix "${tumor_target_cnn}" "${tumor_antitarget_cnn}" "${reference_cnn}" \
        -o "${cnr_file}"
    run_cmd "CNVkit分段" \
        "${cnvkit}" segment ${CNVKIT_SEGMENT_PARAMS:-} "${cnr_file}" -o "${cns_file}"
    run_cmd "CNVkit拷贝数调用" \
        "${cnvkit}" call ${CNVKIT_CALL_PARAMS:-} "${cns_file}" -o "${call_file}"

    run_cmd "CNVkit导出SEG" \
        "${cnvkit}" export seg "${cns_file}" -o "${DIR_CNV}/${SAMPLE_ID}.seg" || \
        log_warn "CNVkit SEG导出失败 (非致命)"
    run_cmd "CNVkit导出VCF" \
        "${cnvkit}" export vcf "${call_file}" -o "${DIR_CNV}/${SAMPLE_ID}.cnv.vcf" || \
        log_warn "CNVkit VCF导出失败 (非致命)"

    write_cnvkit_summary "${call_file}" "${target_bed}" "${reference_cnn}" "${baseline_type}"
}

main() {
    log_step "步骤8: CNV 拷贝数变异分析"

    if [ "${SKIP_CNV:-true}" = true ] || [ "${RUN_CNV:-false}" = false ]; then
        log_warn "跳过CNV分析 (SKIP_CNV=${SKIP_CNV:-true}, RUN_CNV=${RUN_CNV:-false})"
        return 0
    fi

    local input_bam
    input_bam=$(get_final_bam) || return 1
    check_file "输入BAM" "${input_bam}" || return 1

    local target_bed="${CNVKIT_TARGET_BED:-${INTERVAL_FILE:-}}"
    if [ -z "${target_bed}" ] || [ ! -f "${target_bed}" ]; then
        log_error "CNV目标区域BED不存在: ${target_bed}"
        return 1
    fi

    mkdir -p "${DIR_CNV}"
    local method
    method=$(choose_cnv_method)
    if [ "${method}" = "invalid" ]; then
        log_error "未知 CNV_METHOD: ${CNV_METHOD:-}"
        return 1
    fi
    if [ "${CNV_METHOD:-auto}" = "depth" ]; then
        log_warn "CNV_METHOD=depth为兼容旧配置，将按depth_qc运行且不输出正式CNV calls"
    fi
    if [ "${method}" = "depth_qc" ] && [ "${CNV_REQUIRE_REFERENCE:-false}" = true ]; then
        log_error "CNV_REQUIRE_REFERENCE=true，但CNVkit或matched/pooled normal reference不可用"
        return 1
    fi

    log_info "CNV分析方法: ${method}"
    case "${method}" in
        cnvkit) run_cnvkit_cnv "${input_bam}" "${target_bed}" ;;
        depth_qc) run_depth_qc "${input_bam}" "${target_bed}" ;;
    esac

    log_info "CNV模块完成! 摘要: ${DIR_CNV}/${SAMPLE_ID}_cnv_summary.txt"
}

main "$@"
