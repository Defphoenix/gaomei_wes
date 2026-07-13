#!/bin/bash
#===============================================================================
# 08_cnv.sh - CNV 拷贝数变异分析
#
# 优先使用 CNVkit；如果环境中没有 CNVkit，则使用 mosdepth 对目标 BED
# 做轻量 depth-ratio CNV，输出与下游一致的 cnr/cns/call.cns/summary 文件。
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/../config.sh}"
source "${CONFIG_FILE}"
source "${SCRIPT_DIR}/utils.sh"

choose_cnv_method() {
    local requested="${CNV_METHOD:-auto}"
    if [ "${requested}" = "cnvkit" ]; then
        echo "cnvkit"
    elif [ "${requested}" = "depth" ]; then
        echo "depth"
    elif command -v "${TOOL_CNVKIT:-cnvkit.py}" >/dev/null 2>&1; then
        echo "cnvkit"
    else
        echo "depth"
    fi
}

run_depth_cnv() {
    local input_bam="$1"
    local target_bed="$2"

    check_tool "mosdepth" "${TOOL_MOSDEPTH:-mosdepth}" || return 1

    local prefix="${DIR_CNV}/${SAMPLE_ID}.cnv_depth"
    local regions_gz="${prefix}.regions.bed.gz"
    local cnr_file="${DIR_CNV}/${SAMPLE_ID}.cnr"
    local cns_file="${DIR_CNV}/${SAMPLE_ID}.cns"
    local call_file="${DIR_CNV}/${SAMPLE_ID}.call.cns"
    local vcf_file="${DIR_CNV}/${SAMPLE_ID}.cnv.vcf"
    local seg_file="${DIR_CNV}/${SAMPLE_ID}.seg"

    run_cmd "mosdepth 目标区域深度" \
        "${TOOL_MOSDEPTH:-mosdepth}" \
        -t "${MOSDEPTH_THREADS:-2}" \
        -n \
        --by "${target_bed}" \
        "${prefix}" \
        "${input_bam}"

    check_file "mosdepth区域深度" "${regions_gz}" || return 1

    gzip -cd "${regions_gz}" | awk -v sample="${SAMPLE_ID}" '
        BEGIN{OFS="\t"; print "chromosome","start","end","gene","depth","log2","weight"}
        {rows[NR]=$0; depth[NR]=$4; sum+=$4; n++}
        END{
            mean=(n>0 ? sum/n : 0);
            if (mean <= 0) mean=1;
            for (i=1; i<=n; i++) {
                split(rows[i], a, "\t");
                ratio=depth[i]/mean;
                log2=(ratio>0 ? log(ratio)/log(2) : -10);
                gene="bin_"i;
                print a[1],a[2],a[3],gene,depth[i],log2,1;
            }
        }' > "${cnr_file}"

    awk -F'\t' -v OFS='\t' '
        NR==1{print "chromosome","start","end","gene","log2","probes"; next}
        {print $1,$2,$3,$4,$6,1}
    ' "${cnr_file}" > "${cns_file}"

    awk -F'\t' -v OFS='\t' -v low="${CNV_NEUTRAL_LOW:-0.75}" -v high="${CNV_NEUTRAL_HIGH:-1.25}" '
        NR==1{print "chromosome","start","end","gene","log2","probes","depth_ratio","cn","call"; next}
        {
            ratio=2^$5;
            cn=(ratio < low ? 1 : (ratio > high ? 3 : 2));
            call=(cn < 2 ? "loss" : (cn > 2 ? "gain" : "neutral"));
            print $1,$2,$3,$4,$5,$6,ratio,cn,call;
        }
    ' "${cns_file}" > "${call_file}"

    awk -F'\t' -v OFS='\t' -v sample="${SAMPLE_ID}" '
        BEGIN{print "ID","chrom","loc.start","loc.end","num.mark","seg.mean"}
        NR>1{print sample,$1,$2+1,$3,$6,$5}
    ' "${cns_file}" > "${seg_file}"

    {
        echo "##fileformat=VCFv4.2"
        echo "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"CNV type\">"
        echo "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">"
        echo "##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Estimated copy number\">"
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        awk -F'\t' 'NR>1 && $9!="neutral"{printf "%s\t%d\t%s_%d_%d\tN\t<%s>\t.\tPASS\tSVTYPE=%s;END=%d;CN=%d\n",$1,$2+1,$9,$2+1,$3,toupper($9),toupper($9),$3,$8}' "${call_file}"
    } > "${vcf_file}"

    write_cnv_summary "${call_file}" "${target_bed}" "depth"
}

run_cnvkit_cnv() {
    local input_bam="$1"
    local target_bed="$2"

    check_tool "CNVkit" "${TOOL_CNVKIT:-cnvkit.py}" || return 1

    local cnr_file="${DIR_CNV}/${SAMPLE_ID}.cnr"
    local cns_file="${DIR_CNV}/${SAMPLE_ID}.cns"
    local call_file="${DIR_CNV}/${SAMPLE_ID}.call.cns"
    local target_cnn="${DIR_CNV}/${SAMPLE_ID}.targetcoverage.cnn"
    local antitarget_bed="${CNVKIT_ANTITARGET_BED:-}"
    local antitarget_cnn="${DIR_CNV}/${SAMPLE_ID}.antitargetcoverage.cnn"
    local ref_param=""

    if [ -n "${CNVKIT_REFERENCE:-}" ] && [ -f "${CNVKIT_REFERENCE}" ]; then
        ref_param="--reference ${CNVKIT_REFERENCE}"
    fi

    run_cmd "CNVkit target coverage" \
        "${TOOL_CNVKIT:-cnvkit.py}" coverage \
        "${input_bam}" \
        "${target_bed}" \
        -o "${target_cnn}"

    if [ -z "${antitarget_bed}" ]; then
        antitarget_bed="${DIR_CNV}/${SAMPLE_ID}.antitarget.bed"
        run_cmd "CNVkit antitarget" \
            "${TOOL_CNVKIT:-cnvkit.py}" antitarget \
            "${target_bed}" \
            -o "${antitarget_bed}" || antitarget_bed=""
    fi

    if [ -n "${antitarget_bed}" ] && [ -f "${antitarget_bed}" ]; then
        run_cmd "CNVkit antitarget coverage" \
            "${TOOL_CNVKIT:-cnvkit.py}" coverage \
            "${input_bam}" \
            "${antitarget_bed}" \
            -o "${antitarget_cnn}" || antitarget_cnn=""
    fi

    if [ -n "${antitarget_cnn}" ] && [ -f "${antitarget_cnn}" ]; then
        run_cmd "CNVkit fix" \
            "${TOOL_CNVKIT:-cnvkit.py}" fix \
            "${target_cnn}" \
            "${antitarget_cnn}" \
            ${ref_param} \
            -o "${cnr_file}"
    else
        run_cmd "CNVkit fix target-only" \
            "${TOOL_CNVKIT:-cnvkit.py}" fix \
            "${target_cnn}" \
            ${ref_param} \
            -o "${cnr_file}"
    fi

    run_cmd "CNVkit segment" \
        "${TOOL_CNVKIT:-cnvkit.py}" segment \
        ${CNVKIT_SEGMENT_PARAMS:-} \
        "${cnr_file}" \
        -o "${cns_file}"

    run_cmd "CNVkit call" \
        "${TOOL_CNVKIT:-cnvkit.py}" call \
        ${CNVKIT_CALL_PARAMS:-} \
        "${cns_file}" \
        -o "${call_file}"

    run_cmd "CNVkit export SEG" \
        "${TOOL_CNVKIT:-cnvkit.py}" export seg \
        "${cns_file}" \
        -o "${DIR_CNV}/${SAMPLE_ID}.seg" || log_warn "CNVkit SEG导出失败"

    run_cmd "CNVkit export VCF" \
        "${TOOL_CNVKIT:-cnvkit.py}" export vcf \
        "${call_file}" \
        -o "${DIR_CNV}/${SAMPLE_ID}.cnv.vcf" || log_warn "CNVkit VCF导出失败"

    write_cnv_summary "${call_file}" "${target_bed}" "cnvkit"
}

write_cnv_summary() {
    local call_file="$1"
    local target_bed="$2"
    local method="$3"
    local summary="${DIR_CNV}/${SAMPLE_ID}_cnv_summary.txt"

    {
        echo "CNV分析摘要 - ${SAMPLE_ID}"
        echo "================================"
        echo "方法: ${method}"
        echo "目标区域: ${target_bed}"
        echo "调用结果: ${call_file}"
        echo ""
        if [ -f "${call_file}" ]; then
            echo "片段统计:"
            awk -F'\t' 'NR>1{total++; call=$NF; counts[call]++} END{print "  总片段: "total+0; for (k in counts) print "  "k": "counts[k]}' "${call_file}"
        else
            echo "调用结果未生成"
        fi
    } > "${summary}"
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
    log_info "CNV分析方法: ${method}"

    case "${method}" in
        cnvkit) run_cnvkit_cnv "${input_bam}" "${target_bed}" ;;
        depth) run_depth_cnv "${input_bam}" "${target_bed}" ;;
        *) log_error "未知 CNV_METHOD: ${method}"; return 1 ;;
    esac

    log_info "CNV分析完成!"
    log_info "CNR: ${DIR_CNV}/${SAMPLE_ID}.cnr"
    log_info "CNS: ${DIR_CNV}/${SAMPLE_ID}.cns"
    log_info "CALL: ${DIR_CNV}/${SAMPLE_ID}.call.cns"
}

main "$@"
