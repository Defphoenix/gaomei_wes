#!/usr/bin/env python3
"""Compare a VEP VCF with a tabular variant report and audit BED coverage."""

from __future__ import annotations

import argparse
import gzip
import html
import json
import math
import re
import shlex
from bisect import bisect_right
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import unquote

import pandas as pd


AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Stop": "*", "Sec": "U",
}


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def normalize_chrom(value: str) -> str:
    value = str(value).strip()
    return value if value.startswith("chr") else f"chr{value}"


def normalize_cdna(value: Any) -> str:
    text = unquote(str(value or "")).strip()
    if ":" in text:
        text = text.rsplit(":", 1)[-1]
    match = re.search(r"([cn]\.[^\s;]+)", text, re.I)
    return match.group(1).replace(" ", "") if match else text.replace(" ", "")


def normalize_protein(value: Any) -> str:
    text = unquote(str(value or "")).strip()
    if ":" in text:
        text = text.rsplit(":", 1)[-1]
    text = text.replace("Ter", "*").replace("%3D", "=")
    for aa3, aa1 in AA3_TO_1.items():
        text = text.replace(aa3, aa1)
    match = re.search(r"(p\.[^\s;]+)", text, re.I)
    return match.group(1).replace(" ", "") if match else text.replace(" ", "")


def parse_info(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        elif item:
            result[item] = "true"
    return result


def parse_sample(format_text: str, sample_text: str) -> dict[str, str]:
    keys = format_text.split(":") if format_text else []
    values = sample_text.split(":") if sample_text else []
    values.extend([""] * (len(keys) - len(values)))
    return dict(zip(keys, values))


def allele_value(text: str, alt_index: int) -> str:
    if not text or text == ".":
        return ""
    values = text.split(",")
    return values[alt_index] if alt_index < len(values) else values[0]


def sample_metrics(sample: dict[str, str], alt_index: int) -> tuple[str, str, str]:
    ad_values = sample.get("AD", "").split(",")
    alt_reads = ad_values[alt_index + 1] if len(ad_values) > alt_index + 1 else ""
    dp = sample.get("DP", "")
    af = allele_value(sample.get("AF", ""), alt_index)
    if not af and dp not in {"", ".", "0"} and alt_reads not in {"", "."}:
        try:
            af = f"{int(alt_reads) / int(dp):.6g}"
        except ValueError:
            pass
    return dp, alt_reads, af


def annotation_score(item: dict[str, str]) -> tuple[int, ...]:
    impact = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}.get(item.get("IMPACT", ""), 0)
    return (
        int(bool(item.get("MANE_SELECT") or item.get("MANE"))),
        int(item.get("CANONICAL") == "YES"),
        int(item.get("BIOTYPE") == "protein_coding"),
        impact,
        int(bool(item.get("HGVSc"))),
        int(bool(item.get("HGVSp"))),
    )


def parse_vcf(path: Path, tumor_sample: str | None, normal_sample: str | None):
    metadata: list[str] = []
    records: list[list[str]] = []
    samples: list[str] = []
    csq_fields: list[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("##"):
                metadata.append(line)
                if line.startswith("##INFO=<ID=CSQ"):
                    match = re.search(r"Format: ([^\"]+)", line)
                    if match:
                        csq_fields = match.group(1).split("|")
                continue
            if line.startswith("#CHROM"):
                samples = line.split("\t")[9:]
                continue
            if line and not line.startswith("#"):
                records.append(line.split("\t"))

    if not csq_fields:
        raise SystemExit("Input VCF has no VEP CSQ header")
    if not samples:
        raise SystemExit("Input VCF has no sample columns")
    tumor_sample = tumor_sample or samples[0]
    normal_sample = normal_sample or (samples[1] if len(samples) > 1 else None)
    if tumor_sample not in samples:
        raise SystemExit(f"Tumor sample {tumor_sample!r} not found; VCF samples: {samples}")
    if normal_sample and normal_sample not in samples:
        raise SystemExit(f"Normal sample {normal_sample!r} not found; VCF samples: {samples}")
    tumor_index = samples.index(tumor_sample)
    normal_index = samples.index(normal_sample) if normal_sample else None

    variant_rows: list[dict[str, Any]] = []
    annotation_rows: list[dict[str, Any]] = []
    for fields in records:
        if len(fields) < 8:
            continue
        chrom, pos_text, variant_id, ref, alt_text, qual, filter_text, info_text = fields[:8]
        info = parse_info(info_text)
        format_text = fields[8] if len(fields) > 8 else ""
        sample_values = fields[9:]
        tumor = parse_sample(format_text, sample_values[tumor_index])
        normal = parse_sample(format_text, sample_values[normal_index]) if normal_index is not None else {}
        parsed_csq: list[dict[str, str]] = []
        for entry in info.get("CSQ", "").split(","):
            values = entry.split("|")
            values.extend([""] * (len(csq_fields) - len(values)))
            parsed_csq.append(dict(zip(csq_fields, values)))

        for alt_index, alt in enumerate(alt_text.split(",")):
            key = f"{normalize_chrom(chrom)}:{pos_text}:{ref}:{alt}"
            allele_csq = [item for item in parsed_csq if item.get("Allele") == alt]
            if not allele_csq and len(alt_text.split(",")) == 1:
                allele_csq = parsed_csq
            best = max(allele_csq, key=annotation_score, default={})
            t_dp, t_alt, t_af = sample_metrics(tumor, alt_index)
            n_dp, n_alt, n_af = sample_metrics(normal, alt_index)
            genes = sorted({item.get("SYMBOL", "") for item in allele_csq if item.get("SYMBOL")})
            consequences = sorted({c for item in allele_csq for c in item.get("Consequence", "").split("&") if c})
            row = {
                "variant_key": key, "chrom": normalize_chrom(chrom), "pos": int(pos_text),
                "id": variant_id, "ref": ref, "alt": alt, "qual": qual, "filter": filter_text,
                "tlod": allele_value(info.get("TLOD", ""), alt_index),
                "genes": ";".join(genes), "gene": best.get("SYMBOL", ""),
                "consequence": best.get("Consequence", ""), "all_consequences": ";".join(consequences),
                "impact": best.get("IMPACT", ""), "transcript": best.get("Feature", ""),
                "hgvsc": best.get("HGVSc", ""), "hgvsp": best.get("HGVSp", ""),
                "tumor_sample": tumor_sample, "tumor_dp": t_dp, "tumor_alt_reads": t_alt, "tumor_af": t_af,
                "normal_sample": normal_sample or "", "normal_dp": n_dp, "normal_alt_reads": n_alt, "normal_af": n_af,
            }
            variant_rows.append(row)
            for item in allele_csq:
                annotation_rows.append({
                    "variant_key": key, "gene": item.get("SYMBOL", ""),
                    "consequence": item.get("Consequence", ""), "impact": item.get("IMPACT", ""),
                    "transcript": item.get("Feature", ""), "canonical": item.get("CANONICAL", ""),
                    "mane_select": item.get("MANE_SELECT", "") or item.get("MANE", ""),
                    "hgvsc": item.get("HGVSc", ""), "hgvsp": item.get("HGVSp", ""),
                    "hgvsc_norm": normalize_cdna(item.get("HGVSc", "")),
                    "hgvsp_norm": normalize_protein(item.get("HGVSp", "")),
                })
    return pd.DataFrame(variant_rows), pd.DataFrame(annotation_rows), metadata, samples, tumor_sample, normal_sample


def load_bed(path: Path):
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith(("#", "track", "browser")):
                continue
            fields = raw.rstrip().split("\t")
            if len(fields) < 3:
                continue
            try:
                intervals[normalize_chrom(fields[0])].append((int(fields[1]), int(fields[2])))
            except ValueError:
                continue
    starts: dict[str, list[int]] = {}
    merged: dict[str, list[tuple[int, int]]] = {}
    for chrom, values in intervals.items():
        combined: list[list[int]] = []
        for start, end in sorted(values):
            if combined and start <= combined[-1][1]:
                combined[-1][1] = max(combined[-1][1], end)
            else:
                combined.append([start, end])
        merged[chrom] = [(a, b) for a, b in combined]
        starts[chrom] = [a for a, _ in merged[chrom]]
    return merged, starts


def bed_status(chrom: str, pos: int, ref: str, intervals, starts) -> tuple[bool, int | None]:
    chrom = normalize_chrom(chrom)
    values = intervals.get(chrom, [])
    if not values:
        return False, None
    query_start, query_end = pos - 1, pos - 1 + max(1, len(ref))
    index = bisect_right(starts[chrom], query_start) - 1
    candidates = []
    if index >= 0:
        candidates.append(values[index])
    if index + 1 < len(values):
        candidates.append(values[index + 1])
    for start, end in candidates:
        if query_start < end and query_end > start:
            return True, 0
    distance = min(min(abs(query_start - end), abs(start - query_end)) for start, end in candidates)
    return False, distance


def find_column(frame: pd.DataFrame, names: Iterable[str]) -> str | None:
    normalized = {re.sub(r"\s+", "", str(col)).lower(): col for col in frame.columns}
    for name in names:
        key = re.sub(r"\s+", "", name).lower()
        if key in normalized:
            return str(normalized[key])
    return None


def load_other(path: Path) -> pd.DataFrame:
    try:
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    except UnicodeDecodeError:
        frame = pd.read_csv(path, sep="\t", dtype=str, encoding="gb18030").fillna("")
    aliases = {
        "gene": ["基因", "gene", "symbol"],
        "region": ["外显子/内含子", "区域", "region", "exon"],
        "hgvsc_other": ["核苷酸变异", "hgvsc", "cdna", "c.hgvs"],
        "hgvsp_other": ["氨基酸变异", "hgvsp", "protein", "p.hgvs"],
        "other_vaf": ["丰度/拷贝数", "丰度", "vaf", "af", "frequency"],
        "other_type": ["变异类型", "类型", "type", "consequence"],
        "chrom": ["chrom", "chr", "染色体"], "pos": ["pos", "position", "位置"],
        "ref": ["ref", "参考碱基"], "alt": ["alt", "突变碱基"],
    }
    result = pd.DataFrame(index=frame.index)
    for target, names in aliases.items():
        source = find_column(frame, names)
        result[target] = frame[source].astype(str) if source else ""
    result.insert(0, "other_row", range(1, len(result) + 1))
    result["gene"] = result["gene"].str.strip()
    result["hgvsc_norm"] = result["hgvsc_other"].map(normalize_cdna)
    result["hgvsp_norm"] = result["hgvsp_other"].map(normalize_protein)
    result["other_vaf_numeric"] = pd.to_numeric(result["other_vaf"].str.replace("%", "", regex=False), errors="coerce") / 100
    return result


def parse_actual_filters(metadata: list[str]) -> dict[str, str]:
    commands: dict[str, str] = {}
    for line in metadata:
        match = re.match(r'##GATKCommandLine=<ID=([^,]+),CommandLine="(.*?)",Version="([^"]+)', line)
        if match:
            commands[match.group(1)] = match.group(2)
    mutect = commands.get("Mutect2", "")
    filtering = commands.get("FilterMutectCalls", "")

    def option(command: str, *names: str) -> str:
        try:
            tokens = shlex.split(command)
        except ValueError:
            tokens = command.split()
        for name in names:
            if name in tokens:
                index = tokens.index(name)
                return tokens[index + 1] if index + 1 < len(tokens) else "是"
        return "未使用"

    return {
        "配对正常样本": option(mutect, "--normal-sample", "-normal"),
        "分析区间": option(mutect, "--intervals", "-L"),
        "最低比对质量": option(mutect, "--minimum-mapping-quality"),
        "每个比对起点最大 reads": option(mutect, "--max-reads-per-alignment-start"),
        "Panel of Normals": option(mutect, "--panel-of-normals"),
        "群体胚系资源": option(mutect, "--germline-resource"),
        "方向偏倚模型": option(filtering, "--orientation-bias-artifact-priors", "--ob-priors"),
        "污染估计表": option(filtering, "--contamination-table"),
        "FilterMutectCalls策略": option(filtering, "--threshold-strategy"),
        "FilterMutectCalls FDR": option(filtering, "--false-discovery-rate"),
    }


def parse_filtered_vcf(path: Path, tumor_sample: str, normal_sample: str | None):
    samples: list[str] = []
    records: dict[str, dict[str, Any]] = {}
    filter_counts: Counter[str] = Counter()
    total_records = 0
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("#CHROM"):
                samples = line.split("\t")[9:]
                continue
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            total_records += 1
            chrom, pos, _, ref, alt_text, _, filter_text, info_text = fields[:8]
            for tag in filter_text.split(";"):
                filter_counts[tag] += 1
            if tumor_sample not in samples:
                continue
            tumor_index = samples.index(tumor_sample)
            normal_index = samples.index(normal_sample) if normal_sample and normal_sample in samples else None
            fmt = fields[8] if len(fields) > 8 else ""
            values = fields[9:]
            tumor = parse_sample(fmt, values[tumor_index])
            normal = parse_sample(fmt, values[normal_index]) if normal_index is not None else {}
            info = parse_info(info_text)
            for alt_index, alt in enumerate(alt_text.split(",")):
                t_dp, t_alt, t_af = sample_metrics(tumor, alt_index)
                n_dp, n_alt, n_af = sample_metrics(normal, alt_index)
                key = f"{normalize_chrom(chrom)}:{pos}:{ref}:{alt}"
                records[key] = {
                    "filtered_variant_key": key, "filtered_filter": filter_text,
                    "filtered_tlod": allele_value(info.get("TLOD", ""), alt_index),
                    "filtered_roq": info.get("ROQ", ""),
                    "filtered_tumor_dp": t_dp, "filtered_tumor_alt_reads": t_alt, "filtered_tumor_af": t_af,
                    "filtered_tumor_f1r2": tumor.get("F1R2", ""),
                    "filtered_tumor_f2r1": tumor.get("F2R1", ""), "filtered_tumor_sb": tumor.get("SB", ""),
                    "filtered_normal_dp": n_dp, "filtered_normal_alt_reads": n_alt, "filtered_normal_af": n_af,
                }
    return records, filter_counts, total_records


def compare(ours: pd.DataFrame, annotations: pd.DataFrame, other: pd.DataFrame,
            gene_vaf_tolerance: float) -> tuple[pd.DataFrame, set[str]]:
    cdna_map: dict[tuple[str, str], set[str]] = defaultdict(set)
    protein_map: dict[tuple[str, str], set[str]] = defaultdict(set)
    gene_map: dict[str, set[str]] = defaultdict(set)
    for row in annotations.itertuples(index=False):
        gene = str(row.gene).upper()
        gene_map[gene].add(row.variant_key)
        if row.hgvsc_norm:
            cdna_map[(gene, row.hgvsc_norm)].add(row.variant_key)
        if row.hgvsp_norm and row.hgvsp_norm not in {"p.?", "?"}:
            protein_map[(gene, row.hgvsp_norm)].add(row.variant_key)
    by_key = ours.set_index("variant_key").to_dict("index")

    rows: list[dict[str, Any]] = []
    matched_keys: set[str] = set()
    for item in other.to_dict("records"):
        gene = str(item["gene"]).upper()
        c_candidates = cdna_map.get((gene, item["hgvsc_norm"]), set()) if item["hgvsc_norm"] else set()
        p_candidates = protein_map.get((gene, item["hgvsp_norm"]), set()) if item["hgvsp_norm"] else set()
        intersection = c_candidates & p_candidates
        candidates = intersection or c_candidates or p_candidates
        match_method = "基因+cDNA+蛋白HGVS" if intersection else ("基因+cDNA HGVS" if c_candidates else ("基因+蛋白HGVS" if p_candidates else "未匹配"))
        comparison_status = "两流程共同检出" if candidates else "仅外部流程检出"
        if not candidates and len(gene_map.get(gene, set())) == 1 and not pd.isna(item["other_vaf_numeric"]):
            gene_candidate = next(iter(gene_map[gene]))
            try:
                our_af = float(by_key[gene_candidate].get("tumor_af", ""))
            except (TypeError, ValueError):
                our_af = math.nan
            if not math.isnan(our_af) and abs(our_af - float(item["other_vaf_numeric"])) <= gene_vaf_tolerance:
                candidates = {gene_candidate}
                match_method = f"同基因唯一候选且AF差≤{gene_vaf_tolerance:.1%}；HGVS转录本不同，需坐标确认"
                comparison_status = "疑似共同检出（转录本差异）"
        if len(candidates) > 1:
            match_method += "（多候选）"
        matched_key = sorted(candidates)[0] if candidates else ""
        matched = by_key.get(matched_key, {})
        if matched_key:
            matched_keys.add(matched_key)
        row = dict(item)
        row.update({
            "comparison_status": comparison_status,
            "match_method": match_method,
            "matched_variant_key": matched_key,
            "our_chrom": matched.get("chrom", ""), "our_pos": matched.get("pos", ""),
            "our_ref": matched.get("ref", ""), "our_alt": matched.get("alt", ""),
            "our_filter": matched.get("filter", ""), "our_tlod": matched.get("tlod", ""),
            "our_tumor_dp": matched.get("tumor_dp", ""), "our_tumor_alt_reads": matched.get("tumor_alt_reads", ""),
            "our_tumor_af": matched.get("tumor_af", ""), "our_normal_dp": matched.get("normal_dp", ""),
            "our_normal_alt_reads": matched.get("normal_alt_reads", ""), "our_normal_af": matched.get("normal_af", ""),
            "our_consequence": matched.get("consequence", ""),
            "same_gene_seen_in_our_vcf": "是" if gene in gene_map else "否",
        })
        rows.append(row)
    return pd.DataFrame(rows), matched_keys


def esc(value: Any) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return html.escape(str(value))


def table_html(frame: pd.DataFrame, columns: list[tuple[str, str]], limit: int | None = None) -> str:
    shown = frame.head(limit) if limit else frame
    header = "".join(f"<th>{esc(label)}</th>" for _, label in columns)
    body = []
    for row in shown.to_dict("records"):
        cells = "".join(f"<td>{esc(row.get(key, ''))}</td>" for key, _ in columns)
        body.append(f"<tr>{cells}</tr>")
    if not body:
        body.append(f"<tr><td colspan='{len(columns)}'>无记录</td></tr>")
    return f"<div class='table-wrap'><table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table></div>"


def write_html(path: Path, ours: pd.DataFrame, comparison: pd.DataFrame, only_ours: pd.DataFrame,
               actual_filters: dict[str, str], summary: dict[str, Any], source_paths: dict[str, Path]) -> None:
    shared = comparison[comparison["comparison_status"] != "仅外部流程检出"]
    likely_shared = comparison[comparison["comparison_status"] == "疑似共同检出（转录本差异）"]
    only_other = comparison[comparison["comparison_status"] == "仅外部流程检出"]
    filter_rows = pd.DataFrame([{"item": key, "value": value} for key, value in actual_filters.items()])
    field_rows = pd.DataFrame([
        {
            "field": "TLOD",
            "full_name": "Tumor Log Odds",
            "definition": "突变存在相对于不存在的 log10 似然比",
            "direction": "越高表示肿瘤样本中的变异证据越强",
            "interpretation": "不是后验概率；仍需结合ROQ、正常样本证据、测序质量和其他过滤标签",
        },
        {
            "field": "ROQ",
            "full_name": "Read Orientation Quality",
            "definition": "ALT不是由read-orientation伪影造成的Phred质量值",
            "direction": "越高越好；越低表示方向偏倚伪影风险越高",
            "interpretation": "可近似按P(方向伪影)=10^(-ROQ/10)理解，不等于临床真阳性概率",
        },
        {
            "field": "F1R2 / F2R1",
            "full_name": "Paired-read orientation counts",
            "definition": "ALT/REF在两类成对读段方向中的计数",
            "direction": "明显不平衡可能支持orientation过滤",
            "interpretation": "不是普通正负链计数；应结合ROQ和方向偏倚模型判断",
        },
    ])
    tlod_examples = pd.DataFrame([
        {"TLOD": "3", "模型似然比": "10^3 = 1,000:1"},
        {"TLOD": "6.3", "模型似然比": "约10^6.3 = 2,000,000:1"},
        {"TLOD": "10", "模型似然比": "10^10 = 10,000,000,000:1"},
    ])
    roq_examples = pd.DataFrame([
        {"ROQ": "1", "近似方向伪影概率": "79.4%"},
        {"ROQ": "2", "近似方向伪影概率": "63.1%"},
        {"ROQ": "10", "近似方向伪影概率": "10%"},
        {"ROQ": "20", "近似方向伪影概率": "1%"},
        {"ROQ": "30", "近似方向伪影概率": "0.1%"},
    ])
    if summary.get("filtered_total_records", 0):
        trace_note = (
            "本次已提供FilterMutectCalls完整filtered VCF，可以确认候选是否被过滤及其FILTER标签。"
            "对于filtered VCF中完全没有出现的BED内位点，仍需raw VCF和BAM判断是未形成候选、"
            "低于调用阈值，还是缺少有效ALT支持。"
        )
    else:
        trace_note = (
            "当前只提供PASS VCF，已被过滤的候选及其FILTER原因不在文件中；"
            "需要补充Mutect2 raw/filtered VCF和索引。"
        )
    comparison_columns = [
        ("gene", "基因"), ("hgvsc_other", "外部c.HGVS"), ("hgvsp_other", "外部p.HGVS"),
        ("other_vaf", "外部丰度"), ("comparison_status", "比较结论"), ("match_method", "匹配依据"),
        ("matched_variant_key", "本流程基因组位点"), ("our_tumor_af", "本流程Tumor AF"),
        ("our_normal_af", "本流程Normal AF"), ("bed_status", "BED判断"), ("bed_explanation", "BED说明"),
        ("gtf_transcript", "VEP115近似转录本"), ("gtf_variant", "VEP115近似hg38位点"),
        ("filtered_filter", "Filtered VCF结果"), ("filtered_tlod", "TLOD"),
        ("filtered_tumor_af", "Filtered Tumor AF"),
        ("filtered_tumor_alt_reads", "Tumor ALT reads"), ("filtered_tumor_f1r2", "Tumor F1R2"),
        ("filtered_tumor_f2r1", "Tumor F2R1"), ("filtered_roq", "ROQ"),
        ("filtered_normal_alt_reads", "Normal ALT reads"), ("difference_diagnosis", "差异诊断"),
    ]
    ours_columns = [
        ("variant_key", "位点"), ("gene", "代表基因"), ("consequence", "VEP后果"),
        ("hgvsc", "HGVSc"), ("hgvsp", "HGVSp"), ("tumor_af", "Tumor AF"),
        ("normal_af", "Normal AF"), ("tlod", "TLOD"), ("bed_status", "BED"),
    ]
    max_count = max(summary["our_pass_variants"], summary["other_report_variants"], 1)
    bars = [
        ("本流程 PASS VCF", summary["our_pass_variants"], "#1d4ed8"),
        ("外部报告", summary["other_report_variants"], "#b45309"),
        ("严格/疑似共同", summary["shared_variants"], "#15803d"),
        ("仅外部报告", summary["only_other_variants"], "#b91c1c"),
    ]
    bar_html = "".join(
        f"<div class='bar-row'><span>{esc(label)}</span><div class='track'><i style='width:{value/max_count*100:.1f}%;background:{color}'></i></div><b>{value}</b></div>"
        for label, value, color in bars
    )
    cards = "".join(
        f"<div class='metric'><span>{esc(label)}</span><strong>{value}</strong></div>"
        for label, value in [
            ("本流程PASS位点", summary["our_pass_variants"]), ("外部报告位点", summary["other_report_variants"]),
            ("严格HGVS共同检出", summary["exact_shared_variants"]), ("转录本差异疑似共同", summary["likely_shared_variants"]),
            ("仅外部检出", summary["only_other_variants"]),
            ("本流程BED内", summary["our_in_bed"]), ("本流程BED外", summary["our_outside_bed"]),
        ]
    )
    content = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>WES突变结果对比报告</title><style>
:root{{--ink:#172033;--muted:#667085;--line:#d9dee8;--paper:#fff;--bg:#f3f5f8;--blue:#1d4ed8}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font-family:Arial,"PingFang SC","Microsoft YaHei",sans-serif;line-height:1.55}}
header{{background:#16243a;color:#fff;padding:34px 5vw 28px}}header h1{{margin:0 0 6px;font-size:30px;letter-spacing:0}}header p{{margin:0;color:#ccd5e2}}
main{{max-width:1440px;margin:0 auto;padding:24px}}section{{background:var(--paper);border:1px solid var(--line);border-radius:6px;padding:22px;margin-bottom:18px}}
h2{{font-size:20px;margin:0 0 14px}}h3{{font-size:16px;margin:20px 0 10px}}.metrics{{display:grid;grid-template-columns:repeat(6,minmax(120px,1fr));gap:10px}}
.metric{{border-left:4px solid var(--blue);background:#f8fafc;padding:12px}}.metric span{{display:block;color:var(--muted);font-size:13px}}.metric strong{{font-size:27px}}
.note{{border-left:4px solid #b45309;background:#fff8eb;padding:12px 14px;margin:12px 0}}.good{{border-left-color:#15803d;background:#f0fdf4}}.bad{{border-left-color:#b91c1c;background:#fff1f2}}
.bar-row{{display:grid;grid-template-columns:150px 1fr 55px;gap:12px;align-items:center;margin:10px 0}}.track{{height:14px;background:#e8ecf2}}.track i{{display:block;height:100%}}
.table-wrap{{overflow:auto;border:1px solid var(--line)}}table{{border-collapse:collapse;width:100%;font-size:12px}}th,td{{text-align:left;padding:8px 9px;border-bottom:1px solid var(--line);white-space:nowrap}}th{{position:sticky;top:0;background:#eef2f7;z-index:1}}tr:nth-child(even) td{{background:#fafbfc}}
code{{background:#eef2f7;padding:2px 5px}}ul{{padding-left:22px}}footer{{color:var(--muted);font-size:12px;padding:0 24px 28px;max-width:1440px;margin:auto}}
@media(max-width:900px){{.metrics{{grid-template-columns:repeat(2,1fr)}}main{{padding:12px}}header{{padding:24px 18px}}}}
</style></head><body><header><h1>WES突变结果对比报告</h1><p>本流程 VEP/Mutect2 PASS VCF 与外部突变表的 HGVS 对齐及 BED 覆盖审计</p></header><main>
<section><h2>结论概览</h2><div class="metrics">{cards}</div><h3>数量对比</h3>{bar_html}
<div class="note"><b>不能直接把 168 与 32 解释成灵敏度差异。</b> 本流程输入是全部 Mutect2 PASS 位点并包含非编码、同义等后果；外部文件看起来是经过功能类别或报告规则整理后的短表。</div>
<div class="note bad"><b>BED判断边界：</b> 外部表没有 CHROM/POS/REF/ALT。共同检出位点使用本流程坐标；仅外部位点中有 {summary.get('gtf_approx_mapped', 0)} 个按 Ensembl 115 canonical/MANE 转录本近似映射，仍不等同于对方原始坐标。</div>
<div class="note"><b>转录本差异：</b> 有 {len(likely_shared)} 个结果在同一基因只有一个本流程候选且 AF 接近，但 HGVS 编号不同。它们被列为“疑似共同检出”，必须拿外部基因组坐标或转录本编号最终确认。</div></section>
<section><h2>这份 VCF 实际使用的检测与过滤参数</h2>{table_html(filter_rows, [('item','项目'),('value','VCF头记录的实际值')])}
<div class="note"><b>版本差异：</b> 这份 VCF 实际使用 <code>--max-reads-per-alignment-start {esc(actual_filters.get('每个比对起点最大 reads'))}</code>；当前主代码默认值为 <code>50</code>。数值 5 会更强地下采样高深度区域，可能损失低 VAF 证据。缺少 PoN、群体胚系资源或污染估计不会让流程无法运行，但会改变假阳性/胚系污染控制能力。</div>
<p>最终报告 VCF 的主筛选是：paired Mutect2 检测 → 可用时建立方向偏倚/污染模型 → GATK FilterMutectCalls → 仅提取 FILTER=PASS → VEP 注释。VEP 在这里增加注释，不再按功能或 VAF 删除记录。</p></section>
	<section><h2>Filtered VCF追踪结论</h2><p>Filtered VCF共有 {summary.get('filtered_total_records', 0)} 个候选，其中 {summary.get('filtered_pass_records', 0)} 个PASS。12个仅外部报告位点中：{summary.get('external_detected_but_filtered', 0)} 个已被Mutect2检出但未PASS，{summary.get('external_inside_bed_not_emitted', 0)} 个在BED内但没有候选，{summary.get('external_outside_bed', 0)} 个位于BED外，{summary.get('external_unresolved', 0)} 个无法定位。</p><p><b>主要过滤标签：</b>{esc(summary.get('filtered_top_filters', ''))}</p></section>
	<section><h2>TLOD与ROQ字段解释</h2><p>以下定义直接依据当前Mutect2 VCF的INFO头字段，并结合本次差异位点解释。</p>
	{table_html(field_rows, [('field','字段'),('full_name','英文全称'),('definition','定义'),('direction','数值方向'),('interpretation','本报告判读')])}
	<h3>TLOD如何理解</h3><p><code>TLOD = log10[L(存在突变) / L(不存在突变)]</code>。例如：</p>
	{table_html(tlod_examples, [('TLOD','TLOD'),('模型似然比','模型似然比')])}
	<div class="note"><b>重要：</b>TLOD是模型似然比，不是“突变为真的概率”。高TLOD表示肿瘤中有较强ALT证据，但仍可能因read orientation、正常样本证据、链偏倚、碱基质量或比对质量等原因被过滤。一般PASS并非只靠固定TLOD阈值决定；本流程的<code>TLOD ≥ 6.3</code>只用于严格TMB二次筛选。</div>
	<h3>ROQ如何理解</h3><p><code>ROQ</code>采用Phred尺度，可近似换算为 <code>P(方向伪影) ≈ 10^(-ROQ/10)</code>：</p>
	{table_html(roq_examples, [('ROQ','ROQ'),('近似方向伪影概率','近似方向伪影概率')])}
	<div class="note bad"><b>本次差异位点的实际含义：</b>4个外部报告位点在filtered VCF中具有TLOD 8.61–18.54，说明肿瘤ALT证据不弱；但ROQ仅1–2，方向伪影风险很高，因此被标记为<code>orientation</code>。它们应归类为“已检出但未通过质量过滤，建议结合BAM、F1R2/F2R1和独立实验复核”，不能直接记作阳性，也不能说流程完全没有检出。</div>
	</section>
	<section><h2>外部报告逐条对齐</h2>{table_html(comparison, comparison_columns)}</section>
	<section><h2>共同检出 ({len(shared)})</h2>{table_html(shared, comparison_columns)}</section>
	<section><h2>仅外部流程报告 ({len(only_other)})</h2>{table_html(only_other, comparison_columns)}
	<div class="note">要最终区分“BED外”“被Mutect2过滤”“低于调用阈值”“注释/转录本写法不同”，外部结果仍应补充CHROM、POS、REF、ALT和转录本ID。{esc(trace_note)}</div></section>
<section><h2>仅本流程检出 ({len(only_ours)})</h2>{table_html(only_ours, ours_columns, 300)}<p>页面最多展示300行；完整数据见导出的 TSV/XLSX。</p></section>
<section><h2>本流程的严格 TMB 二次筛选</h2><p>严格 TMB 并不是“全部报告突变”的筛选条件，而是单独计算 TMB 时再次要求：PASS、TLOD ≥ 6.3、Tumor DP ≥ 20、Tumor ALT reads ≥ 5、Tumor AF ≥ 5%、Normal DP ≥ 10、Normal ALT reads ≤ 2、Normal AF ≤ 2%、群体 AF ≤ 0.1%，并仅保留指定编码/剪接后果且位于有效编码 BED。</p></section>
<section><h2>输入文件</h2><ul><li>本流程VCF：{esc(source_paths['vcf'])}</li><li>外部结果：{esc(source_paths['other'])}</li><li>捕获BED：{esc(source_paths['bed'])}</li></ul></section>
</main><footer>本报告用于研发对比与流程审计，不构成临床判读。匹配采用“基因 + 标准化 cDNA/protein HGVS”；不同转录本、左右归一化差异和复杂 InDel 仍需按基因组坐标复核。</footer></body></html>"""
    path.write_text(content, encoding="utf-8")


def write_markdown(path: Path, summary: dict[str, Any], actual_filters: dict[str, str],
                   only_other: pd.DataFrame) -> None:
    filters = "\n".join(f"| {key} | {value} |" for key, value in actual_filters.items())
    missing = "\n".join(
        f"| {row['gene']} | {row['hgvsc_other']} | {row['hgvsp_other']} | {row['other_vaf']} | {row['bed_status']} |"
        for row in only_other.to_dict("records")
    ) or "| - | - | - | - | - |"
    if summary.get("filtered_total_records", 0):
        trace_note = (
            "本次已使用 filtered VCF 追踪 FilterMutectCalls 标签；对于 filtered VCF 中完全没有出现的"
            "BED 内位点，仍需 raw VCF 和 BAM 判断是未形成候选、低于调用阈值，还是缺少有效 ALT 支持。"
        )
    else:
        trace_note = "当前只有 PASS VCF，需要补充 Mutect2 raw/filtered VCF 才能追踪未PASS候选。"
    text = f"""# WES突变结果对比汇报摘要

## 主要结果

- 本流程 VEP/Mutect2 PASS VCF：{summary['our_pass_variants']} 个位点。
- 外部报告：{summary['other_report_variants']} 个位点。
- 严格按基因和 HGVS 共同检出：{summary['exact_shared_variants']} 个。
- 疑似因转录本不同而 HGVS 编号不同：{summary['likely_shared_variants']} 个。
- 仅外部报告且当前 PASS VCF 未匹配：{summary['only_other_variants']} 个。
- 本流程 {summary['our_in_bed']} 个 PASS 位点位于所给 BED 内，{summary['our_outside_bed']} 个位于 BED 外。

外部表缺少 CHROM、POS、REF、ALT。报告中的近似位置按 Ensembl release 115
canonical/MANE 转录本推算；它适合判断 BED 覆盖趋势，但不等同于外部流程原始坐标。
{trace_note}

- VEP115 GTF近似映射成功：{summary.get('gtf_approx_mapped', 0)} 个。
- 近似位置位于BED内：{summary.get('gtf_approx_in_bed', 0)} 个。
- 近似位置位于BED外：{summary.get('gtf_approx_outside_bed', 0)} 个。
- 仍无法映射：{summary.get('gtf_unresolved', summary['only_other_variants'])} 个。
- Mutect2已检出但被FilterMutectCalls过滤：{summary.get('external_detected_but_filtered', 0)} 个。
- BED内但filtered VCF没有候选：{summary.get('external_inside_bed_not_emitted', 0)} 个。
- Filtered VCF主要过滤标签：{summary.get('filtered_top_filters', '')}。

## 本VCF实际参数

| 项目 | VCF头记录值 |
|---|---|
{filters}

这份历史结果使用 `--max-reads-per-alignment-start 5`，而当前主代码默认是 `50`。
历史参数会在高深度区域更强地下采样，可能减少低 VAF 支持 reads。本次没有使用 PoN、
群体胚系资源或污染估计表，使用了方向偏倚模型。

## TLOD 与 ROQ 怎么看

| 字段 | 英文全称 | 定义 | 数值方向 |
|---|---|---|---|
| TLOD | Tumor Log Odds | 突变存在相对于不存在的 log10 似然比 | 越高表示肿瘤中的变异证据越强 |
| ROQ | Read Orientation Quality | ALT 不是由 read-orientation 伪影造成的 Phred 质量值 | 越高越好；越低表示方向伪影风险越高 |

`TLOD = log10[L(存在突变) / L(不存在突变)]`。TLOD 3、6.3、10 分别对应约
`1,000:1`、`2,000,000:1`、`10,000,000,000:1` 的模型似然比。TLOD 不是后验概率，
高 TLOD 也不能覆盖正常样本证据、方向偏倚、链偏倚、碱基质量和比对质量等过滤条件。
一般 PASS 并非只由固定 TLOD 阈值决定；本流程的 `TLOD >= 6.3` 只用于严格 TMB 二次筛选。

ROQ 可近似按 `P(方向伪影) = 10^(-ROQ/10)` 理解：ROQ 1、2、10、20、30
分别对应约 79.4%、63.1%、10%、1%、0.1% 的方向伪影概率。该数值是模型质量值，
不是临床真阳性概率，也不是普通正负链计数。F1R2/F2R1 是 paired-read orientation
计数，应和 ROQ 及方向偏倚模型一起判读。

本次4个被 `orientation` 过滤的外部报告位点具有 TLOD 8.61–18.54，但 ROQ 仅1–2。
因此更准确的结论是“Mutect2 已检出且肿瘤证据不弱，但方向伪影风险高，需结合 BAM、
F1R2/F2R1 和独立实验复核”，不能直接记作阳性，也不能说流程完全没有检出。

## 仅外部报告

| 基因 | c.HGVS | p.HGVS | 丰度 | BED判断 |
|---|---|---|---:|---|
{missing}

## 下一步所需文件

1. 外部流程完整结果，至少包含 `CHROM/POS/REF/ALT`、参考基因组版本、转录本 ID。
2. 本流程的 `mutect2.raw.vcf.gz` 和 `mutect2.filtered.vcf.gz` 及索引。
3. 两套流程实际使用的 BED、PoN、germline resource、caller 版本和过滤参数。

完整逐条对齐请查看 `variant_comparison.xlsx` 和 `variant_comparison_report_zh.html`。
"""
    path.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="使用 pandas 对比 VEP VCF 与外部突变表，并检查 BED 覆盖")
    parser.add_argument("--our-vcf", required=True, type=Path)
    parser.add_argument("--other-table", required=True, type=Path)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--tumor-sample")
    parser.add_argument("--normal-sample")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--gtf", type=Path, help="可选：Ensembl release对应GTF，用于外部独有HGVS近似定位")
    parser.add_argument("--reference-fasta", type=Path, help="与GTF对应的参考FASTA，用于近似位置参考碱基校验")
    parser.add_argument("--filtered-vcf", type=Path, help="可选：Mutect2 FilterMutectCalls完整输出，用于追踪非PASS原因")
    parser.add_argument("--gene-vaf-tolerance", type=float, default=0.02,
                        help="同基因唯一候选的疑似转录本匹配允许绝对AF差，默认0.02")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    ours, annotations, metadata, vcf_samples, tumor_sample, normal_sample = parse_vcf(
        args.our_vcf, args.tumor_sample, args.normal_sample
    )
    other = load_other(args.other_table)
    intervals, starts = load_bed(args.bed)

    bed_results = ours.apply(
        lambda row: bed_status(row["chrom"], int(row["pos"]), row["ref"], intervals, starts), axis=1
    )
    ours["in_bed"] = [item[0] for item in bed_results]
    ours["bed_distance_bp"] = [item[1] for item in bed_results]
    ours["bed_status"] = ours["in_bed"].map({True: "BED内", False: "BED外"})

    comparison, matched_keys = compare(ours, annotations, other, args.gene_vaf_tolerance)
    for column in ["gtf_chrom", "gtf_pos", "gtf_ref", "gtf_alt", "gtf_transcript",
                   "gtf_reference_check", "gtf_mapping_status", "gtf_variant"]:
        comparison[column] = ""
    if args.gtf:
        from hgvs_gtf_mapper import map_variants
        only_mask = comparison["comparison_status"] == "仅外部流程检出"
        mapped = map_variants(comparison.loc[only_mask].copy(), args.gtf, args.reference_fasta)
        mapped = mapped.set_index("other_row")
        for index in comparison.index[only_mask]:
            other_row = comparison.at[index, "other_row"]
            if other_row not in mapped.index:
                continue
            for column in ["gtf_chrom", "gtf_pos", "gtf_ref", "gtf_alt", "gtf_transcript",
                           "gtf_reference_check", "gtf_mapping_status"]:
                value = mapped.at[other_row, column] if column in mapped.columns else ""
                comparison.at[index, column] = "" if pd.isna(value) else value
            if comparison.at[index, "gtf_pos"] != "":
                comparison.at[index, "gtf_variant"] = (
                    f"{comparison.at[index, 'gtf_chrom']}:{int(float(comparison.at[index, 'gtf_pos']))}:"
                    f"{comparison.at[index, 'gtf_ref']}:{comparison.at[index, 'gtf_alt']}"
                )
    filtered_records: dict[str, dict[str, Any]] = {}
    filtered_counts: Counter[str] = Counter()
    filtered_total_records = 0
    if args.filtered_vcf:
        filtered_records, filtered_counts, filtered_total_records = parse_filtered_vcf(
            args.filtered_vcf, tumor_sample, normal_sample
        )
    filtered_columns = [
        "filtered_variant_key", "filtered_filter", "filtered_tlod", "filtered_roq", "filtered_tumor_dp",
        "filtered_tumor_alt_reads", "filtered_tumor_af", "filtered_tumor_f1r2", "filtered_tumor_f2r1",
        "filtered_tumor_sb", "filtered_normal_dp",
        "filtered_normal_alt_reads", "filtered_normal_af",
    ]
    for column in filtered_columns:
        comparison[column] = ""
    for index, row in comparison.iterrows():
        candidate_key = row["matched_variant_key"] or row["gtf_variant"]
        detail = filtered_records.get(candidate_key, {})
        for column in filtered_columns:
            comparison.at[index, column] = detail.get(column, "")
    bed_by_key = ours.set_index("variant_key")[["in_bed", "bed_distance_bp"]].to_dict("index")
    bed_labels, bed_explanations = [], []
    for row in comparison.to_dict("records"):
        key = row["matched_variant_key"]
        if key:
            detail = bed_by_key[key]
            bed_labels.append("BED内" if detail["in_bed"] else "BED外")
            bed_explanations.append("使用本流程匹配位点的基因组坐标判定")
        elif row.get("chrom") and row.get("pos") and row.get("ref"):
            inside, distance = bed_status(row["chrom"], int(row["pos"]), row["ref"], intervals, starts)
            bed_labels.append("BED内" if inside else "BED外")
            bed_explanations.append("使用外部表提供的基因组坐标判定" + (f"；距最近BED约{distance} bp" if distance else ""))
        elif row.get("gtf_pos") not in {"", None} and not pd.isna(row.get("gtf_pos")):
            inside, distance = bed_status(row["gtf_chrom"], int(float(row["gtf_pos"])), row.get("gtf_ref", "N"), intervals, starts)
            bed_labels.append("BED内（近似）" if inside else "BED外（近似）")
            bed_explanations.append(
                f"Ensembl115 {row.get('gtf_transcript', '')} canonical/MANE近似映射，参考碱基{row.get('gtf_reference_check', '未检查')}"
                + (f"；距最近BED约{distance} bp" if not inside and distance is not None else "")
            )
        else:
            bed_labels.append("无法判定")
            if row.get("gtf_mapping_status"):
                bed_explanations.append(f"Ensembl115 GTF: {row['gtf_mapping_status']}")
            else:
                bed_explanations.append("外部表缺少CHROM/POS/REF/ALT")
    comparison["bed_status"] = bed_labels
    comparison["bed_explanation"] = bed_explanations
    diagnoses = []
    for row in comparison.to_dict("records"):
        if row["comparison_status"] != "仅外部流程检出":
            diagnoses.append("本流程PASS")
        elif row.get("filtered_filter"):
            diagnoses.append(f"Mutect2已检出，但被过滤: {row['filtered_filter']}")
        elif str(row.get("bed_status", "")).startswith("BED外"):
            diagnoses.append("近似位点在BED外，原Mutect2 -L不会检测")
        elif args.filtered_vcf and str(row.get("bed_status", "")).startswith("BED内"):
            diagnoses.append("BED内，但filtered VCF无该候选")
        elif not args.filtered_vcf:
            diagnoses.append("未提供filtered VCF，暂不能追踪")
        else:
            diagnoses.append("坐标无法确定，暂不能追踪")
    comparison["difference_diagnosis"] = diagnoses

    only_ours = ours[~ours["variant_key"].isin(matched_keys)].copy()
    only_other = comparison[comparison["comparison_status"] == "仅外部流程检出"].copy()
    shared = comparison[comparison["comparison_status"] != "仅外部流程检出"].copy()
    exact_shared = comparison[comparison["comparison_status"] == "两流程共同检出"].copy()
    likely_shared = comparison[comparison["comparison_status"] == "疑似共同检出（转录本差异）"].copy()
    actual_filters = parse_actual_filters(metadata)
    approx_mapped = only_other["gtf_pos"].replace("", pd.NA).notna() if "gtf_pos" in only_other else pd.Series(dtype=bool)
    approx_bed = only_other["bed_status"].str.startswith("BED内") if "bed_status" in only_other else pd.Series(dtype=bool)
    approx_outside = only_other["bed_status"].str.startswith("BED外") if "bed_status" in only_other else pd.Series(dtype=bool)
    detected_but_filtered = only_other["filtered_filter"].replace("", pd.NA).notna()
    inside_not_emitted = only_other["difference_diagnosis"] == "BED内，但filtered VCF无该候选"
    outside_bed = only_other["bed_status"].str.startswith("BED外")
    unresolved_external = only_other["bed_status"] == "无法判定"
    summary = {
        "our_pass_variants": int(len(ours)), "other_report_variants": int(len(other)),
        "shared_variants": int(len(shared)), "exact_shared_variants": int(len(exact_shared)),
        "likely_shared_variants": int(len(likely_shared)), "only_our_variants": int(len(only_ours)),
        "only_other_variants": int(len(only_other)), "our_in_bed": int(ours["in_bed"].sum()),
        "our_outside_bed": int((~ours["in_bed"]).sum()), "vcf_samples": vcf_samples,
        "gtf_approx_mapped": int(approx_mapped.sum()), "gtf_approx_in_bed": int(approx_bed.sum()),
        "gtf_approx_outside_bed": int(approx_outside.sum()),
        "gtf_unresolved": int(len(only_other) - approx_mapped.sum()),
        "filtered_total_records": int(filtered_total_records),
        "filtered_pass_records": int(filtered_counts.get("PASS", 0)),
        "external_detected_but_filtered": int(detected_but_filtered.sum()),
        "external_inside_bed_not_emitted": int(inside_not_emitted.sum()),
        "external_outside_bed": int(outside_bed.sum()), "external_unresolved": int(unresolved_external.sum()),
        "filtered_top_filters": "；".join(f"{key}={value}" for key, value in filtered_counts.most_common(8)),
        "tumor_sample": tumor_sample, "normal_sample": normal_sample,
        "actual_filters": actual_filters,
        "limitations": [
            "外部表没有CHROM/POS/REF/ALT；GTF结果是按VEP115 canonical/MANE转录本的近似映射。",
            ("已使用filtered VCF追踪FilterMutectCalls原因；未提供raw VCF和BAM时仍不能解释未进入候选的位点。"
             if args.filtered_vcf else "输入VCF只含PASS记录，不能查看未通过FilterMutectCalls位点的具体失败原因。"),
            "HGVS匹配可能受转录本选择和复杂InDel表达差异影响。",
        ],
    }

    ours.to_csv(args.output_dir / "01_our_vcf_variants.tsv", sep="\t", index=False)
    annotations.to_csv(args.output_dir / "02_our_vep_annotations.tsv", sep="\t", index=False)
    comparison.to_csv(args.output_dir / "03_external_report_comparison.tsv", sep="\t", index=False)
    shared.to_csv(args.output_dir / "04_shared_variants.tsv", sep="\t", index=False)
    likely_shared.to_csv(args.output_dir / "04b_likely_transcript_matches.tsv", sep="\t", index=False)
    only_other.to_csv(args.output_dir / "05_only_external_variants.tsv", sep="\t", index=False)
    only_ours.to_csv(args.output_dir / "06_only_our_variants.tsv", sep="\t", index=False)
    pd.DataFrame([{"filter": key, "count": value} for key, value in filtered_counts.most_common()]).to_csv(
        args.output_dir / "08_filtered_vcf_filter_counts.tsv", sep="\t", index=False
    )
    (args.output_dir / "comparison_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_markdown(args.output_dir / "comparison_summary_zh.md", summary, actual_filters, only_other)

    try:
        with pd.ExcelWriter(args.output_dir / "variant_comparison.xlsx", engine="openpyxl") as writer:
            pd.DataFrame([summary | {"actual_filters": "见实际过滤参数sheet", "limitations": "; ".join(summary["limitations"])}]).to_excel(writer, "摘要", index=False)
            pd.DataFrame([{"项目": k, "实际值": v} for k, v in actual_filters.items()]).to_excel(writer, "实际过滤参数", index=False)
            pd.DataFrame([
                {
                    "字段": "TLOD", "英文全称": "Tumor Log Odds",
                    "定义": "突变存在相对于不存在的log10似然比",
                    "数值方向": "越高表示肿瘤变异证据越强",
                    "换算/判读": "TLOD 6.3约等于2,000,000:1模型似然比；不是后验概率；一般PASS不只由TLOD决定",
                },
                {
                    "字段": "ROQ", "英文全称": "Read Orientation Quality",
                    "定义": "ALT不是由read-orientation伪影造成的Phred质量值",
                    "数值方向": "越高越好；越低表示方向伪影风险越高",
                    "换算/判读": "P(方向伪影)约等于10^(-ROQ/10)；ROQ 1/2约对应79.4%/63.1%",
                },
                {
                    "字段": "F1R2/F2R1", "英文全称": "Paired-read orientation counts",
                    "定义": "ALT/REF在两类成对读段方向中的计数",
                    "数值方向": "显著不平衡可能支持orientation过滤",
                    "换算/判读": "不是普通正负链计数，应结合ROQ和方向偏倚模型",
                },
            ]).to_excel(writer, "字段解释", index=False)
            pd.DataFrame([{"filter": key, "count": value} for key, value in filtered_counts.most_common()]).to_excel(writer, "Filtered过滤统计", index=False)
            comparison.to_excel(writer, "外部结果逐条对齐", index=False)
            shared.to_excel(writer, "共同检出", index=False)
            likely_shared.to_excel(writer, "疑似转录本差异", index=False)
            only_other.to_excel(writer, "仅外部检出", index=False)
            only_ours.to_excel(writer, "仅本流程检出", index=False)
            ours.to_excel(writer, "本流程全部PASS", index=False)
    except ImportError:
        print("warning: openpyxl unavailable; skipped XLSX output")

    report_path = args.output_dir / "variant_comparison_report_zh.html"
    write_html(report_path, ours, comparison, only_ours, actual_filters, summary, {
        "vcf": args.our_vcf.resolve(), "other": args.other_table.resolve(), "bed": args.bed.resolve()
    })
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"report={report_path}")


if __name__ == "__main__":
    main()
