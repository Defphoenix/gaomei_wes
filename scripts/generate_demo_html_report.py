#!/usr/bin/env python3
"""Build a Bootstrap/ECharts HTML report for synthetic WES demo projects."""

from __future__ import annotations

import argparse
import html
import json
import os
import subprocess
from pathlib import Path

BCFTOOLS = os.environ.get("BCFTOOLS", "bcftools")

ACTIONABILITY = {
    "ERBB2_HER2": {
        "alteration": "HER2/ERBB2 改变或 HER2 阳性相关肿瘤",
        "therapy": "可关注 trastuzumab deruxtecan；HER2 阳性场景下还可能涉及 trastuzumab/pertuzumab/T-DM1 等",
        "evidence": "与癌种、HER2 检测和药品适应证强相关",
        "note": "需要结合 HER2 IHC/ISH 或经验证的 ERBB2 改变解读；仅凭 WES SNV 不能直接指导用药。",
    },
    "ERBB2": {
        "alteration": "HER2/ERBB2 改变或 HER2 阳性相关肿瘤",
        "therapy": "可关注 trastuzumab deruxtecan；HER2 阳性场景下还可能涉及 trastuzumab/pertuzumab/T-DM1 等",
        "evidence": "与癌种、HER2 检测和药品适应证强相关",
        "note": "需要结合 HER2 IHC/ISH 或经验证的 ERBB2 改变解读；仅凭 WES SNV 不能直接指导用药。",
    },
    "BRCA1": {
        "alteration": "BRCA1/2 致病或疑似致病改变、HRD 相关场景",
        "therapy": "可关注 PARP 抑制剂，例如 olaparib；部分场景下铂类敏感性也有参考意义",
        "evidence": "与癌种、胚系/体系来源、致病性和适应证强相关",
        "note": "需要判断变异致病性和胚系/体系来源；BRCA1 VUS 不能直接视为可用药靶点。",
    },
    "TP53": {
        "alteration": "TP53 突变",
        "therapy": "目前缺少广泛批准的 TP53 直接靶向治疗；可关注临床试验和标准治疗背景",
        "evidence": "临床试验/预后或生物学标志物",
        "note": "TP53 通常具有重要生物学意义，但常规临床中一般不是直接可用药靶点。",
    },
}


def run_lines(cmd: list[str]) -> list[str]:
    try:
        out = subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL)
    except Exception:
        return []
    return [line for line in out.splitlines() if line.strip()]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    if not lines:
        return []
    cols = lines[0].split("\t")
    return [dict(zip(cols, line.split("\t"))) for line in lines[1:] if line.strip()]


def count_vcf(vcf: Path) -> int:
    if not vcf.exists():
        return 0
    return len(run_lines([BCFTOOLS, "view", "-H", str(vcf)]))


def query_vcf(vcf: Path) -> list[dict[str, str]]:
    rows = []
    if not vcf.exists():
        return rows
    fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER[\t%AD\t%AF\t%DP]\n"
    for line in run_lines([BCFTOOLS, "query", "-f", fmt, str(vcf)]):
        parts = line.split("\t")
        row = {
            "chrom": parts[0],
            "pos": parts[1],
            "id": parts[2],
            "ref": parts[3],
            "alt": parts[4],
            "filter": parts[5],
            "tumor_ad": parts[-3] if len(parts) >= 9 else "",
            "tumor_af": parts[-2] if len(parts) >= 9 else "",
            "tumor_dp": parts[-1] if len(parts) >= 9 else "",
        }
        rows.append(row)
    return rows


def truth_lookup(truth: list[dict[str, str]]) -> dict[tuple[str, str, str, str], dict[str, str]]:
    lookup = {}
    for row in truth:
        key = (row.get("chrom", ""), row.get("pos", ""), row.get("ref", ""), row.get("alt", ""))
        lookup[key] = row
    return lookup


def annotate_gene(row: dict[str, str], truth_by_variant: dict[tuple[str, str, str, str], dict[str, str]]) -> str:
    key = (row.get("chrom", ""), row.get("pos", ""), row.get("ref", ""), row.get("alt", ""))
    truth_row = truth_by_variant.get(key, {})
    return truth_row.get("gene", "")


def read_neoantigen_manifest(path: Path, pair_id: str) -> list[dict[str, str]]:
    peptide_table = path.with_name(path.name.replace("_neoantigen_manifest.tsv", "_neoantigen_peptides.tsv"))
    source = peptide_table if peptide_table.exists() else path
    rows = read_tsv(source)
    if not rows:
        return []
    rows = []
    for row in read_tsv(source):
        if source == path and row.get("status") != "emitted":
            continue
        rows.append(
            {
                "pair_id": pair_id,
                "gene": row.get("gene", ""),
                "variant": row.get("variant_id", ""),
                "effect": row.get("event_type") or row.get("consequence", ""),
                "mer": row.get("mer", ""),
                "peptide_start": row.get("mutant_peptide_start", ""),
                "peptide_end": row.get("mutant_peptide_end", ""),
                "changed_aa": row.get("amino_acids", ""),
                "wt_peptide": row.get("wildtype_peptide") or row.get("wildtype_window", ""),
                "mt_peptide": row.get("mutant_peptide") or row.get("mutant_window", ""),
                "logic": "来自VEP蛋白位置/氨基酸改变 + 配置的蛋白FASTA；枚举所有覆盖突变位点的指定mer候选肽。",
            }
        )
    return rows


def make_actionability_rows(variant_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    seen = set()
    rows = []
    for row in variant_rows:
        gene = row.get("gene", "")
        if not gene or gene in seen:
            continue
        seen.add(gene)
        item = ACTIONABILITY.get(gene, {})
        if not item:
            rows.append(
                {
                    "gene": gene,
                    "alteration": "检出突变",
                    "therapy": "当前演示知识库暂无条目",
                    "evidence": "NA",
                    "note": "正式版本建议接入 CIViC/OncoKB/CGI 或院内知识库。",
                }
            )
            continue
        rows.append({"gene": gene, **item})
    return rows


def parse_fastp(path: Path) -> dict[str, float]:
    if not path.exists():
        return {}
    data = json.loads(path.read_text())
    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    return {
        "before_reads": before.get("total_reads", 0),
        "after_reads": after.get("total_reads", 0),
        "q30_rate": after.get("q30_rate", 0.0),
    }


def parse_flagstat(path: Path) -> dict[str, str]:
    if not path.exists():
        return {"mapped_rate": "NA", "proper_rate": "NA"}
    text = path.read_text()
    mapped = "NA"
    proper = "NA"
    for line in text.splitlines():
        if " mapped (" in line and "primary" not in line:
            mapped = line.split("(")[1].split(":")[0].strip()
        if " properly paired " in line:
            proper = line.split("(")[1].split(":")[0].strip()
    return {"mapped_rate": mapped, "proper_rate": proper}


def parse_coverage(path: Path) -> dict[str, str]:
    if not path.exists():
        return {"mean_depth": "NA"}
    mean_depth = "NA"
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 4 and fields[0] == "total_region":
            mean_depth = fields[3]
    return {"mean_depth": mean_depth}


def html_table(rows: list[dict[str, str]], cols: list[str | tuple[str, str]]) -> str:
    if not rows:
        return "<p class='text-muted mb-0'>暂无记录。</p>"
    keys = [c[0] if isinstance(c, tuple) else c for c in cols]
    labels = [c[1] if isinstance(c, tuple) else c for c in cols]
    head = "".join(f"<th>{html.escape(label)}</th>" for label in labels)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{html.escape(str(row.get(key, '')))}</td>" for key in keys) + "</tr>")
    return f"<div class='table-responsive'><table class='table table-sm table-striped align-middle'><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table></div>"


def main() -> int:
    global BCFTOOLS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--project", action="append", default=[], help="pair_id=/path/to/project")
    parser.add_argument("--bcftools", default=BCFTOOLS, help="bcftools executable path")
    args = parser.parse_args()
    BCFTOOLS = args.bcftools

    truth = read_tsv(args.truth)
    projects = {}
    for item in args.project:
        pair_id, path = item.split("=", 1)
        projects[pair_id] = Path(path)

    cards = []
    variant_chart = []
    depth_chart = []
    all_variant_rows = []
    all_neo_rows = []

    for pair_id, project in projects.items():
        result_dir = project / "results" / f"{pair_id}"
        variants_dir = result_dir / "variants"
        vcf = variants_dir / f"{pair_id}.mutect2.pass.vcf.gz"
        raw_vcf = variants_dir / f"{pair_id}.mutect2.raw.vcf.gz"
        vep_vcf = result_dir / "annotation" / f"{pair_id}.vep.vcf.gz"
        manual_vcf = result_dir / "annotation" / f"{pair_id}.vep.manual_filtered.vcf.gz"
        pass_count = count_vcf(vcf)
        raw_count = count_vcf(raw_vcf)
        vep_count = count_vcf(vep_vcf)
        manual_count = count_vcf(manual_vcf)
        variant_chart.append({"sample": pair_id, "raw": raw_count, "pass": pass_count, "vep": vep_count, "manual": manual_count})

        rows = query_vcf(vcf)
        truth_by_variant = truth_lookup(truth)
        for row in rows:
            row["pair_id"] = pair_id
            row["gene"] = annotate_gene(row, truth_by_variant)
            all_variant_rows.append(row)

        tumor_id = pair_id.rsplit("_vs_", 1)[0]
        normal_id = pair_id.rsplit("_vs_", 1)[1] if "_vs_" in pair_id else ""
        all_neo_rows.extend(
            read_neoantigen_manifest(
                result_dir / "neoantigen" / f"{pair_id}_neoantigen_manifest.tsv",
                pair_id,
            )
        )
        tumor_fastp = parse_fastp(project / "results" / tumor_id / "trimmed" / f"{tumor_id}_fastp.json")
        normal_fastp = parse_fastp(project / "results" / normal_id / "trimmed" / f"{normal_id}_fastp.json")
        tumor_flag = parse_flagstat(project / "results" / tumor_id / "aligned" / f"{tumor_id}.dedup.flagstat.txt")
        normal_flag = parse_flagstat(project / "results" / normal_id / "aligned" / f"{normal_id}.dedup.flagstat.txt")
        cov = parse_coverage(result_dir / "coverage" / f"{pair_id}.mosdepth.summary.txt")
        depth_chart.append({"sample": pair_id, "depth": cov.get("mean_depth", "0")})

        cards.append(
            f"""
            <div class="col-md-6">
              <div class="card h-100 shadow-sm">
                <div class="card-body">
                  <h5 class="card-title">{html.escape(pair_id)}</h5>
                  <p class="mb-1"><b>原始突变数:</b> {raw_count}</p>
                  <p class="mb-1"><b>PASS突变数:</b> {pass_count}</p>
                  <p class="mb-1"><b>VEP注释记录:</b> {vep_count}</p>
                  <p class="mb-1"><b>人工阈值过滤通过:</b> {manual_count}</p>
                  <p class="mb-1"><b>肿瘤比对率:</b> {html.escape(tumor_flag['mapped_rate'])}; <b>正常比对率:</b> {html.escape(normal_flag['mapped_rate'])}</p>
                  <p class="mb-0"><b>质控后reads:</b> 肿瘤 {tumor_fastp.get('after_reads', 'NA')}；正常 {normal_fastp.get('after_reads', 'NA')}</p>
                </div>
              </div>
            </div>
            """
        )

    truth_cols = [("pair_id", "配对ID"), ("tumor_sample", "肿瘤样本"), ("gene", "基因"), ("chrom", "染色体"), ("pos", "位置"), ("ref", "参考"), ("alt", "突变"), ("expected_af", "设计AF")]
    variant_cols = [("pair_id", "配对ID"), ("gene", "基因"), ("chrom", "染色体"), ("pos", "位置"), ("id", "ID"), ("ref", "参考"), ("alt", "突变"), ("filter", "过滤"), ("tumor_ad", "肿瘤AD"), ("tumor_af", "肿瘤AF"), ("tumor_dp", "肿瘤DP")]
    neo_cols = [("pair_id", "配对ID"), ("gene", "基因"), ("variant", "突变"), ("effect", "突变类型"), ("mer", "mer"), ("peptide_start", "起点"), ("peptide_end", "终点"), ("changed_aa", "氨基酸变化"), ("wt_peptide", "野生型肽段"), ("mt_peptide", "突变型肽段"), ("logic", "生成逻辑")]
    action_cols = [("gene", "基因"), ("alteration", "改变类型"), ("therapy", "潜在药物/方向"), ("evidence", "证据级别提示"), ("note", "说明")]
    actionability_rows = make_actionability_rows(all_variant_rows)

    payload = {
        "variantChart": variant_chart,
        "depthChart": depth_chart,
    }
    html_doc = f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>WES肿瘤-正常配对分析演示报告</title>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
  <script src="https://cdn.jsdelivr.net/npm/echarts@5.5.1/dist/echarts.min.js"></script>
  <style>
    body {{ background: #f6f8fb; }}
    .hero {{ background: #111827; color: white; padding: 28px 0; }}
    .card {{ border-radius: 8px; }}
    .chart {{ height: 320px; }}
    code {{ color: #0f766e; }}
  </style>
</head>
<body>
  <section class="hero">
    <div class="container">
      <h1 class="h3 mb-2">WES肿瘤-正常配对分析演示报告</h1>
      <p class="mb-0">覆盖 TP53 / BRCA1 / ERBB2(HER2) 的模拟靶区数据，用于展示从FASTQ到突变、新抗原候选窗口和用药提示的分析闭环。</p>
    </div>
  </section>
  <main class="container py-4">
    <div class="row g-3 mb-4">{''.join(cards)}</div>
    <div class="row g-3 mb-4">
      <div class="col-lg-7"><div class="card shadow-sm"><div class="card-body"><h2 class="h5">突变数量</h2><div id="variantChart" class="chart"></div></div></div></div>
      <div class="col-lg-5"><div class="card shadow-sm"><div class="card-body"><h2 class="h5">靶区平均深度</h2><div id="depthChart" class="chart"></div></div></div></div>
    </div>
    <div class="card shadow-sm mb-4"><div class="card-body"><h2 class="h5">设计突变真值</h2>{html_table(truth, truth_cols)}</div></div>
    <div class="card shadow-sm mb-4"><div class="card-body"><h2 class="h5">检出的PASS突变</h2>{html_table(all_variant_rows, variant_cols)}</div></div>
    <div class="card shadow-sm mb-4"><div class="card-body"><h2 class="h5">新抗原候选肽</h2><p class="text-muted small">演示逻辑：基于VEP蛋白改变和配置的蛋白FASTA重建突变蛋白，枚举所有覆盖突变位点的8-15mer候选肽。真实生产流程应结合转录本选择、表达量和HLA结合预测综合筛选。</p>{html_table(all_neo_rows, neo_cols)}</div></div>
    <div class="card shadow-sm"><div class="card-body"><h2 class="h5">药物/可干预性提示</h2><p class="text-muted small">仅用于演示，不构成临床用药建议。正式解读需要结合癌种、病理、胚系/体系来源、变异致病性、药品适应证、指南和OncoKB/CIViC/CGI等知识库。</p>{html_table(actionability_rows, action_cols)}</div></div>
  </main>
  <script>
    const payload = {json.dumps(payload, ensure_ascii=False)};
    const vc = echarts.init(document.getElementById('variantChart'));
    vc.setOption({{
      tooltip: {{ trigger: 'axis' }},
      legend: {{ data: ['原始', 'PASS', 'VEP', '人工过滤'] }},
      xAxis: {{ type: 'category', data: payload.variantChart.map(x => x.sample) }},
      yAxis: {{ type: 'value' }},
      series: [
        {{ name: '原始', type: 'bar', data: payload.variantChart.map(x => x.raw) }},
        {{ name: 'PASS', type: 'bar', data: payload.variantChart.map(x => x.pass) }},
        {{ name: 'VEP', type: 'bar', data: payload.variantChart.map(x => x.vep) }},
        {{ name: '人工过滤', type: 'bar', data: payload.variantChart.map(x => x.manual) }}
      ]
    }});
    const dc = echarts.init(document.getElementById('depthChart'));
    dc.setOption({{
      tooltip: {{ trigger: 'axis' }},
      xAxis: {{ type: 'category', data: payload.depthChart.map(x => x.sample) }},
      yAxis: {{ type: 'value', name: '深度' }},
      series: [{{ name: '平均深度', type: 'line', smooth: true, data: payload.depthChart.map(x => Number(x.depth) || 0) }}]
    }});
    window.addEventListener('resize', () => {{ vc.resize(); dc.resize(); }});
  </script>
</body>
</html>
"""
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html_doc)
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
