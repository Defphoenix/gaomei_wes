#!/usr/bin/env python3
"""Create a normalized MSI call table from MSIsensor-pro outputs.

MSIsensor-pro output names and table layouts differ by subcommand/version.  This
parser is intentionally permissive: it scans all text-like files with the given
prefix, extracts total/unstable site counts or an MSI score when present, and
writes one stable TSV for downstream reports.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable, Optional


NUM_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)")


def normalize_key(raw: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", raw.lower()).strip("_")


def first_number(text: str) -> Optional[float]:
    match = NUM_RE.search(text)
    if not match:
        return None
    return float(match.group(0))


def text_files(prefix: Path) -> list[Path]:
    parent = prefix.parent
    stem = prefix.name
    files = sorted(parent.glob(f"{stem}*"))
    # Do not feed our own normalized outputs back into a later parse. Otherwise
    # stale summaries can be mistaken for raw MSIsensor metrics (for example a
    # sample ID such as SH05677 being parsed as an MSI score of 5677).
    generated = {
        f"{stem}_call.tsv",
        f"{stem}_call_summary.txt",
        f"{stem}_summary.txt",
    }
    return [
        path
        for path in files
        if path.is_file()
        and path.name not in generated
        and path.stat().st_size < 20_000_000
    ]


def safe_read_lines(path: Path) -> list[str]:
    try:
        return path.read_text(errors="ignore").splitlines()
    except Exception:
        return []


def parse_key_value_line(line: str) -> tuple[str, Optional[float]]:
    if ":" in line:
        key, value = line.split(":", 1)
        return normalize_key(key), first_number(value)
    fields = re.split(r"\s+", line.strip())
    if len(fields) >= 2 and not NUM_RE.fullmatch(fields[0]):
        return normalize_key(fields[0]), first_number(" ".join(fields[1:]))
    return "", None


def assign_metric(metrics: dict[str, float], key: str, value: Optional[float]) -> None:
    if value is None or not key:
        return
    if key in {"total_number_of_sites", "total_sites", "totalsites", "total"}:
        metrics.setdefault("total_sites", value)
    elif key in {
        "number_of_somatic_sites",
        "number_of_unstable_sites",
        "somatic_sites",
        "unstable_sites",
        "unstable",
        "somatic",
    }:
        metrics.setdefault("unstable_sites", value)
    elif key in {
        "msi",
        "msi_score",
        "msi_score_percent",
        "score",
        "percentage",
        "percent",
        "msi_percent",
    }:
        metrics.setdefault("msi_score", value)


def parse_header_table(lines: Iterable[str], metrics: dict[str, float]) -> None:
    rows = [line.strip() for line in lines if line.strip() and not line.startswith("#")]
    for idx, line in enumerate(rows[:-1]):
        header = re.split(r"\t+|\s{2,}", line)
        values = re.split(r"\t+|\s{2,}", rows[idx + 1])
        if len(header) < 2 or len(header) != len(values):
            continue
        row = {normalize_key(k): v for k, v in zip(header, values)}
        for key, value in row.items():
            assign_metric(metrics, key, first_number(value))


def parse_metrics(files: list[Path]) -> dict[str, float | str]:
    metrics: dict[str, float | str] = {}
    for path in files:
        lines = safe_read_lines(path)
        parse_header_table(lines, metrics)  # type: ignore[arg-type]
        for line in lines:
            key, value = parse_key_value_line(line)
            assign_metric(metrics, key, value)  # type: ignore[arg-type]
    total = metrics.get("total_sites")
    unstable = metrics.get("unstable_sites")
    if "msi_score" not in metrics and isinstance(total, (int, float)) and total > 0 and isinstance(unstable, (int, float)):
        metrics["msi_score"] = unstable / total * 100
    score = metrics.get("msi_score")
    if isinstance(score, (int, float)) and not 0 <= score <= 100:
        metrics.pop("msi_score", None)
    return metrics


def classify(score: Optional[float], low: float, high: float) -> tuple[str, str]:
    if score is None:
        return "NOT_DETERMINED", "no_parseable_msi_score"
    if score >= high:
        return "MSI-H", f"score>={high:g}"
    if score >= low:
        return "MSI-L", f"{low:g}<=score<{high:g}"
    return "MSS", f"score<{low:g}"


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--mode", required=True)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--low-threshold", type=float, default=3.0)
    parser.add_argument("--high-threshold", type=float, default=20.0)
    args = parser.parse_args()

    files = text_files(args.prefix)
    metrics = parse_metrics(files)
    score = metrics.get("msi_score")
    status, reason = classify(score if isinstance(score, (int, float)) else None, args.low_threshold, args.high_threshold)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    result_files = ",".join(str(path) for path in files)
    args.out.write_text(
        "\t".join(
            [
                "sample",
                "mode",
                "msi_status",
                "msi_score_percent",
                "unstable_sites",
                "total_sites",
                "low_threshold_percent",
                "high_threshold_percent",
                "reason",
                "result_files",
            ]
        )
        + "\n"
        + "\t".join(
            [
                args.sample,
                args.mode,
                status,
                fmt(score),
                fmt(metrics.get("unstable_sites")),
                fmt(metrics.get("total_sites")),
                fmt(args.low_threshold),
                fmt(args.high_threshold),
                reason,
                result_files,
            ]
        )
        + "\n"
    )

    args.summary.write_text(
        "\n".join(
            [
                f"MSI判定摘要 - {args.sample}",
                "================================",
                f"模式: {args.mode}",
                f"MSI状态: {status}",
                f"MSI score(%): {fmt(score) or 'NA'}",
                f"不稳定位点数: {fmt(metrics.get('unstable_sites')) or 'NA'}",
                f"总位点数: {fmt(metrics.get('total_sites')) or 'NA'}",
                f"阈值: MSS < {args.low_threshold:g}%, MSI-L {args.low_threshold:g}-{args.high_threshold:g}%, MSI-H >= {args.high_threshold:g}%",
                f"原因: {reason}",
                "",
                "结果文件:",
                *(f"  {path}" for path in files),
                "",
                "说明: WES MSI 判定依赖微卫星位点列表、覆盖度和捕获区域；正式报告建议结合癌种、样本质量和验证集阈值校准。",
            ]
        )
        + "\n"
    )
    print(f"msi_status={status}")
    print(f"msi_call={args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
