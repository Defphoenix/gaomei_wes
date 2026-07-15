#!/usr/bin/env python3
"""Strict, auditable TMB calculation from a VEP-annotated somatic VCF."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def parse_number(value: str | None) -> float | None:
    if value is None or value in {"", ".", "NA", "-"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def allele_value(value: str | None, alt_index: int) -> float | None:
    if not value:
        return None
    values = value.split(",")
    return parse_number(values[alt_index] if alt_index < len(values) else values[0])


def parse_info(text: str) -> dict[str, str]:
    info: dict[str, str] = {}
    for item in text.split(";"):
        key, sep, value = item.partition("=")
        info[key] = value if sep else "true"
    return info


def parse_format(format_text: str, sample_text: str) -> dict[str, str]:
    keys = format_text.split(":")
    values = sample_text.split(":")
    return dict(zip(keys, values))


def sample_metrics(sample: dict[str, str], alt_index: int) -> tuple[int | None, int | None, float | None]:
    dp_value = parse_number(sample.get("DP"))
    ad = sample.get("AD", "").split(",")
    alt_reads = parse_number(ad[alt_index + 1]) if len(ad) > alt_index + 1 else None
    af = allele_value(sample.get("AF"), alt_index)
    if af is None and dp_value and alt_reads is not None:
        af = alt_reads / dp_value
    return (
        int(dp_value) if dp_value is not None else None,
        int(alt_reads) if alt_reads is not None else None,
        af,
    )


def load_bed(path: Path) -> tuple[dict[str, list[tuple[int, int]]], int]:
    raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Invalid BED line: {line.rstrip()}")
            start, end = int(fields[1]), int(fields[2])
            if end > start:
                raw[fields[0]].append((start, end))

    merged: dict[str, list[tuple[int, int]]] = {}
    total = 0
    for chrom, intervals in raw.items():
        current: list[list[int]] = []
        for start, end in sorted(intervals):
            if not current or start > current[-1][1]:
                current.append([start, end])
            else:
                current[-1][1] = max(current[-1][1], end)
        merged[chrom] = [(start, end) for start, end in current]
        total += sum(end - start for start, end in merged[chrom])
    return merged, total


def in_bed(intervals: dict[str, list[tuple[int, int]]], chrom: str, pos: int) -> bool:
    zero_based = pos - 1
    return any(start <= zero_based < end for start, end in intervals.get(chrom, []))


def select_annotation(annotations: list[dict[str, str]], allowed: set[str]) -> tuple[dict[str, str] | None, set[str]]:
    candidates: list[tuple[int, dict[str, str], set[str]]] = []
    for annotation in annotations:
        consequences = set(annotation.get("Consequence", "").split("&"))
        matched = consequences & allowed
        if not matched:
            continue
        priority = 2
        if annotation.get("MANE_SELECT") not in {None, "", "-"}:
            priority = 0
        elif annotation.get("CANONICAL", "").upper() == "YES":
            priority = 1
        candidates.append((priority, annotation, matched))
    if not candidates:
        return None, set()
    candidates.sort(key=lambda item: item[0])
    return candidates[0][1], candidates[0][2]


def max_population_af(annotation: dict[str, str], fields: list[str]) -> float | None:
    values: list[float] = []
    for field in fields:
        for raw in annotation.get(field, "").split("&"):
            value = parse_number(raw)
            if value is not None:
                values.append(value)
    return max(values) if values else None


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--accepted", required=True)
    parser.add_argument("--rejected", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--tumor-sample", required=True)
    parser.add_argument("--normal-sample", default="")
    parser.add_argument("--effective-coding-bed", default="")
    parser.add_argument("--denominator-mb", type=float, default=0.0)
    parser.add_argument("--denominator-validated", choices=("true", "false"), default="false")
    parser.add_argument("--min-qual", type=float, default=0.0)
    parser.add_argument("--min-tlod", type=float, default=6.3)
    parser.add_argument("--min-tumor-dp", type=int, default=20)
    parser.add_argument("--min-tumor-alt-reads", type=int, default=5)
    parser.add_argument("--min-tumor-af", type=float, default=0.05)
    parser.add_argument("--min-normal-dp", type=int, default=10)
    parser.add_argument("--max-normal-alt-reads", type=int, default=2)
    parser.add_argument("--max-normal-af", type=float, default=0.02)
    parser.add_argument("--max-population-af", type=float, default=0.001)
    parser.add_argument("--consequences", required=True)
    parser.add_argument("--population-fields", default="MAX_AF,gnomADe_AF,gnomADg_AF,AF")
    args = parser.parse_args()

    vcf_path = Path(args.vcf)
    csq_fields: list[str] = []
    sample_names: list[str] = []
    records: list[list[str]] = []
    with open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                match = re.search(r"Format: ([^\"]+)", line)
                if match:
                    csq_fields = match.group(1).rstrip(">\n").split("|")
            elif line.startswith("#CHROM"):
                sample_names = line.rstrip("\n").split("\t")[9:]
            elif not line.startswith("#") and line.strip():
                records.append(line.rstrip("\n").split("\t"))

    if not csq_fields:
        raise SystemExit("Input VCF does not contain a parseable VEP CSQ header")
    if args.tumor_sample not in sample_names:
        raise SystemExit(f"Tumor sample '{args.tumor_sample}' not found; VCF samples: {','.join(sample_names)}")
    if args.normal_sample and args.normal_sample not in sample_names:
        raise SystemExit(f"Normal sample '{args.normal_sample}' not found; VCF samples: {','.join(sample_names)}")
    tumor_index = sample_names.index(args.tumor_sample)
    normal_index = sample_names.index(args.normal_sample) if args.normal_sample in sample_names else None

    bed_intervals: dict[str, list[tuple[int, int]]] = {}
    denominator_source = "configured_mb"
    denominator_bp = int(args.denominator_mb * 1_000_000)
    if args.effective_coding_bed:
        bed_intervals, denominator_bp = load_bed(Path(args.effective_coding_bed))
        denominator_source = str(Path(args.effective_coding_bed))
    if denominator_bp <= 0:
        raise SystemExit("TMB denominator is zero; configure an effective coding BED or positive denominator Mb")

    allowed = {item for item in args.consequences.split(",") if item}
    population_fields = [item for item in args.population_fields.split(",") if item]
    accepted: list[dict[str, object]] = []
    rejected: list[dict[str, object]] = []
    seen: set[tuple[str, int, str, str]] = set()

    for fields in records:
        if len(fields) < 8:
            continue
        chrom, pos_text, _, ref, alt_text, qual_text, filter_text, info_text = fields[:8]
        pos = int(pos_text)
        info = parse_info(info_text)
        format_text = fields[8] if len(fields) > 8 else ""
        samples = fields[9:]
        tumor = parse_format(format_text, samples[tumor_index])
        normal = parse_format(format_text, samples[normal_index]) if normal_index is not None else {}
        csq_entries = []
        for entry in info.get("CSQ", "").split(","):
            values = entry.split("|")
            values.extend([""] * (len(csq_fields) - len(values)))
            csq_entries.append(dict(zip(csq_fields, values)))

        for alt_index, alt in enumerate(alt_text.split(",")):
            key = (chrom, pos, ref, alt)
            if key in seen:
                continue
            seen.add(key)
            reasons: list[str] = []
            if filter_text not in {"PASS", "."}:
                reasons.append("FILTER_not_PASS")
            if bed_intervals and not in_bed(bed_intervals, chrom, pos):
                reasons.append("outside_effective_coding_bed")
            qual = parse_number(qual_text)
            if args.min_qual > 0 and (qual is None or qual < args.min_qual):
                reasons.append("low_or_missing_QUAL")
            tlod = allele_value(info.get("TLOD"), alt_index)
            if args.min_tlod > 0 and (tlod is None or tlod < args.min_tlod):
                reasons.append("low_or_missing_TLOD")

            tumor_dp, tumor_alt, tumor_af = sample_metrics(tumor, alt_index)
            if tumor_dp is None or tumor_dp < args.min_tumor_dp:
                reasons.append("low_tumor_DP")
            if tumor_alt is None or tumor_alt < args.min_tumor_alt_reads:
                reasons.append("low_tumor_ALT_reads")
            if tumor_af is None or tumor_af < args.min_tumor_af:
                reasons.append("low_tumor_AF")

            normal_dp = normal_alt = None
            normal_af = None
            if normal_index is not None:
                normal_dp, normal_alt, normal_af = sample_metrics(normal, alt_index)
                if normal_dp is None or normal_dp < args.min_normal_dp:
                    reasons.append("low_normal_DP")
                if normal_alt is None or normal_alt > args.max_normal_alt_reads:
                    reasons.append("high_normal_ALT_reads")
                if normal_af is None or normal_af > args.max_normal_af:
                    reasons.append("high_normal_AF")

            allele_annotations = [item for item in csq_entries if item.get("Allele") == alt]
            if not allele_annotations and len(alt_text.split(",")) == 1:
                allele_annotations = csq_entries
            annotation, matched = select_annotation(allele_annotations, allowed)
            if annotation is None:
                reasons.append("non_eligible_VEP_consequence")
                annotation = allele_annotations[0] if allele_annotations else {}
            pop_af = max_population_af(annotation, population_fields)
            if pop_af is not None and pop_af > args.max_population_af:
                reasons.append("population_AF_above_threshold")

            row: dict[str, object] = {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene": annotation.get("SYMBOL", ""),
                "consequence": "&".join(sorted(matched)) if matched else annotation.get("Consequence", ""),
                "feature": annotation.get("Feature", ""),
                "hgvsc": annotation.get("HGVSc", ""),
                "hgvsp": annotation.get("HGVSp", ""),
                "qual": qual_text,
                "tlod": "" if tlod is None else tlod,
                "tumor_dp": "" if tumor_dp is None else tumor_dp,
                "tumor_alt_reads": "" if tumor_alt is None else tumor_alt,
                "tumor_af": "" if tumor_af is None else f"{tumor_af:.6g}",
                "normal_dp": "" if normal_dp is None else normal_dp,
                "normal_alt_reads": "" if normal_alt is None else normal_alt,
                "normal_af": "" if normal_af is None else f"{normal_af:.6g}",
                "max_population_af": "" if pop_af is None else f"{pop_af:.6g}",
                "status": "accepted" if not reasons else "rejected",
                "reasons": ";".join(reasons),
            }
            (accepted if not reasons else rejected).append(row)

    denominator_mb = denominator_bp / 1_000_000
    tmb = len(accepted) / denominator_mb
    summary = {
        "input_vcf": str(vcf_path),
        "tumor_sample": args.tumor_sample,
        "normal_sample": args.normal_sample if normal_index is not None else "",
        "total_vcf_records": len(records),
        "accepted_variants": len(accepted),
        "rejected_variants": len(rejected),
        "denominator_bp": denominator_bp,
        "denominator_mb": denominator_mb,
        "denominator_source": denominator_source,
        "denominator_validated": args.denominator_validated == "true",
        "tmb_mutations_per_mb": tmb,
        "consequences": sorted(allowed),
        "thresholds": vars(args),
    }
    output_fields = [
        "chrom", "pos", "ref", "alt", "gene", "consequence", "feature", "hgvsc", "hgvsp",
        "qual", "tlod", "tumor_dp", "tumor_alt_reads", "tumor_af", "normal_dp",
        "normal_alt_reads", "normal_af", "max_population_af", "status", "reasons",
    ]
    write_tsv(Path(args.accepted), accepted, output_fields)
    write_tsv(Path(args.rejected), rejected, output_fields)
    Path(args.summary_json).write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with Path(args.summary_tsv).open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key in ("total_vcf_records", "accepted_variants", "rejected_variants", "denominator_bp", "denominator_mb", "denominator_source", "denominator_validated", "tmb_mutations_per_mb"):
            handle.write(f"{key}\t{summary[key]}\n")
    print(f"accepted_variants={len(accepted)}")
    print(f"denominator_mb={denominator_mb:.6f}")
    print(f"tmb_mutations_per_mb={tmb:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
