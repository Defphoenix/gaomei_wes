#!/usr/bin/env python3
"""Apply an auditable, configurable somatic evidence filter to a VEP VCF."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from pathlib import Path
from typing import TextIO


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def parse_number(value: str | None) -> float | None:
    if value is None or value in {"", ".", "NA", "-"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_info(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.split(";"):
        key, sep, value = item.partition("=")
        result[key] = value if sep else "true"
    return result


def parse_format(format_text: str, sample_text: str) -> dict[str, str]:
    return dict(zip(format_text.split(":"), sample_text.split(":")))


def allele_value(value: str | None, alt_index: int) -> float | None:
    if not value:
        return None
    values = value.split(",")
    selected = values[alt_index] if alt_index < len(values) else values[0]
    return parse_number(selected)


def sample_metrics(sample: dict[str, str], alt_index: int) -> tuple[int | None, int | None, float | None]:
    dp_number = parse_number(sample.get("DP"))
    ad_values = sample.get("AD", "").split(",")
    alt_number = parse_number(ad_values[alt_index + 1]) if len(ad_values) > alt_index + 1 else None
    af = allele_value(sample.get("AF"), alt_index)
    if af is None and dp_number and alt_number is not None:
        af = alt_number / dp_number
    return (
        int(dp_number) if dp_number is not None else None,
        int(alt_number) if alt_number is not None else None,
        af,
    )


def vep_allele(ref: str, alt: str) -> str:
    """Return VEP's minimal allele representation for a biallelic variant."""
    ref_work, alt_work = ref, alt
    while ref_work and alt_work and ref_work[0] == alt_work[0]:
        ref_work, alt_work = ref_work[1:], alt_work[1:]
    while ref_work and alt_work and ref_work[-1] == alt_work[-1]:
        ref_work, alt_work = ref_work[:-1], alt_work[:-1]
    return alt_work or "-"


def parse_csq(info: dict[str, str], csq_fields: list[str]) -> list[dict[str, str]]:
    annotations: list[dict[str, str]] = []
    for entry in info.get("CSQ", "").split(","):
        if not entry:
            continue
        values = entry.split("|")
        values.extend([""] * (len(csq_fields) - len(values)))
        annotations.append(dict(zip(csq_fields, values)))
    return annotations


def annotations_for_alt(
    annotations: list[dict[str, str]], ref: str, alt: str
) -> list[dict[str, str]]:
    aliases = {alt, vep_allele(ref, alt)}
    matched = [item for item in annotations if item.get("Allele") in aliases]
    return matched or annotations


def max_population_af(annotations: list[dict[str, str]], fields: list[str]) -> float | None:
    values: list[float] = []
    for annotation in annotations:
        for field in fields:
            for raw in re.split(r"[&,]", annotation.get(field, "")):
                value = parse_number(raw)
                if value is not None:
                    values.append(value)
    return max(values) if values else None


def display_annotation(annotations: list[dict[str, str]]) -> dict[str, str]:
    if not annotations:
        return {}
    return sorted(
        annotations,
        key=lambda item: (
            item.get("MANE_SELECT", "") in {"", "-"},
            item.get("CANONICAL", "").upper() != "YES",
        ),
    )[0]


def criterion(value: float | int | None, threshold: float | int, direction: str) -> bool:
    if value is None:
        return False
    return value >= threshold if direction == "min" else value <= threshold


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--audit-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--tumor-sample", required=True)
    parser.add_argument("--normal-sample", default="")
    parser.add_argument("--min-tlod", type=float, default=6.3)
    parser.add_argument("--min-tumor-dp", type=int, default=20)
    parser.add_argument("--min-tumor-alt-reads", type=int, default=5)
    parser.add_argument("--min-tumor-af", type=float, default=0.02)
    parser.add_argument("--min-normal-dp", type=int, default=20)
    parser.add_argument("--max-normal-alt-reads", type=int, default=2)
    parser.add_argument("--max-normal-af", type=float, default=0.02)
    parser.add_argument("--max-population-af", type=float, default=0.001)
    parser.add_argument("--population-fields", default="MAX_AF,gnomADe_AF,gnomADg_AF,AF")
    args = parser.parse_args()

    input_path = Path(args.input_vcf)
    output_path = Path(args.output_vcf)
    audit_path = Path(args.audit_tsv)
    summary_path = Path(args.summary_json)
    population_fields = [value for value in args.population_fields.split(",") if value]

    csq_fields: list[str] = []
    sample_names: list[str] = []
    header_lines: list[str] = []
    record_lines: list[str] = []
    with open_text(input_path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                match = re.search(r"Format: ([^\"]+)", line)
                if match:
                    csq_fields = match.group(1).rstrip(">\n").split("|")
            if line.startswith("#"):
                header_lines.append(line)
                if line.startswith("#CHROM"):
                    sample_names = line.rstrip("\n").split("\t")[9:]
            elif line.strip():
                record_lines.append(line)

    if record_lines and not csq_fields:
        raise SystemExit("Input VCF does not contain a parseable VEP CSQ header")
    if args.tumor_sample not in sample_names:
        raise SystemExit(f"Tumor sample '{args.tumor_sample}' not found; VCF samples: {','.join(sample_names)}")
    if args.normal_sample and args.normal_sample not in sample_names:
        raise SystemExit(f"Normal sample '{args.normal_sample}' not found; VCF samples: {','.join(sample_names)}")

    tumor_index = sample_names.index(args.tumor_sample)
    normal_index = sample_names.index(args.normal_sample) if args.normal_sample else None
    rows: list[dict[str, object]] = []
    passed_lines: list[str] = []
    reason_counts: dict[str, int] = {}

    for line in record_lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue
        chrom, pos, _, ref, alt_text, _, filter_text, info_text, format_text = fields[:9]
        alts = alt_text.split(",")
        if len(alts) != 1:
            raise SystemExit(
                f"Manual filter requires biallelic records; normalize the VEP VCF first: {chrom}:{pos}"
            )
        alt = alts[0]
        info = parse_info(info_text)
        sample_values = fields[9:]
        tumor = parse_format(format_text, sample_values[tumor_index])
        normal = parse_format(format_text, sample_values[normal_index]) if normal_index is not None else {}
        tumor_dp, tumor_alt, tumor_af = sample_metrics(tumor, 0)
        normal_dp, normal_alt, normal_af = sample_metrics(normal, 0) if normal_index is not None else (None, None, None)
        tlod = allele_value(info.get("TLOD"), 0)
        annotations = annotations_for_alt(parse_csq(info, csq_fields), ref, alt)
        pop_af = max_population_af(annotations, population_fields)
        annotation = display_annotation(annotations)

        checks = {
            "pass_mutect_filter": filter_text in {"PASS", "."},
            "pass_tlod": criterion(tlod, args.min_tlod, "min"),
            "pass_tumor_dp": criterion(tumor_dp, args.min_tumor_dp, "min"),
            "pass_tumor_alt_reads": criterion(tumor_alt, args.min_tumor_alt_reads, "min"),
            "pass_tumor_af": criterion(tumor_af, args.min_tumor_af, "min"),
            "pass_population_af": pop_af is None or pop_af <= args.max_population_af,
        }
        if normal_index is not None:
            checks.update(
                {
                    "pass_normal_dp": criterion(normal_dp, args.min_normal_dp, "min"),
                    "pass_normal_alt_reads": criterion(normal_alt, args.max_normal_alt_reads, "max"),
                    "pass_normal_af": criterion(normal_af, args.max_normal_af, "max"),
                }
            )
        else:
            checks.update(
                {
                    "pass_normal_dp": True,
                    "pass_normal_alt_reads": True,
                    "pass_normal_af": True,
                }
            )

        reasons = [name[5:] if name.startswith("pass_") else name for name, passed in checks.items() if not passed]
        for reason in reasons:
            reason_counts[reason] = reason_counts.get(reason, 0) + 1
        passed = not reasons
        if passed:
            passed_lines.append(line)

        rows.append(
            {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene": annotation.get("SYMBOL", ""),
                "consequence": annotation.get("Consequence", ""),
                "hgvsc": annotation.get("HGVSc", ""),
                "hgvsp": annotation.get("HGVSp", ""),
                "tlod": "" if tlod is None else tlod,
                "tumor_dp": "" if tumor_dp is None else tumor_dp,
                "tumor_alt_reads": "" if tumor_alt is None else tumor_alt,
                "tumor_af": "" if tumor_af is None else f"{tumor_af:.6g}",
                "normal_dp": "" if normal_dp is None else normal_dp,
                "normal_alt_reads": "" if normal_alt is None else normal_alt,
                "normal_af": "" if normal_af is None else f"{normal_af:.6g}",
                "max_population_af": "" if pop_af is None else f"{pop_af:.6g}",
                **checks,
                "overall_pass": passed,
                "failure_reasons": ";".join(reasons),
            }
        )

    threshold_header = (
        "##gaomei_wes_manual_filter="
        f"TLOD>={args.min_tlod},tumor_DP>={args.min_tumor_dp},"
        f"tumor_ALT>={args.min_tumor_alt_reads},tumor_AF>={args.min_tumor_af},"
        f"normal_DP>={args.min_normal_dp},normal_ALT<={args.max_normal_alt_reads},"
        f"normal_AF<={args.max_normal_af},population_AF<={args.max_population_af};"
        "missing_population_AF=PASS\n"
    )
    with output_path.open("w", encoding="utf-8") as handle:
        inserted = False
        for line in header_lines:
            if line.startswith("#CHROM") and not inserted:
                handle.write(threshold_header)
                inserted = True
            handle.write(line)
        handle.writelines(passed_lines)

    fieldnames = [
        "chrom", "pos", "ref", "alt", "gene", "consequence", "hgvsc", "hgvsp",
        "tlod", "tumor_dp", "tumor_alt_reads", "tumor_af", "normal_dp",
        "normal_alt_reads", "normal_af", "max_population_af", "pass_mutect_filter",
        "pass_tlod", "pass_tumor_dp", "pass_tumor_alt_reads", "pass_tumor_af",
        "pass_normal_dp", "pass_normal_alt_reads", "pass_normal_af",
        "pass_population_af", "overall_pass", "failure_reasons",
    ]
    with audit_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "input_vcf": str(input_path),
        "output_vcf": str(output_path),
        "tumor_sample": args.tumor_sample,
        "normal_sample": args.normal_sample,
        "input_variants": len(rows),
        "passed_variants": len(passed_lines),
        "failed_variants": len(rows) - len(passed_lines),
        "failure_reason_counts": reason_counts,
        "thresholds": {
            "min_tlod": args.min_tlod,
            "min_tumor_dp": args.min_tumor_dp,
            "min_tumor_alt_reads": args.min_tumor_alt_reads,
            "min_tumor_af": args.min_tumor_af,
            "min_normal_dp": args.min_normal_dp,
            "max_normal_alt_reads": args.max_normal_alt_reads,
            "max_normal_af": args.max_normal_af,
            "max_population_af": args.max_population_af,
            "missing_population_af_passes": True,
            "population_fields": population_fields,
        },
    }
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"input_variants={len(rows)}")
    print(f"passed_variants={len(passed_lines)}")
    print(f"failed_variants={len(rows) - len(passed_lines)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
