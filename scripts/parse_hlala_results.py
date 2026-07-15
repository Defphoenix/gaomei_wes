#!/usr/bin/env python3
"""Normalize HLA*LA best-guess output and derive class-I binding alleles."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def normalized_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def find_column(columns: list[str], candidates: tuple[str, ...]) -> str | None:
    keyed = {normalized_key(column): column for column in columns}
    for candidate in candidates:
        if normalized_key(candidate) in keyed:
            return keyed[normalized_key(candidate)]
    return None


def clean_allele(value: str, locus: str) -> str:
    allele = value.strip().replace("HLA-", "")
    allele = re.sub(r"[gGpPnN]$", "", allele)
    if allele and "*" not in allele and locus:
        allele = f"{locus}*{allele}"
    return allele


def binding_allele(allele: str, locus: str) -> str:
    if locus not in {"A", "B", "C"} or "*" not in allele:
        return ""
    name, fields_text = allele.split("*", 1)
    fields = fields_text.split(":")
    if len(fields) < 2:
        return ""
    return f"HLA-{name}*{fields[0]}:{fields[1]}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--alleles-output", required=True)
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    with input_path.open(encoding="utf-8") as handle:
        rows = [line.rstrip("\n") for line in handle if line.strip() and not line.startswith("#")]
    if not rows:
        raise SystemExit(f"HLA*LA result is empty: {input_path}")

    dialect = csv.excel_tab if "\t" in rows[0] else csv.excel
    reader = csv.DictReader(rows, dialect=dialect)
    columns = reader.fieldnames or []
    locus_column = find_column(columns, ("Locus", "Gene"))
    allele_column = find_column(columns, ("Allele", "AlleleFull", "HLAType"))
    chromosome_column = find_column(columns, ("Chromosome", "Copy", "Haplotype"))
    q1_column = find_column(columns, ("Q1",))
    q2_column = find_column(columns, ("Q2",))
    average_coverage_column = find_column(columns, ("AverageCoverage",))
    first_decile_column = find_column(columns, ("CoverageFirstDecile",))
    minimum_coverage_column = find_column(columns, ("MinimumCoverage",))
    kmers_column = find_column(columns, ("proportionkMersCovered",))
    perfect_g_column = find_column(columns, ("perfectG",))
    if not locus_column or not allele_column:
        raise SystemExit(f"Unsupported HLA*LA columns: {', '.join(columns)}")

    output_rows: list[dict[str, str]] = []
    binding: list[str] = []
    for index, row in enumerate(reader, start=1):
        locus = row.get(locus_column, "").strip().replace("HLA-", "")
        allele_raw = row.get(allele_column, "").strip()
        allele = clean_allele(allele_raw, locus)
        if not allele or allele in {"NA", "-"}:
            continue
        fields = allele.split("*", 1)[1].split(":") if "*" in allele else []
        bind = binding_allele(allele, locus)
        if bind and bind not in binding:
            binding.append(bind)
        output_rows.append(
            {
                "sample": args.sample,
                "locus": locus,
                "copy": row.get(chromosome_column, str(index)) if chromosome_column else str(index),
                "hla_la_call": allele_raw,
                "g_group": allele_raw,
                "full_allele": allele,
                "field_count": str(len(fields)),
                "legacy_resolution": f"{len(fields) * 2}-digit" if fields else "unknown",
                "resolution_note": "HLA*LA G-group call; not guaranteed full-genomic resolution",
                "binding_allele": bind,
                "quality_1": row.get(q1_column, "") if q1_column else "",
                "quality_2": row.get(q2_column, "") if q2_column else "",
                "average_coverage": row.get(average_coverage_column, "") if average_coverage_column else "",
                "coverage_first_decile": row.get(first_decile_column, "") if first_decile_column else "",
                "minimum_coverage": row.get(minimum_coverage_column, "") if minimum_coverage_column else "",
                "proportion_kmers_covered": row.get(kmers_column, "") if kmers_column else "",
                "perfect_g": row.get(perfect_g_column, "") if perfect_g_column else "",
            }
        )

    if not output_rows:
        raise SystemExit("No HLA alleles could be parsed from HLA*LA output")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(output_rows[0])
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    Path(args.alleles_output).write_text(",".join(binding) + "\n", encoding="utf-8")
    print(f"hla_calls={len(output_rows)}")
    print(f"binding_alleles={len(binding)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
