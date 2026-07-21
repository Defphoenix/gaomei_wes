#!/usr/bin/env python3
"""Summarize read-level evidence at selected SNV loci from paired BAM files."""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def query_index_at(cigar: str, alignment_start: int, target_pos: int) -> int | None:
    ref_pos = alignment_start
    query_pos = 0
    for length_text, operation in CIGAR_RE.findall(cigar):
        length = int(length_text)
        if operation in "M=X":
            if ref_pos <= target_pos < ref_pos + length:
                return query_pos + target_pos - ref_pos
            ref_pos += length
            query_pos += length
        elif operation in "DN":
            if ref_pos <= target_pos < ref_pos + length:
                return None
            ref_pos += length
        elif operation in "IS":
            query_pos += length
    return None


def orientation(flag: int) -> str:
    if flag & 0x40:
        return "F2R1" if flag & 0x10 else "F1R2"
    if flag & 0x80:
        return "F1R2" if flag & 0x10 else "F2R1"
    return "unpaired"


def read_observations(samtools: str, bam: Path, chrom: str, pos: int) -> list[dict[str, Any]]:
    command = [samtools, "view", str(bam), f"{chrom}:{pos}-{pos}"]
    process = subprocess.run(command, check=True, text=True, capture_output=True)
    observations: list[dict[str, Any]] = []
    for line in process.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        flag = int(fields[1])
        if flag & (0x4 | 0x100 | 0x200 | 0x800):
            continue
        query_index = query_index_at(fields[5], int(fields[3]), pos)
        if query_index is None or query_index >= len(fields[9]):
            continue
        base = fields[9][query_index].upper()
        base_quality = ord(fields[10][query_index]) - 33 if fields[10] != "*" else 0
        observations.append({
            "qname": fields[0], "flag": flag, "base": base, "bq": base_quality,
            "mapq": int(fields[4]), "cigar": fields[5], "reverse": bool(flag & 0x10),
            "read_number": 1 if flag & 0x40 else (2 if flag & 0x80 else 0),
            "duplicate": bool(flag & 0x400), "proper_pair": bool(flag & 0x2),
            "orientation": orientation(flag), "query_pos": query_index + 1,
            "distance_to_end": min(query_index, len(fields[9]) - query_index - 1),
            "read_length": len(fields[9]), "soft_clipped": "S" in fields[5],
            "indel_cigar": bool(re.search(r"[ID]", fields[5])),
            "alignment_start": int(fields[3]),
        })
    return observations


def summarize(observations: list[dict[str, Any]], ref: str, alt: str) -> dict[str, Any]:
    usable = [row for row in observations if not row["duplicate"]]
    high_quality = [row for row in usable if row["mapq"] >= 30 and row["bq"] >= 20]
    raw_alt_rows = [row for row in observations if row["base"] == alt]
    raw_ref_rows = [row for row in observations if row["base"] == ref]
    fragments: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in high_quality:
        fragments[row["qname"]].append(row)

    fragment_bases: list[str] = []
    alt_fragment_rows: list[dict[str, Any]] = []
    for rows in fragments.values():
        bases = {row["base"] for row in rows}
        if len(bases) == 1:
            base = next(iter(bases))
            fragment_bases.append(base)
            if base == alt:
                alt_fragment_rows.append(max(rows, key=lambda row: (row["bq"], row["mapq"])))

    raw_counts = Counter(row["base"] for row in observations)
    hq_counts = Counter(row["base"] for row in high_quality)
    fragment_counts = Counter(fragment_bases)
    alt_rows = [row for row in high_quality if row["base"] == alt]
    orientation_counts = Counter(row["orientation"] for row in alt_fragment_rows)
    strand_counts = Counter("reverse" if row["reverse"] else "forward" for row in alt_fragment_rows)
    read_number_counts = Counter(f"R{row['read_number']}" for row in alt_fragment_rows)
    fragment_depth = len(fragment_bases)
    alt_fragments = fragment_counts[alt]
    downsample_groups: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in usable:
        if row["mapq"] >= 30 and row["bq"] >= 10:
            downsample_groups[row["alignment_start"]].append(row)
    cap5_rows = [row for rows in downsample_groups.values() for row in rows[:5]]
    return {
        "raw_read_depth": len(observations),
        "raw_ref_reads": raw_counts[ref], "raw_alt_reads": raw_counts[alt],
        "raw_alt_mapq0_reads": sum(row["mapq"] == 0 for row in raw_alt_rows),
        "raw_alt_low_mapq_reads": sum(row["mapq"] < 30 for row in raw_alt_rows),
        "raw_alt_low_bq_reads": sum(row["bq"] < 20 for row in raw_alt_rows),
        "raw_alt_softclipped_reads": sum(row["soft_clipped"] for row in raw_alt_rows),
        "raw_alt_indel_cigar_reads": sum(row["indel_cigar"] for row in raw_alt_rows),
        "raw_alt_complex_cigar_reads": sum(
            row["soft_clipped"] or row["indel_cigar"] for row in raw_alt_rows
        ),
        "raw_ref_softclipped_reads": sum(row["soft_clipped"] for row in raw_ref_rows),
        "nonduplicate_reads": len(usable), "hq_read_depth": len(high_quality),
        "hq_ref_reads": hq_counts[ref], "hq_alt_reads": hq_counts[alt],
        "hq_other_reads": len(high_quality) - hq_counts[ref] - hq_counts[alt],
        "hq_fragment_depth": fragment_depth, "hq_ref_fragments": fragment_counts[ref],
        "hq_alt_fragments": alt_fragments,
        "hq_alt_vaf": round(alt_fragments / fragment_depth, 5) if fragment_depth else 0.0,
        "alt_forward_fragments": strand_counts["forward"],
        "alt_reverse_fragments": strand_counts["reverse"],
        "alt_f1r2_fragments": orientation_counts["F1R2"],
        "alt_f2r1_fragments": orientation_counts["F2R1"],
        "alt_read1_fragments": read_number_counts["R1"],
        "alt_read2_fragments": read_number_counts["R2"],
        "alt_mean_bq": round(statistics.mean(row["bq"] for row in alt_rows), 2) if alt_rows else 0,
        "alt_mean_mapq": round(statistics.mean(row["mapq"] for row in alt_rows), 2) if alt_rows else 0,
        "alt_median_end_distance": round(statistics.median(row["distance_to_end"] for row in alt_rows), 2) if alt_rows else 0,
        "alt_softclipped_reads": sum(row["soft_clipped"] for row in alt_rows),
        "alt_unique_alignment_starts": len({row["alignment_start"] for row in alt_rows}),
        "duplicate_reads": sum(row["duplicate"] for row in observations),
        "cap5_retained_reads_estimate": len(cap5_rows),
        "cap5_alt_reads_estimate": sum(row["base"] == alt for row in cap5_rows),
    }


def preliminary_reason(tumor: dict[str, Any], normal: dict[str, Any]) -> str:
    alt = tumor["hq_alt_fragments"]
    depth = tumor["hq_fragment_depth"]
    vaf = tumor["hq_alt_vaf"]
    reasons: list[str] = []
    if depth == 0:
        if tumor["raw_read_depth"] and tumor["raw_alt_reads"] and tumor["raw_alt_low_mapq_reads"] == tumor["raw_alt_reads"]:
            return "覆盖与ALT均存在，但ALT全部MAPQ<30，被Mutect2映射质量过滤"
        return "肿瘤无Mutect2可用的高质量覆盖"
    if alt == 0:
        reasons.append("肿瘤无高质量ALT片段")
    elif alt < 3:
        reasons.append(f"肿瘤仅{alt}个高质量ALT片段")
    if 0 < vaf < 0.02:
        reasons.append(f"肿瘤ALT VAF低({vaf:.2%})")
    if alt >= 3 and min(tumor["alt_f1r2_fragments"], tumor["alt_f2r1_fragments"]) == 0:
        reasons.append("ALT仅见于一种paired-read orientation")
    if alt >= 2 and tumor["alt_median_end_distance"] <= 5:
        reasons.append("ALT集中在读段末端")
    if tumor["raw_alt_reads"] and tumor["raw_alt_complex_cigar_reads"] / tumor["raw_alt_reads"] >= 0.4:
        reasons.append("ALT经常与soft-clip/indel复杂CIGAR共现")
    if normal["hq_alt_fragments"] >= 2 or normal["hq_alt_vaf"] >= 0.02:
        reasons.append("正常样本也有ALT支持")
    return "；".join(reasons) if reasons else "存在一定ALT证据，需结合Mutect2局部组装/下采样进一步判断"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tumor-bam", required=True, type=Path)
    parser.add_argument("--normal-bam", required=True, type=Path)
    parser.add_argument("--loci-file", required=True, type=Path)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    details: dict[str, Any] = {}
    with args.loci_file.open(encoding="utf-8", newline="") as handle:
        for locus in csv.DictReader(handle, delimiter="\t"):
            chrom, pos = locus["chrom"], int(locus["pos"])
            key = f"{chrom}:{pos}:{locus['ref']}:{locus['alt']}"
            tumor_obs = read_observations(args.samtools, args.tumor_bam, chrom, pos)
            normal_obs = read_observations(args.samtools, args.normal_bam, chrom, pos)
            tumor = summarize(tumor_obs, locus["ref"], locus["alt"])
            normal = summarize(normal_obs, locus["ref"], locus["alt"])
            row = {"gene": locus.get("gene", ""), "variant": key, **locus}
            row.update({f"tumor_{name}": value for name, value in tumor.items()})
            row.update({f"normal_{name}": value for name, value in normal.items()})
            row["preliminary_reason"] = preliminary_reason(tumor, normal)
            rows.append(row)
            details[key] = {"tumor": tumor_obs, "normal": normal_obs}

    output_tsv = args.output_dir / "bam_locus_evidence.tsv"
    with output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "bam_locus_read_details.json").write_text(
        json.dumps(details, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"loci={len(rows)}")
    print(f"summary={output_tsv}")


if __name__ == "__main__":
    main()
