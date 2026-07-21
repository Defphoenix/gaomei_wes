#!/usr/bin/env python3
"""Approximate coding HGVS to genomic coordinates using an Ensembl GTF."""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd


ATTR_RE = re.compile(r'(\S+) "([^"]*)"')
SNV_RE = re.compile(r'^c\.(\d+)([+-]\d+)?([ACGT])>([ACGT])$', re.I)
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def parse_attributes(text: str) -> tuple[dict[str, str], set[str]]:
    pairs = ATTR_RE.findall(text)
    attrs: dict[str, str] = {}
    tags: set[str] = set()
    for key, value in pairs:
        if key == "tag":
            tags.add(value)
        else:
            attrs[key] = value
    return attrs, tags


def load_transcripts(gtf: Path, genes: set[str]) -> dict[str, list[dict[str, Any]]]:
    wanted = re.compile(r'gene_name "(' + "|".join(re.escape(gene) for gene in sorted(genes)) + r')"')
    transcripts: dict[str, dict[str, Any]] = {}
    with gtf.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not wanted.search(line):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "CDS"}:
                continue
            attrs, tags = parse_attributes(fields[8])
            transcript_id = attrs.get("transcript_id", "")
            gene = attrs.get("gene_name", "")
            if not transcript_id or gene not in genes:
                continue
            record = transcripts.setdefault(transcript_id, {
                "gene": gene, "transcript": transcript_id,
                "transcript_version": attrs.get("transcript_version", ""),
                "chrom": fields[0], "strand": 1 if fields[6] == "+" else -1,
                "biotype": attrs.get("transcript_biotype", ""), "tags": set(), "cds": [],
            })
            record["tags"].update(tags)
            if attrs.get("transcript_version"):
                record["transcript_version"] = attrs["transcript_version"]
            if fields[2] == "CDS":
                record["cds"].append((int(fields[3]), int(fields[4])))
    by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in transcripts.values():
        if record["cds"]:
            by_gene[record["gene"]].append(record)
    return by_gene


def transcript_score(record: dict[str, Any]) -> tuple[int, int, int, int, int]:
    tags = record["tags"]
    cds_length = sum(end - start + 1 for start, end in record["cds"])
    return (
        int("Ensembl_canonical" in tags), int("MANE_Select" in tags),
        int("gencode_primary" in tags), int(record["biotype"] == "protein_coding"), cds_length,
    )


def map_cds_position(record: dict[str, Any], cds_position: int, intron_offset: int) -> int | None:
    segments = sorted(record["cds"], reverse=record["strand"] == -1)
    consumed = 0
    for start, end in segments:
        length = end - start + 1
        if consumed < cds_position <= consumed + length:
            within = cds_position - consumed - 1
            genomic = start + within if record["strand"] == 1 else end - within
            return genomic + intron_offset * record["strand"]
        consumed += length
    return None


class IndexedFasta:
    def __init__(self, path: Path):
        self.path = path
        self.handle = path.open("rb")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        with Path(f"{path}.fai").open() as handle:
            for line in handle:
                name, length, offset, line_bases, line_width, *_ = line.rstrip().split("\t")
                self.index[name] = (int(length), int(offset), int(line_bases), int(line_width))

    def base(self, chrom: str, pos: int) -> str:
        candidates = [chrom, chrom if chrom.startswith("chr") else f"chr{chrom}"]
        name = next((item for item in candidates if item in self.index), "")
        if not name:
            return ""
        length, offset, line_bases, line_width = self.index[name]
        if pos < 1 or pos > length:
            return ""
        zero = pos - 1
        byte_offset = offset + (zero // line_bases) * line_width + (zero % line_bases)
        self.handle.seek(byte_offset)
        return self.handle.read(1).decode("ascii").upper()

    def close(self) -> None:
        self.handle.close()


def map_variants(frame: pd.DataFrame, gtf: Path, fasta: Path | None = None) -> pd.DataFrame:
    genes = {str(value).strip() for value in frame["gene"] if str(value).strip()}
    by_gene = load_transcripts(gtf, genes)
    genome = IndexedFasta(fasta) if fasta else None
    rows: list[dict[str, Any]] = []
    try:
        for item in frame.to_dict("records"):
            gene = str(item.get("gene", "")).strip()
            hgvsc = str(item.get("hgvsc_other", item.get("hgvsc", ""))).strip()
            match = SNV_RE.match(hgvsc)
            transcripts = sorted(by_gene.get(gene, []), key=transcript_score, reverse=True)
            if not match:
                rows.append(dict(item) | {"gtf_mapping_status": "不支持的HGVS写法"})
                continue
            if not transcripts:
                rows.append(dict(item) | {"gtf_mapping_status": "GTF未找到编码转录本"})
                continue
            cds_position = int(match.group(1))
            intron_offset = int(match.group(2) or 0)
            cdna_ref, cdna_alt = match.group(3).upper(), match.group(4).upper()
            selected = transcripts[0]
            genomic_pos = map_cds_position(selected, cds_position, intron_offset)
            if genomic_pos is None:
                rows.append(dict(item) | {
                    "gtf_transcript": selected["transcript"],
                    "gtf_mapping_status": "canonical CDS长度不足，无法映射",
                })
                continue
            genomic_ref = cdna_ref if selected["strand"] == 1 else cdna_ref.translate(COMPLEMENT)
            genomic_alt = cdna_alt if selected["strand"] == 1 else cdna_alt.translate(COMPLEMENT)
            fasta_ref = genome.base(selected["chrom"], genomic_pos) if genome else ""
            ref_check = "通过" if fasta_ref == genomic_ref else (f"不一致(FASTA={fasta_ref})" if fasta_ref else "未检查")
            version = f".{selected['transcript_version']}" if selected["transcript_version"] else ""
            rows.append(dict(item) | {
                "gtf_chrom": selected["chrom"] if str(selected["chrom"]).startswith("chr") else f"chr{selected['chrom']}",
                "gtf_pos": genomic_pos, "gtf_ref": genomic_ref, "gtf_alt": genomic_alt,
                "gtf_transcript": f"{selected['transcript']}{version}",
                "gtf_tags": ";".join(sorted(selected["tags"])),
                "gtf_reference_check": ref_check,
                "gtf_mapping_status": "近似映射成功" if ref_check == "通过" else "近似映射需复核",
            })
    finally:
        if genome:
            genome.close()
    return pd.DataFrame(rows)


def load_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip().split("\t")
            try:
                chrom = fields[0] if fields[0].startswith("chr") else f"chr{fields[0]}"
                intervals[chrom].append((int(fields[1]), int(fields[2])))
            except (IndexError, ValueError):
                continue
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def overlaps_bed(intervals: dict[str, list[tuple[int, int]]], chrom: str, pos: int) -> bool:
    zero = pos - 1
    return any(start <= zero < end for start, end in intervals.get(chrom, []))


def main() -> None:
    parser = argparse.ArgumentParser(description="用Ensembl GTF将简单c.HGVS近似映射到基因组并检查BED")
    parser.add_argument("--input", required=True, type=Path, help="包含gene和hgvsc_other列的TSV")
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--fasta", type=Path)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    frame = pd.read_csv(args.input, sep="\t", dtype=str).fillna("")
    mapped = map_variants(frame, args.gtf, args.fasta)
    intervals = load_bed(args.bed)
    def classify_bed(row: pd.Series) -> str:
        pos = row.get("gtf_pos")
        if pd.isna(pos) or pos == "":
            return "无法判断"
        return "BED内" if overlaps_bed(intervals, str(row.get("gtf_chrom", "")), int(pos)) else "BED外"

    mapped["gtf_bed_status"] = mapped.apply(classify_bed, axis=1)
    mapped.to_csv(args.output, sep="\t", index=False)
    print(mapped[["gene", "hgvsc_other", "gtf_transcript", "gtf_chrom", "gtf_pos", "gtf_ref", "gtf_alt", "gtf_reference_check", "gtf_bed_status"]].to_string(index=False))


if __name__ == "__main__":
    main()
