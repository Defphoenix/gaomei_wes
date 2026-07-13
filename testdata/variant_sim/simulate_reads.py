#!/usr/bin/env python3
"""Simulate 100 paired-end reads with known variants from a local reference."""

from __future__ import annotations

import argparse
import gzip
import subprocess
from pathlib import Path


COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1].upper()


def fetch_reference(reference: Path, region: str) -> str:
    result = subprocess.run(
        ["samtools", "faidx", str(reference), region],
        check=True,
        text=True,
        capture_output=True,
    )
    return "".join(line.strip() for line in result.stdout.splitlines() if not line.startswith(">")).upper()


def build_haplotype(ref_seq: str, variants: list[dict]) -> tuple[str, dict[int, int]]:
    pieces = []
    ref_to_hap = {}
    ref_idx = 0
    hap_idx = 0
    by_offset = {var["offset"]: var for var in variants}

    while ref_idx < len(ref_seq):
        ref_to_hap[ref_idx] = hap_idx
        var = by_offset.get(ref_idx)
        if var is None:
            pieces.append(ref_seq[ref_idx])
            ref_idx += 1
            hap_idx += 1
            continue

        if var["type"] == "snp":
            pieces.append(var["alt"])
            ref_idx += 1
            hap_idx += 1
        elif var["type"] == "del":
            ref_len = len(var["ref"])
            alt = var["alt"]
            pieces.append(alt)
            ref_idx += ref_len
            hap_idx += len(alt)
        else:
            raise ValueError(f"Unsupported variant type: {var['type']}")

    ref_to_hap[len(ref_seq)] = hap_idx
    return "".join(pieces), ref_to_hap


def write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        for name, seq in records:
            handle.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--region", default="chr22:20000000-20001000")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--sample", default="sim100")
    parser.add_argument("--read-count", type=int, default=100)
    parser.add_argument("--read-len", type=int, default=100)
    parser.add_argument("--insert", type=int, default=240)
    args = parser.parse_args()

    chrom, coords = args.region.split(":")
    region_start = int(coords.split("-")[0])
    ref_seq = fetch_reference(args.reference, args.region)
    if "N" in ref_seq:
        raise ValueError(f"Reference region contains N bases: {args.region}")

    variants = [
        {"name": "known_snp_1", "type": "snp", "offset": 300, "ref": ref_seq[300], "alt": "T" if ref_seq[300] != "T" else "A"},
        {"name": "known_snp_2", "type": "snp", "offset": 470, "ref": ref_seq[470], "alt": "G" if ref_seq[470] != "G" else "C"},
        {"name": "known_del_3bp", "type": "del", "offset": 620, "ref": ref_seq[620:624], "alt": ref_seq[620]},
    ]
    hap_seq, ref_to_hap = build_haplotype(ref_seq, variants)

    r1_records = []
    r2_records = []
    max_start = len(ref_seq) - args.insert - 2
    for i in range(args.read_count):
        ref_start = 120 + (i * 7) % min(520, max_start - 120)
        hap_start = ref_to_hap[ref_start]
        fragment = hap_seq[hap_start : hap_start + args.insert]
        if len(fragment) < args.insert:
            raise ValueError("Fragment exceeded haplotype boundary")

        r1 = fragment[: args.read_len]
        r2 = revcomp(fragment[-args.read_len :])
        name = f"{args.sample}_{i + 1:03d}"
        r1_records.append((f"{name}/1", r1))
        r2_records.append((f"{name}/2", r2))

    write_fastq(args.outdir / f"{args.sample}_R1.fastq.gz", r1_records)
    write_fastq(args.outdir / f"{args.sample}_R2.fastq.gz", r2_records)

    truth_path = args.outdir / f"{args.sample}_truth.tsv"
    with truth_path.open("wt") as handle:
        handle.write("name\tchrom\tpos_1based\tref\talt\ttype\n")
        for var in variants:
            pos = region_start + var["offset"]
            handle.write(
                f"{var['name']}\t{chrom}\t{pos}\t{var['ref']}\t{var['alt']}\t{var['type']}\n"
            )

    bed_path = args.outdir / f"{args.sample}_target.bed"
    with bed_path.open("wt") as handle:
        bed_start = region_start - 1
        bed_end = int(coords.split("-")[1])
        handle.write(f"{chrom}\t{bed_start}\t{bed_end}\n")

    print(f"sample={args.sample}")
    print(f"region={args.region}")
    print(f"read_pairs={args.read_count}")
    print(f"r1={args.outdir / f'{args.sample}_R1.fastq.gz'}")
    print(f"r2={args.outdir / f'{args.sample}_R2.fastq.gz'}")
    print(f"truth={truth_path}")
    print(f"bed={bed_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
