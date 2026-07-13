#!/usr/bin/env python3
"""Generate tiny tumor-normal FASTQ demo datasets for WES pipeline testing.

The reads are simulated from real hg38 reference windows around TP53, BRCA1,
and ERBB2/HER2. Tumor reads contain deterministic SNV alleles at configurable
fractions; normal reads remain reference.
"""

from __future__ import annotations

import argparse
import gzip
import random
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Target:
    gene: str
    chrom: str
    start: int
    end: int
    variant_pos: int


@dataclass(frozen=True)
class Sample:
    sample_id: str
    role: str
    pair_id: str
    mutations: dict[str, float]


TARGETS = [
    Target("TP53", "chr17", 7673770, 7674670, 7674220),
    Target("BRCA1", "chr17", 43070500, 43071400, 43070945),
    Target("ERBB2_HER2", "chr17", 39724250, 39725150, 39724701),
]


SAMPLES = [
    Sample("demo_pair01_normal", "normal", "demo_pair01", {}),
    Sample("demo_pair01_tumor", "tumor", "demo_pair01", {"TP53": 0.35, "ERBB2_HER2": 0.28}),
    Sample("demo_pair02_normal", "normal", "demo_pair02", {}),
    Sample("demo_pair02_tumor", "tumor", "demo_pair02", {"BRCA1": 0.40, "TP53": 0.22}),
]


def run_text(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, text=True)


def fetch_ref(samtools: str, fasta: Path, chrom: str, start: int, end: int) -> str:
    text = run_text([samtools, "faidx", str(fasta), f"{chrom}:{start}-{end}"])
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">")).upper()


def alt_base(ref: str) -> str:
    return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(ref.upper(), "A")


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def mutate_fragment(fragment: str, target: Target, fragment_start: int, alt: str) -> str:
    idx = target.variant_pos - fragment_start
    if 0 <= idx < len(fragment):
        fragment = fragment[:idx] + alt + fragment[idx + 1 :]
    return fragment


def write_fastq_record(handle, name: str, seq: str) -> None:
    qual = "I" * len(seq)
    handle.write(f"@{name}\n{seq}\n+\n{qual}\n")


def simulate_sample(
    sample: Sample,
    fasta: Path,
    out_dir: Path,
    targets: list[Target],
    pairs_per_target: int,
    read_len: int,
    insert_size: int,
    seed: int,
    samtools: str,
) -> list[dict[str, str]]:
    rng = random.Random(seed)
    sample_dir = out_dir / sample.sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    r1_path = sample_dir / f"{sample.sample_id}_R1.fastq.gz"
    r2_path = sample_dir / f"{sample.sample_id}_R2.fastq.gz"
    truth_rows: list[dict[str, str]] = []

    target_refs = {t.gene: fetch_ref(samtools, fasta, t.chrom, t.start, t.end) for t in targets}

    with gzip.open(r1_path, "wt") as r1, gzip.open(r2_path, "wt") as r2:
        read_index = 1
        for target in targets:
            ref_window = target_refs[target.gene]
            ref_at_variant = ref_window[target.variant_pos - target.start]
            alt = alt_base(ref_at_variant)
            allele_fraction = sample.mutations.get(target.gene, 0.0)
            alt_pairs = int(round(pairs_per_target * allele_fraction))

            if allele_fraction > 0:
                truth_rows.append(
                    {
                        "pair_id": sample.pair_id,
                        "tumor_sample": sample.sample_id,
                        "gene": target.gene,
                        "chrom": target.chrom,
                        "pos": str(target.variant_pos),
                        "ref": ref_at_variant,
                        "alt": alt,
                        "expected_af": f"{allele_fraction:.2f}",
                    }
                )

            max_offset = len(ref_window) - insert_size
            for pair_i in range(pairs_per_target):
                # Bias fragments toward the variant position while keeping the
                # rest of the target covered enough for coverage plots.
                if pair_i < int(pairs_per_target * 0.75):
                    center = target.variant_pos - target.start - insert_size // 2
                    offset = max(0, min(max_offset, center + rng.randint(-80, 80)))
                else:
                    offset = rng.randint(0, max_offset)

                fragment_start = target.start + offset
                fragment = ref_window[offset : offset + insert_size]
                carries_alt = pair_i < alt_pairs
                if carries_alt:
                    fragment = mutate_fragment(fragment, target, fragment_start, alt)

                r1_seq = fragment[:read_len]
                r2_seq = reverse_complement(fragment[-read_len:])
                read_name = f"{sample.sample_id}:{target.gene}:{read_index}:{'ALT' if carries_alt else 'REF'}"
                write_fastq_record(r1, f"{read_name}/1", r1_seq)
                write_fastq_record(r2, f"{read_name}/2", r2_seq)
                read_index += 1

    return truth_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--pairs-per-target", type=int, default=330)
    parser.add_argument("--read-len", type=int, default=150)
    parser.add_argument("--insert-size", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20260713)
    parser.add_argument("--samtools", default="/Users/mac/anaconda3/envs/big_wes_pipeline_env/bin/samtools")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fastq_dir = args.out_dir / "fastq"
    fastq_dir.mkdir(parents=True, exist_ok=True)

    bed_path = args.out_dir / "demo_targets.tp53_brca1_her2.bed"
    with bed_path.open("w") as bed:
        for target in sorted(TARGETS, key=lambda item: (item.chrom, item.start)):
            bed.write(f"{target.chrom}\t{target.start - 1}\t{target.end}\t{target.gene}\n")

    truth_rows: list[dict[str, str]] = []
    for i, sample in enumerate(SAMPLES):
        truth_rows.extend(
            simulate_sample(
                sample,
                args.reference,
                fastq_dir,
                TARGETS,
                args.pairs_per_target,
                args.read_len,
                args.insert_size,
                args.seed + i,
                args.samtools,
            )
        )

    truth_path = args.out_dir / "demo_truth_mutations.tsv"
    with truth_path.open("w") as truth:
        cols = ["pair_id", "tumor_sample", "gene", "chrom", "pos", "ref", "alt", "expected_af"]
        truth.write("\t".join(cols) + "\n")
        for row in truth_rows:
            truth.write("\t".join(row[c] for c in cols) + "\n")

    readme = args.out_dir / "README_demo_data.md"
    readme.write_text(
        "\n".join(
            [
                "# Synthetic tumor-normal WES demo data",
                "",
                "Generated samples:",
                "- demo_pair01_normal / demo_pair01_tumor",
                "- demo_pair02_normal / demo_pair02_tumor",
                "",
                "Target genes: TP53, BRCA1, ERBB2/HER2.",
                f"BED: {bed_path}",
                f"Truth TSV: {truth_path}",
                "",
                "These FASTQs are synthetic and are intended only for pipeline demonstration.",
            ]
        )
        + "\n"
    )

    print(f"FASTQ directory: {fastq_dir}")
    print(f"BED: {bed_path}")
    print(f"Truth: {truth_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
