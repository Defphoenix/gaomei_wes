#!/usr/bin/env python3
"""Run or emulate peptide-HLA binding prediction for pipeline integration tests."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable, List, Tuple


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name = ""
    chunks: List[str] = []
    with path.open("rt") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records.append((name, "".join(chunks)))
                name = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if name:
        records.append((name, "".join(chunks)))
    return records


def parse_alleles(raw: str) -> List[str]:
    return [item.strip() for item in raw.replace(";", ",").split(",") if item.strip()]


def simple_score(peptide: str, allele: str) -> Tuple[float, float]:
    hydrophobic = sum(1 for aa in peptide if aa in "AILMFWVY")
    aromatic = sum(1 for aa in peptide if aa in "FWY")
    length_penalty = abs(len(peptide) - 9) * 0.15
    allele_seed = (sum(ord(ch) for ch in allele) % 17) / 100.0
    score = max(0.01, min(0.99, 0.12 + hydrophobic / max(len(peptide), 1) * 0.65 + aromatic * 0.03 - length_penalty + allele_seed))
    affinity_nm = max(20.0, 50000.0 / (1.0 + score * 80.0))
    return round(affinity_nm, 3), round((1.0 - score) * 5.0, 4)


def write_simple(records: List[Tuple[str, str]], alleles: List[str], output: Path, threshold: float) -> None:
    with output.open("wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["peptide_id", "peptide", "allele", "prediction_tool", "affinity_nm", "rank", "binder"])
        for peptide_id, peptide in records:
            for allele in alleles:
                affinity, rank = simple_score(peptide, allele)
                writer.writerow([peptide_id, peptide, allele, "simple", affinity, rank, "yes" if affinity <= threshold else "no"])


def run_mhcflurry(records: List[Tuple[str, str]], alleles: List[str], output: Path, tool: str) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        peptide_csv = Path(tmpdir) / "peptides.csv"
        with peptide_csv.open("wt", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["allele", "peptide", "peptide_id"])
            for peptide_id, peptide in records:
                for allele in alleles:
                    writer.writerow([allele, peptide, peptide_id])
        cmd = [tool, "--out", str(output), "--output-delimiter", "\t", "--no-throw", str(peptide_csv)]
        subprocess.run(cmd, check=True)


def run_netmhcpan(records: List[Tuple[str, str]], alleles: List[str], output: Path, tool: str) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta = Path(tmpdir) / "peptides.fa"
        with fasta.open("wt") as handle:
            for peptide_id, peptide in records:
                handle.write(f">{peptide_id}\n{peptide}\n")
        with output.open("wt") as handle:
            subprocess.run([tool, "-f", str(fasta), "-a", ",".join(alleles), "-BA"], check=True, stdout=handle)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peptides", required=True, type=Path)
    parser.add_argument("--alleles", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tool", default="auto", choices=["auto", "netmhcpan", "mhcflurry", "simple"])
    parser.add_argument("--netmhcpan-bin", default="netMHCpan")
    parser.add_argument("--mhcflurry-bin", default="mhcflurry-predict")
    parser.add_argument("--threshold-nm", type=float, default=500.0)
    args = parser.parse_args()

    records = read_fasta(args.peptides)
    alleles = parse_alleles(args.alleles)
    if not records:
        raise ValueError(f"No peptides found in FASTA: {args.peptides}")
    if not alleles:
        raise ValueError("No HLA alleles provided")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    tool = args.tool
    if tool == "auto":
        if shutil.which(args.netmhcpan_bin):
            tool = "netmhcpan"
        elif shutil.which(args.mhcflurry_bin):
            tool = "mhcflurry"
        else:
            tool = "simple"

    if tool == "netmhcpan":
        if not shutil.which(args.netmhcpan_bin):
            raise FileNotFoundError(f"netMHCpan not found: {args.netmhcpan_bin}")
        run_netmhcpan(records, alleles, args.output, args.netmhcpan_bin)
    elif tool == "mhcflurry":
        if not shutil.which(args.mhcflurry_bin):
            raise FileNotFoundError(f"mhcflurry-predict not found: {args.mhcflurry_bin}")
        run_mhcflurry(records, alleles, args.output, args.mhcflurry_bin)
    else:
        write_simple(records, alleles, args.output, args.threshold_nm)

    print(f"binding_records={len(records) * len(alleles)}")
    print(f"binding_output={args.output}")
    print(f"binding_tool={tool}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
