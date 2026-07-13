#!/usr/bin/env python3
"""Generate candidate neoantigen peptides from a VEP-annotated VCF.

The script reads VEP CSQ annotations, maps affected transcripts/proteins to a
wild-type protein FASTA, reconstructs simple amino-acid changes when possible,
and emits mutant peptide FASTA plus a manifest table.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


SUPPORTED = {
    "missense_variant",
    "protein_altering_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
}
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class CsqRecord:
    chrom: str
    pos: str
    ref: str
    alt: str
    var_id: str
    gene: str
    symbol: str
    feature: str
    protein_id: str
    consequence: str
    protein_position: str
    amino_acids: str
    hgvsp: str
    canonical: str
    mutant_protein: str


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def parse_info(info: str) -> Dict[str, str]:
    parsed: Dict[str, str] = {}
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
        else:
            parsed[item] = "true"
    return parsed


def parse_csq_header(vcf: Path) -> List[str]:
    fmt_re = re.compile(r'Format: ([^"]+)')
    with open_text(vcf) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                match = fmt_re.search(line)
                if not match:
                    raise ValueError("CSQ header found, but no Format field was detected")
                return match.group(1).split("|")
            if line.startswith("#CHROM"):
                break
    raise ValueError("No VEP CSQ header found in VCF")


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []

    def flush() -> None:
        if name is None:
            return
        sequence = "".join(chunks).replace("*", "").upper()
        if not sequence:
            return
        aliases = {name, name.split()[0]}
        for token in re.split(r"[|\s]", name):
            if token:
                aliases.add(token)
                if "." in token:
                    aliases.add(token.split(".", 1)[0])
        for alias in aliases:
            records.setdefault(alias, sequence)

    with path.open("rt") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        flush()
    return records


def first_int(value: str) -> Optional[int]:
    match = re.search(r"\d+", value or "")
    if not match:
        return None
    return int(match.group(0))


def get_field(values: Dict[str, str], *names: str) -> str:
    for name in names:
        if name in values and values[name]:
            return values[name]
    return ""


def iter_csq_records(vcf: Path, csq_fields: Sequence[str]) -> Iterable[CsqRecord]:
    with open_text(vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, var_id, ref, alt, _qual, _flt, info = parts[:8]
            info_map = parse_info(info)
            if "CSQ" not in info_map:
                continue
            record_mutant_protein = get_field(info_map, "MUTANT_PROTEIN", "MutantProtein")
            for raw_csq in info_map["CSQ"].split(","):
                values = dict(zip(csq_fields, raw_csq.split("|")))
                consequences = values.get("Consequence", "")
                if not any(item in SUPPORTED for item in consequences.split("&")):
                    continue
                yield CsqRecord(
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    var_id=var_id if var_id != "." else f"{chrom}:{pos}:{ref}>{alt}",
                    gene=values.get("Gene", ""),
                    symbol=values.get("SYMBOL", ""),
                    feature=values.get("Feature", ""),
                    protein_id=get_field(values, "ENSP", "Protein_id", "Protein"),
                    consequence=consequences,
                    protein_position=values.get("Protein_position", ""),
                    amino_acids=values.get("Amino_acids", ""),
                    hgvsp=values.get("HGVSp", ""),
                    canonical=values.get("CANONICAL", ""),
                    mutant_protein=get_field(values, "MUTANT_PROTEIN", "MutantProtein") or record_mutant_protein,
                )


def choose_best(records: Iterable[CsqRecord]) -> List[CsqRecord]:
    grouped: Dict[str, List[CsqRecord]] = {}
    for rec in records:
        grouped.setdefault(rec.var_id, []).append(rec)
    chosen: List[CsqRecord] = []
    for recs in grouped.values():
        recs.sort(
            key=lambda r: (
                0 if r.canonical.upper() == "YES" else 1,
                0 if r.feature else 1,
                0 if r.protein_id else 1,
            )
        )
        chosen.append(recs[0])
    return chosen


def lookup_protein(rec: CsqRecord, proteins: Dict[str, str]) -> Optional[str]:
    for key in (rec.protein_id, rec.feature):
        if key and key in proteins:
            return proteins[key]
    return None


def reconstruct_mutant(rec: CsqRecord, wt: Optional[str]) -> Tuple[Optional[str], str]:
    if rec.mutant_protein:
        return rec.mutant_protein.replace("*", "").upper(), "custom_mutant_protein"
    if not wt:
        return None, "missing_wildtype_protein"

    pos = first_int(rec.protein_position)
    aa_change = rec.amino_acids
    if not pos or pos < 1 or pos > len(wt):
        return None, "invalid_or_missing_protein_position"

    if "/" not in aa_change:
        return None, "missing_amino_acid_change"
    ref_aa, alt_aa = aa_change.split("/", 1)
    alt_aa = alt_aa.replace("*", "")

    consequences = set(rec.consequence.split("&"))
    if "missense_variant" in consequences or "protein_altering_variant" in consequences:
        if len(alt_aa) == 1 and alt_aa in VALID_AA:
            return wt[: pos - 1] + alt_aa + wt[pos:], "simple_aa_substitution"
        return None, "unsupported_missense_amino_acid_change"

    if "inframe_insertion" in consequences and alt_aa:
        return wt[: pos - 1] + alt_aa + wt[pos - 1 :], "simple_inframe_insertion"

    if "inframe_deletion" in consequences:
        delete_len = len(ref_aa.replace("-", ""))
        if delete_len > 0:
            return wt[: pos - 1] + wt[pos - 1 + delete_len :], "simple_inframe_deletion"
        return None, "unsupported_inframe_deletion"

    if "stop_gained" in consequences:
        return wt[: pos - 1], "stop_gained_truncation"

    if "stop_lost" in consequences:
        return None, "stop_lost_requires_mutant_protein_sequence"

    if "start_lost" in consequences:
        return None, "start_lost_requires_alternative_start_model"

    if "frameshift_variant" in consequences:
        return None, "frameshift_requires_mutant_protein_sequence"

    return None, "unsupported_consequence"


def peptide_window(seq: str, center_pos: Optional[int], flank: int) -> str:
    if center_pos is None or not seq:
        return seq
    center = min(max(center_pos - 1, 0), len(seq) - 1)
    start = max(center - flank, 0)
    end = min(center + flank + 1, len(seq))
    return seq[start:end]


def event_type(rec: CsqRecord) -> str:
    consequences = set(rec.consequence.split("&"))
    if "missense_variant" in consequences or "protein_altering_variant" in consequences:
        return "missense"
    if "stop_gained" in consequences:
        return "stop_gained"
    if "stop_lost" in consequences:
        return "stop_lost"
    if "start_lost" in consequences:
        return "start_lost"
    if "frameshift_variant" in consequences:
        return "frameshift"
    if "inframe_insertion" in consequences:
        return "inframe_insertion"
    if "inframe_deletion" in consequences:
        return "inframe_deletion"
    return "other"


def peptides_for_record(
    rec: CsqRecord,
    wt: Optional[str],
    mutant: str,
    flanks: Sequence[int],
    flank: int,
) -> List[Tuple[int, int, str, str]]:
    pos = first_int(rec.protein_position)
    peptides: List[Tuple[int, int, str, str]] = []
    seen = set()
    for window_flank in flanks:
        wt_window = peptide_window(wt, pos, window_flank) if wt else ""
        mut_window = peptide_window(mutant, pos, window_flank)
        key = (window_flank, mut_window)
        if not mut_window or key in seen:
            continue
        if mut_window == wt_window:
            continue
        if set(mut_window) <= VALID_AA:
            peptides.append((window_flank, len(mut_window), mut_window, wt_window))
            seen.add(key)
    return peptides


def write_outputs(
    records: Sequence[CsqRecord],
    proteins: Dict[str, str],
    fasta_out: Path,
    manifest_out: Path,
    sample: str,
    lengths: Sequence[int],
    flank: int,
) -> None:
    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    manifest_out.parent.mkdir(parents=True, exist_ok=True)

    with fasta_out.open("wt") as fasta, manifest_out.open("wt") as manifest:
        manifest.write(
            "\t".join(
                [
                    "sample",
                    "variant_id",
                    "gene",
                    "feature",
                    "protein_id",
                    "consequence",
                    "event_type",
                    "protein_position",
                    "amino_acids",
                    "method",
                    "status",
                    "peptide_id",
                    "window_flank",
                    "window_length",
                    "wildtype_window",
                    "mutant_window",
                    "reason",
                ]
            )
            + "\n"
        )

        peptide_count = 0
        for rec in records:
            wt = lookup_protein(rec, proteins)
            mutant, method = reconstruct_mutant(rec, wt)
            if not mutant:
                manifest.write(
                    "\t".join(
                        [
                            sample,
                            rec.var_id,
                            rec.symbol or rec.gene,
                            rec.feature,
                            rec.protein_id,
                            rec.consequence,
                            event_type(rec),
                            rec.protein_position,
                            rec.amino_acids,
                            method,
                            "skipped",
                            "",
                            "",
                            "",
                            "",
                            "",
                            method,
                        ]
                    )
                    + "\n"
                )
                continue

            peptides = peptides_for_record(rec, wt, mutant, lengths, flank)
            if not peptides:
                manifest.write(
                    "\t".join(
                        [
                            sample,
                            rec.var_id,
                            rec.symbol or rec.gene,
                            rec.feature,
                            rec.protein_id,
                            rec.consequence,
                            event_type(rec),
                            rec.protein_position,
                            rec.amino_acids,
                            method,
                            "no_novel_peptide",
                            "",
                            "",
                            "",
                            "",
                            "",
                            "no_mutant_peptide_not_seen_in_wildtype_window",
                        ]
                    )
                    + "\n"
                )
                continue

            for idx, (window_flank, window_length, pep, wt_window) in enumerate(peptides, start=1):
                peptide_count += 1
                peptide_id = f"{sample}|{rec.var_id}|{rec.symbol or rec.gene}|{rec.feature}|flank{window_flank}|p{idx}"
                fasta.write(f">{peptide_id}\n{pep}\n")
                manifest.write(
                    "\t".join(
                        [
                            sample,
                            rec.var_id,
                            rec.symbol or rec.gene,
                            rec.feature,
                            rec.protein_id,
                            rec.consequence,
                            event_type(rec),
                            rec.protein_position,
                            rec.amino_acids,
                            method,
                            "emitted",
                            peptide_id,
                            str(window_flank),
                            str(window_length),
                            wt_window,
                            pep,
                            "",
                        ]
                    )
                    + "\n"
                )

    print(f"neoantigen_records={len(records)}")
    print(f"neoantigen_peptides={sum(1 for line in fasta_out.open() if line.startswith('>'))}")


def parse_lengths(raw: str) -> List[int]:
    lengths = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        lengths.append(int(item))
    if not lengths:
        raise ValueError("No peptide lengths were provided")
    return lengths


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vep-vcf", required=True, type=Path)
    parser.add_argument("--protein-fasta", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--lengths", default="8,9,10,11", help="Window flanks around the altered amino acid. 8 means 8 aa upstream + mutated aa + 8 aa downstream, usually 17 aa total.")
    parser.add_argument("--flank", default=30, type=int)
    args = parser.parse_args(argv)

    if not args.vep_vcf.exists():
        raise FileNotFoundError(f"VEP VCF not found: {args.vep_vcf}")
    if not args.protein_fasta.exists():
        raise FileNotFoundError(f"Protein FASTA not found: {args.protein_fasta}")

    csq_fields = parse_csq_header(args.vep_vcf)
    proteins = read_fasta(args.protein_fasta)
    records = choose_best(iter_csq_records(args.vep_vcf, csq_fields))
    write_outputs(
        records=records,
        proteins=proteins,
        fasta_out=args.output_fasta,
        manifest_out=args.manifest,
        sample=args.sample,
        lengths=parse_lengths(args.lengths),
        flank=args.flank,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
