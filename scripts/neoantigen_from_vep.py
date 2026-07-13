#!/usr/bin/env python3
"""Generate candidate neoantigen peptides from a VEP-annotated VCF.

The script reads VEP CSQ annotations, maps affected transcripts/proteins to a
wild-type protein FASTA, reconstructs simple amino-acid changes when possible,
and emits mutant peptide FASTA plus a manifest table.
"""

from __future__ import annotations

import argparse
import csv
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
    hgvsc: str
    hgvsp: str
    canonical: str
    mutant_protein: str


@dataclass
class MutantProtein:
    sequence: str
    method: str
    changed_start: int
    changed_end: int
    wt_prefix: str
    wt_changed: str
    wt_suffix: str
    mut_prefix: str
    mut_changed: str
    mut_suffix: str


@dataclass
class PeptideRecord:
    peptide_id: str
    mer: int
    start: int
    end: int
    peptide: str
    wt_peptide: str


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


def variant_key(chrom: str, pos: str, ref: str, alt: str) -> str:
    return f"{chrom}:{pos}:{ref}>{alt}"


def read_annovar_table(path: Optional[Path]) -> Dict[str, Dict[str, str]]:
    """Read an optional ANNOVAR-style table keyed by chr/start/ref/alt.

    The function is intentionally permissive because ANNOVAR table headers vary
    between exonic_variant_function, multianno, and custom post-processed files.
    """
    if not path or not path.exists():
        return {}
    with path.open("rt", newline="") as handle:
        sample = handle.readline()
        if not sample:
            return {}
        delimiter = "\t" if "\t" in sample else ","
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter=delimiter)
        ann: Dict[str, Dict[str, str]] = {}
        for row in reader:
            chrom = get_field(row, "Chr", "CHROM", "chrom", "#CHROM")
            pos = get_field(row, "Start", "POS", "pos")
            ref = get_field(row, "Ref", "REF", "ref")
            alt = get_field(row, "Alt", "ALT", "alt")
            if chrom and pos and ref and alt:
                ann[variant_key(chrom, pos, ref, alt)] = row
        return ann


def annovar_text(row: Dict[str, str]) -> str:
    if not row:
        return ""
    preferred = [
        "Func.refGene",
        "Gene.refGene",
        "ExonicFunc.refGene",
        "AAChange.refGene",
        "Func.ensGene",
        "Gene.ensGene",
        "ExonicFunc.ensGene",
        "AAChange.ensGene",
    ]
    chunks = []
    for key in preferred:
        value = row.get(key, "")
        if value:
            chunks.append(f"{key}={value}")
    if chunks:
        return ";".join(chunks)
    return ";".join(f"{key}={value}" for key, value in row.items() if value)


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
                    hgvsc=values.get("HGVSc", ""),
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


def split_protein(seq: str, start: int, end: int) -> Tuple[str, str, str]:
    start = max(start, 1)
    end = max(end, start)
    return seq[: start - 1], seq[start - 1 : end], seq[end:]


def reconstruct_mutant(rec: CsqRecord, wt: Optional[str]) -> Tuple[Optional[MutantProtein], str]:
    if rec.mutant_protein:
        seq = rec.mutant_protein.replace("*", "").upper()
        pos = first_int(rec.protein_position) or 1
        prefix, changed, suffix = split_protein(seq, pos, pos)
        return MutantProtein(seq, "custom_mutant_protein", pos, pos, "", "", "", prefix, changed, suffix), "custom_mutant_protein"
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
            seq = wt[: pos - 1] + alt_aa + wt[pos:]
            wt_pre, wt_mid, wt_post = split_protein(wt, pos, pos)
            mut_pre, mut_mid, mut_post = split_protein(seq, pos, pos)
            return MutantProtein(seq, "simple_aa_substitution", pos, pos, wt_pre, wt_mid, wt_post, mut_pre, mut_mid, mut_post), "simple_aa_substitution"
        return None, "unsupported_missense_amino_acid_change"

    if "inframe_insertion" in consequences and alt_aa:
        seq = wt[: pos - 1] + alt_aa + wt[pos - 1 :]
        mut_end = pos + len(alt_aa) - 1
        wt_pre, wt_mid, wt_post = split_protein(wt, pos, pos)
        mut_pre, mut_mid, mut_post = split_protein(seq, pos, mut_end)
        return MutantProtein(seq, "simple_inframe_insertion", pos, mut_end, wt_pre, wt_mid, wt_post, mut_pre, mut_mid, mut_post), "simple_inframe_insertion"

    if "inframe_deletion" in consequences:
        delete_len = len(ref_aa.replace("-", ""))
        if delete_len > 0:
            seq = wt[: pos - 1] + wt[pos - 1 + delete_len :]
            mut_pos = min(pos, max(len(seq), 1))
            wt_pre, wt_mid, wt_post = split_protein(wt, pos, pos + delete_len - 1)
            mut_pre, mut_mid, mut_post = split_protein(seq, mut_pos, mut_pos)
            return MutantProtein(seq, "simple_inframe_deletion", mut_pos, mut_pos, wt_pre, wt_mid, wt_post, mut_pre, mut_mid, mut_post), "simple_inframe_deletion"
        return None, "unsupported_inframe_deletion"

    if "stop_gained" in consequences:
        seq = wt[: pos - 1]
        mut_pos = min(pos, max(len(seq), 1))
        wt_pre, wt_mid, wt_post = split_protein(wt, pos, pos)
        mut_pre, mut_mid, mut_post = split_protein(seq, mut_pos, mut_pos)
        return MutantProtein(seq, "stop_gained_truncation", mut_pos, mut_pos, wt_pre, wt_mid, wt_post, mut_pre, mut_mid, mut_post), "stop_gained_truncation"

    if "stop_lost" in consequences:
        return None, "stop_lost_requires_mutant_protein_sequence"

    if "start_lost" in consequences:
        return None, "start_lost_requires_alternative_start_model"

    if "frameshift_variant" in consequences:
        return None, "frameshift_requires_mutant_protein_sequence"

    return None, "unsupported_consequence"


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


def all_mutant_peptides(
    rec: CsqRecord,
    wt: Optional[str],
    mutant: MutantProtein,
    mers: Sequence[int],
) -> List[PeptideRecord]:
    peptides: List[PeptideRecord] = []
    seen = set()
    for mer in mers:
        if mer <= 0 or len(mutant.sequence) < mer:
            continue
        first_start = max(1, mutant.changed_end - mer + 1)
        last_start = min(mutant.changed_start, len(mutant.sequence) - mer + 1)
        for start in range(first_start, last_start + 1):
            end = start + mer - 1
            peptide = mutant.sequence[start - 1 : end]
            if len(peptide) != mer or not set(peptide) <= VALID_AA:
                continue
            wt_peptide = wt[start - 1 : end] if wt and end <= len(wt) else ""
            if wt and peptide in wt:
                continue
            key = (mer, start, peptide)
            if key in seen:
                continue
            seen.add(key)
            peptide_id = ""
            peptides.append(PeptideRecord(peptide_id, mer, start, end, peptide, wt_peptide))
    return peptides


def write_outputs(
    records: Sequence[CsqRecord],
    proteins: Dict[str, str],
    fasta_out: Path,
    manifest_out: Path,
    protein_table_out: Path,
    peptide_table_out: Path,
    fasta_dir: Path,
    sample: str,
    mers: Sequence[int],
    annovar: Dict[str, Dict[str, str]],
) -> None:
    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    manifest_out.parent.mkdir(parents=True, exist_ok=True)
    protein_table_out.parent.mkdir(parents=True, exist_ok=True)
    peptide_table_out.parent.mkdir(parents=True, exist_ok=True)
    fasta_dir.mkdir(parents=True, exist_ok=True)

    fasta_handles = {mer: (fasta_dir / f"{mer}mer.fa").open("wt") for mer in mers}
    try:
        with (
            fasta_out.open("wt") as fasta,
            manifest_out.open("wt") as manifest,
            protein_table_out.open("wt") as protein_table,
            peptide_table_out.open("wt") as peptide_table,
        ):
            protein_table.write(
                "\t".join(
                    [
                        "sample",
                        "variant_id",
                        "annovar_annotation",
                        "vep_consequence",
                        "gene",
                        "feature",
                        "protein_id",
                        "hgvsc",
                        "hgvsp",
                        "protein_position",
                        "amino_acids",
                        "method",
                        "status",
                        "wildtype_protein",
                        "mutant_protein",
                        "wildtype_prefix",
                        "wildtype_changed",
                        "wildtype_suffix",
                        "mutant_prefix",
                        "mutant_changed",
                        "mutant_suffix",
                        "reason",
                    ]
                )
                + "\n"
            )
            peptide_table.write(
                "\t".join(
                    [
                        "sample",
                        "peptide_id",
                        "variant_id",
                        "gene",
                        "feature",
                        "protein_id",
                        "hgvsc",
                        "hgvsp",
                        "amino_acids",
                        "mer",
                        "mutant_peptide_start",
                        "mutant_peptide_end",
                        "wildtype_peptide",
                        "mutant_peptide",
                    ]
                )
                + "\n"
            )
            manifest.write(
            "\t".join(
                [
                    "sample",
                    "variant_id",
                    "annovar_annotation",
                    "gene",
                    "feature",
                    "protein_id",
                    "consequence",
                    "event_type",
                    "protein_position",
                    "amino_acids",
                    "hgvsc",
                    "hgvsp",
                    "method",
                    "status",
                    "peptide_id",
                    "mer",
                    "mutant_peptide_start",
                    "mutant_peptide_end",
                    "wildtype_peptide",
                    "mutant_peptide",
                    "wildtype_protein",
                    "mutant_protein",
                    "mutant_prefix",
                    "mutant_changed",
                    "mutant_suffix",
                    "reason",
                ]
            )
            + "\n"
            )

            peptide_count = 0
            for rec in records:
                wt = lookup_protein(rec, proteins)
                mutant, method = reconstruct_mutant(rec, wt)
                key = variant_key(rec.chrom, rec.pos, rec.ref, rec.alt)
                ann_txt = annovar_text(annovar.get(key, {}))
                if not mutant:
                    protein_table.write(
                        "\t".join(
                            [
                                sample,
                                rec.var_id,
                                ann_txt,
                                rec.consequence,
                                rec.symbol or rec.gene,
                                rec.feature,
                                rec.protein_id,
                                rec.hgvsc,
                                rec.hgvsp,
                                rec.protein_position,
                                rec.amino_acids,
                                method,
                                "skipped",
                                wt or "",
                                "",
                                "",
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
                    manifest.write(
                        "\t".join(
                            [
                                sample,
                                rec.var_id,
                                ann_txt,
                                rec.symbol or rec.gene,
                                rec.feature,
                                rec.protein_id,
                                rec.consequence,
                                event_type(rec),
                                rec.protein_position,
                                rec.amino_acids,
                                rec.hgvsc,
                                rec.hgvsp,
                                method,
                                "skipped",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                wt or "",
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

                protein_table.write(
                    "\t".join(
                        [
                            sample,
                            rec.var_id,
                            ann_txt,
                            rec.consequence,
                            rec.symbol or rec.gene,
                            rec.feature,
                            rec.protein_id,
                            rec.hgvsc,
                            rec.hgvsp,
                            rec.protein_position,
                            rec.amino_acids,
                            mutant.method,
                            "reconstructed",
                            wt or "",
                            mutant.sequence,
                            mutant.wt_prefix,
                            mutant.wt_changed,
                            mutant.wt_suffix,
                            mutant.mut_prefix,
                            mutant.mut_changed,
                            mutant.mut_suffix,
                            "",
                        ]
                    )
                    + "\n"
                )

                peptides = all_mutant_peptides(rec, wt, mutant, mers)
                if not peptides:
                    manifest.write(
                        "\t".join(
                            [
                                sample,
                                rec.var_id,
                                ann_txt,
                                rec.symbol or rec.gene,
                                rec.feature,
                                rec.protein_id,
                                rec.consequence,
                                event_type(rec),
                                rec.protein_position,
                                rec.amino_acids,
                                rec.hgvsc,
                                rec.hgvsp,
                                mutant.method,
                                "no_novel_peptide",
                                "",
                                "",
                                "",
                                "",
                                "",
                                "",
                                wt or "",
                                mutant.sequence,
                                mutant.mut_prefix,
                                mutant.mut_changed,
                                mutant.mut_suffix,
                                "no_mutant_peptide_not_seen_in_wildtype_protein",
                            ]
                        )
                        + "\n"
                    )
                    continue

                for idx, pep in enumerate(peptides, start=1):
                    peptide_count += 1
                    peptide_id = f"{sample}|{rec.var_id}|{rec.symbol or rec.gene}|{rec.feature}|{pep.mer}mer|{pep.start}-{pep.end}|p{idx}"
                    fasta.write(f">{peptide_id}\n{pep.peptide}\n")
                    fasta_handles[pep.mer].write(f">{peptide_id}\n{pep.peptide}\n")
                    peptide_table.write(
                        "\t".join(
                            [
                                sample,
                                peptide_id,
                                rec.var_id,
                                rec.symbol or rec.gene,
                                rec.feature,
                                rec.protein_id,
                                rec.hgvsc,
                                rec.hgvsp,
                                rec.amino_acids,
                                str(pep.mer),
                                str(pep.start),
                                str(pep.end),
                                pep.wt_peptide,
                                pep.peptide,
                            ]
                        )
                        + "\n"
                    )
                manifest.write(
                    "\t".join(
                        [
                            sample,
                            rec.var_id,
                            ann_txt,
                            rec.symbol or rec.gene,
                            rec.feature,
                            rec.protein_id,
                            rec.consequence,
                            event_type(rec),
                            rec.protein_position,
                            rec.amino_acids,
                            rec.hgvsc,
                            rec.hgvsp,
                            mutant.method,
                            "emitted",
                            peptide_id,
                            str(pep.mer),
                            str(pep.start),
                            str(pep.end),
                            pep.wt_peptide,
                            pep.peptide,
                            wt or "",
                            mutant.sequence,
                            mutant.mut_prefix,
                            mutant.mut_changed,
                            mutant.mut_suffix,
                            "",
                        ]
                    )
                    + "\n"
                )
    finally:
        for handle in fasta_handles.values():
            handle.close()

    print(f"neoantigen_records={len(records)}")
    print(f"neoantigen_peptides={sum(1 for line in fasta_out.open() if line.startswith('>'))}")


def parse_lengths(raw: str) -> List[int]:
    lengths = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        if "-" in item:
            start, end = item.split("-", 1)
            lengths.extend(range(int(start), int(end) + 1))
        else:
            lengths.append(int(item))
    if not lengths:
        raise ValueError("No peptide lengths were provided")
    return sorted(set(lengths))


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vep-vcf", required=True, type=Path)
    parser.add_argument("--protein-fasta", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--protein-table", required=True, type=Path)
    parser.add_argument("--peptide-table", required=True, type=Path)
    parser.add_argument("--fasta-dir", required=True, type=Path)
    parser.add_argument("--annovar-txt", type=Path, help="Optional ANNOVAR multianno/exonic table to merge into neoantigen reports.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--lengths", default="8,9,10,11,12,13,14,15", help="Comma-separated peptide mer lengths to enumerate, for example 8,9,10,11 or 8-15.")
    parser.add_argument("--flank", default=30, type=int, help="Reserved for compatibility; all-mer generation does not use this value.")
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
        protein_table_out=args.protein_table,
        peptide_table_out=args.peptide_table,
        fasta_dir=args.fasta_dir,
        sample=args.sample,
        mers=parse_lengths(args.lengths),
        annovar=read_annovar_table(args.annovar_txt),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
