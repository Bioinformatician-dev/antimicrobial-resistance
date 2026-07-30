"""Normalize raw BLAST hits into a consistent Hit schema.

Reference headers across resfinder/card/ncbi/argannot/megares (as
republished by abricate) share the pattern::

    <db>~~~<gene_allele>~~~<accession>~~~<drug_class_list> <gene_name>

This module parses that header and combines it with BLAST coordinates to
produce one Hit per detected gene, filtered by identity/coverage.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Hit:
    contig: str
    gene: str
    database: str
    accession: str
    drugs: list[str]
    identity: float
    coverage: float
    align_length: int
    ref_length: int
    contig_start: int
    contig_end: int
    strand: str
    evalue: float
    bitscore: float

    def to_dict(self) -> dict:
        d = self.__dict__.copy()
        d["drugs"] = ";".join(self.drugs)
        return d


def _parse_subject_header(sseqid: str) -> tuple[str, str, str, list[str]]:
    """Parse a '<db>~~~<gene>~~~<accession>~~~<drugs>' style subject id."""
    parts = sseqid.split("~~~")
    if len(parts) >= 4:
        db, gene, accession, drug_field = parts[0], parts[1], parts[2], parts[3]
        drugs = [d for d in drug_field.split(";") if d]
        return db, gene, accession, drugs
    # Fallback for reference sets that don't follow the convention.
    return "unknown", sseqid, "", []


def parse_hits(
    raw_hits: list[dict],
    database_name: str,
    min_identity: float = 80.0,
    min_coverage: float = 80.0,
) -> list[Hit]:
    hits: list[Hit] = []
    for row in raw_hits:
        pident = float(row["pident"])
        align_len = int(row["length"])
        slen = int(row["slen"])
        coverage = (align_len / slen) * 100 if slen else 0.0

        if pident < min_identity or coverage < min_coverage:
            continue

        _, gene, accession, drugs = _parse_subject_header(row["sseqid"])

        qstart, qend = int(row["qstart"]), int(row["qend"])
        strand = "+" if qstart <= qend else "-"
        contig_start, contig_end = sorted((qstart, qend))

        hits.append(
            Hit(
                contig=row["qseqid"],
                gene=gene,
                database=database_name,
                accession=accession,
                drugs=drugs,
                identity=round(pident, 2),
                coverage=round(coverage, 2),
                align_length=align_len,
                ref_length=slen,
                contig_start=contig_start,
                contig_end=contig_end,
                strand=strand,
                evalue=float(row["evalue"]),
                bitscore=float(row["bitscore"]),
            )
        )
    return hits
