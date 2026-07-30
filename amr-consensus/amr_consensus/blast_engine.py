"""Thin wrapper around NCBI BLAST+ for gene-presence detection.

We use blastn with a permissive task ('blastn', not 'megablast') so that
divergent gene variants (e.g. tet(M) alleles differing by several percent)
are still recovered, then filter downstream by identity/coverage. This
mirrors the approach used by ABRicate and ResFinder.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

# Column order for BLAST outfmt 6, extended with qlen/slen so we can compute
# reference-gene coverage without a second lookup.
BLAST_OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore qlen slen"
)
BLAST_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
]


class BlastError(RuntimeError):
    pass


def run_blast(query_fasta: Path, blast_db_prefix: Path, threads: int = 4) -> list[dict]:
    """Run blastn of ``query_fasta`` against a prebuilt BLAST db.

    Returns a list of raw hit dicts (one per HSP), un-filtered.
    """
    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", str(blast_db_prefix),
        "-outfmt", BLAST_OUTFMT,
        "-num_threads", str(threads),
        "-task", "blastn",
        "-dust", "no",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise BlastError(f"blastn failed against {blast_db_prefix}:\n{result.stderr}")

    hits = []
    for line in result.stdout.strip().splitlines():
        if not line:
            continue
        values = line.split("\t")
        row = dict(zip(BLAST_FIELDS, values))
        hits.append(row)
    return hits
