"""Fetch curated AMR reference databases and build local BLAST indexes."""
from __future__ import annotations

import subprocess
import urllib.request
from pathlib import Path

from .config import DB_SOURCES


class DatabaseError(RuntimeError):
    pass


def fetch_database(name: str, data_dir: Path, force: bool = False) -> Path:
    """Download a reference FASTA for ``name`` and build a nucleotide BLAST db.

    Returns the path to the BLAST database prefix (usable directly with
    ``blastn -db``).
    """
    if name not in DB_SOURCES:
        raise DatabaseError(
            f"Unknown database '{name}'. Available: {', '.join(sorted(DB_SOURCES))}"
        )

    db_dir = data_dir / name
    db_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = db_dir / "sequences.fasta"
    blast_prefix = db_dir / "sequences"

    if not fasta_path.exists() or force:
        url = DB_SOURCES[name]
        try:
            urllib.request.urlretrieve(url, fasta_path)
        except Exception as exc:
            raise DatabaseError(f"Failed to download '{name}' from {url}: {exc}") from exc

    index_marker = db_dir / "sequences.nsq"
    if not index_marker.exists() or force:
        _makeblastdb(fasta_path, blast_prefix)

    return blast_prefix


def build_local_database(fasta_path: Path, db_dir: Path, name: str = "custom") -> Path:
    """Build a BLAST database from a local/offline FASTA file (no download)."""
    db_dir = db_dir / name
    db_dir.mkdir(parents=True, exist_ok=True)
    blast_prefix = db_dir / "sequences"
    _makeblastdb(fasta_path, blast_prefix)
    return blast_prefix


def _makeblastdb(fasta_path: Path, out_prefix: Path) -> None:
    cmd = [
        "makeblastdb",
        "-in", str(fasta_path),
        "-dbtype", "nucl",
        "-out", str(out_prefix),
        "-title", out_prefix.parent.name,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise DatabaseError(
            f"makeblastdb failed for {fasta_path}:\n{result.stderr}"
        )


def available_databases() -> list[str]:
    return sorted(DB_SOURCES)
