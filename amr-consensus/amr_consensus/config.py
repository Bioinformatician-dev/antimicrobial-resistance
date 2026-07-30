"""Central configuration: database sources, default thresholds."""
from dataclasses import dataclass, field
from pathlib import Path

# Curated, actively-maintained AMR gene databases we know how to fetch and
# parse. Sequences are pulled from the abricate project's GitHub mirror,
# which republishes ResFinder, CARD, NCBI AMRFinderPlus, ARG-ANNOT and
# MEGARes reference sets in a single consistent FASTA format.
DB_SOURCES = {
    "resfinder": "https://raw.githubusercontent.com/tseemann/abricate/master/db/resfinder/sequences",
    "card": "https://raw.githubusercontent.com/tseemann/abricate/master/db/card/sequences",
    "ncbi": "https://raw.githubusercontent.com/tseemann/abricate/master/db/ncbi/sequences",
    "argannot": "https://raw.githubusercontent.com/tseemann/abricate/master/db/argannot/sequences",
    "megares": "https://raw.githubusercontent.com/tseemann/abricate/master/db/megares/sequences",
}

DEFAULT_DATA_DIR = Path.home() / ".amr-consensus" / "db"

# Detection thresholds, matching community-standard defaults used by
# ABRicate / ResFinder for BLAST-based gene-presence calls.
DEFAULT_MIN_IDENTITY = 80.0   # percent identity over the aligned region
DEFAULT_MIN_COVERAGE = 80.0   # percent of the reference gene length covered


@dataclass
class RunConfig:
    input_fasta: Path
    databases: list[str]
    outdir: Path
    min_identity: float = DEFAULT_MIN_IDENTITY
    min_coverage: float = DEFAULT_MIN_COVERAGE
    threads: int = 4
    data_dir: Path = field(default_factory=lambda: DEFAULT_DATA_DIR)
