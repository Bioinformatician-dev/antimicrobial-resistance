#!/usr/bin/env bash
# Convenience wrapper around `amr-consensus setup-db`.
#
# Usage:
#   ./scripts/download_db.sh                 # downloads resfinder + card (defaults)
#   ./scripts/download_db.sh resfinder,card,ncbi,argannot,megares
set -euo pipefail

DBS="${1:-resfinder,card}"

echo "Downloading and indexing: ${DBS}"
echo "(Full databases; this can take a minute or two the first time.)"
amr-consensus setup-db --db "${DBS}"

echo ""
echo "Done. Databases are cached in ~/.amr-consensus/db/"
echo "Run 'amr-consensus run --input your_assembly.fasta --db ${DBS}' to scan a genome."
