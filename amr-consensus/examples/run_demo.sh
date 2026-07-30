#!/usr/bin/env bash
# Runs the full pipeline on the bundled demo data and regenerates
# examples/expected_output/ -- a self-contained "data to result" demo that
# needs no network access (uses --custom-db-fasta, no database download).
#
# Usage: ./examples/run_demo.sh
set -euo pipefail

cd "$(dirname "$0")/.."

echo "Input:     examples/demo_assembly.fasta   (6 synthetic contigs)"
echo "Reference: examples/demo_reference_db.fasta (6 real ResFinder genes)"
echo ""

amr-consensus run \
  --input examples/demo_assembly.fasta \
  --custom-db-fasta examples/demo_reference_db.fasta \
  --outdir examples/expected_output

echo ""
echo "Done. Open examples/expected_output/amr_consensus_report.html in a browser to view the result."
