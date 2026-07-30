# Contributing to amr-consensus

Thanks for considering a contribution! This project is small and focused —
here's how to get productive quickly.

## Development setup

```bash
git clone https://github.com/YOUR_USERNAME/amr-consensus.git
cd amr-consensus
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"

# BLAST+ is required and not installed via pip:
sudo apt-get install ncbi-blast+   # Debian/Ubuntu
brew install blast                 # macOS
```

## Running tests

```bash
pytest tests/ -v
```

Tests run entirely offline against the bundled fixtures in `tests/data/`
(a 6-gene subset of ResFinder and a synthetic assembly built from it) — no
network access or full database download needed.

## Code style

We use [ruff](https://docs.astral.sh/ruff/) for linting:

```bash
ruff check amr_consensus/ tests/
```

CI will fail on lint errors, so please run this before opening a PR.

## Adding a new reference database

Reference databases live in `amr_consensus/config.py::DB_SOURCES`. To add
one:

1. Confirm the source publishes a FASTA with headers amr-consensus can
   parse (`db~~~gene~~~accession~~~drug_list` — see `parser.py`), or add a
   new header-parsing branch if the format differs.
2. Add the download URL to `DB_SOURCES`.
3. Add a test in `tests/test_parser.py` covering the new header format.

## Reporting bugs / requesting features

Please open a GitHub issue with:
- The command you ran
- Expected vs. actual output
- Your `amr-consensus --version` and BLAST+ version (`blastn -version`)

## Pull requests

- Keep PRs focused on one change.
- Add or update tests for any behavior change.
- Update the README if you change CLI flags or output format.
