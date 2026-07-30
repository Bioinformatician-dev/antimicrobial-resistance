"""Command-line interface for amr-consensus."""
from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from .blast_engine import run_blast
from .config import DEFAULT_DATA_DIR, DEFAULT_MIN_COVERAGE, DEFAULT_MIN_IDENTITY
from .consensus import build_consensus
from .db_manager import available_databases, build_local_database, fetch_database
from .parser import parse_hits
from .report import write_reports

app = typer.Typer(
    name="amr-consensus",
    help="Detect antimicrobial resistance genes by reconciling calls across multiple reference databases.",
    add_completion=False,
)
console = Console()


@app.command()
def setup_db(
    databases: str = typer.Option(
        "resfinder,card", "--db", help="Comma-separated list of databases to download."
    ),
    data_dir: Path = typer.Option(DEFAULT_DATA_DIR, help="Where to store downloaded databases."),
    force: bool = typer.Option(False, help="Re-download and rebuild even if already present."),
):
    """Download curated reference databases and build local BLAST indexes."""
    db_list = [d.strip() for d in databases.split(",") if d.strip()]
    for db in db_list:
        console.print(f"[bold]Fetching[/bold] {db} ...")
        try:
            prefix = fetch_database(db, data_dir, force=force)
            console.print(f"  [green]ready[/green] -> {prefix}")
        except Exception as exc:  # noqa: BLE001
            console.print(f"  [red]failed[/red]: {exc}")
            raise typer.Exit(1)


@app.command()
def list_dbs():
    """List supported reference databases."""
    for name in available_databases():
        console.print(f"  - {name}")


@app.command()
def run(
    input_fasta: Path = typer.Option(..., "--input", "-i", exists=True, help="Assembly/contigs FASTA to scan."),
    databases: str = typer.Option(
        "resfinder,card", "--db", help="Comma-separated list of databases to query."
    ),
    outdir: Path = typer.Option(Path("amr_results"), "--outdir", "-o", help="Output directory."),
    min_identity: float = typer.Option(DEFAULT_MIN_IDENTITY, help="Minimum percent identity."),
    min_coverage: float = typer.Option(DEFAULT_MIN_COVERAGE, help="Minimum percent reference coverage."),
    threads: int = typer.Option(4, help="Threads for blastn."),
    data_dir: Path = typer.Option(DEFAULT_DATA_DIR, help="Directory holding prebuilt databases."),
    custom_db_fasta: Path = typer.Option(
        None, help="Optional local reference FASTA to use instead of downloading (offline mode)."
    ),
):
    """Scan an assembly/contigs FASTA for AMR genes and produce a consensus report."""
    db_list = [d.strip() for d in databases.split(",") if d.strip()]
    all_hits = []

    if custom_db_fasta:
        console.print(f"[bold]Building local database[/bold] from {custom_db_fasta}")
        prefix = build_local_database(custom_db_fasta, outdir / "_db", name="custom")
        db_list = ["custom"]
        raw = run_blast(input_fasta, prefix, threads=threads)
        all_hits.extend(parse_hits(raw, "custom", min_identity, min_coverage))
    else:
        for db in db_list:
            console.print(f"[bold]Querying[/bold] {db} ...")
            try:
                prefix = fetch_database(db, data_dir)
            except Exception as exc:  # noqa: BLE001
                console.print(f"  [red]skipping {db}[/red]: {exc}")
                continue
            raw = run_blast(input_fasta, prefix, threads=threads)
            hits = parse_hits(raw, db, min_identity, min_coverage)
            console.print(f"  {len(hits)} hits >= {min_identity}% id / {min_coverage}% cov")
            all_hits.extend(hits)

    calls = build_consensus(all_hits, total_databases=len(db_list))
    paths = write_reports(
        calls, outdir, sample_name=input_fasta.stem,
        databases=db_list, min_identity=min_identity, min_coverage=min_coverage,
    )

    _print_summary(calls, paths)


def _print_summary(calls, paths) -> None:
    table = Table(title="AMR Consensus Calls")
    table.add_column("Contig")
    table.add_column("Gene(s)")
    table.add_column("Confidence")
    table.add_column("Identity")
    table.add_column("Coverage")
    table.add_column("Databases")
    table.add_column("Drug classes")

    for c in calls:
        style = {"high": "green", "medium": "yellow", "low": "red"}[c.confidence]
        gene_display = c.primary_gene + (f" (+{c.n_allele_variants - 1})" if c.n_allele_variants > 1 else "")
        table.add_row(
            c.contig,
            gene_display,
            f"[{style}]{c.confidence}[/{style}] ({c.confidence_score})",
            f"{c.best_identity}%",
            f"{c.best_coverage}%",
            ";".join(sorted(set(c.supporting_databases))),
            ";".join(sorted({d for d in c.drugs if d})) or "-",
        )

    console.print(table)
    console.print("\n[bold]Reports written:[/bold]")
    for kind, path in paths.items():
        console.print(f"  {kind.upper():5s} -> {path}")


if __name__ == "__main__":
    app()
