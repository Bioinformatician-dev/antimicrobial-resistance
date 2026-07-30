from typer.testing import CliRunner

from amr_consensus.cli import app

runner = CliRunner()


def test_cli_run_offline_mode(tmp_path, mini_db_fasta, test_assembly_fasta):
    """End-to-end CLI run using --custom-db-fasta, so it needs no network access."""
    outdir = tmp_path / "results"
    result = runner.invoke(
        app,
        [
            "run",
            "--input", str(test_assembly_fasta),
            "--custom-db-fasta", str(mini_db_fasta),
            "--outdir", str(outdir),
        ],
    )
    assert result.exit_code == 0, result.output
    assert (outdir / "amr_consensus_report.tsv").exists()
    assert (outdir / "amr_consensus_report.json").exists()
    assert (outdir / "amr_consensus_report.html").exists()
    assert "mecA" in result.output


def test_cli_list_dbs():
    result = runner.invoke(app, ["list-dbs"])
    assert result.exit_code == 0
    assert "resfinder" in result.output
    assert "card" in result.output
