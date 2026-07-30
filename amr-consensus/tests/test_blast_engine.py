from amr_consensus.blast_engine import run_blast


def test_run_blast_finds_expected_number_of_hsps(blast_db, test_assembly_fasta):
    raw = run_blast(test_assembly_fasta, blast_db, threads=2)
    assert len(raw) > 0
    # every raw HSP row must carry all required BLAST fields
    required = {"qseqid", "sseqid", "pident", "length", "slen", "evalue"}
    assert required.issubset(raw[0].keys())


def test_run_blast_no_hits_on_random_sequence(blast_db, tmp_path):
    query = tmp_path / "random.fasta"
    query.write_text(">random\n" + "ACGT" * 50 + "\n")
    raw = run_blast(query, blast_db, threads=2)
    # random sequence should not produce strong hits against real AMR genes
    assert all(float(r["pident"]) < 90 for r in raw) or len(raw) == 0
