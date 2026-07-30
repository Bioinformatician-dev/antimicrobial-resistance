from amr_consensus.blast_engine import run_blast
from amr_consensus.parser import _parse_subject_header, parse_hits


def test_parse_subject_header_standard_format():
    db, gene, acc, drugs = _parse_subject_header(
        "resfinder~~~blaCTX-M-150_1~~~KF769131~~~Cefotaxime;Ceftazidime"
    )
    assert db == "resfinder"
    assert gene == "blaCTX-M-150_1"
    assert acc == "KF769131"
    assert drugs == ["Cefotaxime", "Ceftazidime"]


def test_parse_subject_header_fallback():
    db, gene, _acc, drugs = _parse_subject_header("some_plain_id")
    assert db == "unknown"
    assert gene == "some_plain_id"
    assert drugs == []


def test_parse_hits_filters_by_identity_and_coverage(blast_db, test_assembly_fasta):
    raw = run_blast(test_assembly_fasta, blast_db, threads=2)

    lenient = parse_hits(raw, "resfinder", min_identity=50, min_coverage=50)
    strict = parse_hits(raw, "resfinder", min_identity=99, min_coverage=99)

    assert len(strict) <= len(lenient)
    # only the near-perfect mecA / blaCTX-M-150 hits should survive a 99% filter
    strict_genes = {h.gene for h in strict}
    assert "mecA_1" in strict_genes


def test_parse_hits_excludes_truncated_and_negative_control(blast_db, test_assembly_fasta):
    raw = run_blast(test_assembly_fasta, blast_db, threads=2)
    hits = parse_hits(raw, "resfinder", min_identity=80, min_coverage=80)
    contigs_with_hits = {h.contig for h in hits}

    # truncated blaTEM (60% coverage) must be filtered out by --mincov
    assert "contig_5_truncated_partial_blaTEM" not in contigs_with_hits
    # negative control contig must have no hits
    assert "contig_6_no_amr_genes" not in contigs_with_hits
    # the four "real" positives must be present
    assert contigs_with_hits == {
        "contig_1_plasmid_ctxm",
        "contig_2_chromosome_mecA",
        "contig_3_integron_sul1",
        "contig_4_borderline_tetM",
    }
