from amr_consensus.consensus import _overlaps, build_consensus
from amr_consensus.parser import Hit


def _make_hit(contig, start, end, gene, db, identity=95.0):
    return Hit(
        contig=contig, gene=gene, database=db, accession="X1", drugs=["Cefotaxime"],
        identity=identity, coverage=100.0, align_length=end - start, ref_length=end - start,
        contig_start=start, contig_end=end, strand="+", evalue=1e-50, bitscore=500.0,
    )


def test_overlaps_true_for_overlapping_ranges():
    a = _make_hit("c1", 100, 500, "geneA", "resfinder")
    b = _make_hit("c1", 120, 480, "geneA_variant", "card")
    assert _overlaps(a, b)


def test_overlaps_false_for_disjoint_ranges():
    a = _make_hit("c1", 100, 200, "geneA", "resfinder")
    b = _make_hit("c1", 5000, 5200, "geneB", "card")
    assert not _overlaps(a, b)


def test_build_consensus_merges_multi_db_agreement():
    hits = [
        _make_hit("c1", 100, 900, "blaCTX-M-15", "resfinder", identity=99.0),
        _make_hit("c1", 105, 895, "CTX-M-15", "card", identity=98.5),
    ]
    calls = build_consensus(hits, total_databases=2)
    assert len(calls) == 1
    call = calls[0]
    assert call.n_databases == 2
    assert call.confidence == "high"


def test_build_consensus_single_db_lower_confidence_than_multi_db():
    single = build_consensus(
        [_make_hit("c1", 100, 900, "geneA", "resfinder", identity=85.0)], total_databases=3
    )
    multi = build_consensus(
        [
            _make_hit("c1", 100, 900, "geneA", "resfinder", identity=85.0),
            _make_hit("c1", 100, 900, "geneA", "card", identity=85.0),
            _make_hit("c1", 100, 900, "geneA", "ncbi", identity=85.0),
        ],
        total_databases=3,
    )
    assert single[0].confidence_score < multi[0].confidence_score


def test_build_consensus_separates_distinct_loci_on_same_contig():
    hits = [
        _make_hit("c1", 100, 900, "geneA", "resfinder"),
        _make_hit("c1", 5000, 5800, "geneB", "resfinder"),
    ]
    calls = build_consensus(hits, total_databases=1)
    assert len(calls) == 2


def test_build_consensus_collapses_many_allele_variants_to_one_primary_gene():
    """Regression test: dozens of near-identical alleles hitting the same
    locus (common with blaCTX-M/mecA-style gene families) must collapse
    into ONE call with a single primary_gene, not flood gene_names into an
    unreadable blob."""
    hits = [
        _make_hit("c1", 100, 900, f"blaCTX-M-{i}", "resfinder", identity=90.0 + (i % 5))
        for i in range(50)
    ]
    calls = build_consensus(hits, total_databases=1)
    assert len(calls) == 1
    call = calls[0]
    assert call.n_allele_variants == 50
    assert call.primary_gene in {h.gene for h in hits}
    # primary_gene must be a single gene name, not a joined blob
    assert ";" not in call.primary_gene
