"""Reconcile hits from multiple AMR databases into consensus calls.

Different databases (ResFinder, CARD, NCBI, ARG-ANNOT, MEGARes) curate
overlapping but non-identical sets of reference alleles. The same gene on
the same contig is frequently detected by more than one database under
slightly different names or coordinates. This module clusters same-locus
hits (by contig + overlapping coordinates) and produces one consensus call
per locus, with a confidence score reflecting cross-database agreement.
"""
from __future__ import annotations

from dataclasses import dataclass

from .parser import Hit


# Confidence tiers, roughly following the logic: more independent databases
# agreeing on a locus, and higher identity, both raise confidence that the
# call reflects a real resistance determinant rather than a spurious
# low-identity match.
def _confidence_label(score: float) -> str:
    if score >= 0.85:
        return "high"
    if score >= 0.6:
        return "medium"
    return "low"


@dataclass
class ConsensusCall:
    contig: str
    locus_start: int
    locus_end: int
    primary_gene: str
    gene_names: list[str]
    n_allele_variants: int
    supporting_databases: list[str]
    n_databases: int
    best_identity: float
    mean_identity: float
    best_coverage: float
    drugs: list[str]
    confidence_score: float
    confidence: str

    def to_dict(self) -> dict:
        d = self.__dict__.copy()
        d["gene_names"] = ";".join(sorted(set(self.gene_names)))
        d["supporting_databases"] = ";".join(sorted(set(self.supporting_databases)))
        d["drugs"] = ";".join(sorted(set(self.drugs)))
        return d


def _overlaps(a: Hit, b: Hit, min_frac: float = 0.5) -> bool:
    lo = max(a.contig_start, b.contig_start)
    hi = min(a.contig_end, b.contig_end)
    overlap = max(0, hi - lo)
    shorter = min(a.contig_end - a.contig_start, b.contig_end - b.contig_start)
    if shorter <= 0:
        return False
    return (overlap / shorter) >= min_frac


def build_consensus(hits: list[Hit], total_databases: int) -> list[ConsensusCall]:
    """Cluster hits into per-locus consensus calls.

    ``total_databases`` is the number of databases actually queried in this
    run, used to normalize the confidence score (agreement across all
    queried databases -> higher confidence than agreement across a subset).
    """
    by_contig: dict[str, list[Hit]] = {}
    for h in hits:
        by_contig.setdefault(h.contig, []).append(h)

    calls: list[ConsensusCall] = []
    for contig, contig_hits in by_contig.items():
        contig_hits.sort(key=lambda h: h.contig_start)
        clusters: list[list[Hit]] = []
        for h in contig_hits:
            placed = False
            for cluster in clusters:
                if any(_overlaps(h, member) for member in cluster):
                    cluster.append(h)
                    placed = True
                    break
            if not placed:
                clusters.append([h])

        for cluster in clusters:
            dbs = {h.database for h in cluster}
            genes = [h.gene for h in cluster]
            drugs = [d for h in cluster for d in h.drugs]
            identities = [h.identity for h in cluster]
            best_hit = max(cluster, key=lambda h: h.bitscore)

            db_agreement = len(dbs) / total_databases if total_databases else 1.0
            identity_factor = best_hit.identity / 100.0
            score = round(0.6 * db_agreement + 0.4 * identity_factor, 3)

            calls.append(
                ConsensusCall(
                    contig=contig,
                    locus_start=min(h.contig_start for h in cluster),
                    locus_end=max(h.contig_end for h in cluster),
                    primary_gene=best_hit.gene,
                    gene_names=genes,
                    n_allele_variants=len(set(genes)),
                    supporting_databases=sorted(dbs),
                    n_databases=len(dbs),
                    best_identity=max(identities),
                    mean_identity=round(sum(identities) / len(identities), 2),
                    best_coverage=best_hit.coverage,
                    drugs=drugs,
                    confidence_score=score,
                    confidence=_confidence_label(score),
                )
            )

    calls.sort(key=lambda c: (c.contig, c.locus_start))
    return calls
