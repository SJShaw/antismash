# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import defaultdict
import json
import logging
from tempfile import NamedTemporaryFile
from typing import Dict, List

from antismash.common import fasta, path, subprocessing
from antismash.common.secmet import (
    CDSFeature,
    Record,
)

from .components import calculate_component_score
from .ordering import calculate_order_score
from .results import ClusterCompareResults

LEGACY_SCORES = {  # 1192
    "BGC0001192.1": (14, 35747.0/35747, 1),
    "BGC0001153.1": (14, 34042.0/35747, 2),
    "BGC0000408.1": (14, 33984.0/35747, 3),
    "BGC0001715.1": (12, 24511.0/35747, 4),
    "BGC0000400.1": (12, 17761.0/35747, 5),
    "BGC0000403.1": (12, 14707.0/35747, 6),
    "BGC0000449.1": (11, 16504.0/35747, 7),
    "BGC0001536.1": (11, 13776.0/35747, 8),
}

KIRRO_SCORES = {   #kirromycin/1070
    "BGC0001070.1": (132, 72791.0/72791.0, 1),
    "BGC0001807.1": (14, 6854.0/72791.0, 2),
    "BGC0001183.1": (10, 5345.0/72791.0, 3),
    "BGC0001230.1": (10, 3036.0/72791.0, 4),
    "BGC0001766.1": (10, 3036.0/72791.0, 5),
    "BGC0001657.1": (10, 1629.0/72791.0, 6),
    "BGC0001448.1": (9, 2491.0/72791.0, 7),
    "BGC0001567.1": (9, 1507.0/72791.0, 8),
    "BGC0001052.1": (8, 3201.0/72791.0, 9),
    "BGC0001975.1": (8, 2338.0/72791.0, 10),
}


class ReferenceScorer:
    def __init__(self, accession, data, hits_by_gene):
        self.accession = accession
        self.data = data
        self._raw_hits = hits_by_gene  # TODO remove
        self.hits_by_gene = trim_to_best_hit(hits_by_gene)
        self._identity = -1.
        self._order = -1.
        self._components = -1.

    @property
    def identity(self) -> float:
        if self._identity < 0:
            self._identity = calculate_identity_score(self.hits_by_gene)
        return self._identity

    @property
    def order(self) -> float:
        if self._order < 0:
            self._order = calculate_order_score(self.hits_by_gene, self.data)
        return self._order

    @property
    def components(self) -> float:
        if self._components < 0:
            self._components = calculate_component_score(self.hits_by_gene, self.data)
        return self._components

    def __repr__(self) -> str:
        return "ReferenceScorer(%s)" % (str(self))

    def __str__(self) -> str:
        return "%s: raw_id=%.2f, order=%.2f, comp=%.2f" % (
            self.accession,
            self.identity,
            self.order,
            self.components,
        )


def run(record: Record):
    if not record.get_regions():
        return ClusterCompareResults(record.id)

    # TODO handle custom databases
    with open(path.get_full_path(__file__, "data", "data.json")) as handle:
        ref_data = json.loads(handle.read())
    _, by_reference = find_diamond_matches(record, path.get_full_path(__file__, "data", "proteins.dmnd"))
    scores = {accession: ReferenceScorer(accession, ref_data[accession], hits) for accession, hits in by_reference.items()}

    max_id = max(1, max(score.identity for score in scores.values()))
    ranked_scores = sorted(scores.items(), key=lambda x: calculate_final_score(x[1], max_id), reverse=True)
    for acc, score in ranked_scores:
        print("%s: %.3f ..%d.. id=%.2f, order=%.2f, comp=%.2f" % (
            acc,
            calculate_final_score(score, max_id),
            len(score.hits_by_gene),
            score.identity / max_id,
            score.order,
            score.components,
        ), end="")
        prev_scores = KIRRO_SCORES if "AM746336" in record.id else LEGACY_SCORES
        if acc in prev_scores:
            print(", legacy=", prev_scores[acc])
        else:
            print()

    # TODO: set scoring/similarity to be by protocluster

    raise NotImplementedError()


def calculate_final_score(score, max_id):
    return ((score.identity / max_id) * score.order * score.components) ** (1/3)


def calculate_identity_score(hits):
    return sum(hit.identity_score for hit in hits.values())


def trim_to_best_hit(hits_by_reference_gene):
    pairs = {}
    for ref, hits in hits_by_reference_gene.items():
        for hit in hits:
            pairs[(ref, hit)] = hit.identity_score

    best_first = (i[0] for i in sorted(pairs.items(), key=lambda x: x[1], reverse=True))

    mapping = {}

    features = set()
    for ref, hit in best_first:
        if hit.cds in features or ref in mapping:
            continue
        mapping[ref] = hit
        features.add(hit.cds)
    assert 1 <= len(mapping) <= len(hits_by_reference_gene)
    return mapping


def find_diamond_matches(record: Record, database: str):
    """ Runs diamond, comparing all features in the given regions to the given database

        Arguments:
            regions: the regions to use features from
            database: the path of the database to compare to

        Returns:
    """  # TODO comment
    logging.info("Comparing regions to reference database")
    extra_args = [
        "--compress", "0",
        "--max-target-seqs", "10000",
        "--evalue", "1e-05",
        "--outfmt", "6",  # 6 is blast tabular format, just as in blastp
    ]
    features = []
    for region in record.get_regions():
        features.extend(region.cds_children)
    assert features

    with NamedTemporaryFile() as temp_file:
        temp_file.write(fasta.get_fasta_from_features(features, numeric_names=True).encode())
        raw = subprocessing.run_diamond_search(temp_file.name, database, mode="blastp", opts=extra_args)

    return blast_parse(raw, dict(enumerate(features)))


class Hit:
    def __init__(self, reference_cluster: str, reference_id: str, cds: CDSFeature,
                 percent_identity: int, blast_score: int, percent_coverage: int, evalue: float) -> None:
        self.reference_cluster = reference_cluster  # TODO actually reference record, need to split into regions/protoclusters
        self.reference_id = reference_id
        self.cds = cds
        self.percent_identity = percent_identity
        self.blast_score = blast_score
        self.percent_coverage = percent_coverage
        self.evalue = evalue

    @property
    def identity_score(self) -> float:
        return self.blast_score * self.percent_coverage

    def __repr__(self) -> str:
        return "%s(pid=%d,cov%d)" % (self.reference_id, self.percent_identity, self.percent_coverage)


def blast_parse(diamond_output: str, inputs_to_features,
                min_seq_coverage: float = 30., min_perc_identity: float = 25.):
    """ Parses blast output into a usable form, limiting to a single best hit
        for every query. Results can be further trimmed by minimum thresholds of
        both coverage and percent identity.

        Arguments:
            diamond_output: the output from diamond in blast format
            record: used to get all gene ids in the cluster, and used as a
                    backup to fetch sequence length if missing from seqlengths
            min_seq_coverage: the exclusive lower bound of sequence coverage for a match
            min_perc_identity: the exclusive lower bound of identity similarity for a match

        Returns:
            a tuple of
                a dictionary mapping query id to Query instance
                a dictionary mapping cluster number to
                    a list of Query instances from that cluster
    """  # TODO update args
    hits_by_cds = defaultdict(lambda: defaultdict(list))  # type: Dict[CDSFeature, Dict[str, List[Hit]]]
    hits_by_reference = defaultdict(lambda: defaultdict(list))
    for line in diamond_output.splitlines():
        hit = parse_hit(line.split("\t"), inputs_to_features)
        if not (hit.percent_identity >= min_perc_identity and hit.percent_coverage >= min_seq_coverage):
            continue
        hits_by_cds[hit.cds][hit.reference_cluster].append(hit)
        hits_by_reference[hit.reference_cluster][hit.reference_id].append(hit)
    return hits_by_cds, hits_by_reference


def parse_hit(parts: List[str], feature_mapping: Dict[str, CDSFeature]) -> None:  # TODO

    # 0. 	 qseqid 	 query (e.g., gene) sequence id
    # 1. 	 sseqid 	 subject (e.g., reference genome) sequence id
    # 2. 	 pident 	 percentage of identical matches
    # 3. 	 length 	 alignment length
    # 4. 	 mismatch 	 number of mismatches
    # 5. 	 gapopen 	 number of gap openings
    # 6. 	 qstart 	 start of alignment in query
    # 7. 	 qend 	 end of alignment in query
    # 8. 	 sstart 	 start of alignment in subject
    # 9. 	 send 	 end of alignment in subject
    # 10. 	 evalue 	 expect value
    # 11. 	 bitscore 	 bit score

    assert len(parts) == 12
    cds = feature_mapping[int(parts[0])]
    reference_cluster, reference_id = parts[1].split("|")
    perc_ident = int(float(parts[2]) + 0.5)  # the ceiling is enough precision
    evalue = float(parts[10])
    blastscore = int(float(parts[11]) + 0.5)  # the ceiling is enough precision
    perc_coverage = (float(parts[3]) / len(cds.translation)) * 100
    return Hit(reference_cluster, reference_id, cds, perc_ident, blastscore, perc_coverage, evalue)
