# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import defaultdict
import logging
import math
import os
from tempfile import NamedTemporaryFile
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Set,
    Tuple,
    TypeVar,
)

from antismash.common import fasta, path, subprocessing
from antismash.common.secmet import (
    CDSFeature,
    Record,
)
from antismash.common.secmet.features import CDSCollection, Region, Protocluster
from antismash.config import get_config

from .components import calculate_component_score, calculate_component_score_ref_in_query, calculate_component_score_query_in_ref, gather_query_components, Components
from .data_structures import load_data, Hit, ReferenceRecord, ReferenceRegion, ScoresByRegion, ScoresByProtocluster, ReferenceScorer, ReferenceProtocluster, ReferenceArea, HitsByReference
from .ordering import calculate_order_score, calculate_order_score_ref_in_query, calculate_order_score_query_in_ref
from .results import ClusterCompareResults, VariantResults

HitsByCDS = Dict[CDSFeature, Dict[str, List[Hit]]]
HitsByReferenceName = Dict[str, Dict[str, List[Hit]]]  # RefRegion -> RefCDSName -> Hits
Ranking = List[ReferenceScorer]
T = TypeVar("T")

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

KIRRO_SCORES = {   # kirromycin/1070
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


def filter_by_query_area(area: CDSCollection, hits_by_reference: HitsByReference) -> HitsByReference:
    hits_for_area = defaultdict(lambda: defaultdict(list))  # type: HitsByReference
    for ref_area, hits_by_ref_area in hits_by_reference.items():
        for ref_id, hits in hits_by_ref_area.items():
            for hit in hits:
                if hit.cds in area:
                    hits_for_area[ref_area][ref_id].append(hit)
    return hits_for_area


def filter_by_reference_protocluster(area: ReferenceProtocluster, hits_by_reference: HitsByReference) -> HitsByReference:
    result = defaultdict(dict)  # type: HitsByReference
    for ref_area, hits_by_ref_area in hits_by_reference.items():
        for ref_cds, hits in hits_by_ref_area.items():
            if ref_cds in area.cdses:
                result[ref_area][ref_cds] = hits
    return result


def rank_scores(scores: Ranking) -> Ranking:
    if not scores:
        return []
    max_id = max(1, max(score.raw_identity for score in scores))
    relevant = []
    for score in scores:
        assert score.raw_identity <= max_id
        score.calc_identity(max_id)
        if score.identity < 0.05:  # TODO: make customisable
            continue
        relevant.append(score)

    return sorted(relevant, reverse=True)


def score_query_area(query_area: CDSCollection, hits_by_reference: HitsByReference, query_components: Dict[CDSCollection, Components], order: Callable, component: Callable) -> Ranking:
    local_hits = filter_by_query_area(query_area, hits_by_reference)

    if not local_hits:
        return []

    best_hits = {ref: trim_to_best_hit(hits) for ref, hits in local_hits.items()}

    scores = []  # type: Ranking
    for reference_area, hits in best_hits.items():
        scorer = ReferenceScorer(hits, reference_area, query_area.cds_children, query_components[query_area],
                                 calculate_identity_score, order, component)
        if scorer.order < 0.05 or scorer.components < 0.05:  # TODO: make customisable
            continue
        scores.append(scorer)
    return rank_scores(scores)


def score_as_protoclusters(label: str, region: Region, hits_by_reference: HitsByReference, query_components: Dict[CDSCollection, Components], order: Callable, component: Callable) -> VariantResults:
    local_hits = filter_by_query_area(region, hits_by_reference)

    # rank the results for the full region independently, as they're the sum of protocluster scores
    reference_total_scores = defaultdict(float)  # type: Dict[ReferenceRegion, float]

    scores_by_protocluster = defaultdict(dict)  # type: Dict[int, Dict[ReferenceRegion, ReferenceScorer]]
    for protocluster in region.get_unique_protoclusters():
        for scorer in score_query_area(protocluster, local_hits, query_components, order, component):
            reference_total_scores[scorer.reference] += scorer.final_score
            scores_by_protocluster[protocluster.get_protocluster_number()][scorer.reference] = scorer

    region_ranking = sorted(reference_total_scores.items(), key=lambda x: x[1], reverse=True)
    region_ranking, scores_by_protocluster = apply_limits_to_rankings(region_ranking, scores_by_protocluster)
    return VariantResults(label, region_ranking, scores_by_protocluster, local_hits)


def score_as_region(label: str, region: Region, hits_by_reference: HitsByReference, query_components: Dict[CDSCollection, Components], order: Callable, component: Callable) -> VariantResults:
    local_hits = filter_by_query_area(region, hits_by_reference)
    ranking = score_query_area(region, local_hits, query_components, order, component)[:50]  # TODO: make customisable
    region_ranking = sorted(((scorer.reference, scorer.final_score) for scorer in ranking), key=lambda x: x[1])
    return VariantResults(label, region_ranking, ranking, local_hits)


def score_against_protoclusters(label: str, region: Region, hits_by_reference: HitsByReference, query_components: Dict[CDSCollection, Components], order: Callable, component: Callable) -> VariantResults:
    score_matrix = defaultdict(lambda: defaultdict(dict))  # type: Dict[int, Dict[ReferenceRegion, Dict[ReferenceProtocluster, ReferenceScorer]]]
    reference_best_scores = defaultdict(lambda: defaultdict(float))  # type: Dict[Protocluster, Dict[ReferenceRegion, float]]
    local_hits = filter_by_query_area(region, hits_by_reference)
    for ref_region in local_hits:
        hits_for_ref_region = {ref_region: local_hits[ref_region]}
        for ref_protocluster in ref_region.protoclusters:
            hits = filter_by_reference_protocluster(ref_protocluster, hits_for_ref_region)
            for protocluster in region.get_unique_protoclusters():
                for scorer in score_query_area(protocluster, hits, query_components, order, component):
                    assert scorer.reference is ref_region, "%s  %s" % (scorer.reference, ref_region)
                    reference_best_scores[protocluster][ref_region] = max(scorer.final_score, reference_best_scores[protocluster][ref_region])  # TODO: this might be misrepresentative
                    score_matrix[protocluster.get_protocluster_number()][ref_region][ref_protocluster] = scorer

    reference_total_scores = defaultdict(float)  # type: Dict[ReferenceRegion, float]
    for ref_region_to_score in reference_best_scores.values():
        for ref_region, score in ref_region_to_score.items():
            reference_total_scores[ref_region] += score
    region_ranking = sorted(reference_total_scores.items(), key=lambda x: x[1], reverse=True)  # TODO: misrepresentative continued
    region_ranking, score_matrix = apply_limits_to_rankings(region_ranking, score_matrix)
    return VariantResults(label, region_ranking, score_matrix, local_hits)  # TODO: hit building isn't any good here


ValuesByReference = Dict[int, Dict[ReferenceRegion, Any]]


def apply_limits_to_rankings(region_ranking: ScoresByRegion, score_matrix: ValuesByReference,
                             max_reference_regions: int = 50) -> Tuple[ScoresByRegion, ValuesByReference]:  # TODO use a const/config option
    region_ranking = region_ranking[:max_reference_regions]
    regions = {rank[0] for rank in region_ranking}
    limited_matrix = {}  # type: ValuesByReference
    for key, pairs_by_ref_region in score_matrix.items():
        limited_matrix[key] = {ref: value for ref, value in pairs_by_ref_region.items() if ref in regions}
    return region_ranking, limited_matrix


def run(record: Record) -> ClusterCompareResults:
    if not record.get_regions():
        return ClusterCompareResults(record.id, {})
    logging.debug("loading reference database")
    data_dir = os.path.join(get_config().database_dir, "cluster_compare")
    references = load_data(os.path.join(data_dir, "data.json"))  # TODO: delay and restrict to only those with hits
    logging.debug("reference database loaded")
    # TODO keep hits_by_cds for performance
    _, hits_by_name = find_diamond_matches(record, os.path.join(data_dir, "proteins.dmnd"))
    hits = convert_to_references(hits_by_name, references)
    import time
    fastest = 1000.
    slowest = 0.

    results_per_region = {}  # type: Dict[int, Dict[str, VariantResults]]
    for region in record.get_regions():
        query_components = {}  # type: Dict[CDSCollection, Components]
        query_components[region] = gather_query_components(region.cds_children)
        for protocluster in region.get_unique_protoclusters():
            query_components[protocluster] = gather_query_components(protocluster.cds_children)
        start = time.time()
        logging.debug("analysing %s region: %s", record, region)
        # TODO: keep only results that are non-zero
        region_results = {}
        region_results["ALL_TO_ALL_best"] = score_as_region("ALL_TO_ALL_best", region, hits, query_components, calculate_order_score, calculate_component_score)
        region_results["ALL_TO_ALL_QiR"] = score_as_region("ALL_TO_ALL_QiR", region, hits, query_components, calculate_order_score_query_in_ref, calculate_component_score_query_in_ref)
        region_results["ALL_TO_ALL_RiQ"] = score_as_region("ALL_TO_ALL_RiQ", region, hits, query_components, calculate_order_score_ref_in_query, calculate_component_score_ref_in_query)

        region_results["PC_TO_ALL_best"] = score_as_protoclusters("PC_TO_ALL_best", region, hits, query_components, calculate_order_score, calculate_component_score)
        region_results["PC_TO_ALL_QiR"] = score_as_protoclusters("PC_TO_ALL_QiR", region, hits, query_components, calculate_order_score_query_in_ref, calculate_component_score_query_in_ref)
        region_results["PC_TO_ALL_RiQ"] = score_as_protoclusters("PC_TO_ALL_RiQ", region, hits, query_components, calculate_order_score_ref_in_query, calculate_component_score_ref_in_query)

        region_results["PC_TO_PC_best"] = score_against_protoclusters("PC_TO_ALL_best", region, hits, query_components, calculate_order_score, calculate_component_score)
        region_results["PC_TO_PC_QiR"] = score_against_protoclusters("PC_TO_ALL_QiR", region, hits, query_components, calculate_order_score_query_in_ref, calculate_component_score_query_in_ref)
        region_results["PC_TO_PC_RiQ"] = score_against_protoclusters("PC_TO_ALL_RiQ", region, hits, query_components, calculate_order_score_ref_in_query, calculate_component_score_ref_in_query)

        results_per_region[region.get_region_number()] = region_results
        total = time.time() - start
        fastest = min(fastest, total)
        slowest = max(slowest, total)

    logging.critical("fastest region time: %s, slowest: %s", fastest, slowest)
    return ClusterCompareResults(record.id, results_per_region)


def calculate_identity_score(score: float, max_score: float) -> float:
    normalised = score / max_score
    # avoid division by 0
    if normalised > 1. - 1e-8:
        return normalised
    # logit function, shifted from x range of (-6, 6) to (0, 1)
    # tiny values can result in small negatives
    # very high normalised values can result in just over 1 (e.g. .997 -> 1.003)
    result = min(1, max(0, math.log(normalised/(1-normalised)) / 12 + 0.5))
    assert 0 <= result <= 1, normalised
    return result


def trim_to_best_hit(hits_by_reference_gene: Dict[str, List[Hit]]) -> Dict[str, Hit]:
    pairs = {}
    for ref, hits in hits_by_reference_gene.items():
        for hit in hits:
            pairs[(ref, hit)] = hit.identity_score

    best_first = (i[0] for i in sorted(pairs.items(), key=lambda x: x[1], reverse=True))

    mapping = {}  # type: Dict[str, Hit]

    features = set()  # type: Set[CDSFeature]
    for ref, hit in best_first:
        if hit.cds in features or ref in mapping:
            continue
        mapping[ref] = hit
        assert hit.cds not in features
        features.add(hit.cds)
    assert 1 <= len(mapping) <= len(hits_by_reference_gene)
    return mapping


def convert_to_references(hits_by_name: HitsByReferenceName, references: Dict[str, ReferenceRecord]) -> HitsByReference:
    results = defaultdict(dict)  # type: HitsByReference
    for name, hits in hits_by_name.items():
        reference_record = references[name]
        for cds_id, cds_hits in hits.items():
            cds_name = reference_record.cds_mapping[str(cds_id)]
            for region in reference_record.regions:
                assert region.accession == reference_record.accession
                if cds_name in region.cdses:
                    results[region][cds_name] = cds_hits
    return results


def find_diamond_matches(record: Record, database: str) -> Tuple[HitsByCDS, HitsByReferenceName]:
    """ Runs diamond, comparing all features in the given regions to the given database

        Arguments:
            regions: the regions to use features from
            database: the path of the database to compare to

        Returns:
    """  # TODO comment
    # TODO gather by region in question, not record
    logging.info("Comparing regions to reference database")
    extra_args = [
        "--compress", "0",
        "--max-target-seqs", "10000",
        "--evalue", "1e-05",
        "--outfmt", "6",  # 6 is blast tabular format, just as in blastp
    ]
    features = record.get_cds_features_within_regions()

    with NamedTemporaryFile() as temp_file:
        temp_file.write(fasta.get_fasta_from_features(features, numeric_names=True).encode())
        raw = subprocessing.run_diamond_search(temp_file.name, database, mode="blastp", opts=extra_args)

    return blast_parse(raw, {i: feature for i, feature in enumerate(features)})


def blast_parse(diamond_output: str, inputs_to_features: Dict[int, CDSFeature],
                min_seq_coverage: float = 30., min_perc_identity: float = 25.) -> Tuple[HitsByCDS, HitsByReferenceName]:
    """ Parses blast output into a usable form, limiting to a single best hit
        for every query. Results can be further trimmed by minimum thresholds of
        both coverage and percent identity.

        Arguments:
            diamond_output: the output from diamond in blast format
            inputs_to_features: a mapping of integer used in fasta supplied to diamond to original CDS feature
            min_seq_coverage: the exclusive lower bound of sequence coverage for a match
            min_perc_identity: the exclusive lower bound of identity similarity for a match

        Returns:
            a tuple of
                a dictionary mapping query id to Query instance
                a dictionary mapping cluster number to
                    a list of Query instances from that cluster
    """  # TODO update args
    hits_by_cds = defaultdict(lambda: defaultdict(list))  # type: HitsByCDS
    hits_by_reference = defaultdict(lambda: defaultdict(list))  # type: HitsByReferenceName
    for line in diamond_output.splitlines():
        hit = parse_hit(line.split("\t"), inputs_to_features)
        if not (hit.percent_identity >= min_perc_identity and hit.percent_coverage >= min_seq_coverage):
            continue
        hits_by_cds[hit.cds][hit.reference_record].append(hit)
        hits_by_reference[hit.reference_record][hit.reference_id].append(hit)
    return hits_by_cds, hits_by_reference


def parse_hit(parts: List[str], feature_mapping: Dict[int, CDSFeature]) -> Hit:

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
    reference_record, reference_id = parts[1].split("|")
    perc_ident = int(float(parts[2]) + 0.5)  # the ceiling is enough precision
    evalue = float(parts[10])
    blastscore = int(float(parts[11]) + 0.5)  # the ceiling is enough precision
    perc_coverage = (float(parts[3]) / len(cds.translation)) * 100
    return Hit(reference_record, reference_id, cds, perc_ident, blastscore, perc_coverage, evalue)
