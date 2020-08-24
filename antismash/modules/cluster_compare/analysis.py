# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import defaultdict
import json
import logging
import math
from tempfile import NamedTemporaryFile
from typing import (
    Any,
    Dict,
    List,
    Set,
    Tuple,
)

from antismash.common import fasta, path, subprocessing
from antismash.common.secmet import (
    CDSFeature,
    Record,
)
from antismash.common.secmet.features import CDSCollection

from .components import calculate_component_score
from .data_structures import load_data, Hit, RawRecord, RawRegion, ScoresByRegion, ScoresByProtocluster, ReferenceScorer
from .ordering import calculate_order_score
from .results import ClusterCompareResults, VariantResults

HitsByCDS = Dict[CDSFeature, Dict[str, List[Hit]]]
HitsByReference = Dict[str, Dict[str, List[Hit]]]

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


def filter_by_area(area: CDSCollection, hits_by_reference: Dict[str, Dict[str, List[Hit]]]) -> Dict[str, Dict[str, List[Hit]]]:
    region_cds_names = {cds.get_name() for cds in area.cds_children}
    hits_for_area = defaultdict(lambda: defaultdict(list))  # type: Dict[str, Dict[str, List[Hit]]]
    for cluster, hits_by_cluster in hits_by_reference.items():
        for ref_id, hits in hits_by_cluster.items():
            for hit in hits:
                if hit.cds.get_name() in region_cds_names:
                    hits_for_area[cluster][ref_id].append(hit)
    return hits_for_area


def run_PC_TO_ALL(record: Record, ref_data, processed, by_reference) -> VariantResults:
    scores_by_region = {}  # type: ScoresByRegion
    scores_by_protocluster = {}
    hits_by_region = {}

    for region in record.get_regions():
        hits_for_region = filter_by_area(region, by_reference)
        best_hits_for_region = {ref: trim_to_best_hit(hits) for ref, hits in hits_for_region.items()}
        if not best_hits_for_region:
            continue
        scores_within_region = defaultdict(list)  # type: ScoresByProtocluster
        for protocluster in region.get_unique_protoclusters():
            hits_for_protocluster = filter_by_area(protocluster, hits_for_region)
            scores = {accession: ReferenceScorer(accession, ref_data[accession], trim_to_best_hit(hits), processed[accession], protocluster.cds_children,
                                                 calculate_identity_score, calculate_order_score, calculate_component_score)
                      for accession, hits in hits_for_protocluster.items()}
            if not scores:
                continue
            max_id = max(1, max(score.raw_identity for score in scores.values()))
            for score in scores.values():
                assert score.raw_identity <= max_id
                score.calc_identity(max_id)
            scores_by_protocluster[protocluster.get_protocluster_number()] = scores
            scores_within_region[region.get_region_number()].extend(scores.items())

        # rank the results for a region differently
        region_ranking = defaultdict(float)  # type: Dict[str, float]
        for ref_accession, score in scores_within_region[region.get_region_number()]:
            region_ranking[ref_accession] += score.final_score
        ranking = sorted(region_ranking.items(), key=lambda x: x[1], reverse=True)
        scores_by_region[region.get_region_number()] = [(acc, (score, processed[acc].regions[0])) for acc, score in ranking]  # TODO handle multiple regions in a ref record
        hits_by_region[region.get_region_number()] = best_hits_for_region

    return VariantResults("PC_TO_ALL", scores_by_region, scores_by_protocluster, hits_by_region)


def run(record: Record) -> ClusterCompareResults:
    results = ClusterCompareResults(record.id)
    if not record.get_regions():
        return results

    # TODO handle custom databases
    with open(path.get_full_path(__file__, "data", "data.json")) as handle:
        ref_data = json.loads(handle.read())
    processed = load_data(path.get_full_path(__file__, "data", "data.json"))
    _, by_reference = find_diamond_matches(record, path.get_full_path(__file__, "data", "proteins.dmnd"))

    results.add_variant_results(run_PC_TO_ALL(record, ref_data, processed, by_reference))
    return results


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


def find_diamond_matches(record: Record, database: str) -> Tuple[Dict[CDSFeature, Dict[str, List[Hit]]],
                                                                 Dict[str, Dict[str, List[Hit]]]]:
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
    features = []  # type: List[CDSFeature]
    for region in record.get_regions():
        features.extend(region.cds_children)
    assert features

    with NamedTemporaryFile() as temp_file:
        temp_file.write(fasta.get_fasta_from_features(features, numeric_names=True).encode())
        raw = subprocessing.run_diamond_search(temp_file.name, database, mode="blastp", opts=extra_args)

    return blast_parse(raw, {i: feature for i, feature in enumerate(features)})


def blast_parse(diamond_output: str, inputs_to_features: Dict[int, CDSFeature],
                min_seq_coverage: float = 30., min_perc_identity: float = 25.) -> Tuple[HitsByCDS, HitsByReference]:
    """ Parses blast output into a usable form, limiting to a single best hit
        for every query. Results can be further trimmed by minimum thresholds of
        both coverage and percent identity.

        Arguments:
            diamond_output: the output from diamond in blast format
            min_seq_coverage: the exclusive lower bound of sequence coverage for a match
            min_perc_identity: the exclusive lower bound of identity similarity for a match

        Returns:
            a tuple of
                a dictionary mapping query id to Query instance
                a dictionary mapping cluster number to
                    a list of Query instances from that cluster
    """  # TODO update args
    hits_by_cds = defaultdict(lambda: defaultdict(list))  # type: HitsByCDS
    hits_by_reference = defaultdict(lambda: defaultdict(list))  # type: HitsByReference
    for line in diamond_output.splitlines():
        hit = parse_hit(line.split("\t"), inputs_to_features)
        if not (hit.percent_identity >= min_perc_identity and hit.percent_coverage >= min_seq_coverage):
            continue
        hits_by_cds[hit.cds][hit.reference_cluster].append(hit)
        hits_by_reference[hit.reference_cluster][hit.reference_id].append(hit)
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
    reference_cluster, reference_id = parts[1].split("|")
    perc_ident = int(float(parts[2]) + 0.5)  # the ceiling is enough precision
    evalue = float(parts[10])
    blastscore = int(float(parts[11]) + 0.5)  # the ceiling is enough precision
    perc_coverage = (float(parts[3]) / len(cds.translation)) * 100
    return Hit(reference_cluster, reference_id, cds, perc_ident, blastscore, perc_coverage, evalue)
