# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from typing import (
    Any,
    Dict,
    List,
    Sequence,
    Tuple,
)

from antismash.common.secmet import CDSFeature

from .data_structures import Hit, ReferenceArea, ReferenceCDS

Pairing = Tuple[int, int, bool]


NORMALISATION_TARGET = (0.8, 0.8)  # (desired score, required percentage contiguous matches)
EXTRA_SEGMENT_PENALTY = 0.75  # for multiple disjoint segments
REVERSED_SEGMENT_PENALTY = 0.9  # for when the strand of a segment doesn't match


def calculate_order_score_ref_in_query(area_features: Sequence[CDSFeature], hits: Dict[str, Hit], reference: ReferenceArea) -> float:
    return calculate_order_score(area_features, hits, reference, ref_in_query=True)


def calculate_order_score_query_in_ref(area_features: Sequence[CDSFeature], hits: Dict[str, Hit], reference: ReferenceArea) -> float:
    return calculate_order_score(area_features, hits, reference, query_in_ref=True)


def calculate_order_score(area_features: Sequence[CDSFeature], hits: Dict[str, Hit], reference: ReferenceArea, ref_in_query: bool = False, query_in_ref: bool = False) -> float:
    if not hits:
        return 0.
    # build lists of reference features and features ordered by location
    reference_features = reference.cdses

    segments = find_segments(hits, area_features, reference_features, reference.cds_mapping)
    assert segments

    assert not (ref_in_query and query_in_ref)
    if ref_in_query:
        max_possible = len(area_features)
    elif query_in_ref:
        max_possible = len(reference_features)
    else:
        max_possible = min(len(area_features), len(reference_features))

    return score_segments(segments, max_possible)


def _build_segments_from_pairings(pairings: Sequence[Pairing]) -> List[List[Pairing]]:
    segments = [[pairings[0]]]
    for pairing in pairings[1:]:
        prev_cds, prev_ref, prev_strand_match = segments[-1][-1]
        cds, ref, strand_match = pairing
        if any([
            prev_strand_match != strand_match,  # strands no longer match
            cds != prev_cds + 1,  # CDSes aren't contiguous
            abs(ref - prev_ref) != 1,   # references aren't contiguous
        ]):
            segments.append([])
        segments[-1].append(pairing)
    return segments


def find_segments(hits: Dict[str, Hit], features: Sequence[CDSFeature], reference_features: Dict[str, ReferenceCDS], mapping: Dict[int, str]) -> List[List[Pairing]]:

    # assumes features are sorted

    pairings = []
    for ref_cds_name, hit in hits.items():
        cds_index = features.index(hit.cds) + 1  # TODO: performance
        reference_index = int(hit.reference_id)
        reference_strand = reference_features[ref_cds_name].location.strand
        pairings.append((cds_index, reference_index, hit.cds.location.strand == reference_strand))

    pairings.sort()

    return _build_segments_from_pairings(pairings)


def score_segments(segments: List[List[Pairing]], max_possible: int) -> float:
    # calculate the base of the exponential scoring function
    base = (1/NORMALISATION_TARGET[1])**(1/((1-NORMALISATION_TARGET[0])*max_possible))

    result = 0.
    # score each segment
    for segment in segments:
        result += base**len(segment)

    # penalise for each extra segment
    result *= EXTRA_SEGMENT_PENALTY**(len(segments)-1)

    # penalise for having incorrect strands for a segment
    reversed_segments = []  # type: List[List[Pairing]]
    # TODO: reversed is always empty

    result *= REVERSED_SEGMENT_PENALTY**len(reversed_segments)

    # normalise by theoretical perfect score
    result = result / base**max_possible

    return result
