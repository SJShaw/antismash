# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

NORMALISATION_TARGET = (0.8, 0.8)  # (desired score, required percentage contiguous matches)
EXTRA_SEGMENT_PENALTY = 0.75  # for multiple disjoint segments
REVERSED_SEGMENT_PENALTY = 0.9  # for when the strand of a segment doesn't match


def calculate_order_score(area_features, hits, ref_data):
    if not hits:
        return 0.
    # build lists of reference features and features ordered by location
    ref_names = list(ref_data["regions"][0]["cdses"])  # TODO handle protoclusters instead of regions
    reference_features = ref_data["regions"][0]["cdses"]
#    reference_feature_order = {v: k for k, v in ref_data["regions"][0]["cds_mapping"]}

    segments = find_segments(hits, area_features, reference_features, ref_data["cds_mapping"])
    assert segments, hits[0].reference_id
    max_possible = min(len(area_features), len(reference_features))

    return score_segments(segments, max_possible)


def _build_segments_from_pairings(pairings, loud=False):
    segments = [[pairings[0]]]
    if loud:
        print("prev_pair, pairing, [strandcont, cds_contig, ref_contig]")
    for pairing in pairings[1:]:
        prev_cds, prev_ref, prev_strand_match = segments[-1][-1]
        cds, ref, strand_match = pairing
        if any([
            prev_strand_match != strand_match,  # strands no longer match
            cds != prev_cds + 1,  # CDSes aren't contiguous
            abs(ref - prev_ref) != 1,   # references aren't contiguous
        ]):
            if loud:
                print(segments[-1][-1], pairing, [
                    prev_strand_match != strand_match,  # strands no longer match
                    cds != prev_cds + 1,  # CDSes aren't contiguous
                    abs(ref - prev_ref) != 1,   # references aren't contiguous
                ])
            segments.append([])
        segments[-1].append(pairing)
    if loud:
        print("total segment counts", [len(i) for i in segments])
    return segments


def find_segments(hits, features, reference_features, mapping):

    # features should always be sorted
    assert sorted(features) == features

    pairings = []
    for hit in hits.values():
        cds_index = features.index(hit.cds) + 1
        reference_index = int(hit.reference_id)
        reference_strand = 1 if "+" in reference_features[mapping[hit.reference_id]]["location"] else -1
        pairings.append((cds_index, reference_index, hit.cds.location.strand == reference_strand))

    pairings.sort()

    return _build_segments_from_pairings(pairings)


def score_segments(segments, max_possible):
    # calculate the base of the exponential scoring function
    base = (1/NORMALISATION_TARGET[1])**(1/((1-NORMALISATION_TARGET[0])*max_possible))

    result = 0.
    # score each segment
    for segment in segments:
        result += base**len(segment)

    # penalise for each extra segment
    result *= EXTRA_SEGMENT_PENALTY**(len(segments)-1)

    # penalise for having incorrect strands for a segment
    reversed_segments = []
    result *= REVERSED_SEGMENT_PENALTY**len(reversed_segments)

    # normalise by theoretical perfect score
    result = result / base**max_possible

    return result
