# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML generation for the clusterblast variants """

from typing import Any, Dict, List

from antismash.common import path
from antismash.common.html_renderer import HTMLSections, FileTemplate, Markup
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer
from antismash.common.secmet import Record, Region, locations

from .results import ClusterCompareResults


def will_handle(_products: List[str]) -> bool:
    """ Clusterblast is relevant to every region, so return True for every
        product """
    return True

def scores_within_limit(scores):
    return [score for score in scores[:10] if score[1].hits_by_gene]


def generate_html(region_layer: RegionLayer, results: ClusterCompareResults,
                  record_layer: RecordLayer, options_layer: OptionsLayer
                  ) -> HTMLSections:
    """ Generates the HTML sections for all variants of clusterblast
    """

    html = HTMLSections("cluster-compare")
    base_tooltip = ("Shows %s that are similar to the current region. Genes marked with the "
                    "same colour are interrelated. White genes have no relationship.<br>"
                    "Click on reference genes to show details of similarities to "
                    "genes within the current region.")

    # TODO  variants
    tooltip = base_tooltip % "clusters from the MiBIG database"
    tooltip += "<br>Click on an accession to open that entry in the MiBIG database."
    scores = scores_within_limit(results.scores_by_region.get(region_layer.get_region_number(), []))
    div = generate_div(region_layer, record_layer, options_layer, "mibig", tooltip, scores)
    html.add_detail_section("known-cluster-compare", div, "known-cluster-compare")

    return html


def generate_div(region_layer: RegionLayer, record_layer: RecordLayer,
                 options_layer: OptionsLayer, search_type: str,
                 tooltip: str, results: str) -> Markup:  # TODO fix typing
    """ Generates the specific HTML section of the body for a given variant of
        clusterblast
    """
    template = FileTemplate(path.get_full_path(__file__, "templates", "%s.html" % search_type))
    return template.render(record=record_layer, region=region_layer, options=options_layer, tooltip=tooltip, results=results)


def generate_javascript_data(record: Record, region: Region, results: ClusterCompareResults) -> Dict[str, Any]:
    data = {
        "reference_clusters": {},
    }
    # TODO variants
    known = scores_within_limit(results.scores_by_region.get(region.get_region_number(), []))
    if not known:
        return {}
    for accession, score in known:
        ref = score.reference
        genes = ref.get_cds_json()
        for name, gene in genes.items():
            gene["locus_tag"] = name
            location = locations.location_from_string(gene["location"])
            gene["start"] = location.start
            gene["end"] = location.end
            gene["strand"] = location.strand
            gene["linked"] = {}
        ref_entry = {  # TODO: merge these for different protoclusters
            "links": [],
            "start": ref.start,
            "end": ref.end,
        }
        data["reference_clusters"][accession] = ref_entry

        mismatching_strands = 0
        for ref_cds_id, hit in score.hits_by_gene.items():
            assert locations.locations_overlap(hit.cds.location, region.location)
            query_cds = hit.cds
            query_point = query_cds.location.start + (query_cds.location.end - query_cds.location.start) // 2
            ref_cds = score.reference.cdses[score.reference.cds_mapping[ref_cds_id]]
            ref_location = locations.location_from_string(ref_cds.location)
            subject_point = ref_location.start + (ref_location.end - ref_location.start) // 2
            if query_cds.location.strand != ref_location.strand:
                mismatching_strands += 1
            genes[ref_cds.name]["linked"][region.get_region_number()] = query_cds.get_name()
            ref_entry["links"].append({
                "query": query_cds.get_name(),
                "subject": ref_cds.name,
                "query_loc": query_point,
                "subject_loc": subject_point,
            })
        ref_entry["reverse"] = mismatching_strands > len(ref_entry["links"]) / 2
        ref_entry["genes"] = list(genes.values())
    return data
