# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML generation for the clusterblast variants """

from typing import Any, Dict, List

from antismash.common import path
from antismash.common.html_renderer import HTMLSections, FileTemplate, Markup
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer
from antismash.common.secmet import Record

from .results import ClusterCompareResults


def will_handle(_products: List[str]) -> bool:
    """ Clusterblast is relevant to every region, so return True for every
        product """
    return True


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
    div = generate_div(region_layer, record_layer, options_layer, "mibig", tooltip, results.scores)
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


def generate_javascript_data(record: Record, results: ClusterCompareResults) -> Dict[str, Any]:
    data = {
        "reference_clusters": {},
    }
    known = results.scores  # TODO
    assert known
    for _, score in known[:10]:  # TODO: use actual set
        assert ref.protein_details, ref
        genes = ref.cdses
        ref_entry = {  # TODO: merge these for different protoclusters
            "links": [],
        }
        data["reference_clusters"][ref.get_name()] = ref_entry

        mismatching_strands = 0
        for query, subject in score.scored_pairings:
            query_cds = record.get_cds_by_name(query.id)
            query_point = query_cds.location.start + (query_cds.location.end - query_cds.location.start) // 2
            subject_point = subject.start + (subject.end - subject.start) // 2
            if query_cds.location.strand != (-1 if subject.strand == "-" else 1):
                mismatching_strands += 1
            genes.get(subject.locus_tag, genes.get(subject.name))["linked"] = True
            ref_entry["links"].append({
                "query": query.id,
                "subject": subject.locus_tag,
                "query_loc": query_point,
                "subject_loc": subject_point,
            })
        ref_entry["reverse"] = mismatching_strands > len(ref_entry["links"]) / 2
        ref_entry["genes"] = list(genes.values())
        ref_entry["start"] = min(x["start"] for x in ref_entry["genes"])
        ref_entry["end"] = max(x["end"] for x in ref_entry["genes"])
        print(ref, len(ref_entry["links"]))
    return data
