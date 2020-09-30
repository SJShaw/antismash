# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML generation for the clusterblast variants """

from typing import Any, Dict, List, Union

from antismash.common import path
from antismash.common.html_renderer import HTMLSections, FileTemplate, Markup
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer
from antismash.common.secmet import Record, Region, locations

from .analysis import trim_to_best_hit
from .data_structures import ScoresByProtocluster, ReferenceScorer
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

    # TODO adjust tooltip by variant
    tooltip = base_tooltip % "clusters from the MiBIG database"
    tooltip += "<br>Click on an accession to open that entry in the MiBIG database."
    for variant, result in sorted(results.results_by_region[region_layer.get_region_number()].items()):
        scores = result.scores_by_region[:10]  # TODO: make limit customisable
        tag = "%s-cluster-compare" % (variant.replace(" ", "-"))
        search_type = "mibig"  # TODO: fix this for only known variants
        if "PC_TO_PC" in variant:
            search_type = "matrix"
        elif "ALL_TO_ALL" in variant:
            search_type = "single"  # TODO: limit results to 10/custom option
        div = generate_div(tag, region_layer, record_layer, options_layer, search_type, tooltip, scores, result.scores_by_protocluster)
        html.add_detail_section(tag, div, tag)

    return html


class Row:
    class SubRow:
        def __init__(self, ref_product: str, scores_per_proto: ScoresByProtocluster) -> None:
            self.product = ref_product
            self.scores_per_proto = scores_per_proto

    def __init__(self, accession: str, reg, pro, scorer, proto_results, region) -> None:
        self.accession = accession
        self.subrows = []
        assert not scorer, type(scorer)


def generate_div(tag: str, region_layer: RegionLayer, record_layer: RecordLayer,
                 options_layer: OptionsLayer, search_type: str,
                 tooltip: str, results: List[Any], proto_results: Union[ScoresByProtocluster, List[ReferenceScorer]]) -> Markup:  # TODO fix typing
    """ Generates the specific HTML section of the body for a given variant of
        clusterblast
    """
    template = FileTemplate(path.get_full_path(__file__, "templates", "%s.html" % search_type))
    return template.render(tag=tag, record=record_layer, region=region_layer, options=options_layer, tooltip=tooltip, results=results, proto_results=proto_results)


def generate_javascript_data(record: Record, region: Region, results: ClusterCompareResults) -> Dict[str, Any]:
    data = {
    }
    import logging; logging.critical("cluster compare: skipping JS output components")
    return data
    for variant, result in results.results_by_region[region.get_region_number()].items():
        variant_data = {
            "reference_clusters": {}
        }  # type: Dict[str, Dict[str, Any]]

        scores = result.scores_by_region[:10]  # TODO: use matching custom option as above
        if not scores:
            continue

        data[variant] = variant_data

        for reference, score in scores:
            ref_entry = {  # TODO: merge these for different protoclusters
                "start": reference.start,
                "end": reference.end,
                "links": [],  # added to later
                "genes": [],  # added to later
                "reverse": False,  # potentially changed later
            }
            genes = {}
            for name, cds in reference.cdses.items():
                gene_json = cds.get_minimal_json()
                gene_json["linked"] = {}
                genes[cds.name] = gene_json
            # TODO: match reference names/labels for both rows in HTML and javascript for drawing
            accession = "%s: %s-%s (%s)" % (reference.accession, reference.start, reference.end, ", ".join(reference.products))
            variant_data["reference_clusters"][accession] = ref_entry

            mismatching_strands = 0
            for ref_cds_id, hit in trim_to_best_hit(result.hits_by_region.get(reference, {})).items():
                assert locations.locations_overlap(hit.cds.location, region.location)
                query_cds = hit.cds
                query_point = query_cds.location.start + (query_cds.location.end - query_cds.location.start) // 2
                ref_cds = reference.cdses[ref_cds_id]
                subject_point = ref_cds.location.start + (ref_cds.location.end - ref_cds.location.start) // 2
                if query_cds.location.strand != ref_cds.location.strand:
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
