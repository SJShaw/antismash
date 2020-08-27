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
    for variant, result in sorted(results.get_all_variants().items()):
        scores = result.scores_by_region.get(region_layer.get_region_number(), [])[:10]  # TODO: make limit customisable
        tag = "%s-cluster-compare" % (variant.replace(" ", "-"))
        search_type = "mibig"  # TODO: fix this for only known variants
        if "PC_TO_PC" in variant:
            search_type = "matrix"
        div = generate_div(tag, region_layer, record_layer, options_layer, search_type, tooltip, scores, result.scores_by_protocluster)
        html.add_detail_section(tag, div, tag)

    return html


class Row:
    class SubRow:
        def __init__(self, ref_product, scores_per_proto):
            self.product = ref_product
            self.scores_per_proto = scores_per_proto

    def __init__(self, accession, reg, pro, scorer, proto_results, region):
        self.accession = accession
        self.subrows = []
        assert not scorer, type(scorer)


def generate_div(tag: str, region_layer: RegionLayer, record_layer: RecordLayer,
                 options_layer: OptionsLayer, search_type: str,
                 tooltip: str, results: List[Any], proto_results: Dict[int, Dict[Any, Dict[Any, Any]]]) -> Markup:  # TODO fix typing
    """ Generates the specific HTML section of the body for a given variant of
        clusterblast
    """
    large_rows = []
    if search_type == "matrix":
        for proto_number, results in proto_results.items():
            for acc_reg_pro, scorer_dict in results.items():
                import logging; logging.critical(acc_reg_pro); logging.critical(" == "); logging.critical(scorer_dict)
                acc, reg, pro = acc_reg_pro
                import logging; logging.critical(scorer_dict)
                large_rows.append(Row(acc, reg, pro, scorer_dict, proto_results, region_layer))
    template = FileTemplate(path.get_full_path(__file__, "templates", "%s.html" % search_type))
    return template.render(tag=tag, record=record_layer, region=region_layer, options=options_layer, tooltip=tooltip, results=results, proto_results=proto_results, large_rows=large_rows)


def generate_javascript_data(record: Record, region: Region, results: ClusterCompareResults) -> Dict[str, Any]:
    data = {
    }
    for variant, result in results.get_all_variants().items():
        variant_data = {
            "reference_clusters": {}
        }  # type: Dict[str, Dict[str, Any]]

        scores = result.scores_by_region.get(region.get_region_number(), [])[:10]
        if not scores:
            continue

        data[variant] = variant_data

        for accession, totalscore_ref in scores:
            _, ref = totalscore_ref
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
            if not isinstance(accession, str):
                accession = "%s: %s-%s (%s)" % (accession[0], accession[1].start, accession[1].end, accession[2].product)
            variant_data["reference_clusters"][accession] = ref_entry

            mismatching_strands = 0
            for ref_cds_id, hit in result.hits_by_region.get(region.get_region_number(), {}).get(accession, {}).items():
                assert locations.locations_overlap(hit.cds.location, region.location)
                query_cds = hit.cds
                query_point = query_cds.location.start + (query_cds.location.end - query_cds.location.start) // 2
                ref_cds = ref.cdses[ref.cds_mapping[ref_cds_id]]
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
