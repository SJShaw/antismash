# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML output for the Type II PKS module """

from typing import List

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSections, docs_link
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .results import ClusterCompareResults


def will_handle(_products: List[str]) -> bool:
    return True


def generate_html(region_layer: RegionLayer, results: ClusterCompareResults,
                  _record_layer: RecordLayer, _options_layer: OptionsLayer) -> HTMLSections:
    """ Generate HTML with results """
    html = HTMLSections("cluster-compare")

    tooltip_content = ("TODO")  # TODO

    return html
