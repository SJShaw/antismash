# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the nrps_pks module """

from typing import Any, Dict, Optional

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record

from .data_structures import ScoresByRegion, ScoresByProtocluster


class VariantResults:
    def __init__(self, variant_name: str, scores_by_region: ScoresByRegion, scores_by_protocluster: ScoresByProtocluster, hits_by_region) -> None:
        assert variant_name
        self.variant_name = variant_name
        self.scores_by_region = scores_by_region
        self.scores_by_protocluster = scores_by_protocluster
        self.hits_by_region = hits_by_region


class ClusterCompareResults(ModuleResults):
    """ The results of cluster comparison """
    _schema_version = 1

    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)
        self._variants = {}  # type: Dict[str, VariantResults]

    def add_variant_results(self, result: VariantResults) -> None:
        variant = result.variant_name
        if variant in self._variants:
            raise ValueError("Results for %s already exist" % variant)
        self._variants[variant] = result

    def get_variant_results(self, variant: str) -> VariantResults:
        return self._variants[variant]

    def get_all_variants(self) -> Dict[str, VariantResults]:
        return dict(self._variants)  # TODO: make more resilient

    def to_json(self) -> Dict[str, Any]:
        return {}  # TODO

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["ClusterCompareResults"]:
        return None   # TODO

    def add_to_record(self, record: Record) -> None:
        return  # TODO
