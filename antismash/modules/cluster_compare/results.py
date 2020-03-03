# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the nrps_pks module """

from typing import Any, Dict, Optional

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record


class ClusterCompareResults(ModuleResults):
    """ The results of cluster comparison """
    _schema_version = 1

    def __init__(self, record_id, scores) -> None:
        super().__init__(record_id)
        self.scores = scores

    def to_json(self) -> Dict[str, Any]:
        return {}  # TODO

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["ClusterCompareResults"]:
        return None   # TODO

    def add_to_record(self, record: Record) -> None:
        return  # TODO
