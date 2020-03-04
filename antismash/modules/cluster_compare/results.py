# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the nrps_pks module """

from typing import Any, Dict, List, Optional, Tuple

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.common.secmet.features import CDSCollection

from .data_structures import ReferenceArea

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


#class AreaResult:
#    """ Stores results for a specific area in a record, for a particular
#        flavour of clusterblast.
#    """
#    __slots__ = ["area", "ranking", "total_hits", "prefix"]

#    def __init__(self, area: CDSCollection, ranking: List[Tuple[ReferenceArea, ReferenceScorer]],
#                 reference_proteins: Dict[str, Protein], prefix: str) -> None:
#        """ Arguments:
#                cluster: the cluster feature
#                ranking: a list of tuples in the form (ReferenceArea, Score)
#                reference_proteins: used to generate details for SVG output,
#                                    only relevant portions are stored
#                prefix: an identifier for use in marking SVGs such that
#                        javascript on the results page can differentiate between
#                        types of clusterblast
#        """
#        self.area = area
#        self.ranking = ranking[:10]  # TODO
#        self.total_hits = len(ranking)
#        self.prefix = prefix

#    def jsonify(self) -> Dict[str, Any]:
#        """ Convert the object into a simple dictionary for use in storing
#            results.

#            The function from_json() should reconstruct a new and equal
#            AreaResult from the results of this function.

#            Returns:
#                a dict containing the object data in basic types
#        """
#        ranking = []
#        for cluster, score in self.ranking:
#            scoring = {key: getattr(score, key) for key in score.__slots__ if key != "scored_pairings"}
#            json_cluster = {key: getattr(cluster, key) for key in cluster.__slots__}
#            scoring["pairings"] = [(query.entry, query.index, vars(subject))
#                                   for query, subject in score.scored_pairings]
#            ranking.append((json_cluster, scoring))

#        if isinstance(self.area, Region):
#            area_type = "region"
#            area_number = self.area.get_region_number()
#        elif isinstance(self.area, Protocluster):
#            area_type = "protocluster"
#            area_number = self.area.get_protocluster_number()
#        else:
#            raise ValueError("unhandled area type: %s" % type(self.area))

#        result = {
#            "area_number": area_number,
#            "area_type": area_type,
#            "total_hits": self.total_hits,
#            "ranking": ranking,
#            "prefix": self.prefix,
#        }

#        return result

#    @staticmethod
#    def from_json(json: Dict[str, Any], record: Record,
#                  reference_proteins: Dict[str, Protein]) -> "AreaResult":
#        """ Convert a simple dictionary into a new AreaResult object.

#            The function AreaResult.jsonify() should reconstruct the data
#            provided here.

#            Arguments:
#                json: the dict of data to construct with
#                record: the record used to create the data
#                reference_proteins: a dict mapping protein name to Protein,
#                                    used instead of duplicated storing of
#                                    many Protiens

#            Returns:
#                a dict containing the object data in basic types
#        """
#        ranking = []
#        for cluster, details in json["ranking"]:
#            ref_cluster = ReferenceArea(cluster["accession"], cluster["cluster_label"],
#                                           cluster["proteins"], cluster["description"],
#                                           cluster["cluster_type"], cluster["tags"])
#            score = Score()
#            pairings = details["pairings"]
#            for key, val in details.items():
#                if key == "pairings":
#                    continue
#                setattr(score, key, val)
#            for pairing in pairings:  # entry, index, dict of subject
#                query = Query(pairing[0], pairing[1])
#                subject = Subject.from_dict(pairing[2])
#                score.scored_pairings.append((query, subject))
#            ranking.append((ref_cluster, score))

#        region = record.get_region(json["region_number"])
#        result = AreaResult(region, ranking, reference_proteins, json["prefix"])
#        result.total_hits = json["total_hits"]
#        return result
