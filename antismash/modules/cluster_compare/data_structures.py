# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of data structures for use in the rest of the module"""

import json

from typing import Any, Dict, List


class RawCDS:
    def __init__(self, name, function, components, location) -> None:
        self.name = name
        self.function = function
        self.components = components
        self.location = location

    @classmethod
    def from_json(cls, name: str, data: Dict[str, Any]) -> "RawCDS":
        return cls(name, data["function"], data["components"], data["location"])


class RawProtocluster:
    def __init__(self, cores: List[RawCDS], product: str, location: str) -> None:
        self.cores = cores
        self.product = product
        self.location = location

    @classmethod
    def from_json(cls, data: Dict[str, Any], cdses):
        cores = [cdses[core] for core in data["core_cdses"]]
        return cls(cores, data["product"], data["location"])

class RawRegion:
    def __init__(self, protoclusters, cdses, products, cds_mapping, raw_cdses, start, end):
        self.protoclusters = protoclusters
        self.cdses = cdses
        self.products = products
        self.cds_mapping = cds_mapping
        self.raw_cdses = raw_cdses
        self.start = start
        self.end = end

    @classmethod
    def from_json(cls, data: Dict[str, Any], cds_mapping) -> "RawRegion":
        cdses = {name: RawCDS.from_json(name, cds) for name, cds in data["cdses"].items()}
    
        return cls([RawProtocluster.from_json(proto, cdses) for proto in data["protoclusters"]],
                   cdses,
                   data["products"],
                   cds_mapping,
                   data["cdses"],
                   data["start"],
                   data["end"],
                  )

    @property
    def product_string(self) -> str:
        return ", ".join(self.products)

    def get_cds_json(self) -> str:
        return self.raw_cdses


class RawRecord:
    def __init__(self, accession, regions, cds_mapping):
        self.accession = accession
        self.regions = regions
        self.cds_mapping = cds_mapping

    @classmethod
    def from_json(cls, accession: str, data: Dict[str, Any]) -> "RawRegion":
        regions = [RawRegion.from_json(region, data["cds_mapping"]) for region in data["regions"]]
        return RawRecord(accession, regions, data["cds_mapping"])


def load_data(filename: str):
    with open(filename) as handle:
        raw = json.loads(handle.read())
    return {accession: RawRecord.from_json(accession, record) for accession, record in raw.items()}
        



class ReferenceGene:
    """ A reference gene referred to/contained by ReferenceAreas """
    __slots__ = ["name", "location"]
    def __init__(self, name: str, location: str) -> None:
        self.name = name
        self.location = location


class ReferenceArea:
    """ A reference cluster container, as read from a database of
        antismash-predicted areas.
    """
    __slots__ = ["accession", "location", "gene_names", "description", "kind"]

    def __init__(self, accession: str, location: str, gene_names: List[str],
                 description: str, kind: str) -> None:
        self.accession = accession
        self.location = location
        self.gene_names = gene_names
        self.description = description
        self.kind = kind

    def get_name(self) -> str:
        """ Returns a short label for the cluster that """
        return "%s_%s" % (self.accession, self.location)
