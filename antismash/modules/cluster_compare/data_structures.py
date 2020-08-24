# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of data structures for use in the rest of the module"""

import json
from typing import Any, Callable, Dict, List, Tuple

from antismash.common.secmet import CDSFeature


class Hit:
    def __init__(self, reference_cluster: str, reference_id: str, cds: CDSFeature,
                 percent_identity: int, blast_score: int, percent_coverage: float, evalue: float) -> None:
        self.reference_cluster = reference_cluster  # TODO actually reference record, need to split into regions/protoclusters
        self.reference_id = reference_id
        self.cds = cds
        self.percent_identity = percent_identity
        self.blast_score = blast_score
        self.percent_coverage = percent_coverage
        self.evalue = evalue

    @property
    def identity_score(self) -> float:
        return self.blast_score * self.percent_coverage

    def __repr__(self) -> str:
        return "%s(pid=%d,cov%d)" % (self.reference_id, self.percent_identity, self.percent_coverage)


class RawCDS:
    def __init__(self, name: str, function: str, components: Dict[str, List[Any]], location: str) -> None:
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
    def from_json(cls, data: Dict[str, Any], cdses: Dict[str, RawCDS]) -> "RawProtocluster":
        cores = [cdses[core] for core in data["core_cdses"]]
        return cls(cores, data["product"], data["location"])


class RawRegion:
    def __init__(self, protoclusters: List[RawProtocluster], cdses: Dict[str, RawCDS],
                 products: List[str], cds_mapping: Dict[int, str], raw_cdses: Dict[str, Any], start: int, end: int) -> None:
        self.protoclusters = protoclusters
        self.cdses = cdses
        self.products = products
        self.cds_mapping = cds_mapping
        self.raw_cdses = raw_cdses
        self.start = start
        self.end = end

    @classmethod
    def from_json(cls, data: Dict[str, Any], cds_mapping: Dict[int, str]) -> "RawRegion":
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

    def get_cds_json(self) -> Dict[str, Any]:
        return self.raw_cdses


class RawRecord:
    def __init__(self, accession: str, regions: List[RawRegion], cds_mapping: Dict[int, str]) -> None:
        self.accession = accession
        self.regions = regions
        self.cds_mapping = cds_mapping

    @classmethod
    def from_json(cls, accession: str, data: Dict[str, Any]) -> "RawRecord":
        regions = [RawRegion.from_json(region, data["cds_mapping"]) for region in data["regions"]]
        return RawRecord(accession, regions, data["cds_mapping"])


def load_data(filename: str) -> Dict[str, RawRecord]:
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
        return "%s_%s" % (self.accession, self.location)


class ReferenceScorer:
    def __init__(self, accession: str, data: Dict[str, Any], best_hits: Dict[str, Hit],
                 reference: RawRecord, area_features: Tuple[CDSFeature, ...],
                 ident_calculator: Callable, order_calculator: Callable, component_calculator: Callable
                 ) -> None:  # TODO parse data to object
        self.accession = accession
        self.data = data
        self.hits_by_gene = best_hits
        self._identity = -1.
        self._order = -1.
        self._components = -1.
        self.reference = reference.regions[0]  # TODO
        self._raw_identity = -1.
        self._max_id = -1.
        self._area_features = area_features
        self._ident_calculator = ident_calculator
        self._order_calculator = order_calculator
        self._component_calculator = component_calculator

    def calc_identity(self, max_id: float) -> float:
        self._max_id = max_id
        return self.identity

    @property
    def final_score(self) -> float:
        return (self.identity * self.order * self.components) ** (1/3)

    @property
    def raw_identity(self) -> float:
        if self._raw_identity < 0:
            self._raw_identity = sum(hit.identity_score for hit in self.hits_by_gene.values())
        return self._raw_identity

    @property
    def identity(self) -> float:
        assert self._max_id > 0
        if self._identity < 0:
            self._identity = self._ident_calculator(self.raw_identity, self._max_id)
        assert 0 <= self._identity <= 1, (self._identity, self.raw_identity, self._max_id)
        return self._identity

    @property
    def order(self) -> float:
        if self._order < 0:
            self._order = self._order_calculator(sorted(self._area_features), self.hits_by_gene, self.data)
        return self._order

    @property
    def components(self) -> float:
        if self._components < 0:
            self._components = self._component_calculator(self._area_features, self.data)
        return self._components

    def __repr__(self) -> str:
        return "ReferenceScorer(%s)" % (str(self))

    def __str__(self) -> str:
        return "%s: raw_id=%.2f, order=%.2f, comp=%.2f" % (
            self.accession,
            self.identity,
            self.order,
            self.components,
        )

    def table_string(self) -> str:
        return "%.2f (id:%.2f, order:%.2f, components:%.2f)" % (
            self.final_score,
            self.identity,
            self.order,
            self.components,
        )


ScoresByRegion = Dict[int, List[Tuple[str, Tuple[float, RawRegion]]]]
ScoresByProtocluster = Dict[int, List[Tuple[str, ReferenceScorer]]]
