# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of data structures for use in the rest of the module"""

import json
from typing import Any, Callable, Dict, List, Optional, Tuple

from antismash.common.secmet import CDSFeature
from antismash.common.secmet.locations import location_from_string, Location


class Hit:
    def __init__(self, reference_record: str, reference_id: str, cds: CDSFeature,
                 percent_identity: int, blast_score: int, percent_coverage: float, evalue: float) -> None:
        self.reference_record = reference_record
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
        return "%s->%s(pid=%d,cov%d)" % (self.cds.get_name(), self.reference_id, self.percent_identity, self.percent_coverage)


class ReferenceCDS:
    def __init__(self, name: str, function: str, components: Dict[str, List[Any]], location: Location) -> None:
        self.name = name
        self.function = function
        self.components = components
        self.location = location

    def overlaps_with(self, area: "ReferenceArea") -> bool:
        return self.location.end > area.start and self.location.start < area.end

    @classmethod
    def from_json(cls, name: str, data: Dict[str, Any]) -> "ReferenceCDS":
        return cls(name, data["function"], data["components"], location_from_string(data["location"]))

    def __str__(self) -> str:
        return "ReferenceCDS(%s, %s)" % (self.name, self.location)

    def get_minimal_json(self) -> Dict[str, Any]:
        # to match IOrf in antismash-js for the purposes of drawing
        return {
            "locus_tag": self.name,
            "start": self.location.start,
            "end": self.location.end,
            "strand": self.location.strand,
            "function": self.function,
        }


class ReferenceArea:
    def __init__(self, accession: str, start: int, end: int, cds_mapping: Dict[str, str], cdses: Dict[str, ReferenceCDS], products: List[str]) -> None:
        self.accession = accession
        self.start = start
        self.end = end
        self.cds_mapping = cds_mapping
        self.cdses = {name: cds for name, cds in cdses.items() if cds.overlaps_with(self)}
        self.products = products

    def get_product_string(self) -> str:
        return ", ".join(self.products)

    def get_identifier(self) -> str:
        return "%s (%s-%s)" % (self.accession, self.start, self.end)


class ReferenceProtocluster(ReferenceArea):
    def __init__(self, accession: str, start: int, end: int, cds_mapping: Dict[str, str], cdses: Dict[str, ReferenceCDS], cores: List[ReferenceCDS], product: str) -> None:
        super().__init__(accession, start, end, cds_mapping, cdses, [product])
        self.cores = cores
        self.product = product

    @classmethod
    def from_json(cls, accession: str, data: Dict[str, Any], cdses: Dict[str, ReferenceCDS], cds_mapping: Dict[str, str]) -> "ReferenceProtocluster":
        cores = [cdses[core] for core in data["core_cdses"]]
        location = location_from_string(data["location"])  # TODO: this is silly
        return cls(accession, location.start, location.end, cds_mapping, cdses, cores, data["product"])

    def __str__(self) -> str:
        return "ReferenceProtocluster(%s, %s, %s...)" % (self.start, self.end, self.product)

    def __repr__(self) -> str:
        return str(self)


class ReferenceRegion(ReferenceArea):
    def __init__(self, accession: str, start: int, end: int, protoclusters: List[ReferenceProtocluster], cdses: Dict[str, ReferenceCDS],
                 products: List[str], cds_mapping: Dict[str, str]) -> None:
        super().__init__(accession, start, end, cds_mapping, cdses, products)
        self.protoclusters = protoclusters

    @classmethod
    def from_json(cls, accession: str, data: Dict[str, Any], cds_mapping: Dict[str, str]) -> "ReferenceRegion":
        cdses = {name: ReferenceCDS.from_json(name, cds) for name, cds in data["cdses"].items()}

        return cls(
            accession,
            data["start"],
            data["end"],
            [ReferenceProtocluster.from_json(accession, proto, cdses, cds_mapping) for proto in data["protoclusters"]],
            cdses,
            data["products"],
            cds_mapping,
        )

    @property
    def product_string(self) -> str:
        return ", ".join(self.products)

    def __str__(self) -> str:
        return "ReferenceRegion(%s, %s, %s, %s...)" % (self.accession, self.start, self.end, self.products)

    def __repr__(self) -> str:
        return str(self)


class ReferenceRecord:
    def __init__(self, accession: str, regions: List[ReferenceRegion], cds_mapping: Dict[str, str]) -> None:
        self.accession = accession
        self.regions = regions
        self.cds_mapping = cds_mapping

    @classmethod
    def from_json(cls, accession: str, data: Dict[str, Any]) -> "ReferenceRecord":
        regions = [ReferenceRegion.from_json(accession, region, data["cds_mapping"]) for region in data["regions"]]
        return ReferenceRecord(accession, regions, data["cds_mapping"])


def load_data(filename: str) -> Dict[str, ReferenceRecord]:
    with open(filename) as handle:
        raw = json.loads(handle.read())
    return {accession: ReferenceRecord.from_json(accession, record) for accession, record in raw.items()}


class ReferenceScorer:
    def __init__(self, best_hits: Dict[str, Hit],
                 reference: ReferenceArea, query_features: Tuple[CDSFeature, ...], query_components,
                 ident_calculator: Callable, order_calculator: Callable, component_calculator: Callable) -> None:  # TODO parse data to object
        self.accession = reference.accession
        self.hits_by_gene = best_hits
        for hit in best_hits.values():
            assert hit.reference_record == reference.accession, "%s != %s" % (hit.reference_record, reference.accession)
        self._identity = -1.
        self._order = -1.
        self._components = -1.
        self.reference = reference
        self._raw_identity = -1.
        self._max_id = -1.
        self._query_features = query_features
        self._ident_calculator = ident_calculator
        self._order_calculator = order_calculator
        self._component_calculator = component_calculator
        self._query_components = query_components

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
            self._order = self._order_calculator(sorted(self._query_features), self.hits_by_gene, self.reference)
        return self._order

    @property
    def components(self) -> float:
        if self._components < 0:
            self._components = self._component_calculator(self._query_components, self.reference)
        return self._components

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return "ReferenceScorer(%s: raw_id=%.2f, order=%.2f, comp=%.2f)" % (
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

    def __lt__(self, other: Any) -> bool:
        if not isinstance(other, ReferenceScorer):
            raise TypeError("cannot compare ReferenceScorer to %s" % type(other))
        return self.final_score < other.final_score


HitsByReference = Dict[ReferenceRegion, Dict[str, List[Hit]]]
ScoresByRegion = List[Tuple[ReferenceArea, float]]
ScoresByProtocluster = Dict[int, Dict[ReferenceArea, ReferenceScorer]]
