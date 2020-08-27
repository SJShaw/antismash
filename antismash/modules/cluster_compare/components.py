# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import defaultdict
from typing import (
    Any,  # TODO preferably remove this
    Dict,
    Sequence,
    Tuple,
)

from antismash.common.secmet import CDSFeature
from antismash.common.secmet.locations import location_from_string
from antismash.common.secmet.qualifiers.gene_functions import GeneFunction

SubComponents = Dict[Any, int]


class Components:
    def __init__(self, nrps: SubComponents, pks: SubComponents, secmet: SubComponents, functions: SubComponents) -> None:
        self.nrps = nrps
        self.pks = pks
        self.secmet = secmet
        self.functions = functions

    def __str__(self) -> str:
        import json
        things = {  # TODO: better naming
            "nrps": self.nrps,
            "pks": self.pks,
            "secmet": self.secmet,
            "functions": self.functions,
        }
        parts = []
        for key, val in things.items():
            parts.append("%s  %s" % (key, json.dumps({str(k): v for k, v in val.items()}, indent=1)))
        return "\n".join(parts)


def calculate_component_score_ref_in_query(area_features: Sequence[CDSFeature], ref_data: Dict[str, Any], limit_to_area: Tuple[int, int] = None) -> float:
    return calculate_component_score(area_features, ref_data, ref_in_query=True, limit_to_area=limit_to_area)


def calculate_component_score_query_in_ref(area_features: Sequence[CDSFeature], ref_data: Dict[str, Any], limit_to_area: Tuple[int, int] = None) -> float:
    return calculate_component_score(area_features, ref_data, query_in_ref=True, limit_to_area=limit_to_area)


def calculate_component_score(area_features: Sequence[CDSFeature], ref_data: Dict[str, Any], ref_in_query: bool = False, query_in_ref: bool = False, limit_to_area: Tuple[int, int] = None) -> float:
    assert not (ref_in_query and query_in_ref)
    if limit_to_area:
        assert limit_to_area[0] < limit_to_area[1]
    # TODO properly handle multiple regions
    ref = gather_reference_components(ref_data["regions"][0], limit_to_area=limit_to_area)  # TODO should be handled further up
    # TODO don't repeat the query gather here, do it once per area
    query = gather_query_components(area_features)
    return compare(ref, query, ref_in_query, query_in_ref)


def compare(ref: Components, query: Components, ref_in_query: bool = False, query_in_ref: bool = False) -> float:
    nrps = compare_modules(ref.nrps, query.nrps, ref_in_query, query_in_ref)
    pks = compare_modules(ref.pks, query.pks, ref_in_query, query_in_ref)
    secmet = compare_combos(ref.secmet, query.secmet, ref_in_query, query_in_ref)
    functions = compare_combos(ref.functions, query.functions, ref_in_query, query_in_ref)  # TODO: skip if a minimal run? smcogs missing will cause scores to plummet
    max_modules = sum(ref.pks.values()) + sum(ref.nrps.values())
    modules = 1.
    if max_modules:
        nrps_weighting = sum(ref.nrps.values()) / max_modules
        if 0.001 <= nrps_weighting <= 0.999:
            modules = (nrps * nrps_weighting + pks * (1-nrps_weighting)) / 2
        elif nrps_weighting >= 0.999:
            modules = nrps
        elif nrps_weighting <= 0.001:
            modules = pks
    assert 0 <= modules <= 1, modules
    assert 0 <= secmet <= 1, secmet
    assert 0 <= functions <= 1, functions
    if max_modules:
        return sum([modules, secmet, functions]) / 3
    return (secmet + functions) / 2


def compare_combos(ref: SubComponents, query: SubComponents, ref_in_query: bool = False, query_in_ref: bool = False) -> float:
    if not ref:
        return 1.
    if not query:
        return 0.  # TODO bidirectional, maybe
    ref_combos = set(ref)
    query_combos = set(query)
    ref_extra = ref_combos.difference(query_combos)
    query_extra = query_combos.difference(ref_combos)

    if ref_in_query:
        max_possible = sum(query.values())
    elif query_in_ref:
        max_possible = sum(ref.values())
    else:
        max_possible = min(sum(query.values()), sum(ref.values()))  # TODO directionality would be good to have

    assert max_possible

    found = 0
    for combo in ref_combos.intersection(query_combos):
        found += query[combo]
    return min(found / max_possible, 1.)


def compare_modules(ref: SubComponents, query: SubComponents, ref_in_query: bool = False, query_in_ref: bool = False) -> float:
    if not ref:
        if ref_in_query:
            return 0.
        return 1.
    if not query:
        if ref_in_query:
            return 1.
        return 0.  # TODO bidirectional, maybe
    ref_combos = set(ref)
    query_combos = set(query)
    ref_extra = ref_combos.difference(query_combos)
    query_extra = query_combos.difference(ref_combos)

    # TODO: consider if the modules are different, but still complete, then count them as 0.5
    # since having the same number of modules is still useful even if the
    # modification domains don't exactly match
    # TODO: handle trans-AT specifically, as it's very different

    found = 0.
    if ref_in_query:
        max_possible = sum(ref.values())
        for combo, count in ref.items():
            found += query.get(combo, 0)
    elif query_in_ref:
        max_possible = sum(query.values())
        for combo, count in query.items():
            found += ref.get(combo, 0)
    else:
        max_possible = min(sum(query.values()), sum(ref.values()))
        for combo in ref_combos.intersection(query_combos):
            found += query[combo]

    assert max_possible
    found = min(found, max_possible)
    result = found / max_possible

    assert 0 <= result <= 1, result
    return result


def gather_reference_components(ref_data: Dict[str, Any], limit_to_area: Tuple[int, int] = None) -> Components:
    nrps = defaultdict(int)  # type: SubComponents
    pks = defaultdict(int)  # type: SubComponents
    secmet = defaultdict(int)  # type: SubComponents
    functions = defaultdict(int)  # type: SubComponents

    for cds in ref_data["cdses"].values():
        if limit_to_area:
            # TODO calculate elsewhere and make performant
            loc = location_from_string(cds["location"])
            if not (loc.end > limit_to_area[0] and loc.start < limit_to_area[1]):
                continue

        if cds["function"] != "other":
            functions[cds["function"]] += 1

        if cds["components"]["secmet"]:
            for domain in cds["components"]["secmet"]:
                secmet[domain] += 1

        for module in cds["components"]["modules"]:
            if not module["complete"]:  # TODO check if incomplete is useful
                continue
            if module["type"] == "pks":
                target = pks
            elif module["type"] == "nrps":
                target = nrps
            else:
                continue  # TODO possibly handle unknowns
            target[tuple(module["domains"])] += 1

    return Components(nrps, pks, secmet, functions)


def gather_query_components(area_features: Sequence[CDSFeature]) -> Components:
    nrps = defaultdict(int)  # type: SubComponents
    pks = defaultdict(int)  # type: SubComponents
    secmet = defaultdict(int)  # type: SubComponents
    functions = defaultdict(int)  # type: SubComponents

    for cds in area_features:
        if cds.gene_function != GeneFunction.OTHER:
            functions[str(cds.gene_function)] += 1

        if cds.sec_met:
            for domain in cds.sec_met.domain_ids:
                secmet[domain] += 1

        for module in cds.modules:
            if not module.is_complete():
                continue
            if str(module.module_type) == "pks":
                target = pks
            elif str(module.module_type) == "nrps":  # TODO use actual moduletype enum
                target = nrps
            else:
                assert False, module.module_type  # TODO possibly handle unknowns
            target[tuple(domain.domain.split("_")[0] if domain.domain.startswith("Condensation_") else domain.domain for domain in module.domains if domain.domain)] += 1
    return Components(nrps, pks, secmet, functions)
