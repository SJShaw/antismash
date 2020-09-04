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
from antismash.common.secmet.qualifiers.gene_functions import GeneFunction

from .data_structures import ReferenceArea, Components, SubComponents


def calculate_component_score_ref_in_query(query_components: Components, reference: ReferenceArea) -> float:
    return calculate_component_score(query_components, reference, ref_in_query=True)


def calculate_component_score_query_in_ref(query_components: Components, reference: ReferenceArea) -> float:
    return calculate_component_score(query_components, reference, query_in_ref=True)


def calculate_component_score(query_components: Components, reference: ReferenceArea, ref_in_query: bool = False, query_in_ref: bool = False) -> float:
    assert not (ref_in_query and query_in_ref)
    ref = gather_reference_components(reference)  # TODO should be handled further up
    return compare(ref, query_components, ref_in_query, query_in_ref)


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


def gather_reference_components(reference: ReferenceArea) -> Components:
    existing = reference.get_component_data()
    if existing is not None:
        return existing
    nrps = defaultdict(int)  # type: SubComponents
    pks = defaultdict(int)  # type: SubComponents
    secmet = defaultdict(int)  # type: SubComponents
    functions = defaultdict(int)  # type: SubComponents

    for cds in reference.cdses.values():
        if cds.function != "other":
            functions[cds.function] += 1

        if cds.components["secmet"]:
            for domain in cds.components["secmet"]:
                secmet[domain] += 1

        for module in cds.components["modules"]:
            if not module["complete"]:  # TODO check if incomplete is useful
                continue
            if module["type"] == "pks":
                target = pks
            elif module["type"] == "nrps":
                target = nrps
            else:
                continue  # TODO possibly handle unknowns
            target[tuple(module["domains"])] += 1

    results = Components(nrps, pks, secmet, functions)
    reference.set_component_data(results)
    return results


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
