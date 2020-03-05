# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import defaultdict

from antismash.common.secmet.qualifiers.gene_functions import GeneFunction


class Components:
    def __init__(self, nrps, pks, secmet, functions):
        self.nrps = nrps
        self.pks = pks
        self.secmet = secmet
        self.functions = functions

    def __str__(self):
        import json
        things = {
            "nrps": self.nrps,
            "pks": self.pks,
            "secmet": self.secmet,
            "functions": self.functions,
        }
        parts = []
        for key, val in things.items():
            parts.append("%s  %s" % (key, json.dumps({str(k): v for k, v in val.items()}, indent=1)))
        return "\n".join(parts)


def calculate_component_score(area_features, hits, ref_data, loud=False):
    ref = gather_reference_components(ref_data["regions"][0])  # TODO should be handled further up
    # TODO don't repeat the query gather here, do it once per area
    query = gather_query_components(area_features)
    return compare(ref, query, loud=loud)


def compare(ref, query, loud=False):
    nrps = compare_modules(ref.nrps, query.nrps)
    if loud:
        print("nrps", nrps)
    pks = compare_modules(ref.pks, query.pks, loud)
    if loud:
        print("pks", pks)
    secmet = compare_combos(ref.secmet, query.secmet)
    if loud:
        print("secmet", secmet)
    functions = compare_combos(ref.functions, query.functions)  # TODO: skip if a minimal run? smcogs missing will cause scores to plummet
    if loud:
        print("functions", functions)
    max_modules = sum(ref.pks.values()) + sum(ref.nrps.values())
    modules = 1
    if max_modules:
        nrps_weighting = sum(ref.nrps.values()) / max_modules
        if loud: print(nrps, "*", nrps_weighting, "+", pks, "*", (1-nrps_weighting))
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
        return (modules * secmet * functions) ** (1/3)#sum([modules, secmet, functions]) / 3
    else:
        return (secmet * functions) ** .5


def compare_combos(ref, query, loud=False):
    if not ref:
        return 1.
    if not query:
        return 0.  # TODO bidirectional, maybe
    if loud:
        print("ref", ref)
        print("query", query)
    ref_combos = set(ref)
    query_combos = set(query)
    ref_extra = ref_combos.difference(query_combos)
    query_extra = query_combos.difference(ref_combos)

    max_possible = min(sum(query.values()), sum(ref.values()))  # TODO directionality would be good to have

    assert max_possible

    found = 0
    for combo in ref_combos.intersection(query_combos):
        found += query[combo]
    return min(found / max_possible, 1.)


def compare_modules(ref, query, loud=False):
    if not ref:
        return 1.
    if not query:
        return 0.  # TODO bidirectional, maybe
    if loud:
        print("ref", ref)
        print("query", query)
    ref_combos = set(ref)
    query_combos = set(query)
    ref_extra = ref_combos.difference(query_combos)
    query_extra = query_combos.difference(ref_combos)

    max_possible = max(sum(query.values()), sum(ref.values()))  # TODO directionality would be good to have

    assert max_possible

    found = 0
    for combo in ref_combos.intersection(query_combos):
        found += query[combo]

    # if they're different, but still complete, then count them as half
    # as having the same number of modules is still useful even if the
    # modification domains don't exactly match
    # TODO: handle trans-AT specifically, as it's very different
    found += min(0, len(ref_extra - query_extra)) / 2

    return found / max_possible


def gather_reference_components(ref_data):
    nrps = defaultdict(int)
    pks = defaultdict(int)
    secmet = defaultdict(int)
    functions = defaultdict(int)

    for cds in ref_data["cdses"].values():
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


def gather_query_components(area_features):
    nrps = defaultdict(int)
    pks = defaultdict(int)
    secmet = defaultdict(int)
    functions = defaultdict(int)

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
            target[tuple(domain.domain.split("_")[0] if domain.domain.startswith("Condensation_") else domain.domain for domain in module.domains)] += 1
    return Components(nrps, pks, secmet, functions)
