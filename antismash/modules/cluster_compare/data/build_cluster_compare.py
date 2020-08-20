#!/usr/bin/env python3
import glob
import json
import os
import sys
from typing import (
    Any,
    Dict,
    IO,
)

from antismash.common import secmet
from antismash.modules import cluster_compare


class Counter:
    def __init__(self, start: int = 0) -> None:
        self.current = start

    def next(self) -> int:
        self.current += 1
        return self.current


def convert_module(module: secmet.Module) -> Dict[str, Any]:
    return {
        "domains": [domain.domain for domain in module.domains],
        "type": str(module.module_type),
        "complete": module.is_complete(),
    }


def convert_cds(cds: secmet.CDSFeature) -> Dict[str, Any]:
    result = {
        "location": str(cds.location),
        "function": str(cds.gene_function),
        "components": {
            "secmet": cds.sec_met.domain_ids,
            "modules": [convert_module(module) for module in cds.modules],
        },
    }
    return result


def convert_protocluster(protocluster: secmet.Protocluster) -> Dict[str, Any]:
    result = {
        "product": protocluster.product,
        "core_cdses": [cds.get_name() for cds in protocluster.definition_cdses],
        "location": str(protocluster.location),
    }
    return result


def convert_region(region: secmet.Region, cds_mapping: Dict[int, str], cds_index: Counter, fasta: IO) -> Dict[str, Any]:
    result = {
        "products": region.products,
        "protoclusters": [convert_protocluster(pc) for pc in region.get_unique_protoclusters()],
        "cdses": {cds.get_name(): convert_cds(cds) for cds in region.cds_children},
        "start": min(cds.location.start for cds in region.cds_children),  # trim any intergenic areas
        "end": max(cds.location.end for cds in region.cds_children),
    }
    for cds in region.cds_children:
        index = cds_index.next()
        fasta.write(">%s|%d\n%s\n" % (region.parent_record.id, index, cds.translation))
        cds_mapping[index] = cds.get_name()
    return result


def convert_record(record: secmet.Record, fasta: IO) -> Dict[str, Any]:
    result = {
        "regions": [],
        "cds_mapping": {},
    }  # type: Dict[str, Any]
    cds_index = Counter()
    for region in record.get_regions():
        result["regions"].append(convert_region(region, result["cds_mapping"], cds_index, fasta))
    return result


def convert_all(input_dir: str, output_dir: str) -> None:
    files = [name for name in glob.glob(os.path.join(input_dir, "*", "*.gbk")) if "region" not in name]
    assert files
    result = {}
    with open(os.path.join(output_dir, "proteins.fasta"), "w") as fasta:
        for filename in files:
            print(filename)
            for record in secmet.Record.from_genbank(filename):
                result[record.id] = convert_record(record, fasta)
    with open(os.path.join(output_dir, "data.json"), "w") as handle:
        handle.write(json.dumps(result, indent=1))


if __name__ == "__main__":
    output = os.path.join(os.path.dirname(cluster_compare.__file__), "data")
    if len(sys.argv) != 2:
        print("Usage: %s directory_containing_antismash_output_dirs\nOutputs to %s" % (sys.argv[0], output))
        sys.exit(1)
    convert_all(os.path.abspath(sys.argv[1]), output)
