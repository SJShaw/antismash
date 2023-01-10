# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A visualition module for profile hits from rule detection """

from collections import defaultdict
import logging
from typing import Any, Dict, List, Optional, Union

from antismash.common import hmmer
from antismash.common.hmm_rule_parser import cluster_prediction
from antismash.common.secmet import Record
from antismash.common.secmet.locations import FeatureLocation, convert_protein_position_to_dna
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from antismash.detection import DetectionStage, hmm_detection

from .rulevis_results import RuleVisResults


NAME = "rulevis"
SHORT_DESCRIPTION = "Rule domain visualisation"
DETECTION_STAGE = DetectionStage.PER_AREA

MIN_SCORE = 0.
MAX_EVALUE = 0.01


def get_arguments() -> ModuleArgs:
    """ Builds the module args """
    args = ModuleArgs('RuleVis options', 'rulevis')
    return args


def is_enabled(options: ConfigType) -> bool:
    """  Uses the supplied options to determine if the module should be run """
    return True


def check_prereqs(options: ConfigType) -> List[str]:
    """ Ensure at least one database exists and is valid """
    return []


def check_options(_options: ConfigType) -> List[str]:
    """ Check the requested PFAM database exists """
    # No options to check
    return []


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[RuleVisResults]:
    """ Rebuild previous results """
    if not previous:
        return None
    results = RuleVisResults.from_json(previous, record)
    if not results:
        return None
    if results.score > MIN_SCORE or results.evalue < MAX_EVALUE:
        # new values too lenient, discard resuts
        return None
    return results.refilter(MAX_EVALUE, MIN_SCORE)


def run_on_record(record: Record, results: Optional[RuleVisResults],
                  options: ConfigType) -> RuleVisResults:
    """ Run hmmsearch against RuleVis for all CDS features within the record """

    logging.info('Running RuleVis search')

    if results:
        return results

    ruleset = hmm_detection.get_ruleset(options)

    all_hits = []

    # get the HMMer profile results
    results_by_id = cluster_prediction.find_hmmer_hits(record, ruleset.hmm_profiles, ruleset.database_file,
                                                       ruleset.get_equivalence_groups())
    for hits in results_by_id.values():
        all_hits.extend(hits)

    dynamic_results = cluster_prediction.find_dynamic_hits(record, list(ruleset.dynamic_profiles.values()), results_by_id)
    for dynamic_hits in dynamic_results.values():
        all_hits.extend(dynamic_hits)

    domain_counters: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

    def convert(hit: Union[cluster_prediction.HMMerHit, cluster_prediction.ProfileHit]) -> hmmer.HmmerHit:
        cds = record.get_cds_by_name(hit.hit_id)
        if isinstance(hit, cluster_prediction.HMMerHit):
            query_start = hit.query_start
            query_end = hit.query_end
            location = convert_protein_position_to_dna(query_start, query_end, cds.location)
        else:
            query_start = 0
            query_end = len(cds.translation)
            location = cds.location
        label = f"{cds.get_name()}_{hit.query_id}.{domain_counters[cds.get_name()][hit.query_id]}"
        domain_counters[cds.get_name()][hit.query_id] += 1
        return hmmer.HmmerHit(
            location=str(location),
            label=label,
            locus_tag=cds.get_name(),
            domain=hit.query_id,
            evalue=hit.evalue,
            score=hit.bitscore,
            identifier=hit.query_id,
            description="",
            protein_start=query_start,
            protein_end=query_end,
            translation=cds.translation[query_start:query_end],
        )

    converted_hits = sorted((convert(hit) for hit in all_hits),
                            key=lambda x: (record.get_cds_by_name(x.locus_tag), -len(x.translation)))
    return RuleVisResults(record.id, MAX_EVALUE, MIN_SCORE, converted_hits)
