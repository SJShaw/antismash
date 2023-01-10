# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" TIGRFam anotation for only clusters """

import logging
import os
from typing import Any, Dict, List, Optional

from antismash.common import hmmer, path
from antismash.common.secmet import Record
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

    features = []
    for region in record.get_regions():
        features.extend(list(region.cds_children))
    db_path = path.get_full_path(hmm_detection.__file__, "data", "bgc_seeds.hmm")
    hmmer_results = hmmer.run_hmmer(record, features, MAX_EVALUE, MIN_SCORE, db_path,
                                    "rulevis", filter_overlapping=False, use_cut_tc=False)
    return RuleVisResults.from_hmmer_results(hmmer_results)
