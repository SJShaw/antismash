# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import json
import unittest
from unittest.mock import patch

from antismash.common import (
    fasta,
    path,
    subprocessing,
)
from antismash.common.secmet.test.helpers import DummyCDS
from antismash.common.signature import HmmSignature
from antismash.common.test.helpers import FakeHSPHit, FakeHit
from antismash.detection.genefunctions.halogenases import (
    HalogenaseHmmResult,
    FlavinDependentHalogenase as FDH,
    Match,
)
from antismash.detection.genefunctions.halogenases.flavin_dependent import substrate_analysis
from antismash.detection.genefunctions.halogenases.flavin_dependent.substrate_analysis import (
    run_halogenase_phmms,
    categorize_on_substrate_level,
)


def test_categorisation_with_no_hits():
    with patch.object(substrate_analysis, "retrieve_fdh_signature_residues", return_value={}):
        result = categorize_on_substrate_level(DummyCDS(), FDH("dummy"), [])
    assert result is None


def test_conversion_methods():
    fdh = FDH("dummy")
    fdh.add_potential_matches([
        Match(
            "prof_A",
            "flavin",
            "flavin-dependent",
            confidence=0.5,
            consensus_residues="MAGIC",
            substrate="some_sub",
            target_positions=(3, 4),
            number_of_decorations="deco",
        ),
        Match(
            "prof_B",
            "flavin",
            "flavin-dependent",
            confidence=0.1,
            consensus_residues="BEANS",
            substrate="other_sub",
            target_positions=(1,),
        ),
    ])
    original = fdh.to_json()
    assert isinstance(original, dict)
    raw = json.loads(json.dumps(original))

    rebuilt = FDH.from_json(raw)
    assert isinstance(rebuilt, FDH)
    assert rebuilt.to_json() == original


@patch.object(subprocessing, "run_hmmsearch",
              return_value=[FakeHit(1, 2, 1000, "foo")])
def test_run_halogenase_phmms(_patched):
    signature = HmmSignature("foo", "description", 300, "dummy_path")
    for value in _patched.return_value:
        value.hsps = [FakeHSPHit("foo", "query", bitscore=250)]

    negative_test_halogenase_hmms_by_id = run_halogenase_phmms("", [signature])
    assert not negative_test_halogenase_hmms_by_id

    for value in _patched.return_value:
        value.hsps = [FakeHSPHit("foo", "query", bitscore=1000)]

    positive_test_halogenase_hmms_by_id = run_halogenase_phmms("", [signature])
    assert positive_test_halogenase_hmms_by_id
    for hit in positive_test_halogenase_hmms_by_id["foo"]:
        assert isinstance(hit, HalogenaseHmmResult)


class TestGetBest(unittest.TestCase):
    def setUp(self):
        self.high_confidence_match = Match(
            "prof_A",
            "flavin",
            "flavin-dependent",
            confidence=0.9,
            consensus_residues="MAGIC",
            substrate="some_sub",
            target_positions=(3, 4),
            number_of_decorations="deco",
        )
        self.low_confidence_match = Match(
            "prof_B",
            "flavin",
            "flavin-dependent",
            confidence=0.1,
            consensus_residues="BEANS",
            substrate="other_sub",
            target_positions=(1,),
        )
        self.matches = [self.high_confidence_match, self.low_confidence_match]

    def test_no_matches(self):
        assert not FDH("dummy").get_best_matches()

    def test_winner(self):
        best = FDH("dummy", potential_matches=self.matches).get_best_matches()
        assert len(best) == 1
        assert best[0] is self.high_confidence_match

    def test_multiple_equal(self):
        matches = [self.high_confidence_match, self.high_confidence_match]
        fdh = FDH("dummy", potential_matches=matches)
        assert fdh.potential_matches == matches

    def test_single(self):
        matches = [self.low_confidence_match]
        best = FDH("test_enzyme", potential_matches=matches).get_best_matches()
        assert best == matches
