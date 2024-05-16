# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import json
import unittest
from unittest.mock import Mock, patch

from antismash.common import secmet, subprocessing, utils
from antismash.common.secmet.test.helpers import DummyCDS, DummyRecord
from antismash.common.signature import HmmSignature
from antismash.common.test.helpers import FakeHSPHit, FakeHit
from antismash.detection.genefunctions.halogenases.data_structures import (
    HalogenaseHmmResult,
    FlavinDependentHalogenase as _FDH,
    Match,
    MotifDetails,
)
from antismash.detection.genefunctions.halogenases.flavin_dependent import substrate_analysis
from antismash.detection.genefunctions.halogenases.flavin_dependent.substrate_analysis import (
    categorize_on_substrate_level,
    extract_residues,
    fdh_specific_analysis,
    run_halogenase_phmms,
)


class FDH(_FDH):
    def __init__(self, name="dummy", conventionality_residues={"dummy": "ABCDEF"}, potential_matches=None):
        super().__init__(name, conventionality_residues, potential_matches or [])


class BlockHmmer(unittest.TestCase):
    def setUp(self):
        self.hmmpfam = patch.object(subprocessing, "run_hmmpfam2", side_effect=ValueError("hmmpfam2 should not run"))
        self.hmmpfam.start()
        self.hmmscan = patch.object(subprocessing.hmmscan, "run_hmmscan", side_effect=ValueError("hmmscan should not run"))
        self.hmmscan.start()
        self.hmmsearch = patch.object(subprocessing.hmmsearch, "run_hmmsearch", side_effect=ValueError("hmmsearch should not run"))
        self.hmmsearch.start()

    def tearDown(self):
        self.hmmpfam.stop()
        self.hmmscan.stop()
        self.hmmsearch.stop()


class TestComponents(BlockHmmer):
    def test_categorisation_with_no_hits(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues", return_value={}):
            result = categorize_on_substrate_level(DummyCDS(), FDH("dummy"), [])
        assert result is None

    def test_conversion_methods(self):
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
    def test_run_halogenase_phmms(self, _patched):
        signature = HmmSignature("foo", "description", 300, "dummy_path")
        for value in _patched.return_value:
            value.hsps = [FakeHSPHit("foo", "query", bitscore=250)]

        negative_test_halogenase_hmms_by_id = run_halogenase_phmms("", [signature])
        assert not negative_test_halogenase_hmms_by_id

        for value in _patched.return_value:
            value.hsps = [FakeHSPHit("foo", "query_name", bitscore=1000)]

        results = run_halogenase_phmms("", [signature])["query_name"]
        assert results
        for hit in results:
            assert isinstance(hit, HalogenaseHmmResult)

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="PREFIXWIWVIRYGMIGDAASVIDAYYSQGVSLALVT")
    def test_conventional_non_specific(self, _patched_extract):
        name = "CtoA"
        result = substrate_analysis.categorize_on_consensus_level(
            DummyCDS(locus_tag=name),
            [],
            [HalogenaseHmmResult(name, 200, "all_conventional_FDH", "flavin-dependent")],
        )
        assert result.conventionality_residues == {"W.W.I.": "WIWVIR"}
        assert not result.potential_matches


class TestGetBest(BlockHmmer):
    def setUp(self):
        super().setUp()
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


class TestExtract(BlockHmmer):
    def test_no_alignment(self):
        with patch.object(utils, "extract_from_alignment", side_affect=RuntimeError("should not have been called")):
            with patch.object(substrate_analysis, "get_alignment_against_profile", return_value=None):
                assert extract_residues("MAGIC", [1, 5, 6], Mock()) is None

    def test_alignment_passed(self):
        positions = [1, 5, 6]
        alignment = Mock()
        with patch.object(utils, "extract_from_alignment", return_value="dummy") as patched_util:
            with patch.object(substrate_analysis, "get_alignment_against_profile", return_value=alignment):
                result = extract_residues("MAGIC", positions, Mock())
                assert result == "dummy"
            patched_util.assert_called_once_with(alignment, positions)

    def test_invalid_positions(self):
        with self.assertRaisesRegex(ValueError, "without positions"):
            extract_residues("MAGIC", [], Mock())


class TestSpecificAnalysis(BlockHmmer):
    def specific_analysis_test(self, name, *, bitscore=200, profile_id="dummy_profile"):
        cds = DummyCDS(locus_tag=name)
        record = DummyRecord(features=[cds])

        scan_result = [Mock(id=name, hits=True)]
        with patch.object(subprocessing.hmmscan, "run_hmmscan", return_value=scan_result):
            with patch.object(secmet.Record, "get_cds_features_within_regions", return_value=[cds]):
                hmm_result = HalogenaseHmmResult(name, bitscore, profile_id, "dummy_path")
                with patch.object(substrate_analysis, "run_halogenase_phmms", return_value={name: [hmm_result]}):
                    return fdh_specific_analysis(record)

    @patch.object(subprocessing.hmmscan, "run_hmmscan", return_value=[])
    def test_no_hits(self, _patched_hmmscan):
        cds = DummyCDS()
        record = DummyRecord(features=[cds])
        with patch.object(secmet.Record, "get_cds_features_within_regions", return_value=[cds]):
            results = fdh_specific_analysis(record)

        assert results == []

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="VALAMIVALAMI")
    def test_unconventional(self, _patched_extract_residues):
        name = "VatD"
        result = self.specific_analysis_test(name, bitscore=200, profile_id="unconventional_FDH")
        assert len(result) == 1
        assert result[0].conventionality_residues == {}
        assert not result[0].is_conventional()
        assert not result[0].potential_matches

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="WIWVIRYGMIGDAASVIDAYYSQGVSLALVT")
    def test_conventional(self, _patched_extract_residues):
        name = "CtoA"
        result = self.specific_analysis_test(name, bitscore=200, profile_id="all_general_FDH")
        assert len(result) == 1
        assert result[0].conventionality_residues == {"W.W.I.": "WIWVIR"}


class TestDataStructures(unittest.TestCase):
    def new_via_json(self, item):
        data = item.to_json()
        # ensure it's proper json by converting to a string and back
        data = json.loads(json.dumps(data))
        new = type(item).from_json(data)
        assert new.to_json() == data
        return new

    def test_match_conversion(self):
        match = Match(profile="dummy_profile", cofactor="some_cofactor", family="a family", confidence=0.2, consensus_residues="MAGIC",
                      substrate="a substrate", target_positions=(5,7), number_of_decorations="two_or_three_maybe")

        rebuilt = self.new_via_json(match)
        assert rebuilt == match

    def test_match_equality(self):
        values = {
            "profile": "dummy_profile",
            "cofactor": "some_cofactor",
            "family": "a family",
            "confidence": 0.2,
            "consensus_residues": "MAGIC",
        }

        match = Match(**values)
        assert Match(**values) == match

        # and with any modification, equality should fail
        for key, val in values.items():
            copy = dict(values)
            copy[key] = val * 2  # arbitrary, any difference will do
            assert Match(**copy) != match

    def test_motif_string_equality(self):
        motif = MotifDetails(name="motif name", positions=(1, 6), residues="MT", substrate="something", decorations="dec")
        assert motif == "MT" == motif.residues

    def test_motif_class_equality(self):
        # only the residues are meaningful for comparison
        residues = "HEY"
        first = MotifDetails(name="motif name", positions=(1, 2, 6), residues=residues, substrate="something", decorations="dec")
        different = MotifDetails(name="motif name", positions=(1, 6), residues="ME", substrate="something", decorations="dec")
        assert first.residues != different.residues
        assert first != different
        same = MotifDetails(name="new name", positions=(5, 6, 12), residues=residues, substrate="a", decorations="other")
        assert first.residues == same.residues
        assert first == same
