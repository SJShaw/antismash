# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

""" KtzQ: Trp-7
    KtzR: Trp-6
    mibH: Trp-5
    ChlB4: Orsellinic-like
    End30: Hpg
    BhaA: Tyrosine
"""

import json
import unittest
from unittest.mock import patch

from antismash.common import (
    fasta,
    path,
    subprocessing,
    utils,
)
from antismash.common.secmet.test.helpers import (
    DummyFeature,
    DummyCDS,
)
from antismash.common.signature import HmmSignature
from antismash.common.test.helpers import FakeHSPHit, FakeHit
from antismash.detection.genefunctions.halogenases.flavin_dependent.subgroups import (
    indolic,
)

from antismash.detection.genefunctions.halogenases import (
    HalogenaseHmmResult,
    FlavinDependentHalogenase as FDH,
    Match,
)

from antismash.detection.genefunctions.halogenases.flavin_dependent import substrate_analysis
from antismash.detection.genefunctions.halogenases.flavin_dependent.substrate_analysis import (
    run_halogenase_phmms,
    extract_residues,
    categorize_on_substrate_level,
)

TRANSLATIONS = fasta.read_fasta(path.get_full_path(__file__, "data", "translations.fasta"))
TRANSLATIONS = {key.rsplit("|", 1)[-1]: value for key, value in TRANSLATIONS.items()}


def create_motif_residue_mapping(profile):
    return {motif.name: motif.residues for motif in profile.motifs}


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
        value.hsps = [FakeHSPHit("foo", "foo", bitscore=250)]

    negative_test_halogenase_hmms_by_id = run_halogenase_phmms("", [signature])
    assert not negative_test_halogenase_hmms_by_id

    for value in _patched.return_value:
        value.hsps = [FakeHSPHit("foo", "foo", bitscore=1000)]

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


def test_search_signature_residues():
    target_positions = indolic.TRP_6_MOTIF.positions
    profile = indolic.TRP_6

    hmm_result = HalogenaseHmmResult(
        hit_id="query",
        bitscore=1000,
        query_id=profile.profile_name,
        profile="dummy_path",
    )

    with patch.object(subprocessing.hmmpfam, "run_hmmpfam2") as run_hmmpfam2:
        run_hmmpfam2.return_value = []
        signature_residues = extract_residues(TRANSLATIONS["ktzR"],
                                              target_positions, hmm_result)
        assert signature_residues is None

        run_hmmpfam2.return_value = [FakeHit(1, 2, 1000, "foo")]
        # checking if hit_id != query_id it breaks
        for result in run_hmmpfam2.return_value:
            result.hsps = [FakeHSPHit("foo", hit_id="foo")]
            for hit in result.hsps:
                hit_query = DummyFeature()
                hit_profile = DummyFeature()
                hit_query.seq = "ghjvghjkbln"
                hit_profile.seq = "xfdhcgkbjlnkml"
                hit.aln = [hit_profile, hit_query]

        signature_residues = extract_residues(list(TRANSLATIONS.values())[0],
                                              target_positions, hmm_result)
        assert signature_residues is None

        # checking if hit_id == query_id it runs til the reference target_positions extraction
        for result in run_hmmpfam2.return_value:
            result.hsps = [FakeHSPHit("trp_6_7_FDH", hit_id="trp_6_7_FDH")]
            for hit in result.hsps:
                hit_query = DummyFeature()
                hit_profile = DummyFeature()
                hit_query.seq = TRANSLATIONS["ktzR"]
                hit_profile.seq = "xfdhcgkbjlnkml"
                hit.aln = [hit_profile, hit_query]

        with patch.object(utils, "extract_by_reference_positions", return_value="dummy") as patched:
            signature_residues = extract_residues(TRANSLATIONS["ktzR"], target_positions, hmm_result)
            assert signature_residues == "dummy"
            patched.assert_called_once_with("xfdhcgkbjlnkml", TRANSLATIONS["ktzR"], list(target_positions))
