# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from antismash.common import (
    secmet,
    subprocessing,
)
from antismash.common.test.helpers import (
    DummyCDS,
    DummyRecord,
    FakeHSPHit,
    FakeHit,
)
from antismash.detection.genefunctions.halogenases import (
    HalogenaseHmmResult,
    FlavinDependentHalogenase as FDH,
    Match,
    specific_analysis,
)
from antismash.detection.genefunctions.halogenases.flavin_dependent import (
    substrate_analysis,
)
from antismash.detection.genefunctions.halogenases.flavin_dependent.substrate_analysis import (
    categorize_on_substrate_level,
    fdh_specific_analysis,
)
from antismash.detection.genefunctions.halogenases.flavin_dependent.subgroups import (
    indolic,
)


def create_flavin_match(profile, confidence=0., consensus_residues="", substrate=None,
                        target_positions=None, number_of_decorations=None,
                        ):
    return Match(
        profile,
        "flavin",
        "flavin-dependent",
        confidence=confidence,
        consensus_residues=consensus_residues,
        substrate=substrate,
        target_positions=target_positions,
        number_of_decorations=number_of_decorations,
    )


def create_motif_residue_mapping(profile):
    return {motif.name: motif.residues for motif in profile.motifs}


class IndolicBase(unittest.TestCase):
    def setUp(self):
        self.test_trp_5_match = create_flavin_match(
            "trp_5_FDH", confidence=1, substrate="tryptophan", target_positions=5,
            number_of_decorations="mono",
        )
        self.test_trp_6_7_match = create_flavin_match(
            "trp_6_7_FDH", confidence=1, substrate="tryptophan",
            target_positions=6, number_of_decorations="mono",
        )

        self.trp_5_hmm_result = HalogenaseHmmResult(
            hit_id="query",
            bitscore=1000,
            query_id="trp_5_FDH",
            profile=indolic.SPECIFIC_PROFILES[0].path,
        )
        self.trp_6_7_hmm_result = HalogenaseHmmResult(
            hit_id="query",
            bitscore=1000,
            query_id="trp_6_7_FDH",
            profile=indolic.SPECIFIC_PROFILES[1].path,
        )

        tryptophan_single_matches = [self.test_trp_6_7_match]
        tryptophan_matches = [self.test_trp_5_match, self.test_trp_6_7_match]

        # Trp-5 halogenase
        self.trp_5_enzyme_with_matches = FDH("mibH", potential_matches=tryptophan_matches)
        # Trp-6 halogenase
        self.trp_enzyme_with_matches = FDH("ktzR", potential_matches=tryptophan_single_matches)
        # Trp-7 halogenase
        self.trp_with_no_matches = FDH("")

    def specific_analysis_test(self, name, fake_hit,
                               fake_hsp, categorize_return_value,
                               best_matches,
                               analysis_function):
        fake_hsp.id = name
        fake_hsp.hits = fake_hit

        for item in fake_hit:
            item.hsps = [fake_hsp]  # TODO not sure why this is necessary

        with patch.object(subprocessing.hmmscan, "run_hmmscan", return_value=[fake_hsp]):
            with patch.object(subprocessing, "run_hmmsearch", return_value=fake_hit):
                with patch.object(substrate_analysis, "categorize_on_substrate_level",
                                  return_value=categorize_return_value):
                    with patch.object(FDH, "get_best_matches", return_value=best_matches):
                        cds = DummyCDS(locus_tag=name)
                        record = DummyRecord(features=[cds])
                        with patch.object(secmet.Record, "get_cds_features_within_regions", return_value=[cds]):
                            return analysis_function(record)


class TestIndolic(IndolicBase):
    def test_strong_trp_5(self):
        cds = DummyCDS()
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={"trp_5_FDH": indolic.TRP_5_MOTIF.residues}) as patched:
            categorize_on_substrate_level(DummyCDS(), self.trp_with_no_matches, [self.trp_5_hmm_result])
            patched.assert_called_once_with(
                cds.translation, self.trp_5_hmm_result, (indolic.TRP_5_MOTIF,),
            )

        match = self.trp_with_no_matches.potential_matches[0]
        assert match.profile == "trp_5_FDH"
        assert match.confidence == 1
        assert match.substrate == "tryptophan"
        assert match.number_of_decorations == "mono"
        assert match.target_positions == (5,)

    def test_weak_trp_5(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={motif.name: motif.residues for motif in indolic.TRP_5.motifs}):
            low_quality_hit = HalogenaseHmmResult("trp_5_FDH", 380, "trp_5_FDH", "trp_5_FDH")
            categorize_on_substrate_level(DummyCDS(), self.trp_with_no_matches, [low_quality_hit])
        match = self.trp_with_no_matches.potential_matches[0]
        assert match.profile == "trp_5_FDH"
        assert match.confidence == 0.5
        assert match.substrate == "tryptophan"
        assert match.number_of_decorations == "mono"
        assert match.target_positions == (5,)

    def test_trp_6(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=create_motif_residue_mapping(indolic.TRP_6)):
            categorize_on_substrate_level(DummyCDS(), self.trp_with_no_matches,
                                          [self.trp_6_7_hmm_result])
        match = self.trp_with_no_matches.potential_matches[0]
        assert match.profile == "trp_6_7_FDH"
        assert match.confidence == 1.0
        assert match.substrate == "tryptophan"
        assert match.number_of_decorations == "mono"
        assert match.target_positions == (6,)

    def test_trp_7(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=create_motif_residue_mapping(indolic.TRP_7)):
            categorize_on_substrate_level(DummyCDS(), self.trp_with_no_matches, [self.trp_6_7_hmm_result])
        match = self.trp_with_no_matches.potential_matches[0]
        assert match.profile == "trp_6_7_FDH"
        assert match.confidence == 1.
        assert match.target_positions == (7,)
        assert match.substrate == "tryptophan"
        assert match.number_of_decorations == "mono"

    @patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                  return_value=create_motif_residue_mapping(indolic.TRP_5))
    def test_categorise_substrate_no_match(self, _patched_retrieve_fdh_signature_residues):
        result = categorize_on_substrate_level(DummyCDS(), self.trp_with_no_matches, [])
        assert result is None

    @patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                  return_value=create_motif_residue_mapping(indolic.TRP_5))
    def test_categorise_substrate_good_match(self, _patched_consensus_sig):
        cds = DummyCDS()
        result = categorize_on_substrate_level(cds, FDH(cds.get_name()), [self.trp_5_hmm_result])
        assert isinstance(result, FDH)
        assert result.cds_name == cds.get_name()
        assert result.potential_matches == [
            create_flavin_match(profile="trp_5_FDH", confidence=1.0,
                                consensus_residues=indolic.TRP_5_MOTIF.residues,
                                substrate="tryptophan", target_positions=(5,),
                                number_of_decorations="mono")
        ]


class TestSpecificAnalysis(IndolicBase):
    def test_one_best_match(self):
        result = self.specific_analysis_test("ktzR",
                                             FakeHit(1, 2, 900, "ktzR"),
                                             FakeHSPHit("trp_6_7_FDH", "ktzR", bitscore=1500),
                                             self.trp_enzyme_with_matches,
                                             [self.test_trp_6_7_match],
                                             fdh_specific_analysis)
        assert len(result) == 1
        assert result[0].family == FDH.family
        assert result[0].cds_name == "ktzR"
        assert result[0].cofactor == FDH.cofactor
        assert result[0].potential_matches
        best = result[0].get_best_matches()
        assert best[0].substrate == "tryptophan"
        assert best[0].number_of_decorations == "mono"
        assert best[0].target_positions == 6

    def test_more_best_match(self):
        result = self.specific_analysis_test("mibH", FakeHit(1, 2, 500, "mibH"),
                                             FakeHSPHit("trp_5_FDH", "mibH", bitscore=800),
                                             self.trp_5_enzyme_with_matches,
                                             [self.test_trp_5_match, self.test_trp_6_7_match],
                                             fdh_specific_analysis)
        assert len(result) == 1
        assert result[0].potential_matches
        assert result[0].consensus_residues is None
        assert not result[0].number_of_decorations
        assert result[0].target_positions is None

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
        result = self.specific_analysis_test(name, FakeHit(1, 2, 200, "unconventional_FDH"),
                                             FakeHSPHit("unconventional_FDH", name, bitscore=200),
                                             FDH(name), [],
                                             fdh_specific_analysis)
        assert len(result) == 1
        assert not result[0].consensus_residues
        assert not result[0].potential_matches

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="WIWVIRYGMIGDAASVIDAYYSQGVSLALVT")
    def test_conventional(self, _patched_extract_residues):
        name = "CtoA"
        result = self.specific_analysis_test(name,
                                             FakeHit(1, 2, 200, "all_general_FDH"),
                                             FakeHSPHit("all_general_FDH", name, bitscore=200),
                                             FDH(name), [],
                                             fdh_specific_analysis)
        assert len(result) == 1

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="WIWVIRYGMIGDAASVIDAYYSQGVSLALVT")
    def test_good_match(self, _patched_extract_residues):
        name = "CtoA"
        result = self.specific_analysis_test(name, FakeHit(1, 2, 200, "all_general_FDH"),
                                             FakeHSPHit("all_general_FDH", name, bitscore=200),
                                             FDH(name), [],
                                             specific_analysis)
        assert len(result) == 1
        assert result[0].consensus_residues == {"W.W.I.": "WIWVIR"}

    @patch.object(substrate_analysis, "extract_residues", return_value="VALAMI")
    def test_hits_but_no_motifs(self, _patched_extract_residues):
        name = "VatD"
        result = self.specific_analysis_test(name, FakeHit(1, 2, 200, "unconventional_FDH"),
                                             FakeHSPHit("unconventional_FDH", name, bitscore=200),
                                             FDH(name), [],
                                             specific_analysis)
        assert len(result) == 1
        assert result[0].consensus_residues == {}

    @patch.object(substrate_analysis, "extract_residues", return_value="")
    def test_no_matches(self, _patched_extract_residues):
        name = "VatD"
        result = self.specific_analysis_test(name,
                                             FakeHit(1, 2, 200, "unconventional_FDH"),
                                             FakeHSPHit("unconventional_FDH", name, bitscore=200),
                                             FDH(name), [], fdh_specific_analysis)
        assert len(result) == 1
        assert not result[0].potential_matches
