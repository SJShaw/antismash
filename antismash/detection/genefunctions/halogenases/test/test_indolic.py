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
)
from antismash.detection.genefunctions.halogenases import (
    HalogenaseHmmResult,
    FlavinDependentHalogenase as _FDH,
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


class FDH(_FDH):
    def __init__(self, name="dummy", conventionality_residues="ABCDEF", potential_matches=None):
        super().__init__(name, conventionality_residues, potential_matches or [])


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
