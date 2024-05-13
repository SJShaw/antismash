# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from antismash.common import fasta, path, subprocessing
from antismash.common.test.helpers import (
    DummyCDS,
    DummyFeature,
    FakeHSPHit,
    FakeHit,
)
from antismash.detection.genefunctions.halogenases import (
    HalogenaseHmmResult,
    FlavinDependentHalogenase as FDH,
)
from antismash.detection.genefunctions.halogenases.flavin_dependent import (
    substrate_analysis,
)
from antismash.detection.genefunctions.halogenases.flavin_dependent.substrate_analysis import (
    categorize_on_substrate_level,
)
from antismash.detection.genefunctions.halogenases.flavin_dependent.subgroups import (
    phenolic,
)

TRANSLATIONS = fasta.read_fasta(path.get_full_path(__file__, "data", "translations.fasta"))
TRANSLATIONS = {key.rsplit("|", 1)[-1]: value for key, value in TRANSLATIONS.items()}


def create_motif_residue_mapping(profile):
    return {motif.name: motif.residues for motif in profile.motifs}


class TestPhenolic(unittest.TestCase):
    def setUp(self):
        self.tyrosine_hmm_result = HalogenaseHmmResult(
            hit_id="tyrosine-like_hpg_FDH",
            bitscore=1000,
            query_id="tyrosine-like_hpg_FDH",
            profile=phenolic.SPECIFIC_PROFILES[0].path,
        )
        self.less_confident_tyrosine_hmm_result = HalogenaseHmmResult(
            hit_id="tyrosine-like_hpg_FDH",
            bitscore=310,
            query_id="tyrosine-like_hpg_FDH",
            profile=phenolic.SPECIFIC_PROFILES[0].path
        )
        self.hpg_hmm_result = HalogenaseHmmResult(
            hit_id="tyrosine-like_hpg_FDH",
            bitscore=600,
            query_id="tyrosine-like_hpg_FDH",
            profile=phenolic.SPECIFIC_PROFILES[0].path,
        )
        self.cycline_orsellinic_hmm_result = HalogenaseHmmResult(
            hit_id="cycline_orsellinic_FDH",
            bitscore=600,
            query_id="cycline_orsellinic_FDH",
            profile=phenolic.SPECIFIC_PROFILES[1].path,
        )

        # tyrosine
        self.tyr_enzyme = FDH("BhaA")

        # Hpg
        self.hpg_enzyme = FDH("End30")

        # other phenolic-substrate halogenase (orsellinic-like)
        self.orsellinic_with_no_matches = FDH("ChlB4")

        self.tyrosine_match_not_hpg = {
            "Tyr": "GFQRLGDAGLSGVPSYGADPSGLYW",
            "Hpg": "VALAMI",
        }

    def test_categorising_orsellinic(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=create_motif_residue_mapping(phenolic.ORSELLINIC)):
            categorize_on_substrate_level(DummyCDS(), self.orsellinic_with_no_matches,
                                          [self.cycline_orsellinic_hmm_result])
        match = self.orsellinic_with_no_matches.potential_matches[0]
        assert match.profile == "cycline_orsellinic_FDH"
        assert match.confidence == 1
        assert match.substrate == "cycline_orsellinic-like"
        assert match.target_positions == (6, 8)

    def test_categorising_hpg(self):
        with patch.object(substrate_analysis, "extract_residues",
                          return_value=""):
            assert not categorize_on_substrate_level(DummyCDS(), self.orsellinic_with_no_matches,
                                                     [self.cycline_orsellinic_hmm_result])

        # Tyr and Hpg halogenases
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=create_motif_residue_mapping(phenolic.TYR_HPG)):
            categorize_on_substrate_level(DummyCDS(translation=TRANSLATIONS["End30"]), self.hpg_enzyme,
                                          [self.hpg_hmm_result])
        assert len(self.hpg_enzyme.potential_matches) == 2

        hpg_match = self.hpg_enzyme.potential_matches[0]
        assert hpg_match.profile == "tyrosine-like_hpg_FDH"
        assert hpg_match.confidence == 1.
        assert hpg_match.substrate == "Hpg"
        assert hpg_match.target_positions == (6, 8)

        tyr_match = self.hpg_enzyme.potential_matches[1]
        assert tyr_match.profile == "tyrosine-like_hpg_FDH"
        assert tyr_match.confidence == 0.8
        assert tyr_match.substrate == "Tyr"
        assert tyr_match.target_positions == (6, 8)

    def test_categorising_tyrosine(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=self.tyrosine_match_not_hpg):
            assert categorize_on_substrate_level(DummyCDS(translation=TRANSLATIONS["BhaA"]), self.tyr_enzyme,
                                                 [self.less_confident_tyrosine_hmm_result])
        match = self.tyr_enzyme.potential_matches[0]
        assert match.profile == "tyrosine-like_hpg_FDH"
        assert match.confidence == 0.5
        assert match.substrate == "Tyr"
        assert match.target_positions == (6, 8)

        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={}):
            assert categorize_on_substrate_level(DummyCDS(), self.tyr_enzyme,
                                                 [self.tyrosine_hmm_result]) is None

    def test_retrieve_fdh_signature_residues(self):
        with patch.object(subprocessing.hmmpfam, "run_hmmpfam2") as run_hmmpfam2:
            run_hmmpfam2.return_value = [FakeHit(1, 2, 1000, "foo")]
            # checking if hit_id != query_id it breaks
            for result in run_hmmpfam2.return_value:
                result.hsps = [FakeHSPHit("BhaA", hit_id="tyrosine-like_hpg_FDH")]
                for hit in result.hsps:
                    hit_query = DummyFeature()
                    hit_profile = DummyFeature()
                    hit_query.seq = TRANSLATIONS["BhaA"]
                    hit_profile.seq = TRANSLATIONS["BhaA"]
                    hit.aln = [hit_profile, hit_query]
            cds = DummyCDS(locus_tag="BhaA", translation=TRANSLATIONS["BhaA"])
            assert " " not in cds.translation
            residues = substrate_analysis.retrieve_fdh_signature_residues(
                cds.translation, self.hpg_hmm_result, [phenolic.TYROSINE_LIKE_MOTIF, phenolic.HPG_MOTIF],
            )
        assert isinstance(residues, dict)
        assert residues["Tyr"] is not None and residues["Hpg"] is not None
