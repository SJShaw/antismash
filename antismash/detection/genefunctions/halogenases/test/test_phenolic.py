# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import Mock, patch

from antismash.common import fasta, path, subprocessing
from antismash.common.test.helpers import (
    DummyCDS,
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
        self.hpg_hmm_result = HalogenaseHmmResult(
            hit_id="query_name",
            bitscore=600,
            query_id="tyrosine-like_hpg_FDH",
            profile=phenolic.SPECIFIC_PROFILES[0].path,
        )
        self.cycline_orsellinic_hmm_result = HalogenaseHmmResult(
            hit_id="query_name",
            bitscore=600,
            query_id="cycline_orsellinic_FDH",
            profile=phenolic.SPECIFIC_PROFILES[1].path,
        )

    def test_categorising_orsellinic(self):
        fdh = FDH("")
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=create_motif_residue_mapping(phenolic.ORSELLINIC)):
            categorize_on_substrate_level(DummyCDS(), fdh, [self.cycline_orsellinic_hmm_result])
        assert len(fdh.potential_matches) == 1
        match = fdh.potential_matches[0]
        assert match.profile == "cycline_orsellinic_FDH"
        assert match.confidence == 1
        assert match.substrate == "cycline_orsellinic-like"
        assert match.target_positions == (6, 8)

    def test_no_result(self):
        fdh = FDH("")
        with patch.object(substrate_analysis, "extract_residues",
                          return_value=""):
            result = categorize_on_substrate_level(DummyCDS(), fdh, [self.cycline_orsellinic_hmm_result])
        assert result is None

    def test_both_tyr_and_hpg(self):
        fdh = FDH("")
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=create_motif_residue_mapping(phenolic.TYR_HPG)):
            categorize_on_substrate_level(DummyCDS(), fdh, [self.hpg_hmm_result])

        hpg_match, tyr_match = fdh.potential_matches

        assert hpg_match.profile == "tyrosine-like_hpg_FDH"
        assert hpg_match.confidence == 1.
        assert hpg_match.substrate == "Hpg"
        assert hpg_match.target_positions == (6, 8)

        assert tyr_match.profile == "tyrosine-like_hpg_FDH"
        assert tyr_match.confidence == 0.8
        assert tyr_match.substrate == "Tyr"
        assert tyr_match.target_positions == (6, 8)

    def test_weak_hit_good_motif(self):
        hit = HalogenaseHmmResult(
            hit_id="query",
            bitscore=310,  # weak score
            query_id="tyrosine-like_hpg_FDH",
            profile="dummy_path",
        )
        fdh = FDH("")
        motifs_tyrosine_match_not_hpg = {
            "Tyr": "GFQRLGDAGLSGVPSYGADPSGLYW",
            "Hpg": "VALAMI",
        }
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value=motifs_tyrosine_match_not_hpg):
            assert categorize_on_substrate_level(DummyCDS(), fdh, [hit])
        assert len(fdh.potential_matches) == 1
        match = fdh.potential_matches[0]
        assert match.profile == "tyrosine-like_hpg_FDH"
        assert match.confidence == 0.5
        assert match.substrate == "Tyr"
        assert match.target_positions == (6, 8)

    def test_good_hit_bad_motif(self):
        hit = HalogenaseHmmResult(
            hit_id="tyrosine-like_hpg_FDH",
            bitscore=1000,
            query_id="tyrosine-like_hpg_FDH",
            profile=phenolic.SPECIFIC_PROFILES[0].path,
        )
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={}):
            result = categorize_on_substrate_level(DummyCDS(), FDH(""), [hit])
        assert result is None

    def test_retrieve_fdh_signature_residues(self):
        alignment = FakeHSPHit("dummy", hit_id="tyrosine-like_hpg_FDH")
        alignment.aln = [Mock(seq=TRANSLATIONS["BhaA"]), Mock(seq=TRANSLATIONS["BhaA"])]
        hit = FakeHit(1, 2, 1000, "foo")
        hit.hsps = [alignment]

        with patch.object(subprocessing.hmmpfam, "run_hmmpfam2", return_value=[hit]):
            residues = substrate_analysis.retrieve_fdh_signature_residues(
                "MAGIC", self.hpg_hmm_result, [phenolic.TYROSINE_LIKE_MOTIF, phenolic.HPG_MOTIF],
            )
        assert residues == {
            "Tyr": "MRFGGVTDDNFSWQADVKQQYRANV",
            "Hpg": "MARFGGVTDDVFKNFSWQGCADVKQQYRANV",
        }
