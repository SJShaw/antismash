# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from antismash.common import (
    fasta,
    path,
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

TRANSLATIONS = fasta.read_fasta(path.get_full_path(__file__, "data", "translations.fasta"))
TRANSLATIONS = {key.rsplit("|", 1)[-1]: value for key, value in TRANSLATIONS.items()}


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
        self.trp_with_no_matches = FDH("ktzQ")

    def specific_analysis_test(self, name, fake_hit,
                               fake_hsps, categorize_return_value,
                               get_best_matches_return_value,
                               analysis_function):
        with patch.object(subprocessing.hmmscan, "run_hmmscan",
                          return_value=[fake_hsps]) as _patched_hmmscan:
            for value in _patched_hmmscan.return_value:
                value.id = name
                value.hits = fake_hit
            with patch.object(subprocessing, "run_hmmsearch",
                              return_value=fake_hit) as _patched_run_hmmsearch:
                for value in _patched_run_hmmsearch.return_value:
                    value.hsps = [fake_hsps]
                with patch.object(substrate_analysis, "categorize_on_substrate_level",
                                  return_value=categorize_return_value):
                    with patch.object(FDH, "get_best_matches",
                                      return_value=get_best_matches_return_value):
                        cds = DummyCDS(locus_tag=name,
                                       translation=TRANSLATIONS[name])
                        record = DummyRecord(features=[cds])
                        with patch.object(secmet.Record, "get_cds_features_within_regions",
                                          return_value=[cds]):
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
        cds = DummyCDS(locus_tag="mibH", translation=TRANSLATIONS["mibH"])
        result = categorize_on_substrate_level(cds, FDH("mibH"), [self.trp_5_hmm_result])
        assert isinstance(result, FDH)
        assert result.cds_name == "mibH"
        assert result.potential_matches == [
            create_flavin_match(profile="trp_5_FDH", confidence=1.0,
                                consensus_residues=indolic.TRP_5_MOTIF.residues,
                                substrate="tryptophan", target_positions=(5,),
                                number_of_decorations="mono")
        ]


class TestSpecificAnalysis(IndolicBase):
    def test_one_best_match(self):
        positive_test = self.specific_analysis_test("ktzR",
                                                    FakeHit(1, 2, 900, "ktzR"),
                                                    FakeHSPHit("trp_6_7_FDH", "ktzR", bitscore=1500),
                                                    self.trp_enzyme_with_matches,
                                                    [self.test_trp_6_7_match],
                                                    fdh_specific_analysis)
        assert positive_test[0].family == FDH.family
        assert positive_test[0].cds_name == "ktzR"
        assert positive_test[0].cofactor == FDH.cofactor
        assert positive_test[0].substrates == ("tryptophan",)
        assert positive_test[0].number_of_decorations == "mono"
        assert positive_test[0].potential_matches
        assert positive_test[0].target_positions == 6

    def test_more_best_match(self):
        positive_test = self.specific_analysis_test("mibH", FakeHit(1, 2, 500, "mibH"),
                                                    FakeHSPHit("trp_5_FDH", "mibH", bitscore=800),
                                                    self.trp_5_enzyme_with_matches,
                                                    [self.test_trp_5_match, self.test_trp_6_7_match],
                                                    fdh_specific_analysis)
        assert len(positive_test) == 1
        assert positive_test[0].potential_matches
        assert positive_test[0].consensus_residues is None
        assert not positive_test[0].number_of_decorations
        assert positive_test[0].target_positions is None

    @patch.object(subprocessing.hmmscan, "run_hmmscan", return_value=[])
    def test_no_hits(self, _patched_hmmscan):
        cds = DummyCDS(locus_tag="mibH",
                       translation=TRANSLATIONS["mibH"])
        record = DummyRecord(features=[cds])
        with patch.object(secmet.Record, "get_cds_features_within_regions", return_value=[cds]):
            results = fdh_specific_analysis(record)

        assert results == []

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="VALAMIVALAMI")
    def test_unconventional_fdh_specific_analysis(self, _patched_extract_residues):
        positive_test = self.specific_analysis_test("VatD", FakeHit(1, 2, 200, "unconventional_FDH"),
                                                    FakeHSPHit("unconventional_FDH", "VatD", bitscore=200),
                                                    FDH("VatD"), [],
                                                    fdh_specific_analysis)
        assert isinstance(positive_test[0], FDH)
        assert positive_test[0].substrates is None
        assert not positive_test[0].consensus_residues
        assert not positive_test[0].potential_matches

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="WIWVIRYGMIGDAASVIDAYYSQGVSLALVT")
    def test_conventional_fdh_specific_analysis(self, _patched_extract_residues):
        name = "CtoA"
        positive_test = self.specific_analysis_test(name,
                                                    FakeHit(1, 2, 200, "all_general_FDH"),
                                                    FakeHSPHit("all_general_FDH", name, bitscore=200),
                                                    FDH(name), [],
                                                    fdh_specific_analysis)

        assert not positive_test[0].substrates
        assert not positive_test[0].target_positions

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="WIWVIRYGMIGDAASVIDAYYSQGVSLALVT")
    def test_random(self, _patched_extract_residues):
        name = "CtoA"
        result = substrate_analysis.categorize_on_consensus_level(
            DummyCDS(locus_tag=name, translation=TRANSLATIONS[name]),
            {},
            [HalogenaseHmmResult(name, 200, "all_conventional_FDH", "flavin-dependent")],
        )
        assert result.consensus_residues == {"W.W.I.": "WIWVIR"}
        assert result.substrates is None
        assert not result.potential_matches

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="WIWVIRYGMIGDAASVIDAYYSQGVSLALVT")
    def test_specific_analysis(self, _patched_extract_residues):
        name = "CtoA"
        categorized_halogenase = self.specific_analysis_test(name, FakeHit(1, 2, 200, "all_general_FDH"),
                                                             FakeHSPHit("all_general_FDH", name, bitscore=200),
                                                             FDH(name), [],
                                                             specific_analysis)
        assert categorized_halogenase[0].consensus_residues == {"W.W.I.": "WIWVIR"}

    @patch.object(substrate_analysis, "extract_residues", return_value="VALAMI")
    def test_negative_get_conserved_motifs(self, _patched_extract_residues):
        name = "VatD"
        categorized_halogenase = self.specific_analysis_test(name, FakeHit(1, 2, 200, "unconventional_FDH"),
                                                             FakeHSPHit("unconventional_FDH", name, bitscore=200),
                                                             FDH(name), [],
                                                             specific_analysis)
        assert categorized_halogenase
        assert categorized_halogenase[0].consensus_residues == {}

    @patch.object(substrate_analysis, "extract_residues", return_value="")
    def test_negative_enzyme_hits_fasta(self, _patched_extract_residues):
        name = "VatD"
        no_potential_matches = self.specific_analysis_test(name,
                                                           FakeHit(1, 2, 200, "unconventional_FDH"),
                                                           FakeHSPHit("unconventional_FDH", name, bitscore=200),
                                                           FDH(name), [], fdh_specific_analysis)
        assert not no_potential_matches[0].potential_matches
