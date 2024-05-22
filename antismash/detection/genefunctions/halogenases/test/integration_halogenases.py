# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

import antismash
from antismash.common import fasta, path, secmet
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config, get_config, update_config
from antismash.detection.genefunctions import prepare_data
from antismash.detection.genefunctions.halogenases import specific_analysis as halo_analysis
from antismash.detection.genefunctions.halogenases.data_structures import (
    SPECIFIC_BASE_CONFIDENCE,
    NON_SPECIFIC_BASE_CONFIDENCE,
    UNCONVENTIONAL_BASE_CONFIDENCE,
)

TRANSLATIONS = fasta.read_fasta(path.get_full_path(__file__, "data", "translations.fasta"))
# convert to just gene name, instead of the full info in each FASTA header
TRANSLATIONS = {key.rsplit("|", 1)[-1]: value for key, value in TRANSLATIONS.items()}


class TestHalongenases(unittest.TestCase):
    def setUp(self):
        options = build_config(self.get_args(), isolated=True,
                               modules=antismash.get_all_modules())
        self.old_config = get_config().__dict__
        self.options = update_config(options)
        prepare_data()

    def tearDown(self):
        destroy_config()
        update_config(self.old_config)

    def get_args(self):
        return ["--minimal", "--enable-genefunctions"]

    def generate_record(self, gene_name):
        features = [helpers.DummyCDS(locus_tag=gene_name, translation=TRANSLATIONS[gene_name])]
        record = helpers.DummyRecord(features=features)
        record.get_cds_features_within_regions = lambda: features
        return record

    def get_single_result(self, gene_name):
        record = self.generate_record(gene_name)
        results = halo_analysis(record)
        assert len(results) == 1
        assert results[0].cds_name == gene_name
        return results[0]

    def test_real(self):
        record = secmet.Record.from_genbank(helpers.get_path_to_balhymicin_genbank())[0]
        # the real halogenase had best be present
        assert record.get_cds_by_name("bhaA")
        # and add a dummy second one, just to be sure multiple in a record will function correctly
        dummy_loc = secmet.locations.FeatureLocation(10, 10 + len(TRANSLATIONS["ChlB4"]) * 3, 1)
        dummy_halogenase = secmet.CDSFeature(dummy_loc, locus_tag="ChlB4", translation=TRANSLATIONS["ChlB4"])
        record.add_cds_feature(dummy_halogenase)
        record.add_subregion(secmet.SubRegion(secmet.locations.FeatureLocation(0, len(record)), tool="dummy"))
        record.create_regions()
        assert dummy_halogenase in record.get_cds_features_within_regions()
        results = halo_analysis(record)

        assert len(results) == 2
        first, second = results

        assert first.cds_name == "bhaA"
        assert second.cds_name == "ChlB4"

        assert len(first.potential_matches) == 1
        assert first.potential_matches[0].profile == "tyrosine-like_hpg_FDH"
        assert first.potential_matches[0].confidence == 1.0
        assert first.potential_matches[0].substrate == "Tyr"

        assert len(second.potential_matches) == 2
        matches = sorted(second.potential_matches, key=lambda x: x.profile)
        assert matches[0].profile == "cycline_orsellinic_FDH"
        assert matches[0].confidence == 1.0
        assert matches[1].profile == "pyrrole_FDH"
        assert matches[1].confidence == 1.0

    def test_unconventional(self):
        result = self.get_single_result("VatD")
        assert not result.potential_matches
        assert not result.conventionality_residues
        assert not result.is_conventional()
        assert result.confidence == UNCONVENTIONAL_BASE_CONFIDENCE

    def test_generic_conventional(self):
        result = self.get_single_result("CtoA")
        assert not result.potential_matches
        assert result.conventionality_residues == {"W.W.I.": "WIWVIR"}
        assert result.is_conventional()
        assert result.confidence == NON_SPECIFIC_BASE_CONFIDENCE

    def test_tryptophan_7(self):
        result = self.get_single_result("ktzQ")

        self.fail("need some kind of fallback profile system, otherwise it gets two hits")

        print("\n".join(str(match) for match in result.potential_matches))
        assert len(result.potential_matches) == 1
        assert result.is_conventional()

        match = result.potential_matches[0]
        assert match.profile == "trp_6_7_FDH"
        assert match.confidence == 1.0  # TODO: but should it be?
        assert match.substrate == "tryptophan"
        assert match.target_positions == (7,)
        assert match.number_of_decorations == "mono"
        assert result.confidence == 1.0

    def test_tryptophan_6(self):
        result = self.get_single_result("ktzR")

        assert result.confidence == 1.0
        assert len(result.potential_matches) == 1
        assert result.is_conventional()

        match = result.potential_matches[0]
        assert match.profile == "trp_6_7_FDH"
        assert match.confidence == 1.0
        assert match.substrate == "tryptophan"
        assert match.target_positions == (6,)
        assert match.number_of_decorations == "mono"

    def test_tryptophan_5(self):
        result = self.get_single_result("mibH")

        assert result.confidence == 0.5
        assert len(result.potential_matches) == 1
        assert result.is_conventional()

        match = result.potential_matches[0]
        assert match.profile == "trp_5_FDH"
        assert match.confidence == 0.5
        assert match.substrate == "tryptophan"
        assert match.target_positions == (5,)
        assert match.number_of_decorations == "mono"
        self.fail("this is only a weak hit, where's a strong one?")

    def test_pyrrole_tetra(self):
        result = self.get_single_result("bmp2")

        assert result.confidence == 1.0
        assert len(result.potential_matches) == 1
        assert result.is_conventional()

        match = result.potential_matches[0]
        assert match.profile == "pyrrole_FDH"
        assert match.confidence == 1.0
        assert match.substrate == "pyrrole"
        assert match.target_positions == (5,)
        assert match.number_of_decorations == "tetra"

    def test_pyrrole_mono_di(self):
        result = self.get_single_result("HrmQ")

        assert result.confidence == 1.0
        assert len(result.potential_matches) == 1
        assert result.is_conventional()

        match = result.potential_matches[0]
        assert match.profile == "pyrrole_FDH"
        assert match.confidence == 1.0
        assert match.substrate == "pyrrole"
        assert match.target_positions == (5,)
        assert match.number_of_decorations == "mono_di"

    def test_orsellic(self):
        result = self.get_single_result("XanH")

        assert result.confidence == 1.0
        assert result.is_conventional()

        assert len(result.potential_matches) == 1

        first = result.potential_matches[0]
        assert first.profile == "cycline_orsellinic_FDH"
        assert first.confidence == 1.0
        assert first.substrate == "cycline_orsellinic-like"
        assert first.target_positions == (6, 8)
        assert not first.number_of_decorations

    def test_ambiguous_orsellic_and_pyrrole(self):
        result = self.get_single_result("ChlB4")

        assert result.confidence == 1.0
        assert result.is_conventional()

        assert len(result.potential_matches) == 2

        first, second = sorted(result.potential_matches, key=lambda x: x.profile)
        assert first.profile == "cycline_orsellinic_FDH"
        assert first.confidence == 1.0
        assert first.substrate == "cycline_orsellinic-like"
        assert first.target_positions == (6, 8)
        assert not first.number_of_decorations

        assert second.profile == "pyrrole_FDH"
        assert second.confidence == 1.0
        assert second.substrate == "pyrrole"
        assert first.target_positions == (6, 8)
        assert second.number_of_decorations == "mono_di"

    def check_hpg(self):
        result = self.get_single_result("End30")

        assert result.confidence == 1.0
        assert len(result.potential_matches) == 1
        assert result.is_conventional()

        match = result.potential_matches[0]
        assert match.profile == "tyr-like"
        assert match.confidence == 1.0
        assert match.substrate == "Tyr"
        assert not match.target_positions
        assert not match.number_of_decorations

    def check_tyr(self):
        result = self.get_single_result("bhaA")
        assert len(result.potential_matches) == 1
        assert result.is_conventional()
        assert result.conventionality_residues
        assert result.confidence == 1.0

        match = result.potential_matches[0]
        assert match.profile == "tyrosine-like_hpg_FDH"
        assert match.potential_matches[0].confidence == 1.0
        assert match.potential_matches[0].substrate == "Tyr"
