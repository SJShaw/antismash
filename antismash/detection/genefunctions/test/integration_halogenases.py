# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import json
from tempfile import NamedTemporaryFile
import unittest

from helperlibs.bio import seqio

import antismash
from antismash.common import secmet
from antismash.common.secmet.qualifiers import GeneFunction
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config, get_config, update_config
from antismash.detection.genefunctions import prepare_data
from antismash.detection.genefunctions.halogenases import specific_analysis as halo_analysis

# this particular one is ChlB4 from... TODO?
CHLB4 = """
MQPDFDAAIVGGGPAGSAMASYLAEAGLSVAVFESEMFPRPHIGESLVPATMPVLDEIGV
MPDIEAAGFPKKYGAAWTSAESRDVPHNGFTGLDHDFKAAEVMFVERDQPGVHRDYTFHV
DRGKFDLILLKHAESRGAQVFQKTRVLKADFDTDPDLVTLNCRLGPRTLDFTTRMVIDAS
GRQTMLGNQLKVKVPDPVFNQYAIHAWFEGLDRTAMALDPAKRDYIYVHFLPLEDTWMWQ
IPITDTITSVGVVTQKHRFKAASADREKFFWDIVSSRKDIYDALQKAERIRPFKAEGDYS
YAMRQICGDRFLLIGDAARFVDPIFSSGVSVALNSARLAAKDVIAAHRAGDFRKESFATY
EEKLRRAVRNWYEFISVYYRLNILFTAFVQDPRYRIDVLKMLQGDFYDGEEPKALKAMRD
LVTKVENDPEHLWHPYLGTLRAPSAAPTF
""".replace("\n", "")


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

    def test_multiple(self):
        record = secmet.Record.from_genbank(helpers.get_path_to_balhymicin_genbank())[0]
        # the real halogenase had best be present
        assert record.get_cds_by_name("bhaA")
        # and add a dummy second one, just to be sure multiple in a record will function correctly
        dummy_loc = secmet.locations.FeatureLocation(10, 10 + len(CHLB4) * 3, 1)
        dummy_halogenase = secmet.CDSFeature(dummy_loc, locus_tag="ChlB4", translation=CHLB4)
        record.add_cds_feature(dummy_halogenase)
        record.add_subregion(secmet.SubRegion(secmet.locations.FeatureLocation(0, len(record)), tool="dummy"))
        record.create_regions()
        assert dummy_halogenase in record.get_cds_features_within_regions()
        results = halo_analysis(record)
        assert len(results) == 2
        assert len(results[0].potential_matches) == 1

        bhaA = results[0]
        assert len(bhaA.potential_matches) == 1
        assert bhaA.potential_matches[0].profile == "tyrosine-like_hpg_FDH"
        assert bhaA.potential_matches[0].confidence == 1.0
        assert bhaA.potential_matches[0].substrates == ["Tyr"]

        chlB4 = results[1]
        assert len(chlB4.potential_matches) == 2
        matches = sorted(chlB4.potential_matches, key=lambda x: x.profile)
        assert matches[0].profile == "cycline_orsellinic_FDH"
        assert matches[0].confidence == 1.0
        assert matches[1].profile == "pyrrole_FDH"
        assert matches[1].confidence == 1.0
