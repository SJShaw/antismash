# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from antismash.common import fasta, path
from antismash.common.test.helpers import DummyCDS
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
    pyrrolic,
)


class TestPyrrolic(unittest.TestCase):
    def setUp(self):
        self.pyrrole_enzyme = FDH("bmp2")

        self.strong_pyrrole_hit = HalogenaseHmmResult(
            hit_id="query",
            bitscore=1000,
            query_id="pyrrole_FDH",
            profile=pyrrolic.SPECIFIC_PROFILES[0].path
        )

        self.bad_pyrrole_hit = HalogenaseHmmResult(
            hit_id="query",
            bitscore=1,
            query_id="pyrrole_FDH",
            profile=pyrrolic.SPECIFIC_PROFILES[0].path
        )

    def test_score_below_cutoff(self):
        cds = DummyCDS()
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={}) as patched:
            result = categorize_on_substrate_level(cds, self.pyrrole_enzyme, [self.bad_pyrrole_hit])
            patched.assert_called_once_with(cds.translation, self.bad_pyrrole_hit, pyrrolic.PYRROLE.motifs)
        assert not result

    def test_no_motif(self):
        cds = DummyCDS()
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues", return_value={}) as patched:
            result = categorize_on_substrate_level(cds, self.pyrrole_enzyme, [self.strong_pyrrole_hit])
            patched.assert_called_once_with(cds.translation, self.strong_pyrrole_hit, pyrrolic.PYRROLE.motifs)
        assert not result

    def test_conventional_mono_di(self):
        cds = DummyCDS()
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={"mono_di": "DRSVFW"}) as patched:
            assert categorize_on_substrate_level(cds, self.pyrrole_enzyme,
                                                 [self.strong_pyrrole_hit])
            match = self.pyrrole_enzyme.potential_matches[0]
            assert match.confidence == 1
            assert match.number_of_decorations == "mono_di"
            assert match.profile == "pyrrole_FDH"
            assert match.substrate == "pyrrole"
            patched.assert_called_once_with(cds.translation, self.strong_pyrrole_hit, pyrrolic.PYRROLE.motifs)

    def test_unconventional_mono_di(self):
        cds = DummyCDS()
        name = "unconv_mono_di"
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={name: "YRRNFN"}) as patched:
            assert categorize_on_substrate_level(cds, self.pyrrole_enzyme,
                                                 [self.strong_pyrrole_hit])
            match = self.pyrrole_enzyme.potential_matches[0]
            assert match.confidence == 1
            assert match.number_of_decorations == name
            assert match.profile == "pyrrole_FDH"
            assert match.substrate == "pyrrole"
            patched.assert_called_once_with(cds.translation, self.strong_pyrrole_hit, pyrrolic.PYRROLE.motifs)

    def test_tetra(self):
        cds = DummyCDS()
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={"tetra_mono_di": "RRYFFA"}) as patched:
            assert categorize_on_substrate_level(cds, self.pyrrole_enzyme, [self.strong_pyrrole_hit])
            match = self.pyrrole_enzyme.potential_matches[0]
            assert match.confidence == 1
            assert match.number_of_decorations == "tetra"
            assert match.profile == "pyrrole_FDH"
            assert match.substrate == "pyrrole"
            patched.assert_called_once_with(cds.translation, self.strong_pyrrole_hit, pyrrolic.PYRROLE.motifs)

    def test_residue_mismatch(self):
        residues = {"pyrrole_FDH": "ABCDEF"}
        assert not pyrrolic.PYRROLE.get_matches_from_hit(residues, self.strong_pyrrole_hit)
