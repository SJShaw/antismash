# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.modules.nrps_pks import parsers
from antismash.modules.nrps_pks.data_structures import Prediction, SubstratePrediction
from antismash.modules.nrps_pks.name_mappings import KNOWN_SUBSTRATES, get_substrate_by_name


class DummySubstratePrediction(SubstratePrediction):
    def __init__(self, method, name):
        super().__init__(method)
        self.substrate = get_substrate_by_name(name)

    def get_substrate_classification(self):
        return [self.substrate]


class TestNrpsConsensus(unittest.TestCase):
    def setUp(self) -> None:
        self.go = parsers.generate_nrps_consensus

    def test_single(self):
        data = [
            ("method_a", "ala")
        ]
        assert "Ala" == self.go(_generate_predictions(data))

    def test_majority(self):
        data = [
            ("method_a", "ala"),
            ("method_b", "gly"),
            ("method_c", "ala"),
        ]
        assert "Ala" == self.go(_generate_predictions(data))

    def test_tie(self):
        data = [
            ("method_a", "ala"),
            ("method_b", "gly"),
        ]
        assert "X" == self.go(_generate_predictions(data))

    def test_tie_with_extra_predictions(self):
        data = [
            ("method_a", "ala"),
            ("method_b", "gly"),
            ("method_c", "ala"),
            ("method_d", "gly"),
            ("method_e", "tyr"),
        ]
        assert "X" == self.go(_generate_predictions(data))

    def test_valid_characters(self):
        for substrate in KNOWN_SUBSTRATES:
            assert set(substrate.norine).issubset(parsers.ALLOWABLE_PREDICTION_CHARACTERS)


def _generate_predictions(data: list[tuple[str, str, str]]) -> dict[str, Prediction]:
    return {
        method: DummySubstratePrediction(method, prediction) for method, prediction in data
    }
