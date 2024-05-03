# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Identification and classification of halogenases """

from dataclasses import dataclass, field

from .data_structures import (
    HalogenaseHmmResult,
    FlavinDependentHalogenase,
    Match,
)
from .flavin_dependent.substrate_analysis import fdh_specific_analysis as specific_analysis
