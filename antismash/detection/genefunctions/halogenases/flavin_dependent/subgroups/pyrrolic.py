# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from pathlib import Path

from antismash.common.path import get_full_path
from antismash.detection.genefunctions.halogenases.data_structures import (
    FlavinDependentHalogenase,
    HalogenaseHmmResult,
    MotifDetails,
    Profile,
)

MODIFICATION_COUNT_POSITIONS = (110, 111, 318, 322, 348, 362)

PYRROLE = Profile(
    description="Pyrrole halogenase",
    profile_name="pyrrole_FDH",
    filename=get_full_path(str(Path(__file__).parents[1]), "data", "pyrrole_FDH.hmm"),
    cutoffs=(400,),
    motifs=(
        MotifDetails(
            name="mono_di", positions=MODIFICATION_COUNT_POSITIONS,
            residues="DRSVFW", decorations="mono_di", substrate="pyrrole",
        ),
        MotifDetails(
            name="unconv_mono_di", positions=MODIFICATION_COUNT_POSITIONS, residues="YRRNFN",
            decorations="unconv_mono_di", substrate="pyrrole",
        ),
        MotifDetails(
            name="tetra_mono_di", positions=MODIFICATION_COUNT_POSITIONS, residues="RRYFFA",
            decorations="tetra", substrate="pyrrole",
        ),
    ),
    modification_positions=(5,),
)

VARIANTS = [PYRROLE]

SPECIFIC_PROFILES = [variant.profile for variant in VARIANTS]


def get_matching_profiles(hit: HalogenaseHmmResult) -> list[Profile]:
    return [variant for variant in VARIANTS if variant.profile_name == hit.query_id]
