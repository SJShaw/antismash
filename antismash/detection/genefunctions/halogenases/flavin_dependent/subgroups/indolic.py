# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from pathlib import Path

from antismash.common.path import get_full_path
from antismash.detection.genefunctions.halogenases.data_structures import (
    HalogenaseHmmResult,
    MotifDetails,
    Profile,
)


TRP_5 = Profile(
    description="Tryptophan-5 halogenase",
    profile_name="trp_5_FDH",
    cutoffs=(850, 350),
    filename=get_full_path(str(Path(__file__).parents[1]), "data", "trp_5_FDH.hmm"),
    motifs=(
        MotifDetails(
            name="trp_5_FDH",
            positions=(33, 35, 41, 75, 77, 78, 80, 92, 101, 128, 165, 186, 187, 195, 302,
                       310, 342, 350, 400, 446, 450, 452, 482),
            residues="VSILIREPGLPRGVPRAVLPGEA",
            substrate="tryptophan",
            decorations="mono",
        ),
    ),
    modification_positions=(5,),
)

TRP_6 = Profile(
    description="Tryptophan-6 halogenase",
    profile_name="trp_6_7_FDH",
    cutoffs=(770,),
    filename=get_full_path(str(Path(__file__).parents[1]), "data", "trp_6_7_FDH.hmm"),
    motifs=(
        MotifDetails(
            name="trp_6_7_FDH",
            positions=(19, 37, 45, 73, 75, 90, 129, 130, 142, 157, 181, 192, 194, 219, 221,
                       225, 227, 237, 287, 306, 337, 339, 350, 353, 356, 411, 462, 505),
            residues="TEGCAGFDAYHDRFGNADYGLSIIAKIL",
            substrate="tryptophan",
            decorations="mono",
        ),
    ),
    modification_positions=(6,),
)

TRP_7 = Profile(
    description="Tryptophan-7 halogenase",
    profile_name="trp_6_7_FDH",
    cutoffs=(770,),
    filename=get_full_path(str(Path(__file__).parents[1]), "data", "trp_6_7_FDH.hmm"),
    motifs=(
        MotifDetails(
            name="trp_6_7_FDH",
            positions=tuple(),
            residues="",
            substrate="tryptophan",
            decorations="mono",
        ),
    ),
    modification_positions=(7,),
)

VARIANTS = [TRP_5, TRP_6, TRP_7]

SPECIFIC_PROFILES = [variant.profile for variant in VARIANTS]


def get_matching_profiles(hit: HalogenaseHmmResult) -> list[Profile]:
    return [variant for variant in VARIANTS if variant.profile_name == hit.query_id]
