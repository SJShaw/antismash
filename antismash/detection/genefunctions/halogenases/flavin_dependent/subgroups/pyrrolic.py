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


def update_match(retrieved_residues: dict[str, str],
                 halogenase: FlavinDependentHalogenase,
                 hit: HalogenaseHmmResult) -> None:
    """ Looks whether there are hmm hits that meet the requirement for the categorization
        as a pyrrole halogenase doing mono/di- or tetra-halogenation

        Arguments:
            retrieved_residues: residues of the protein sequence
                                in the place of the signature residues
            halogenase: initiated flavin-dependent halogenase
            hit: details of the hit (e.g. bitscore, name of the profile, etc.)

        Returns:
            if the categorization as Tyr/Hpg/other-phenolic-halogenase could be done
            it instanciates the match including the profile name, cofactor, family,
            position, confidence, signature and substrate,
            otherwise, it doesn't return anything and doesn't instanciate anything
    """
    present = False
    for variant in VARIANTS:
        if hit.query_id == variant.profile_name:
            present = True
            matches = variant.get_matches_from_hit(retrieved_residues, hit)
            halogenase.add_potential_matches(matches)
    if not present:
        raise ValueError(f"unhandled profile: {hit.hit_id}")


def get_matching_profiles(hit: HalogenaseHmmResult) -> list[Profile]:
    return [variant for variant in VARIANTS if variant.profile_name == hit.query_id]
