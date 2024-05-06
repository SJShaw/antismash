# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from pathlib import Path

from antismash.common.secmet import CDSFeature
from antismash.common.path import get_full_path
from antismash.detection.genefunctions.halogenases.data_structures import (
    FlavinDependentHalogenase,
    HalogenaseHmmResult,
    MotifDetails,
    Profile,
)
from antismash.detection.genefunctions.halogenases.flavin_dependent import substrate_analysis

MODIFICATION_COUNT_POSITIONS = (110, 111, 318, 322, 348, 362)

PYRROLE = Profile(
    description="Pyrrole halogenase",
    profile_name="pyrrole_FDH",
    profile_cutoff=400,
    filename=get_full_path(str(Path(__file__).parents[1]), "data", "pyrrole_FDH.hmm"),
    motifs={
        "mono_di": MotifDetails(name="mono_di", positions=MODIFICATION_COUNT_POSITIONS, residues="DRSVFW", decorations="mono_di", substrates=("pyrrole",)),
        "unconv_mono_di": MotifDetails(name="unconv_mono_di", positions=MODIFICATION_COUNT_POSITIONS, residues="YRRNFN", decorations="unconv_mono_di", substrates=("pyrrole",)),
        "tetra": MotifDetails(name="tetra_mono_di", positions=MODIFICATION_COUNT_POSITIONS, residues="RRYFFA", decorations="tetra", substrates=("pyrrole",)),
    },
    modification_positions=[5],
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


def get_consensus_signature(cds: CDSFeature, hit: HalogenaseHmmResult
                            ) -> dict[str, dict[str, str]]:
    """ Retrieves the residues from the substrate-specific,
        pHMMs that are in the positions of the signature residues

        Arguments:
            cds: gene/CDS and its properties
            hit: details of the hit (e.g. bitscore, name of the profile, etc.)

        Returns:
            if the name of the pHMM doesn't match the substrate-specific one's,
            it returns an empty dictionary,
            otherwise, it returns the residues, that are in the same positions as
            the sugnature residues
    """

    signatures = {}
    for variant in VARIANTS:
        if variant.motifs and hit.query_id == variant.profile_name:
            residues = substrate_analysis.retrieve_fdh_signature_residues(cds.translation, hit, variant.motif_positions, variant.motif_names)
            if residues:
                signatures[variant.profile_name] = residues
    return signatures
