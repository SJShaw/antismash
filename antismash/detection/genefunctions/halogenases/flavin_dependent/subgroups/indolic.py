# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from dataclasses import dataclass
from pathlib import Path

from antismash.common.secmet import CDSFeature
from antismash.common.path import get_full_path
from antismash.detection.genefunctions.halogenases.data_structures import (
    FlavinDependentHalogenase,
    HalogenaseHmmResult,
    Match,
    MotifDetails,
    Profile,
)
from antismash.detection.genefunctions.halogenases.flavin_dependent import substrate_analysis


@dataclass(frozen=True, kw_only=True)
class TryptophanProfile(Profile):
    def get_matches_from_hit(self, retrieved_residues: dict[str, str], hit: HalogenaseHmmResult,
                             confidence: float = 1., check_residues: bool = True) -> list[Match]:
        # if there are no motifs, then any found will break the equality check later
        if not self.motif_residues:
            retrieved_residues = {}
        modifier = 1.
        matches = []
        for cutoff in self.cutoffs:
            if hit.bitscore < cutoff:
                modifier = .5
                continue

            match = Match(hit.query_id, "flavin", "FDH",
                          confidence * modifier, consensus_residues="",
                          target_positions=self.modification_positions,
                          number_of_decorations="mono",
                          substrates=("tryptophan",))
            if not check_residues:
                matches.append(match)
                break

            for name, residues in self.motif_residues.items():
                if retrieved_residues[name] != residues:
                    continue
                match.consensus_residues = residues
                matches.append(match)
            break  # lower cutoffs are irrelevant if a higher is satisfied
        return matches


TRP_5_MOTIF = MotifDetails(
    name="trp_5_FDH",
    positions=(33, 35, 41, 75, 77, 78, 80, 92, 101, 128, 165, 186, 187, 195, 302,
               310, 342, 350, 400, 446, 450, 452, 482),
    residues="VSILIREPGLPRGVPRAVLPGEA",
)

TRP_6_MOTIF = MotifDetails(
    name="trp_6_7_FDH",
    positions=(19, 37, 45, 73, 75, 90, 129, 130, 142, 157, 181, 192, 194, 219, 221,
               225, 227, 237, 287, 306, 337, 339, 350, 353, 356, 411, 462, 505),
    residues="TEGCAGFDAYHDRFGNADYGLSIIAKIL",
)

TRP_5 = TryptophanProfile(
    description="Tryptophan-5 halogenase",
    profile_name="trp_5_FDH",
    profile_cutoff=350,
    alternate_cutoffs=[850],
    filename=get_full_path(str(Path(__file__).parents[1]), "data", "trp_5_FDH.hmm"),
    motifs={"trp_5_FDH": TRP_5_MOTIF},
    modification_positions=[5],
)

TRP_6 = TryptophanProfile(
    description="Tryptophan-6 or -7 halogenase",
    profile_name="trp_6_7_FDH",
    profile_cutoff=770,
    filename=get_full_path(str(Path(__file__).parents[1]), "data", "trp_6_7_FDH.hmm"),
    motifs={"trp_6_7_FDH": TRP_6_MOTIF},
    modification_positions=[6],
)

TRP_7 = TryptophanProfile(
    description="Tryptophan-6 or -7 halogenase",
    profile_name="trp_6_7_FDH",
    profile_cutoff=770,
    filename=get_full_path(str(Path(__file__).parents[1]), "data", "trp_6_7_FDH.hmm"),
    motifs={},
    modification_positions=[7],
)

VARIANTS = [TRP_5, TRP_6, TRP_7]

SPECIFIC_PROFILES = [variant.profile for variant in VARIANTS]


def update_match(retrieved_residues: dict[str, str], halogenase: FlavinDependentHalogenase,
                 hit: HalogenaseHmmResult) -> None:
    """ Looks whether there are hmm hits that meet the requirement for the categorization
        as Trp-5, Trp-6, or Trp-7 halogenase

        Arguments:
            retrieved_residues: residues of the protein sequence
                                in the place of the signature residues
            halogenase: initiated flavin-dependent halogenase
            hit: details of the hit (e.g. bitscore, name of the profile, etc.)

        Returns:
            if the categorization as Trp-5/6/7-halogenase could be done, it instanciates the match
            including the profile name, cofactor, family, position,
            confidence, signature and substrate,
            otherwise, it doesn't return anything and doesn't instanciate anything
    """
    matches = []
    if hit.hit_id == "trp_5_FDH":
        matches = TRP_5.get_matches_from_hit(retrieved_residues, hit)
        assert len(matches) == 1, matches
    elif hit.hit_id == "trp_6_7_FDH":
        matches = TRP_6.get_matches_from_hit(retrieved_residues, hit)
        if not matches:
            matches = TRP_7.get_matches_from_hit(retrieved_residues, hit, check_residues=False)

    halogenase.add_potential_matches(matches)


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
            the signature residues
    """
    residues = {}
    for variant in VARIANTS:
        if variant.motifs and hit.query_id == variant.profile_name:
            residues[variant.profile_name] = substrate_analysis.retrieve_fdh_signature_residues(cds.translation,
                                                                                                hit, variant.motif_positions, variant.motif_names)
    return residues
