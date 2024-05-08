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
class TyrosineLikeProfile(Profile):
    def get_matches_from_hit(self, retrieved_residues: dict[str, str], hit: HalogenaseHmmResult,
                             confidence: float = 1.,
                             ) -> list[Match]:
        assert set(retrieved_residues).issubset(set(self.motifs))
        matches = []
        modifier = 1.
        for cutoff in self.cutoffs:
            if hit.bitscore < cutoff:
                modifier *= 0.5
                continue
            matches_hpg = retrieved_residues.get("Hpg") == self.motifs["Hpg"]
            if matches_hpg:
                assert retrieved_residues["Tyr"] == self.motifs["Tyr"]
                matches.append(self.create_match(confidence * modifier, retrieved_residues["Hpg"], self.motifs["Hpg"]))
                # Hpg will also match Tyr, but Hpg is more specific
                confidence = max(confidence - 0.2, 0.)

            if retrieved_residues.get("Tyr") == self.motifs["Tyr"]:
                matches.append(self.create_match(confidence * modifier, retrieved_residues["Tyr"], self.motifs["Tyr"]))

            if matches:
                break
        return matches


TYROSINE_LIKE_MOTIF = MotifDetails.from_dict("Tyr", {
    58: "G", 74: "F", 89: "Q", 92: "R", 99: "L", 107: "G", 149: "D", 150: "A", 152: "G",
    209: "L", 215: "S", 217: "G", 219: "V", 245: "P", 267: "S", 268: "Y", 282: "G",
    284: "A", 289: "D", 290: "P", 293: "S", 295: "G", 305: "L", 331: "Y", 357: "W",
}, substrates=("Tyr",))

assert TYROSINE_LIKE_MOTIF.residues == "GFQRLGDAGLSGVPSYGADPSGLYW"

# Hpg is a subset of Tyrosine-like, with more specific positions required
HPG_MOTIF = MotifDetails.from_other("Hpg", TYROSINE_LIKE_MOTIF, {
    66: "S", 158: "H", 196: "C", 200: "G", 246: "M", 259: "Q",
}, substrates=("Hpg",))

TYR_HPG = TyrosineLikeProfile(
    description="Tyrosine-like substrate halogenase",
    profile_name="tyrosine-like_hpg_FDH",
    profile_cutoff=300,
    filename=get_full_path(str(Path(__file__).parents[1]), "data", "tyrosine-like_hpg_FDH.hmm"),
    motifs={
        "Tyr": TYROSINE_LIKE_MOTIF,
        "Hpg": HPG_MOTIF,
    },
    modification_positions=[6, 8],
    alternate_cutoffs=[390.],
)

ORSELLINIC = Profile(
    description="Orsellinic acid-like or other phenolic substrate halogenase",
    profile_name="cycline_orsellinic_FDH",
    profile_cutoff=500,
    filename=get_full_path(str(Path(__file__).parents[1]), "data", "cycline_orsellinic_FDH.hmm"),
    motifs={
        "cycline_orsellinic-like": MotifDetails.from_dict("cycline_orsellinic-like", {
            23: "L", 27: "G", 39: "P", 40: "R", 59: "G", 74: "G", 109: "R", 113: "D",
            120: "A", 124: "G", 133: "V", 165: "D", 166: "A", 168: "G", 233: "G", 284: "Y",
            291: "G", 303: "F", 305: "D", 306: "P", 309: "S", 311: "G",
        }, substrates=("cycline_orsellinic-like",)),
    },
    modification_positions=[6, 8],
)

VARIANTS = [TYR_HPG, ORSELLINIC]

SPECIFIC_PROFILES = [variant.profile for variant in VARIANTS]


def update_match(retrieved_residues: dict[str, str],
                 halogenase: FlavinDependentHalogenase,
                 hit: HalogenaseHmmResult,
                 ) -> None:
    """ Looks whether there are hmm hits that meet the requirement for the categorization
        as Tyr, Hpg, or cycline/orsellinic-like halogenase

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
    for variant in VARIANTS:
        if variant.profile_name == hit.query_id:
            matches = variant.get_matches_from_hit(retrieved_residues, hit)
            halogenase.add_potential_matches(matches)


def get_consensus_signature(cds: CDSFeature, hit: HalogenaseHmmResult,
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
    for variant in VARIANTS:
        if hit.query_id == variant.profile_name:
            residues = substrate_analysis.retrieve_fdh_signature_residues(
                cds.translation, hit, variant.motif_positions, enzyme_substrates=variant.motif_names,
            )
            return {hit.query_id: residues}

    raise ValueError(f"unhandled profile identifier: {hit.query_id}")
