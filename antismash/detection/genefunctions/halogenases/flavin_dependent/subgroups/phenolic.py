# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from dataclasses import dataclass, field
from functools import cached_property
import logging
from pathlib import Path
from typing import Union

from antismash.common.secmet import CDSFeature
from antismash.common.path import get_full_path
from antismash.common.signature import HmmSignature
from antismash.detection.genefunctions.halogenases.halogenases import (
    Match,
    HalogenaseHmmResult,
    FlavinDependentHalogenase,
)
from antismash.detection.genefunctions.halogenases.flavin_dependent import substrate_analysis


SPECIFIC_PROFILES = [
    HmmSignature(
        "tyrosine-like_hpg_FDH",
        "Tyrosine-like or Hpg substrate halogenase",
        300,
        get_full_path(str(Path(__file__).parents[1]), "data", "tyrosine-like_hpg_FDH.hmm"),
    ),
    HmmSignature(
        "cycline_orsellinic_FDH",
        "Orsellinic acid-like or other phenolic substrate halogenase",
        500,
        get_full_path(str(Path(__file__).parents[1]), "data", "cycline_orsellinic_FDH.hmm")
    ),
]

TYROSINE_LIKE = {
    58: "G", 74: "F", 89: "Q", 92: "R", 99: "L", 107: "G", 149: "D", 150: "A", 152: "G",
    209: "L", 215: "S", 217: "G", 219: "V", 245: "P", 267: "S", 268: "Y", 282: "G",
    284: "A", 289: "D", 290: "P", 293: "S", 295: "G", 305: "L", 331: "Y", 357: "W",
}
TYROSINE_LIKE_POSITIONS = sorted(TYROSINE_LIKE)

# Hpg is a subset of Tyrosine-like, with more specific positions required
HPG = {
    66: "S", 158: "H", 196: "C", 200: "G", 246: "M", 259: "Q",
}
HPG.update(TYROSINE_LIKE)
HPG_POSITIONS = sorted(HPG)

TYR_HPG_RESIDUES = {
    "Tyr": "".join(amino for _, amino in sorted(TYROSINE_LIKE.items())),
    "Hpg": "".join(amino for _, amino in sorted(HPG.items())),
}

OTHER_PHENOLIC = {
    23: "L", 27: "G", 39: "P", 40: "R", 59: "G", 74: "G", 109: "R", 113: "D",
    120: "A", 124: "G", 133: "V", 165: "D", 166: "A", 168: "G", 233: "G", 284: "Y",
    291: "G", 303: "F", 305: "D", 306: "P", 309: "S", 311: "G",
}
OTHER_PHENOLIC_RESIDUES = "".join(amino for _, amino in sorted(OTHER_PHENOLIC.items()))


@dataclass(frozen=True, kw_only=True)
class Variant:
    description: str
    profile_name: str
    profile_cutoff: int
    filename: str

    signature_positions: list[int]
    expected_residues: str

    modification_positions: list[int]

    default_penalty: float = 0.8  # additive, not multiplicative
    alternate_cutoffs: list[tuple[float, float]] = field(default_factory=list)

    def __post_init__(self) -> None:
        assert len(self.signature_positions) == len(self.expected_residues)

    @cached_property
    def profile(self) -> HmmSignature:
        return HmmSignature(self.profile_name, self.description, self.profile_cutoff, self.filename)

    def get_matches(self, retrieved_residues: str, halogenase: FlavinDependentHalogenase, hit: HalogenaseHmmResult) -> list[Match]:
        pass


TYR = Variant(
    description="Tyrosine-like substrate halogenase",
    profile_name="tyrosine-like_hpg_FDH",
    profile_cutoff=300,
    filename=get_full_path(str(Path(__file__).parents[1]), "data", "tyrosine-like_hpg_FDH.hmm"),
    signature_positions=TYROSINE_LIKE_POSITIONS,
    modification_positions=[6, 8],
    expected_residues=TYR_HPG_RESIDUES["Tyr"],
    alternate_cutoffs=[(0., 0.5)],
)

HPG = Variant(
    description="Hpg substrate halogenase",
    profile_name="tyrosine-like_hpg_FDH",
    profile_cutoff=300,
    filename=get_full_path(str(Path(__file__).parents[1]), "data", "tyrosine-like_hpg_FDH.hmm"),
    signature_positions=HPG_POSITIONS,
    modification_positions=[6, 8],
    expected_residues=TYR_HPG_RESIDUES["Hpg"],
    default_penalty=0.2,
    alternate_cutoffs=[(0., 0.5)],
)

ORSELLIC = Variant(
    description="Orsellinic acid-like or other phenolic substrate halogenase",
    profile_name="cycline_orsellinic_FDH",
    profile_cutoff=500,
    filename=get_full_path(str(Path(__file__).parents[1]), "data", "cycline_orsellinic_FDH.hmm"),
    signature_positions=OTHER_PHENOLIC,
    modification_positions=[6, 8],
    expected_residues=OTHER_PHENOLIC_RESIDUES,
)


def search_for_match(
    retrieved_residues: dict[str, str],
    halogenase: FlavinDependentHalogenase,
    hit: HalogenaseHmmResult, positions: list[int],
    cutoffs: list[int], *,
    expected_residues: dict[str, str],
    confidence: float = 1.,
) -> bool:
    """ Looks whether there are hmm hits that meet the requirement for the categorization

        Arguments:
            retrieved_residues: residues of the protein sequence
                                in the place of the signature residues
            halogenase: initiated flavin-dependent halogenase
            hit: details of the hit (e.g. bitscore, name of the profile, etc.)
            positions: the decoration positions
            cutoffs: threshold(s) for the pHMM
            expected_residues: expected, substrate-specific signature residues
            confidence: reliability of the categorization

        Returns:
            True if the hit is one of the phenolic-specific pHMMs, otherwise False
    """

    # check for halogenases with Tyr or Hpg substrates
    modifier = 1.
    if (isinstance(expected_residues, dict) and isinstance(cutoffs, list)
        and isinstance(retrieved_residues, dict)):
        cutoffs.sort(reverse=True)
        substrate_counter = 0
        for subs, sig_res in retrieved_residues.items():
            if sig_res == expected_residues[subs]:
                substrate_counter += 1

        for cutoff in cutoffs:
            # matches the residues for Tyrosine and Hpg as well
            if hit.bitscore < cutoff:
                modifier = .5
                continue
            if substrate_counter == 2:
                halogenase.add_potential_match(Match(hit.query_id, "flavin", "FDH",
                                                    confidence * modifier,
                                                    retrieved_residues["Hpg"],
                                                    target_positions=positions, substrates=["Hpg"]))
                halogenase.add_potential_match(Match(hit.query_id, "flavin", "FDH",
                                                       (confidence * modifier)-0.2,
                                                       retrieved_residues["Tyr"],
                                                       target_positions=positions, substrates=["Tyr"]))
                return True

            if retrieved_residues["Tyr"] == expected_residues["Tyr"]:
                halogenase.add_potential_match(Match(hit.query_id, "flavin", "FDH",
                                                    confidence * modifier, retrieved_residues,
                                                    target_positions=positions, substrates=["Tyr"]))
                halogenase.add_potential_match(Match(hit.query_id, "flavin", "FDH",
                                                    confidence * modifier * .5,
                                                    retrieved_residues["Hpg"],
                                                    target_positions=positions, substrates=["Hpg"]))
                return True

            # TODO: can we never hit Hpg by itself?
            if retrieved_residues["Hpg"] == expected_residues["Hpg"]:
                logging.critical("this isn't covered")
                raise NotImplementedError()

        return False
    if isinstance(cutoffs, int):
        if retrieved_residues != expected_residues or hit.bitscore < cutoffs:
            return False
        halogenase.add_potential_match(Match(hit.query_id, "flavin", "FDH",
                                                    confidence * modifier, retrieved_residues,
                                                    target_positions=positions,
                                                    substrates=["cycline_orsellinic-like"]))
        return True
    return False


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

    if hit.hit_id == "tyrosine-like_hpg_FDH":
        search_for_match(retrieved_residues, halogenase, hit, [6, 8],
                         cutoffs=[SPECIFIC_PROFILES[0].cutoff, 390],
                         expected_residues=TYR_HPG_RESIDUES)
    elif hit.hit_id == "cycline_orsellinic_FDH":
        search_for_match(retrieved_residues, halogenase, hit, [6, 8],
                         cutoffs=SPECIFIC_PROFILES[1].cutoff,
                         expected_residues=OTHER_PHENOLIC_RESIDUES)


def get_consensus_signature(cds: CDSFeature, hit: HalogenaseHmmResult,
                            ) -> Union[dict, dict[str, str]]:
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

    residues = {}
    if hit.query_id == "tyrosine-like_hpg_FDH":
        residues = substrate_analysis.retrieve_fdh_signature_residues(cds.translation, hit,
                                                                      [TYROSINE_LIKE_POSITIONS,
                                                                       HPG_POSITIONS],
                                                                      enzyme_substrates=["Tyr",
                                                                                         "Hpg"])
        return {"tyrosine-like_hpg_FDH": residues}

    if hit.query_id == "cycline_orsellinic_FDH":
        residues = substrate_analysis.retrieve_fdh_signature_residues(cds.translation, hit,
                                                                      OTHER_PHENOLIC)
    return residues
