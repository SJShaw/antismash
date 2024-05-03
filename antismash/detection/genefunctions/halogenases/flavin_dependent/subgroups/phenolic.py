# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from pathlib import Path
from typing import Union

from antismash.common.secmet import CDSFeature
from antismash.common.path import get_full_path
from antismash.common.signature import HmmSignature
from antismash.detection.genefunctions.halogenases.halogenases import (
    Match,
    HalogenaseHmmResult,
    FlavinDependentHalogenases)
from antismash.detection.genefunctions.halogenases.flavin_dependent import substrate_analysis

SPECIFIC_PROFILES = [HmmSignature("tyrosine-like_hpg_FDH",
                                  "Tyrosine-like or Hpg substrate halogenase",
                                  300, get_full_path(str(Path(__file__).parents[1]),
                                                     "data", "tyrosine-like_hpg_FDH.hmm")),
                     HmmSignature("cycline_orsellinic_FDH",
                                  "Orsellinic acid-like or other phenolic substrate halogenase",
                                  500, get_full_path(str(Path(__file__).parents[1]),
                                                     "data", "cycline_orsellinic_FDH.hmm"))]

TYROSINE_LIKE_SIGNATURE = [58, 74, 89, 92, 99, 107, 149, 150,
                            152, 209, 215, 217, 219, 245, 267,
                            268, 282, 284, 289, 290, 293, 295, 305, 331, 357]
HPG_SIGNATURE =  [66, 158, 196, 200, 246, 259]

OTHER_PHENOLIC_SIGNATURE = [23, 27, 39, 40, 59, 74, 109, 113, 120, 124, 133, 165,
                            166, 168, 233, 284, 291, 303, 305, 306, 309, 311]

TYR_HPG_SIGNATURE_RESIDUES = {"Tyr": "GFQRLGDAGLSGVPSYGADPSGLYW",
                              "Hpg": "SHCGMQ"}

OTHER_PHENOLIC_SIGNATURE_RESIDUES = "LGPRGGRDAGVDAGGYGFDPSG"

MODIFICATION_POSITIONS = [6, 8]

def find_tyr_hpg_matches(retrieved_residues: dict[str, str],
                         hit: HalogenaseHmmResult,
                         *,
                         confidence: float = 1.,
                         ) -> list[Match]:
    cutoffs=[SPECIFIC_PROFILES[0].cutoff, 390]
    matches = []
    modifier = 1.
    cutoffs.sort(reverse=True)
    substrate_counter = 0
    for subs, sig_res in retrieved_residues.items():
        if sig_res == TYR_HPG_SIGNATURE_RESIDUES[subs]:
            substrate_counter += 1

    for cutoff in cutoffs:
        # matches the residues for Tyrosine and Hpg as well
        if hit.bitscore < cutoff:
            modifier = .5
            continue
        if substrate_counter == 2:
            matches.append(Match(hit.query_id, "flavin", "FDH", confidence * modifier,
                                 retrieved_residues["Hpg"], target_positions=MODIFICATION_POSITIONS, substrates="Hpg"))
            matches.append(Match(hit.query_id, "flavin", "FDH", (confidence * modifier) - 0.2,
                                 retrieved_residues["Tyr"], target_positions=MODIFICATION_POSITIONS, substrates="Tyr"))
            return matches

        if retrieved_residues["Tyr"] == TYR_HPG_SIGNATURE_RESIDUES["Tyr"]:
            matches.append(Match(hit.query_id, "flavin", "FDH",
                           confidence * modifier, retrieved_residues,
                           target_positions=MODIFICATION_POSITIONS, substrates="Tyr"))
            return matches
    return matches


def find_orsellinic_matches(retrieved_residues: str, hit: HalogenaseHmmResult,
                     *, confidence: float = 1.) -> list[Match]:
    """ Looks whether there are hmm hits that meet the requirement for the categorization

        Arguments:
            retrieved_residues: residues of the protein sequence
                                in the place of the signature residues
            hit: details of the hit (e.g. bitscore, name of the profile, etc.)
            confidence: reliability of the categorization

        Returns:
            if the hit is one of the phenolic-specific pHMMs,
            then it adds the match, without returning anything,
            otherwise, it returns False
    """
    modifier = 1.
    if retrieved_residues != OTHER_PHENOLIC_SIGNATURE_RESIDUES or hit.bitscore < SPECIFIC_PROFILES[1].cutoff:
        return []
    return [
        Match(hit.query_id, "flavin", "FDH",
              confidence * modifier, retrieved_residues,
              target_positions=MODIFICATION_POSITIONS,
              substrates="cycline_orsellinic-like"),
    ]


def update_match(name: str, retrieved_residues: Union[dict[str, str], str],
                 halogenase: FlavinDependentHalogenases,
                 hit: HalogenaseHmmResult) -> None:
    """ Looks whether there are hmm hits that meet the requirement for the categorization
        as Tyr, Hpg, or cycline/orsellinic-like halogenase

        Arguments:
            name: name of the substrate-specific pHMM
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
    matches = []
    if name == "tyrosine-like_hpg_FDH":
        assert isinstance(retrieved_residues, dict)
        matches = find_tyr_hpg_matches(retrieved_residues, hit)
    elif name == "cycline_orsellinic_FDH":
        assert isinstance(retrieved_residues, str)
        matches = find_orsellinic_matches(retrieved_residues, hit)
    else:
        raise ValueError(f"unknown profile name: {name}")

    for match in matches:
        halogenase.add_potential_matches(match)


def get_tyr_hpg_residues(translation: str, hmm_result: HalogenaseHmmResult) -> dict[str, str]:
    """ Get signature residues for an enzyme from each pHMM

        Arguments:
            sequence: protein sequence
            hmm_result: instance of HmmResult class,
                        which contains information about the hit in a pHMM
            signatures: list of the positions that defines the signature residues in a pHMM

        Returns:
            signature residues which were retrieved from a certain pHMM
    """
    signature_residues: dict[str, str] = {}
    substrates_signatures = dict(zip(["Tyr", "Hpg"], [TYROSINE_LIKE_SIGNATURE, HPG_SIGNATURE]))
    for substrate, signature in substrates_signatures.items():
        residues = substrate_analysis.search_residues(translation, signature, hmm_result)
        if residues:
            signature_residues[substrate] = residues
    return signature_residues


def get_consensus_signature(cds: CDSFeature, hit: HalogenaseHmmResult
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
                                                                      [TYROSINE_LIKE_SIGNATURE,
                                                                       HPG_SIGNATURE],
                                                                      enzyme_substrates=["Tyr",
                                                                                         "Hpg"])
        return {"tyrosine-like_hpg_FDH": residues}

    if hit.query_id == "cycline_orsellinic_FDH":
        residues = substrate_analysis.retrieve_fdh_signature_residues(cds.translation, hit,
                                                                      OTHER_PHENOLIC_SIGNATURE)
    return residues
