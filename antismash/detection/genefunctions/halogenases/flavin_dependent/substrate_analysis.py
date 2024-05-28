# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import re
from collections import defaultdict
from typing import Iterable, Optional

from antismash.common.secmet import (
    CDSFeature,
    Record,
)

from antismash.common import (
    fasta,
    subprocessing,
    utils,
)
from antismash.common.subprocessing.hmmpfam import get_alignment_against_profile
from antismash.detection.genefunctions.halogenases.data_structures import (
    FlavinDependentHalogenase,
    HalogenaseHmmResult,
    MotifDetails,
)

from antismash.common.signature import HmmSignature
from antismash.common.subprocessing import hmmscan
from antismash.detection.genefunctions.halogenases.flavin_dependent import substrates
from antismash.detection.genefunctions.halogenases.flavin_dependent.subgroups import (
    indolic,
    phenolic,
    pyrrolic
)

SUBGROUPS = [indolic, phenolic, pyrrolic]


def _get_substrate_specific_profiles() -> list[HmmSignature]:
    """ Collects unique substrate-specific pHMM profiles from the substrate-specific submodules"""
    profiles: dict[str, HmmSignature] = {}
    for submodule in SUBGROUPS:
        for profile in submodule.SPECIFIC_PROFILES:
            existing = profiles.get(profile.name)
            if existing and profile.cutoff >= existing.cutoff:
                continue
            profiles[profile.name] = profile
    return list(profiles.values())


def retrieve_fdh_signature_residues(translation: str, hmm_result: HalogenaseHmmResult,
                                    motifs: Iterable[MotifDetails],
                                    ) -> dict[str, str]:
    """ Extracts residues for each of the given motifs from an HMM hit

        Arguments:
            sequence: protein sequence
            hmm_result: instance of HmmResult class,
                        which contains information about the hit in a pHMM
            motifs: the motifs for which to extract signatures

        Returns:
            signature residues which were retrieved from a certain pHMM
    """
    signature_residues: dict[str, str] = {}
    for motif in motifs:
        if not motif.positions:  # then it's always present, just 'empty'
            signature_residues[motif.name] = ""
            continue
        residues = extract_residues(translation, motif.positions, hmm_result)
        if residues:
            signature_residues[motif.name] = residues
    return signature_residues


def extract_residues(sequence: str, positions: Iterable[int],
                     hmm_result: HalogenaseHmmResult,
                     max_evalue: float = 0.1) -> Optional[str]:
    """ Get the signature residues from the pHMM for the searched protein sequence

        Arguments:
            sequence: protein sequence
            positions: list of position numbers in the pHMM
            hmm_result: properties of the halogenase pHMM hit
            max_evalue: maximum e-value

        Returns:
            residues that are present in the given positions
    """
    if not positions:
        raise ValueError("cannot extract residues without positions")
    alignment = get_alignment_against_profile(sequence, hmm_result.profile,
                                              hmm_result.query_id, max_evalue=max_evalue)
    if not alignment:
        return None
    return utils.extract_from_alignment(alignment, positions)


def search_conserved_motif(cds: CDSFeature, motif_positions: list[int],
                           hmm_result: HalogenaseHmmResult,
                           motif_pattern: str) -> str:
    """ Looks for WxWxIP and Fx.Px.Sx.G conserved motifs, characteristic to FDHs

        Arguments:
            cds: gene/CDS and its properties
            motif_positions: positions of the conserved motifs in the pHMM
            hmm_result: details of the hit (e.g. bitscore, name of the profile, etc.)
            motif_pattern: pattern of the motif in regex

        Returns:
            the conserved motifs, if present, otherwise an empty string
    """
    category = ""

    signature_residues = extract_residues(cds.translation, motif_positions, hmm_result)

    if not signature_residues:
        return category

    motif = re.search(motif_pattern, signature_residues)
    if not motif:
        return category

    return motif[0]


def run_halogenase_phmms(cluster_fasta: str, profiles: list[HmmSignature],
                         ) -> dict[str, list[HalogenaseHmmResult]]:
    """ Check if protein sequences hit any pHMM

        Arguments:
            cluster_fasta: string of protein sequences in a fasta format
            profiles: a list of HMM profiles to run

        Returns:
            a dictionary mapping query names to a list of hits for that query, one per profile
    """
    halogenase_hmms_by_id: dict = defaultdict(list)
    files_run = set()  # so that two info sets sharing the same HMM profile aren't reused
    for sig in profiles:
        if sig.hmm_file in files_run:
            continue
        files_run.add(sig.hmm_file)
        run_results = subprocessing.run_hmmsearch(sig.hmm_file, cluster_fasta)
        for runresult in run_results:
            for hsp in runresult.hsps:
                if hsp.bitscore > sig.cutoff:
                    cds_name = hsp.hit_id
                    hit = HalogenaseHmmResult(cds_name, hsp.bitscore, hsp.query_id, sig.path)
                    halogenase_hmms_by_id[cds_name].append(hit)
    return halogenase_hmms_by_id


def categorize_on_substrate_level(cds: CDSFeature, halogenase: FlavinDependentHalogenase,
                                  hmm_results: list[HalogenaseHmmResult],
                                  ) -> Optional[FlavinDependentHalogenase]:
    """ Check if protein could be categorized as a Flavin-dependent enzyme

        Modifies the HalogenasesResults by adding the label of the enzyme family,
        and possibly the signature and position of halogenation determined.
        If it doesn't meet the requirements, it doesn't make changes on the original instance.

        Arguments:
            cds: the CDS feature to categorize
            halogenase: the halogenase instance for the given CDS
            hmm_results: a list of HMM hits for the CDS

        Returns:
            the given halogenase, if it was modified, otherwise None
    """
    if not hmm_results:
        return None

    matched = False

    for hit in hmm_results:
        for subgroup in SUBGROUPS:
            for profile in subgroup.get_matching_profiles(hit):
                residues = retrieve_fdh_signature_residues(cds.translation, hit, profile.motifs)
                matches = profile.get_matches_from_hit(residues, hit)
                print("profile", profile, "\n", "matches", matches)
                halogenase.add_potential_matches(matches)
                if matches:
                    matched = True
                    break
    return halogenase if matched else None


def categorize_on_consensus_level(cds: CDSFeature, specific_hmm_hits: list[HalogenaseHmmResult],
                                  general_hmm_hits: list[HalogenaseHmmResult],
                                  ) -> FlavinDependentHalogenase:
    assert len(general_hmm_hits) == 1, general_hmm_hits
    # determine conventionality
    conserved_motifs: dict[str, str] = {}
    for hit in general_hmm_hits:
        for motif, positions in substrates.GENERAL_FDH_MOTIFS.items():
            conserved_motif = search_conserved_motif(cds, positions, hit, motif)
            if conserved_motif:
                conserved_motifs[motif] = conserved_motif

    enzyme = FlavinDependentHalogenase(cds.get_name(), conventionality_residues=conserved_motifs)

    if specific_hmm_hits:
        return categorize_on_substrate_level(cds, enzyme, specific_hmm_hits) or enzyme

    return enzyme


def fdh_specific_analysis(record: Record) -> list[FlavinDependentHalogenase]:
    """ Categorization of enzyme, categorizes any halogenase in a cds in regions

        Arguments:
            record: the record instance to analyse

        Returns:
            a list of HalogenasesResults instances representing halogenase enzymes,
                 if there is a clear best match for a given enzyme, then the information
                 about the position of halogenation, the confidence of the categorization,
                 and the characteristic residues is provided.
                 If the enzyme can be categorized into several groups, with the same confidence,
                 then the above mentioned informations are not defined,
                 and the information about the catogries is in the potential_enzymes attribute.
    """

    potential_enzymes = []
    enzymes_with_hits = []
    hmmsearch_fasta = fasta.get_fasta_from_features(record.get_cds_features_within_regions())

    hits = hmmscan.run_hmmscan(substrates.ALL_FDH_PROFILES, hmmsearch_fasta)
    import logging; logging.critical("discarding good info")  # TODO
    for query_result in hits:
        if query_result.hits:
            enzymes_with_hits.append(record.get_cds_by_name(query_result.id))
    if not enzymes_with_hits:
        return []

    hit_enzyme_fasta = fasta.get_fasta_from_features(enzymes_with_hits)
    general_hmm_hits = run_halogenase_phmms(hit_enzyme_fasta, substrates.GENERAL_FDH_PROFILES)
    if general_hmm_hits:
        specific_profiles = _get_substrate_specific_profiles()
        specific_hmm_hits = run_halogenase_phmms(hit_enzyme_fasta, specific_profiles)

    for protein in general_hmm_hits:
        cds = record.get_cds_by_name(protein)
        potential_enzymes.append(categorize_on_consensus_level(cds, specific_hmm_hits[protein],
                                                               general_hmm_hits[protein]))

    return potential_enzymes
