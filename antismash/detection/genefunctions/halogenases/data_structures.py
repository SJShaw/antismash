# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from dataclasses import dataclass, field
from functools import cached_property
from typing import Any, ClassVar, Iterable, Iterator, Optional, Self, Union

from antismash.common.hmmscan_refinement import HMMResult
from antismash.common.signature import HmmSignature


@dataclass
class Match:
    """ Match of the enzyme categorized by check_for_fdh,
        with details about which pHMM (profile) was hit, what position
        the halogenation occurs, what is the confidence of the categorization,
        and what are the signature residues of the protein sequence"""
    profile: str
    cofactor: str
    family: str
    confidence: float
    consensus_residues: str
    substrates: Optional[tuple[str, ...]] = None
    target_positions: Optional[list[int]] = None
    number_of_decorations: str = ""

    def to_json(self) -> dict[str, Any]:
        return vars(self)

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "Match":
        # JSON doesn't have a tuple type, so convert those first
        data["substrates"] = tuple(data["substrates"])
        return cls(**data)


@dataclass
class FlavinDependentHalogenase:
    cds_name: str
    confidence: float = 0
    consensus_residues: Optional[Union[str, dict[str, str]]] = None
    substrates: Optional[tuple[str, ...]] = None
    target_positions: Optional[list[int]] = None
    number_of_decorations: str = ""
    potential_matches: list[Match] = field(default_factory=list)

    cofactor: ClassVar[str] = "flavin"
    family: ClassVar[str] = "FDH"

    def add_potential_match(self, match: Match) -> None:
        """ Adds the features of an enzyme group"""
        self.potential_matches.append(match)

    def add_potential_matches(self, matches: Iterable[Match]) -> None:
        self.potential_matches.extend(matches)

    def get_best_matches(self) -> list[Match]:
        """ If an enzyme meets the requirements for several groups,
            it compares the confidences of the categorizations and
            returns the one with the highest confidence or list of matches.
            If there are more groups with the same confidence, it returns the list of those."""
        best_match = []

        if self.potential_matches:
            if len(self.potential_matches) == 1:
                return [self.potential_matches[0]]

            highest_confidence = max(profile.confidence for profile in self.potential_matches)
            for profile in self.potential_matches:
                if abs(profile.confidence - highest_confidence) <= 0.005:
                    best_match.append(profile)

        return best_match

    def finalize_enzyme(self) -> None:
        """ If there is a best match among the matches based on confidence,
            get that one match and define position, confidence and signature
            in the enzyme instance based on that.
            If there is no one best match, it doesn't change anything."""
        best_matches = self.get_best_matches()
        assert isinstance(best_matches, list), best_matches
        if not best_matches:
            return

        if len(best_matches) == 1:
            best_match = best_matches[0]
            assert best_match.cofactor == self.cofactor
            assert best_match.family == self.family
            self.target_positions = best_match.target_positions
            self.consensus_residues = best_match.consensus_residues
            self.confidence = best_match.confidence
            self.number_of_decorations = best_match.number_of_decorations
            self.substrates = best_match.substrates or tuple()

    def to_json(self) -> dict[str, Any]:
        """ Constructs a JSON representation of this instance """
        potential_matches_json = [match.to_json() for match in self.potential_matches]

        return {
            "cds_name": self.cds_name,
            "family": self.family,
            "cofactor": self.cofactor,
            "substrates": self.substrates if self.substrates else None,
            "target_positions": self.target_positions,
            "number_of_decorations": self.number_of_decorations,
            "consensus_residues": self.consensus_residues,
            "confidence": self.confidence,
            "potential_matches": potential_matches_json
        }

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        """ Constructs an instance from a JSON representation """

        assert data.pop("cofactor") == cls.cofactor
        assert data.pop("family") == cls.family

        cds_name = data["cds_name"]
        substrates = tuple(data["substrates"] or [])
        target_positions = data["target_positions"]
        number_of_decorations = data["number_of_decorations"]
        consensus_residues = data["consensus_residues"]
        confidence = data["confidence"]
        potential_matches = [Match.from_json(profile) for profile in data["potential_matches"]]
        enzyme = cls(cds_name, confidence,
                     consensus_residues, substrates, target_positions,
                     number_of_decorations, potential_matches)
        return enzyme

    def __repr__(self) -> str:
        return f"FlavinDependentHalogenase({self.cds_name=}, {self.confidence=}, {self.potential_matches=})"


class HalogenaseHmmResult(HMMResult):
    """ Enzymes identified as a halogenase

        hit_id: name of the matching profile
        start: start position within the query's translation
        end: end position within the query's translation
        evalue: e-value of the hit
        bitscore: bitscore of the hit
        query_id: name of the profile
        enzyme_type: type of halogenase (e.g. Flavin-dependent, SAM-dependent)
        profile: path to the pHMM file
        internal_hits: any hits contained by this hit
    """
    def __init__(self, hit_id: str, bitscore: float, query_id: str, enzyme_type: str,
                 profile: str, start: int = 0, end: int = 0, evalue: float = 0.0,
                 internal_hits: Iterable[HMMResult] = None) -> None:
        super().__init__(hit_id, start, end, evalue, bitscore, internal_hits=internal_hits)
        self.query_id = query_id
        self.enzyme_type = enzyme_type
        self.profile = profile


@dataclass(frozen=True, kw_only=True)
class MotifDetails:
    name: str
    positions: tuple[int, ...]
    residues: str
    substrates: tuple[str, ...] = field(default_factory=tuple)
    decorations: str = ""

    def __post_init__(self) -> None:
        assert not self.positions or isinstance(self.positions[0], int)
        assert len(list(self.positions)) == len(self.residues)
        assert tuple(sorted(self.positions)) == self.positions

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, MotifDetails):
            return other.residues == self.residues
        if isinstance(other, str):
            return other == self.residues
        return False

    def __iter__(self) -> Iterator[tuple[int, str]]:
        return zip(self.positions, self.residues)

    @classmethod
    def from_dict(cls, name: str, data: dict[int, str], substrates: tuple[str, ...] = None, decorations: str = "") -> Self:
        positions, residues = zip(*sorted(data.items()))
        return cls(name=name, positions=positions, residues="".join(residues), substrates=substrates or tuple(), decorations=decorations)

    @classmethod
    def from_other(cls, name: str, other: Self, additions: dict[int, str], **kwargs: Any) -> Self:
        base = dict(other)
        base.update(additions)
        return cls.from_dict(name, base, **kwargs)


@dataclass(frozen=True, kw_only=True)
class Profile:
    description: str
    profile_name: str
    profile_cutoff: int
    filename: str

    motifs: dict[str, MotifDetails]
    modification_positions: list[int]

    default_penalty: float = 0.8  # additive, not multiplicative
    alternate_cutoffs: list[float] = field(default_factory=list)

    @cached_property
    def profile(self) -> HmmSignature:
        return HmmSignature(self.profile_name, self.description, self.profile_cutoff, self.filename)

    def get_matches_from_hit(self, retrieved_residues: dict[str, str], hit: HalogenaseHmmResult, confidence: float = 1.) -> list[Match]:
        matches = []

        if hit.bitscore < self.profile_cutoff:
            return []

        for name, motif in self.motifs.items():
            if retrieved_residues.get(name) == motif:
                matches.append(self.create_match(confidence, retrieved_residues[name], motif))

        return matches

    def create_match(self, confidence: float, residues: str, motif: MotifDetails) -> Match:
        return Match(self.profile_name, "flavin", "FDH", confidence, residues,
                     target_positions=self.modification_positions, substrates=motif.substrates,
                     number_of_decorations=motif.decorations)

    @cached_property
    def motif_names(self) -> list[str]:
        return list(self.motifs)

    @cached_property
    def motif_positions(self) -> tuple[tuple[int, ...], ...]:
        return tuple(motif.positions for motif in self.motifs.values())

    @cached_property
    def motif_residues(self) -> dict[str, str]:
        return {name: motif.residues for name, motif in self.motifs.items()}

    @cached_property
    def cutoffs(self) -> tuple[float, ...]:
        return tuple(sorted([self.profile_cutoff] + self.alternate_cutoffs, reverse=True))
