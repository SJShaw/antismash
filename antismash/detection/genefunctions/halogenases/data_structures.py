# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from dataclasses import dataclass, field
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
    substrate: Optional[str] = None
    target_positions: Optional[tuple[int, ...]] = None
    number_of_decorations: str = ""

    def to_json(self) -> dict[str, Any]:
        data = dict(vars(self))
        data["target_positions"] = list(data["target_positions"])
        return data

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "Match":
        # JSON doesn't have a tuple type, so convert those first
        data["target_positions"] = tuple(data["target_positions"])
        return cls(**data)


@dataclass
class FlavinDependentHalogenase:
    cds_name: str
    confidence: float = 0
    consensus_residues: Optional[Union[str, dict[str, str]]] = None
    substrates: Optional[tuple[str, ...]] = None
    target_positions: Optional[tuple[int, ...]] = None
    number_of_decorations: str = ""
    potential_matches: list[Match] = field(default_factory=list)

    cofactor: ClassVar[str] = "flavin"
    family: ClassVar[str] = "flavin-dependent"

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
        profile: path to the pHMM file
        internal_hits: any hits contained by this hit
    """
    def __init__(self, hit_id: str, bitscore: float, query_id: str,
                 profile: str, start: int = 0, end: int = 0, evalue: float = 0.0,
                 internal_hits: Iterable[HMMResult] = None) -> None:
        super().__init__(hit_id, start, end, evalue, bitscore, internal_hits=internal_hits)
        self.query_id = query_id
        self.profile = profile


@dataclass(frozen=True, kw_only=True)
class MotifDetails:
    """ A class for holding details about a specific motif.

        Attributes:
            name: the name of the motif
            positions: the positions of the motif's residues within a reference sequence
            residues: the residues at each of the provided positions
            substrates: the substrates to which this motif applies, if any
            decorations: a description of the decoration type or count for the motif
    """
    name: str
    positions: tuple[int, ...] = field(repr=False)
    residues: str
    substrate: str = ""
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
    def from_dict(cls, name: str, data: dict[int, str], **kwargs: Any) -> Self:
        """ Constructs an instance from a dictionary of position->residue, along
            with normal keyword arguments.

            Arguments:
                name: the name for the motif
                data: a dictionary of position and residue pairs
                **: any keyword arguments normally supplied when creating a new instance

            Returns:
                the new instance
        """
        positions, residues = zip(*sorted(data.items()))
        return cls(name=name, positions=positions, residues="".join(residues), **kwargs)

    @classmethod
    def from_other(cls, name: str, other: Self, additions: dict[int, str], **kwargs: Any) -> Self:
        """ Constructs a new instance from an existing instance, using all of the existing
            instances positions, with extra additions and replacing other values if specified.

            Arguments:
                name: the name for the new motif
                additions: a dictionary of position and residue pairs to insert into the existing set
                **: any keyword arguments with which to use for the new motif instead of the old values

            Returns:
                the new instance
        """
        base = dict(other)
        base.update(additions)
        return cls.from_dict(name, base, **kwargs)


@dataclass(frozen=True, kw_only=True)
class Profile:
    """ A class for holding details about a particular halogenase type.

        Attributes:
            description: a description of the halogenase type
            profile_name: the name of the matching pHMM
            filename: the name of the pHMM file itself
            cutoffs: a list of bitscore cutoffs ordered by increasing confidence
            motifs: details for each motif related to the halogenase type
            modification_positions: the positions at which the halogenase will modify a target
    """
    description: str
    profile_name: str
    filename: str = field(repr=False)
    cutoffs: tuple[int, ...]

    motifs: tuple[MotifDetails, ...]
    modification_positions: tuple[int, ...]

    # some pre-cached values, since functools.cached_property ruins some documentation
    _motif_mapping: dict[str, MotifDetails] = field(repr=False, default_factory=dict)
    _hmm_profile: Optional[HmmSignature] = field(repr=False, default=None)

    def __post_init__(self) -> None:
        # some workarounds for pre-caching derived properties within a frozen dataclass
        if sorted(self.cutoffs, reverse=True) != self.cutoffs:
            object.__setattr__(self, "cutoffs", sorted(self.cutoffs, reverse=True))

        object.__setattr__(self, "_hmm_profile",
                           HmmSignature(self.profile_name, self.description, self.cutoffs[-1], self.filename)
                           )

        self._motif_mapping.update({motif.name: motif for motif in self.motifs})

        # value checks
        if len(self._motif_mapping) != len(self.motifs):
            raise ValueError("provided motifs are not uniquely named")
        if not self.cutoffs:
            raise ValueError("at least one cutoff is required for the HMM profile")

    @property
    def profile(self) -> HmmSignature:
        """ The HMM profile for this halogenase type """
        assert self._hmm_profile is not None
        return self._hmm_profile

    def get_matches_from_hit(self, retrieved_residues: dict[str, str], hit: HalogenaseHmmResult,
                             confidence: float = 1.,
                             ) -> list[Match]:
        """ Searches for good motif matches for the halogenase type.

            Arguments:
                retrieved_residues: a dictionary mapping motif name to residues extracted for that motif
                hit: the HMMer hit from the halogenase type's profile
                confidence: a default confidence to use for all matches found

            Returns:
                a list of matches found
        """
        matches = []

        if hit.bitscore < self.cutoffs[-1]:
            return []

        for motif in self.motifs:
            if retrieved_residues.get(motif.name) == motif:
                matches.append(self.create_match(confidence, retrieved_residues[motif.name], motif))

        return matches

    def create_match(self, confidence: float, residues: str, motif: MotifDetails) -> Match:
        """ Creates a match instance for the halogenase type.

            Arguments:
                confidence: the confidence of the match
                residues: the residues of the match, as extracted for the motif
                motif: the details of the motif that triggered the match

            Returns:
                the new match
        """
        return Match(self.profile_name, FlavinDependentHalogenase.cofactor, FlavinDependentHalogenase.family,
                     confidence, residues,
                     target_positions=self.modification_positions, substrate=motif.substrate,
                     number_of_decorations=motif.decorations)

    @property
    def motif_names(self) -> tuple[str, ...]:
        """ The names of the motifs within this halogenase type """
        return tuple(motif.name for motif in self.motifs)
