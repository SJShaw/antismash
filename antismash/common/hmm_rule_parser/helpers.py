# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from typing import Any

import enum
from collections import defaultdict
from operator import xor
from typing import Iterator, Self, Union

class State(enum.StrEnum):
    NEGATED = enum.auto()
    NORMAL = enum.auto()
    BOTH = enum.auto()


class Presence:
    def __init__(self, *, negated: bool | State = False, present: bool = False) -> None:
        if isinstance(negated, State):
            self.state = negated
        elif negated:
            self.state = State.NEGATED
        elif not negated:
            self.state = State.NORMAL
        else:
            raise TypeError(f"bad type: {type(negated)=}")
            assert isinstance(negated, type(self))
            self.state = negated.state
        self.present = present
        if self.state == State.BOTH:
            self.present = True

    def copy(self, negate: bool = False) -> Self:
        if negate:
            return self.flipped()
        return type(self)(negated=self.state, present=self.present)

    def __bool__(self) -> bool:
        return self.state != State.NEGATED and self.present

    def flipped(self) -> Self:
        new = type(self)(negated=State.BOTH, present=True)
        if self.state is State.NORMAL:
            new = type(self)(negated=State.NEGATED, present=self.present)
        elif self.state is State.NEGATED:
            new = type(self)(negated=State.NORMAL, present=self.present)
        return new

    def update(self, other: Self | None, flip: bool = False) -> Self:
        if other is None:
            return self.flipped() if flip else self
        if other.state != self.state:
            self.state = State.BOTH
        if flip:
            return self.flipped()
        return self

    def __repr__(self) -> str:
        return f"Negation({self.state}, {'present' if self.present else 'absent'})"


class PresenceAbsence:
    def __init__(self, states: dict[str, Presence] | None = None) -> None:
        self.states = states or {}

    def __contains__(self, key: Any) -> bool:
        return key in self.states

    def __iter__(self) -> Any:  # TODO
        for key in self.states:
            yield key

    def __getitem__(self, key: Any) -> Presence:
        val = self.states.get(key)
        if val is None:
            raise KeyError(str(key))
        return val

    def __len__(self) -> int:
        return len(self.states)

    def copy(self, *, negate: bool = False) -> Self:
        return type(self)({key: val.copy(negate=negate) for key, val in self.states.items()})

    def items(self) -> Iterator[tuple[str, Presence]]:
        for key, val in self.states.items():
            yield key, val

    def merge(self, other: Self) -> Self:
        result = self.copy()
        for key, val in other.items():
            if key in result:
                result.states[key].update(val)
            else:
                result.states[key] = val.copy()
        return result

    def get_positives(self) -> set[str]:
        return {key for key, val in self.items() if val}

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        parts = []
        for key, val in self.states.items():
            parts.append(f"  {key}: {val}\n")
        return "".join(parts)


class ConditionSatisfaction:
    """ A container for tracking whether a condition was satisfied along with
        what specific subsections of the condition were matched
    """
    def __init__(self, satisfied: bool, anchor: str,
                 matches_in_anchor: PresenceAbsence = None,
                 *,
                 matches_in_neighbours: dict[str, PresenceAbsence] | None = None,
                 negated: bool = False,
                 ) -> None:
        self.satisfied = satisfied
        self.anchor = anchor
        self.negated = negated
        self.matches_in_anchor = matches_in_anchor if matches_in_anchor else PresenceAbsence()
        self.matches_in_neighbours = matches_in_neighbours or {}

    @property
    def positive_matches_in_anchor(self) -> set[str]:
        return self.matches_in_anchor.get_positives()

    def merge(self, other: Self) -> Self:
        assert self.anchor == other.anchor

        anchor_matches = self.matches_in_anchor.merge(other.matches_in_anchor)
        neighbour_matches = {key: val.merge(other.matches_in_neighbours[key]) for key, val in self.matches_in_neighbours.items()}
        neighbour_matches.update({key: val.copy() for key, val in other.matches_in_neighbours.items() if key not in neighbour_matches})

        return type(self)(
            self.satisfied and other.satisfied, self.anchor,
                matches_in_anchor=anchor_matches,
                matches_in_neighbours=neighbour_matches,
                negated=self.negated,
            )

    def copy(self, negate: bool = False) -> Self:
        return type(self)(
            xor(self.satisfied, negate),
            anchor=self.anchor,
            matches_in_anchor=self.matches_in_anchor.copy(negate=negate),
            matches_in_neighbours={key: val.copy(negate=negate) for key, val in self.matches_in_neighbours.items()},
            negated=xor(self.negated, negate),
        )

    @property
    def matches(self) -> set[str]:
        return set(name for name, val in self.matches_in_anchor.items() if val)

    @property
    def ancillary_hits(self) -> dict[str, set[str]]:
        result = {}
        for neighbour, matches in self.matches_in_neighbours.items():
            positives = matches.get_positives()
            if not positives:
                continue
            result[neighbour] = positives
        return result

    # purely for backwards compatibility
    @property
    def met(self) -> bool:
        return bool(self)

    def __bool__(self) -> bool:
        return self.satisfied and bool(self.positive_matches_in_anchor)

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        parts = [(
            f"\n satisfied={self.satisfied} negated={self.negated}"
        )]
        parts.append("  anchor_matches:")
        for key, val in self.matches_in_anchor.items():
            parts.append(f"    {key}: {val}")
        parts.append("  neighbour_matches:")
        for name, matches in self.matches_in_neighbours.items():
            parts.append(f"    {name}")
            for key, val in matches.items():
                parts.append(f"      {key}: {val}")
        return "\n".join(parts)
            


#    def __init__(self, met: bool, matches: Union[set[str], "ConditionMet"] = None,
#                 ancillary_hits: dict[str, set[str]] = None,
#                 *,
#                 negated: bool = False,
#                 presence_and_negation: dict[str, NegationState] = None,
#                 ) -> None:
#        self.presence_and_negation = {key: (value.flipped() if negated else value) for key, value in (presence_and_negation or {}).items()}
#        assert isinstance(met, bool)
#        self.met = met
#        self._matches: set[str] = set()
#        self.ancillary_hits = ancillary_hits or {}
#        if isinstance(matches, ConditionMet) and presence_and_negation is None:
#            self.presence_and_negation = {key: val.copy() for key, val in matches.presence_and_negation.items()}
##        print(f"  ConditionMet: {self.met=} {self.ancillary_hits=}, presence={self.presence_and_negation}")

#    @property
#    def matches(self) -> dict[str, NegationState]:
#        return {key: val.copy() for key, val in self.presence_and_negation.items()}

#    def merge(self, other: Self, *, negate: bool = False, require_all: bool = True, condition: Any = None) -> Self:
#        if require_all:
#            good = self.met and other.met
#        else:
#            good = self.met or other.met
#        presence = self.presence_and_negation.copy()
#        for name, negation in other.presence_and_negation.items():
#            print(f"   in merge   {other=}, {condition=}")
#            if name in presence:
#                presence[name].update(negation)
#            else:
#                presence[name] = negation
#        res = type(self)(good, set(), self.ancillary_hits, negated=negate, presence_and_negation=presence)
#        print("post-merge", res)
#        return res

#    def as_negated(self) -> Self:
#        return type(self)(not self.met, self, self.ancillary_hits, presence_and_negation={key: val.flipped() for key, val in self.presence_and_negation.items()})

#    def __bool__(self) -> bool:
#        return self.met

#    def __repr__(self) -> str:
#        return str(self)

#    def __str__(self) -> str:
#        return (
#            f"satisfied: {self.met}"
#            f", with hits: {self.matches or 'none'}"
#            f" and with others: {self.ancillary_hits or 'none'}"
#        )

