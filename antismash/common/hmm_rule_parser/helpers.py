# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from typing import Any

from dataclasses import dataclass
from operator import xor
from typing import Iterator, Self

#class State(enum.StrEnum):
#    NEGATED = enum.auto()
#    NORMAL = enum.auto()
#    BOTH = enum.auto()


#class Presence:
#    def __init__(self, *, negated: bool | State = False, present: bool = False) -> None:
#        if isinstance(negated, State):
#            self.state = negated
#        elif negated:
#            self.state = State.NEGATED
#        elif not negated:
#            self.state = State.NORMAL
#        else:
#            raise TypeError(f"bad type: {type(negated)=}")
#            assert isinstance(negated, type(self))
#            self.state = negated.state
#        self.present = present
#        if self.state == State.BOTH:
#            self.present = True

#    def copy(self, negate: bool = False) -> Self:
#        if negate:
#            return self.flipped()
#        return type(self)(negated=self.state, present=self.present)

#    def __bool__(self) -> bool:
#        return self.state != State.NEGATED and self.present

#    def flipped(self) -> Self:
#        new = type(self)(negated=State.BOTH, present=True)
#        if self.state is State.NORMAL:
#            new = type(self)(negated=State.NEGATED, present=self.present)
#        elif self.state is State.NEGATED:
#            new = type(self)(negated=State.NORMAL, present=self.present)
#        return new

#    def update(self, other: Self | None, flip: bool = False) -> Self:
#        if other is None:
#            return self.flipped() if flip else self
#        if other.state != self.state:
#            self.state = State.BOTH
#        if flip:
#            return self.flipped()
#        return self

#    def __repr__(self) -> str:
#        return f"Negation({self.state}, {'present' if self.present else 'absent'})"


@dataclass(kw_only=True, frozen=True)
class Presence:
    present: bool = True
    negated: bool = False
    multivalue: bool = False

    def blocks(self) -> bool:
        return self.present and self.negated

    def __bool__(self) -> bool:
        return xor(self.present, self.negated)


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
        return type(self)(self.states.copy())

    def items(self) -> Iterator[tuple[str, Presence]]:
        for key, val in self.states.items():
            yield key, val

    def values(self) -> Iterator[Presence]:
        for val in self.states.values():
            yield val

    def merge(self, other: Self) -> Self:
        result = self.states.copy()
        for name, presence in other.items():
            if name in result and result[name].negated != presence.negated:
                presence = Presence(present=presence.present, negated=True, multivalue=True)
            result[name] = presence
        return type(self)(result)

    def get_positives(self) -> set[str]:
        return {key for key, val in self.items() if val.present and not val.blocks()}

    def get_negatives(self) -> set[str]:
        return {key for key, val in self.items() if val.blocks()}

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
    def __init__(self, anchor: str,
                 *,
                 child_states: list["ConditionSatisfaction"] = None,
                 matches_in_anchor: PresenceAbsence = None,
                 matches_in_neighbours: dict[str, PresenceAbsence] | None = None,
                 negated: bool = False,
                 ) -> None:
        self.child_states = child_states or []
        self.anchor = anchor
        self.negated = negated
        self._matches_in_anchor = matches_in_anchor if matches_in_anchor else PresenceAbsence()
        self._matches_in_neighbours = matches_in_neighbours or {}

    @property
    def matches_in_anchor(self) -> PresenceAbsence:
        if not self.child_states:
            return self._matches_in_anchor
        result = self.child_states[0].matches_in_anchor
        for child in self.child_states[1:]:
            result = result.merge(child.matches_in_anchor)
        return result

    @property
    def matches_in_neighbours(self) -> dict[str, PresenceAbsence]:
        if not self.child_states:
            return self._matches_in_neighbours
        result = {key: val.copy() for key, val in self.child_states[0].matches_in_neighbours.items()}
        for child in self.child_states[1:]:
            for name, presence_absence in child.matches_in_neighbours.items():
                if name in result:
                    result[name] = result[name].merge(presence_absence)
                else:
                    result[name] = presence_absence
        return result

    @property
    def satisfied(self) -> bool:
        return bool(self.positive_matches)

    def blocks(self) -> bool:
        if any(presence.blocks() for presence in self.matches_in_anchor.values()):
            return True
        for neighbour in self.matches_in_neighbours.values():
            if any(presence.blocks() for presence in neighbour.values()):
                return True
        return False

    @property
    def positive_matches_in_anchor(self) -> set[str]:
        print("pos in anc", [(n, bool(m)) for n, m in self.matches_in_anchor.items()])
        return self.matches_in_anchor.get_positives()

    @property
    def positive_matches_in_neighbours(self) -> set[str]:
        matches = set()
        for val in self.ancillary_hits.values():
            matches.update(val)
        return matches

    @property
    def positive_matches(self) -> set[str]:
        print("pos matches input", self.positive_matches_in_anchor, self.positive_matches_in_neighbours)
        print("pos matches output", self.positive_matches_in_anchor.union(self.positive_matches_in_neighbours))
        return self.positive_matches_in_anchor.union(self.positive_matches_in_neighbours)

    def has_any_negative_matches(self) -> bool:
        if self.matches_in_anchor.get_negatives():
            return True
        if any(val.get_negatives() for val in self.matches_in_neighbours.values()):
            return True
        return False

    def copy(self, negate: bool = False) -> Self:
        return type(self)(
            anchor=self.anchor,
            matches_in_anchor=self.matches_in_anchor.copy(negate=negate),
            matches_in_neighbours={key: val.copy(negate=negate) for key, val in self.matches_in_neighbours.items()},
            negated=xor(self.negated, negate),
        )

    @property
    def matches(self) -> set[str]:
        return self.matches_in_anchor.get_positives()

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
        return self.satisfied

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


class AndConditionSatisfied(ConditionSatisfaction):
    @property
    def satisfied(self) -> bool:
        print("and cond sat:", all(child.satisfied for child in self.child_states) , "and" ,any(child.matches_in_anchor.get_positives() for child in self.child_states), "xor", self.negated)
        for i, child in enumerate(self.child_states):
            print("child", i, ":", child, bool(child.matches_in_anchor))
        return xor(
            all(child.satisfied for child in self.child_states) and any(child.matches_in_anchor for child in self.child_states),
            self.negated
        )


class OrConditionSatisfied(ConditionSatisfaction):
    @property
    def satisfied(self) -> bool:
        return xor(
            (
                any(child.satisfied and child.matches_in_anchor for child in self.child_states)
                and not any(child.blocks() for child in self.child_states)
            ),
            self.negated
        )

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

