# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import enum
from collections import defaultdict
from typing import Self, Union

class State(enum.StrEnum):
    NEGATED = enum.auto()
    NORMAL = enum.auto()
    BOTH = enum.auto()

class NegationState:
    def __init__(self, negated: bool | State = False, present: bool = False) -> None:
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

    def copy(self) -> Self:
        return type(self)(self.state, self.present)

    def __bool__(self) -> bool:
        return self.state != State.NEGATED and self.present

    def flipped(self) -> Self:
        if self.state is State.NORMAL:
            return type(self)(State.NEGATED, present=self.present)
#            self.state = State.NEGATED
        elif self.state is State.NEGATED:
            return type(self)(State.NORMAL, present=self.present)
#            self.state = State.NORMAL
        return type(self)(State.BOTH, present=True)

    def update(self, other: Self | None, flip: bool = False) -> Self:
        if other is None:
            return self.flipped() if flip else self
        if other.state != self.state:
            self.state = State.BOTH
        return self.flipped() if flip else self

    def __repr__(self) -> str:
        return f"Negation({self.state}, {'present' if self.present else 'absent'})"


class ConditionMet:
    """ A container for tracking whether a condition was satisfied along with
        what specific subsections of the condition were matched
    """
    def __init__(self, met: bool, matches: Union[set[str], "ConditionMet"] = None,
                 ancillary_hits: dict[str, set[str]] = None,
                 *,
                 negated: bool = False,
                 presence_and_negation: dict[str, NegationState] = None,
                 ) -> None:
        self.presence_and_negation = {key: (value.flipped() if negated else value) for key, value in (presence_and_negation or {}).items()}
        self.negated = negated
        assert isinstance(met, bool)
        self.met = met
        self._matches: set[str] = set()
        self.ancillary_hits = ancillary_hits or {}
        if isinstance(matches, ConditionMet) and presence_and_negation is not None:
            self.presence_and_negation = matches.presence_and_negation.copy()
        print(f"  ConditionMet: {self.met=} {self.matches=}, {self.ancillary_hits=}, presence={self.presence_and_negation}")

    @property
    def matches(self) -> dict[str, NegationState]:
#        if not self.met:  # this condition works for the case of "A or (X and not B) or (B and not Y)" and "A or X and not B",
#            return set()  # but not for nesting like "A or not (X and not B)"
        return {key: val.copy() for key, val in self.presence_and_negation.items()}
        if self.negated:
            return {(name, negated.flipped()) for name, negated in self.presence_and_negation.items()}
        return {(name, negated) for name, negated in self.presence_and_negation.items()}

    def merge(self, other: Self, *, negate: bool = False, require_all: bool = True) -> Self:
        if require_all:
            good = self.met and other.met
        else:
            good = self.met or other.met
        presence = self.presence_and_negation.copy()
        for name, negation in other.presence_and_negation.items():
            if other.negated:
                negation = negation.flipped()
            if name in presence:
                presence[name].update(negation)
            else:
                presence[name] = negation
        return type(self)(good, set(), self.ancillary_hits, negated=negate, presence_and_negation=presence)

    def as_negated(self) -> Self:
        return type(self)(not self.met, self, self.ancillary_hits, negated=self.negated, presence_and_negation=self.presence_and_negation)

    def __bool__(self) -> bool:
        return self.met

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return (
            f"satisfied: {self.met}"
            f", with hits: {self.matches or 'none'}"
            f" and with others: {self.ancillary_hits or 'none'}"
        )

