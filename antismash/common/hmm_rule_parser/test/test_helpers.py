# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,too-many-public-methods

import unittest

from antismash.common.hmm_rule_parser.helpers import (
    AndConditionSatisfied as AndSatisfied,
    ConditionSatisfaction as BaseSatisfied,
    OrConditionSatisfied as OrSatisfied,
    Presence,
    PresenceAbsence as States,
)


class TestPresence(unittest.TestCase):
    def test_blocker(self):
        assert Presence(present=True, negated=True).blocks()
        for pres, neg in [(True, False), (False, True), (False, False)]:
            assert not Presence(present=pres, negated=neg).blocks()

    def test_falsey(self):
        assert not Presence(present=False, negated=False)
        assert not Presence(present=True, negated=True)

    def test_truthy(self):
        assert Presence(present=False, negated=True)
        assert Presence(present=True, negated=False)


class TestSatisfied(unittest.TestCase):
    def setUp(self):
        self.good_present = Presence(present=True, negated=False)
        self.good_missing = Presence(present=False, negated=True)
        self.blocker = Presence(present=True, negated=True)
        self.bad = Presence(present=False, negated=False)

    def test_simple(self):
        def check(sat):
            return sat.satisfied and not sat.blocks()
        assert check(BaseSatisfied(anchor="test", matches_in_anchor=States({"a": self.good_present})))
        # no matches in anchor
        assert not check(BaseSatisfied(anchor="test", matches_in_neighbours={"other": States({"a": self.good_present})}))

        for presence in [self.bad]:
            check(BaseSatisfied(anchor="test", negated=True, matches_in_anchor=States({"a": presence})))
            check(BaseSatisfied(anchor="test", negated=True, matches_in_neighbours={"other": States({"a": presence})}))

#    def test_simple_negated(self):
#        # the case where the whole condition is negated
#        def check(sat):
#            assert not sat.satisfied
#            assert sat.blocks()

#        for presence in [self.good_present, self.good_missing]:
#            check(BaseSatisfied(anchor="test", negated=True, matches_in_anchor=States({"a": presence})))
#            check(BaseSatisfied(anchor="test", negated=True, matches_in_neighbours={"other": States({"a": presence})}))


    def test_children(self):
        a = States({"a": self.good_present})
        b = States({"b": self.good_present})
        anchor_a = BaseSatisfied(anchor="test", matches_in_anchor=a)
        anchor_b = BaseSatisfied(anchor="test", matches_in_anchor=b)
        anchors = BaseSatisfied(anchor="test", child_states=[anchor_a, anchor_b])

        print(anchors)
        assert anchors.positive_matches_in_anchor == {"a", "b"}
        assert not anchors.positive_matches_in_neighbours
        assert anchors.satisfied

        neighbour_b = BaseSatisfied(anchor="test", matches_in_neighbours=States({"x": b}))
        assert not neighbour_b.matches_in_anchor
        mixed = BaseSatisfied(anchor="test", child_states=[anchor_a, neighbour_b])

        print(mixed)
        assert mixed.positive_matches_in_anchor == {"a"}
        assert mixed.positive_matches_in_neighbours == {"b"}
        assert mixed.satisfied

        neighbour_a = BaseSatisfied(anchor="test", matches_in_neighbours=States({"y": a}))
        neighbours = BaseSatisfied(anchor="test", child_states=[neighbour_a, neighbour_b])
        assert not neighbour_a.matches_in_anchor
        assert not neighbour_b.matches_in_anchor
        assert not neighbours.positive_matches_in_anchor
        assert neighbours.positive_matches_in_neighbours == {"a", "b"}
        print(neighbours)
        assert not neighbours.satisfied  # the anchor requires matches
