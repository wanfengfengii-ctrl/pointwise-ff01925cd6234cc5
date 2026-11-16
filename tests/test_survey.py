import io
import json
import random
import unittest
from contextlib import redirect_stdout
from dataclasses import FrozenInstanceError
from fractions import Fraction

from survey import (
    DeleteGroups, Ledger, LedgerSnapshot, ReplicateSnapshot, Reporter, Row,
)
from survey.__main__ import main
from survey.sample import household


class LedgerTests(unittest.TestCase):
    def test_atomic_batch_and_stable_identity(self):
        ledger, _, _ = household()
        old = ledger.snapshot()
        new = ledger.apply(
            [Row("h2", 3, (1, 2), 0, "south"),
             Row("h5", 1, (1, 2), 1, "north")],
            deletes=["h4"],
        )
        self.assertEqual(new.row_ids, ("h1", "h2", "h3", "h5"))
        self.assertEqual(new.revision, 1)
        self.assertEqual(new.rows[1].weight, 3)
        self.assertEqual(old.rows[1].weight, 2)
        self.assertEqual(old.row_ids, ("h1", "h2", "h3", "h4"))
        ledger.apply(order=reversed(new.row_ids))
        self.assertEqual(ledger.snapshot().row_ids, tuple(reversed(new.row_ids)))
        self.assertEqual(ledger.revision, 2)
        ledger.apply()
        self.assertEqual(ledger.revision, 3)

    def test_rejected_batches_leave_same_snapshot(self):
        ledger, _, _ = household()
        saved = ledger.snapshot()
        row = Row("h1", 1, (1, 1), 1, "north")
        cases = [
            dict(upserts=[row, row]),
            dict(upserts=[row], deletes=["h1"]),
            dict(deletes=["missing"]),
            dict(deletes=["h1", "h1"]),
            dict(order=["h1", "h2", "h3"]),
            dict(order=["h1", "h2", "h3", "missing"]),
            dict(order=["h1", "h2", "h3", "h3"]),
            dict(upserts=[Row("bad", 1, (1,), 0, "north")]),
            dict(upserts=[object()]),
        ]
        for kwargs in cases:
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):
                    ledger.apply(**kwargs)
                self.assertIs(ledger.snapshot(), saved)

    def test_source_failure_propagates_without_partial_update(self):
        ledger, _, _ = household()
        saved = ledger.snapshot()
        failure = RuntimeError("row source interrupted")

        def source():
            yield Row("new", 1, (1, 1), 0, "north")
            raise failure

        with self.assertRaises(RuntimeError) as caught:
            ledger.apply(source())
        self.assertIs(caught.exception, failure)
        self.assertIs(ledger.snapshot(), saved)

    def test_snapshot_and_inputs_have_no_mutable_aliases(self):
        controls = [1, 2]
        row = Row("a", 1, controls, 0, "x")
        rows = [row]
        names = ["n", "adults"]
        ledger = Ledger(names, rows)
        controls[0] = 90
        names[0] = "changed"
        rows.clear()
        snap = ledger.snapshot()
        self.assertEqual(snap.control_names, ("n", "adults"))
        self.assertEqual(snap.rows[0].controls, (1, 2))
        with self.assertRaises(FrozenInstanceError):
            snap.rows[0].value = 99
        with self.assertRaises(FrozenInstanceError):
            snap.revision = 99

    def test_empty_ledger_and_complete_replacement(self):
        ledger = Ledger(("n",))
        self.assertEqual(ledger.snapshot().rows, ())
        ledger.apply([Row("a", 1, (1,), 0, "x")], order=["a"])
        ledger.apply(deletes=["a"], order=[])
        self.assertEqual(ledger.snapshot().rows, ())

    def test_invalid_rows_and_metadata(self):
        for weight in [True, 0, -1, 101, float("nan"), float("inf"), 10**400]:
            with self.subTest(weight=weight):
                with self.assertRaises(ValueError):
                    Row("a", weight, (1,), 0, "x")
        for controls in [(), (1, 2, 3, 4, 5), (True,), (101,), (float("nan"),)]:
            with self.assertRaises(ValueError):
                Row("a", 1, controls, 0, "x")
        for value in [True, -101, float("inf")]:
            with self.assertRaises(ValueError):
                Row("a", 1, (1,), value, "x")
        for row_id, group in [("", "x"), ("a", ""), (1, "x")]:
            with self.assertRaises(ValueError):
                Row(row_id, 1, (1,), 0, group)
        row = Row("a", 1, (1,), 0, "x")
        for names, rows in [((), []), (("n", "n"), []), (("n",), [row, row]),
                            (("n", "m"), [row]), ("n", []), (("n",), [1])]:
            with self.assertRaises(ValueError):
                Ledger(names, rows)
        for revision in [-1, True, 1.5]:
            with self.assertRaises(ValueError):
                LedgerSnapshot(revision, ("n",), ())


class DesignTests(unittest.TestCase):
    def test_delete_groups_weights_and_coefficients(self):
        ledger, design, _ = household()
        source = ledger.snapshot()
        snap = design.build(source)
        self.assertIs(snap.source, source)
        self.assertEqual(snap.replicate_ids, ("north", "south"))
        self.assertEqual(snap.replicate_weights, ((0, 4, 0, 1), (2, 0, 3, 0)))
        self.assertEqual(snap.coefficients, (0.5, 0.5))
        self.assertEqual(snap.weights, (1, 2, 1.5, 0.5))
        self.assertEqual(snap.control_names, ("population", "adults"))

    def test_reorder_binds_every_replica_by_row_id(self):
        ledger, design, _ = household()
        old = design.build(ledger.snapshot())
        ledger.apply(order=["h3", "h1", "h4", "h2"])
        new = design.build(ledger.snapshot())
        for before, after in zip(old.replicate_weights, new.replicate_weights):
            self.assertEqual(dict(zip(old.row_ids, before)), dict(zip(new.row_ids, after)))
        self.assertEqual(old.revision, 0)
        self.assertEqual(new.revision, 1)

    def test_invalid_designs_and_missing_groups(self):
        for groups in [[], ["x"], ["x", "x"], list("abcdefghi"), "xy"]:
            with self.assertRaises(ValueError):
                DeleteGroups(groups)
        design = DeleteGroups(("x", "y"))
        for rows in [[], [Row("a", 1, (1,), 0, "x")],
                     [Row("a", 1, (1,), 0, "x"), Row("b", 1, (1,), 0, "z")]]:
            with self.assertRaises(ValueError):
                design.build(Ledger(("n",), rows).snapshot())
        with self.assertRaises(ValueError):
            design.build(None)

    def test_public_replica_snapshot_validation_and_copying(self):
        source = Ledger(("n",), [Row("a", 1, (1,), 1, "x")]).snapshot()
        weights, ids, coefficients = [[0]], ["custom"], [0]
        snap = ReplicateSnapshot(source, ids, weights, coefficients)
        weights[0][0] = 5
        ids[0] = "changed"
        coefficients[0] = 5
        self.assertEqual(snap.replicate_weights, ((0,),))
        self.assertEqual(snap.replicate_ids, ("custom",))
        self.assertEqual(snap.coefficients, (0,))
        with self.assertRaises(FrozenInstanceError):
            snap.coefficients = (1,)
        for bad_ids, bad_weights, bad_coeffs in [
            ([], [], []), (["r"], [], [1]), (["r"], [[1, 2]], [1]),
            (["r"], [[-1]], [1]), (["r"], [[float("inf")]], [1]),
            (["r"], [[1]], [-1]), (["r"], [[1]], [float("nan")]),
            (["r", "r"], [[1], [1]], [1, 1]),
            (["r"], [[True]], [1]), (["r"], [[1]], []),
        ]:
            with self.assertRaises(ValueError):
                ReplicateSnapshot(source, bad_ids, bad_weights, bad_coeffs)


class ReportTests(unittest.TestCase):
    def test_hand_computed_domains_and_empty_domain(self):
        ledger, design, domains = household()
        report = Reporter.raw_ratios(design.build(ledger.snapshot()), domains)
        self.assertEqual(report.revision, 0)
        all_rows, adults, empty = report.estimates
        self.assertEqual((all_rows.count, all_rows.numerator, all_rows.denominator),
                         (4, 3.5, 5))
        self.assertAlmostEqual(all_rows.ratio, 0.7)
        self.assertAlmostEqual(adults.ratio, 0.6)
        self.assertEqual((empty.count, empty.numerator, empty.denominator, empty.ratio),
                         (0, 0, 0, None))

    def test_report_is_frozen_across_edit_and_reorder(self):
        ledger, design, domains = household()
        old = design.build(ledger.snapshot())
        before = Reporter.raw_ratios(old, domains)
        ledger.apply([Row("h1", 1, (1, 1), 1, "south")],
                     order=["h4", "h3", "h2", "h1"])
        current = design.build(ledger.snapshot())
        self.assertEqual(Reporter.raw_ratios(old, domains), before)
        self.assertAlmostEqual(Reporter.raw_ratios(current, domains).estimates[0].ratio, 0.9)
        domains["all"].clear()
        self.assertEqual(before.estimates[0].count, 4)
        with self.assertRaises(FrozenInstanceError):
            before.estimates[0].ratio = 0

    def test_domain_set_semantics_order_and_validation(self):
        ledger, design, _ = household()
        snap = design.build(ledger.snapshot())
        report = Reporter.raw_ratios(snap, {"z": ["h2", "h2"], "a": []})
        self.assertEqual(tuple(r.name for r in report.estimates), ("z", "a"))
        self.assertEqual(report.estimates[0].count, 1)
        self.assertEqual(report.estimates[0].ratio, 1)
        self.assertEqual(Reporter.raw_ratios(snap, {}).estimates, ())
        for domains in [{"ok": ["h1"], "bad": ["missing"]}, {"": []},
                        {"x": "h1"}, {"x": [1]}, []]:
            with self.assertRaises(ValueError):
                Reporter.raw_ratios(snap, domains)
        with self.assertRaises(ValueError):
            Reporter.raw_ratios(ledger.snapshot(), {})

    def test_point_report_uses_only_base_weights(self):
        ledger, design, domains = household()
        snap = design.build(ledger.snapshot())
        other = ReplicateSnapshot(snap.source, ("custom",), ((0, 0, 0, 0),), (0,))
        self.assertEqual(Reporter.raw_ratios(snap, domains),
                         Reporter.raw_ratios(other, domains))

    def test_randomized_independent_fraction_reference(self):
        rng = random.Random(28731)
        for case in range(40):
            group_count = rng.randrange(2, 9)
            n = rng.randrange(group_count, 61)
            groups = tuple("g" + str(i) for i in range(group_count))
            rows = [
                Row("r" + str(i), rng.randrange(1, 13) / 4,
                    (1, rng.randrange(-4, 5), i % 3, (i % 3) * 2),
                    rng.randrange(-10, 11), groups[i % group_count])
                for i in range(n)
            ]
            ledger = Ledger(("n", "x", "z", "twice_z"), rows)
            snap = DeleteGroups(groups).build(ledger.snapshot())
            domains = {
                "d" + str(j): {r.row_id for r in rows if rng.randrange(3) == 0}
                for j in range(7)
            }
            report = Reporter.raw_ratios(snap, domains)
            for result, members in zip(report.estimates, domains.values()):
                selected = [r for r in rows if r.row_id in members]
                den = sum((Fraction(r.weight) for r in selected), Fraction())
                num = sum((Fraction(r.weight) * Fraction(r.value) for r in selected),
                          Fraction())
                self.assertEqual(result.denominator, float(den))
                self.assertEqual(result.numerator, float(num))
                if den:
                    self.assertAlmostEqual(result.ratio, float(num / den), places=12)
                else:
                    self.assertIsNone(result.ratio)
            for group, column in zip(groups, snap.replicate_weights):
                for row, actual in zip(rows, column):
                    expected = (Fraction() if row.group == group else
                                Fraction(row.weight) * group_count / (group_count - 1))
                    self.assertAlmostEqual(actual, float(expected), places=12)
            ledger.apply(order=tuple(reversed(snap.row_ids)))
            reversed_report = Reporter.raw_ratios(
                DeleteGroups(groups).build(ledger.snapshot()), domains)
            self.assertEqual(report.estimates, reversed_report.estimates)

    def test_sixty_rows_four_controls_eight_replicates(self):
        rows = [
            Row(str(i), 1 + i % 5, (1, i % 7, i % 3, 2 * (i % 3)),
                i % 2, str(i % 8))
            for i in range(60)
        ]
        snap = DeleteGroups(tuple(map(str, range(8)))).build(
            Ledger(("n", "x", "z", "twice_z"), rows).snapshot())
        report = Reporter.raw_ratios(snap, {"all": snap.row_ids, "zero": []})
        self.assertEqual(len(snap.replicate_weights), 8)
        self.assertTrue(all(len(column) == 60 for column in snap.replicate_weights))
        self.assertEqual(report.estimates[0].numerator, 90)
        self.assertEqual(report.estimates[0].denominator, 180)
        self.assertEqual(report.estimates[0].ratio, 0.5)
        self.assertIsNone(report.estimates[1].ratio)

    def test_demo_runs_real_revision_and_replication(self):
        output = io.StringIO()
        with redirect_stdout(output):
            main(["demo"])
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["before"], payload["old_snapshot_after_edit"])
        self.assertEqual(payload["before"]["revision"], 0)
        self.assertEqual(payload["current"]["revision"], 1)
        self.assertAlmostEqual(payload["current"]["estimates"][0]["ratio"], 0.9)
        self.assertEqual(payload["current_row_ids"], ["h4", "h3", "h2", "h1"])
        self.assertEqual(payload["replicate_weights"], [[1, 0, 4, 0], [0, 3, 0, 2]])


if __name__ == "__main__":
    unittest.main()
