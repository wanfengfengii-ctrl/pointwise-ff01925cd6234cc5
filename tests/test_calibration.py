"""Independent optimality and calibrated-variance verification.

The production solver in survey.qp is a floating-point active-set method.
These checks derive every small-case optimum independently with exact
Fraction arithmetic: each bounded-variable partition is enumerated, the
free variables solve an exact equality-constrained KKT system, and the
unique feasible stationary point with the smallest objective is the
reference optimum.  No production solver is used to certify itself.
"""
import random
import unittest
from fractions import Fraction
from itertools import product

from survey import (
    Calibrator, CalibrationInfeasible, DeleteGroups, Ledger, Reporter,
    ReplicateSnapshot, Row,
)


def _exact_solve(matrix, rhs):
    """Gaussian elimination over Fraction; returns None when inconsistent."""
    rows = [list(map(Fraction, matrix[i])) + [Fraction(rhs[i])]
            for i in range(len(matrix))]
    width = len(matrix[0]) if matrix else 0
    rank = 0
    for column in range(width):
        pivot = next((r for r in range(rank, len(rows))
                      if rows[r][column] != 0), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_row = rows[rank]
        factor = pivot_row[column]
        for j in range(column, width + 1):
            pivot_row[j] /= factor
        for r in range(len(rows)):
            if r == rank or rows[r][column] == 0:
                continue
            scale = rows[r][column]
            for j in range(column, width + 1):
                rows[r][j] -= scale * pivot_row[j]
        rank += 1
    for row in rows:
        if all(v == 0 for v in row[:width]) and row[width] != 0:
            return None
    answer = [Fraction(0)] * width
    for row in rows:
        leading = next((c for c in range(width) if row[c] != 0), None)
        if leading is not None:
            answer[leading] = row[width]
    return answer


def _safe_solve(matrix, rhs):
    """Exact square solve returning None on singular or inconsistent input."""
    if len(matrix) != len(matrix[0]):
        return None
    try:
        return _exact_solve(matrix, rhs)
    except (IndexError, ZeroDivisionError):
        return None


def _independent_rows(rows):
    """Indices of a maximal independent subset, computed over exact values."""
    if not rows:
        return []
    width = len(rows[0])
    # Carry the original row index alongside each data row.
    work = [list(row) for row in rows]
    tags = list(range(len(rows)))
    rank = 0
    chosen = []
    for column in range(width):
        pivot = next((r for r in range(rank, len(work))
                      if work[r][column] != 0), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        tags[rank], tags[pivot] = tags[pivot], tags[rank]
        pivot_row = work[rank]
        factor = pivot_row[column]
        for r in range(len(work)):
            if r == rank or work[r][column] == 0:
                continue
            scale = work[r][column] / factor
            for j in range(column, width):
                work[r][j] -= scale * pivot_row[j]
        chosen.append(tags[rank])
        rank += 1
    return chosen


def exact_reference(base, matrix, targets, lower, upper):
    """Exact bounded QP optimum, or None when infeasible.

    Minimizes sum (w_i-d_i)^2/d_i subject to matrix*w = targets and
    lower*d_i <= w_i <= upper*d_i on the positive support.
    """
    base = [Fraction(v) for v in base]
    matrix = [[Fraction(v) for v in row] for row in matrix]
    targets = [Fraction(v) for v in targets]
    support = [i for i, d in enumerate(base) if d > 0]
    if not support:
        return [Fraction(0)] * len(base) if all(t == 0 for t in targets) else None
    k = len(support)
    sub = [[row[i] for i in support] for row in matrix]
    lo = [lower * base[i] for i in support]
    inv_d = [1 / base[i] for i in support]
    best = None
    best_objective = None
    for partition in product((-1, 0, 1), repeat=k):
        fixed = {p: (lo[p] if sign == -1 else upper * base[support[p]])
                 for p, sign in enumerate(partition) if sign != 0}
        free = [p for p, sign in enumerate(partition) if sign == 0]
        point = [Fraction(0)] * k
        for p, value in fixed.items():
            point[p] = value
        adjusted = [
            targets[j] - sum(sub[j][p] * point[p] for p in fixed)
            for j in range(len(targets))
        ]
        if free:
            free_count = len(free)
            free_block = [[sub[j][p] for p in free]
                          for j in range(len(targets))]
            independent = _independent_rows(free_block)
            rank = len(independent)
            if rank > free_count:
                # Over-determined: a free-variable solution exists only if
                # the equalities are mutually consistent. When the
                # independent row count equals the free width the system
                # is square; otherwise no unique free solution exists.
                if rank != free_count:
                    continue
                values = _safe_solve(
                    [free_block[r] for r in independent],
                    [adjusted[r] for r in independent])
                if values is None:
                    continue
            else:
                kkt = [[Fraction(0)] * (free_count + rank)
                       for _ in range(free_count + rank)]
                for p in range(free_count):
                    kkt[p][p] = inv_d[free[p]]
                for t, row_index in enumerate(independent):
                    for p in range(free_count):
                        value = free_block[row_index][p]
                        # Canonical target KKT: D z* + A'mu = 1, A z* = b.
                        kkt[p][free_count + t] = value
                        kkt[free_count + t][p] = value
                rhs = [Fraction(1)] * free_count
                rhs.extend(adjusted[row_index] for row_index in independent)
                solution = _safe_solve(kkt, rhs)
                if solution is None:
                    continue
                values = solution[:free_count]
            for p, value in zip(free, values):
                point[p] = value
        else:
            if any(adjusted[j] != 0 for j in range(len(targets))):
                continue
        if any(point[p] < lo[p] or point[p] > upper * base[support[p]]
               for p in range(k)):
            continue
        # Enforce every original equality, including rows that were
        # linearly dependent within the free-variable sub-matrix. The
        # inputs arrive as floats, so compare on a float-scale tolerance
        # rather than demanding exact Fraction equality.
        if any(abs(float(sum(sub[j][p] * point[p] for p in range(k)))
                   - float(targets[j])) > 1e-9
               for j in range(len(targets))):
            continue
        objective = sum(((point[p] - base[support[p]]) ** 2) * inv_d[p]
                        for p in range(k))
        if best_objective is None or objective < best_objective:
            best_objective = objective
            best = point
    if best is None:
        return None
    weights = [Fraction(0)] * len(base)
    for p, i in enumerate(support):
        weights[i] = best[p]
    return weights


class ExactReferenceTests(unittest.TestCase):
    def _case(self, base, controls, targets, lower, upper):
        reference = exact_reference(base, controls, targets, lower, upper)
        rows = [Row("r" + str(i), float(base[i]),
                    tuple(float(controls[j][i]) for j in range(len(targets))),
                    0, "g" + str(i)) for i in range(len(base))]
        names = tuple("c" + str(j) for j in range(len(targets)))
        # An external single-replicate design reuses the same positive
        # base weights, so population and replicate share one exact optimum.
        snapshot = ReplicateSnapshot(
            Ledger(names, rows).snapshot(), ("rep",),
            (tuple(float(v) for v in base),), (1.0,))
        if reference is None:
            with self.assertRaises(CalibrationInfeasible):
                Calibrator.fit(snapshot, names, targets, lower, upper)
            return
        bundle = Calibrator.fit(snapshot, names, targets, lower, upper)
        for label, fitted, d in (
                ("population", bundle.weights, snapshot.weights),
                ("rep", bundle.replicate_weights[0],
                 snapshot.replicate_weights[0])):
            ref = exact_reference(d, controls, targets, lower, upper)
            self.assertIsNotNone(ref, label + " unexpectedly infeasible")
            for actual, expected in zip(fitted, ref):
                self.assertLess(abs(actual - float(expected)),
                                1e-7 * max(1.0, abs(float(expected))))
            objective = sum((fitted[i] - d[i]) ** 2 / d[i]
                            for i in range(len(d)) if d[i] > 0)
            ref_objective = sum(
                ((ref[i] - Fraction(d[i])) ** 2) / Fraction(d[i])
                for i in range(len(d)) if d[i] > 0)
            self.assertLess(
                abs(objective - float(ref_objective)),
                1e-7 * max(1.0, abs(float(ref_objective))))
            for j, target in enumerate(targets):
                residual = sum(controls[j][i] * fitted[i]
                               for i in range(len(d))) - target
                denominator = max(
                    1.0, abs(target),
                    sum(abs(controls[j][i]) * d[i] for i in range(len(d))))
                self.assertLessEqual(abs(residual) / denominator, 1e-8)

    def test_redundant_controls_and_zero_columns(self):
        self._case([1, 1, 1], [[1, 1, 1], [2, 2, 2], [0, 0, 0]],
                   [3, 6, 0], 0.5, 2)
        self._case([2, 1, 1], [[1, -1, 0], [0, 0, 0]],
                   [1, 0], 0.5, 3)

    def test_boundary_optima_lower_and_upper(self):
        self._case([1, 1], [[1, 2]], [0.5], 0.5, 2)
        self._case([1, 1], [[1, 2]], [6], 0.5, 2)
        self._case([1, 2, 1], [[1, 1, 1], [0, 1, 2]],
                   [0.5, 0.25], 0.5, 2)

    def test_zero_support_and_infeasible_replicates(self):
        # Population feasible; deleting one group leaves one row whose
        # maximum control total cannot reach the target.
        base = [1, 1]
        controls = [[1, 2]]
        rows = [Row("a", 1, (1,), 0, "A"), Row("b", 1, (2,), 0, "B")]
        snapshot = DeleteGroups(("A", "B")).build(
            Ledger(("n",), rows).snapshot())
        self.assertIsNotNone(exact_reference(base, controls, [5], 0.5, 2))
        with self.assertRaises(CalibrationInfeasible):
            Calibrator.fit(snapshot, ("n",), (5,), 0.5, 2)
        # Empty support with a nonzero target is infeasible.
        self.assertIsNone(exact_reference([0, 0], [[1, 1]], [1], 0.5, 2))
        self.assertEqual(exact_reference([0, 0], [[1, 1]], [0], 0.5, 2),
                         [Fraction(0), Fraction(0)])

    def test_randomized_exact_small_cases(self):
        rng = random.Random(99173)
        for _ in range(25):
            n = rng.randrange(2, 7)
            m = rng.randrange(1, 3)
            base = [rng.randrange(1, 6) for _ in range(n)]
            controls = [[rng.randrange(-3, 4) for _ in range(n)]
                        for _ in range(m)]
            raw_targets = [
                sum(controls[j][i] * base[i] for i in range(n))
                for j in range(m)
            ]
            scale = rng.choice((0.6, 0.9, 1.0, 1.4, 3.0))
            targets = [round(value * scale) for value in raw_targets]
            lower = rng.choice((0.25, 0.5, 0.8))
            upper = rng.choice((1.5, 2.0, 4.0))
            self._case(base, controls, targets, lower, upper)


class CalibrationWorkflowTests(unittest.TestCase):
    def test_spec_worked_example_values(self):
        ledger = Ledger(("population",), [
            Row("a", 1, (1,), 0, "A"),
            Row("b", 1, (1,), 1, "B"),
        ])
        snapshot = DeleteGroups(("A", "B")).build(ledger.snapshot())
        bundle = Calibrator.fit(snapshot, ("population",), (3,), 0.5, 2)
        self.assertEqual(bundle.revision, 0)
        for actual, expected in zip(bundle.weights, (1.5, 1.5)):
            self.assertAlmostEqual(actual, expected)
        for column, expected in zip(bundle.replicate_weights,
                                    ((0, 3), (3, 0))):
            for actual, value in zip(column, expected):
                self.assertAlmostEqual(actual, value)
        report = Reporter.ratios(
            bundle, {"all": {"a", "b"}, "a-only": {"a"}, "empty": set()})
        all_ratio, a_only, empty = report.estimates
        self.assertAlmostEqual(all_ratio.ratio, 0.5)
        self.assertEqual(len(all_ratio.replicate_ratios), 2)
        self.assertAlmostEqual(all_ratio.replicate_ratios[0], 1)
        self.assertAlmostEqual(all_ratio.replicate_ratios[1], 0)
        self.assertAlmostEqual(all_ratio.variance, 0.25)
        self.assertEqual(all_ratio.missing_replicate_ids, ())
        self.assertEqual(a_only.ratio, 0)
        self.assertEqual(a_only.replicate_ratios, (None, 0))
        self.assertIsNone(a_only.variance)
        self.assertEqual(a_only.missing_replicate_ids, ("A",))
        self.assertIsNone(empty.ratio)
        self.assertIsNone(empty.variance)
        self.assertEqual(empty.replicate_ratios, (None, None))
        self.assertEqual(empty.missing_replicate_ids, ("A", "B"))
        self.assertEqual(empty.count, 0)

    def test_joint_calibration_distinct_from_population_only(self):
        rows = [
            Row("a", 1, (1, 0), 0, "A"),
            Row("b", 1, (1, 1), 0, "B"),
            Row("c", 1, (1, 2), 0, "A"),
            Row("d", 1, (1, 3), 1, "B"),
        ]
        snapshot = DeleteGroups(("A", "B")).build(
            Ledger(("n", "x"), rows).snapshot())
        bundle = Calibrator.fit(snapshot, ("n", "x"), (4, 6), 0.5, 2)
        for column, expected in zip(
                bundle.replicate_weights, ((0, 3, 0, 1), (1, 0, 3, 0))):
            for actual, value in zip(column, expected):
                self.assertAlmostEqual(actual, value)
        estimate = Reporter.ratios(
            bundle, {"all": {"a", "b", "c", "d"}}).estimates[0]
        self.assertAlmostEqual(estimate.ratio, 0.25)
        self.assertAlmostEqual(estimate.replicate_ratios[0], 0.25)
        self.assertAlmostEqual(estimate.replicate_ratios[1], 0)
        self.assertAlmostEqual(estimate.variance, 1 / 32)

    def test_external_snapshot_redundant_bounds(self):
        source = Ledger(
            ("n", "x"),
            [Row("a", 1, (1, 2), 0, "g"), Row("b", 1, (1, 2), 0, "g")],
        ).snapshot()
        snapshot = ReplicateSnapshot(
            source, ("r",), ((2, 2),), (1,))
        # Redundant controls are accepted; target (4, 8) pins the population
        # at its own upper bound 2 and gives the replicate weight 2.
        bundle = Calibrator.fit(snapshot, ("n", "x"), (4, 8), 0.5, 2)
        for actual in bundle.weights:
            self.assertAlmostEqual(actual, 2)
        for actual in bundle.replicate_weights[0]:
            self.assertAlmostEqual(actual, 2)
        # The inconsistent target breaks the duplicated control equation.
        with self.assertRaises(CalibrationInfeasible):
            Calibrator.fit(snapshot, ("n", "x"), (4, 10), 0.5, 2)

    def test_bundle_survives_ledger_edit_and_reorder(self):
        ledger = Ledger(("population",), [
            Row("a", 1, (1,), 0, "A"),
            Row("b", 1, (1,), 1, "B"),
        ])
        snapshot = DeleteGroups(("A", "B")).build(ledger.snapshot())
        bundle = Calibrator.fit(snapshot, ("population",), (3,), 0.5, 2)
        before = Reporter.ratios(
            bundle, {"all": {"a", "b"}}).estimates[0]
        ledger.apply([Row("a", 1, (1,), 1, "A")], order=["b", "a"])
        after = Reporter.ratios(
            bundle, {"all": {"a", "b"}}).estimates[0]
        self.assertEqual(before, after)
        self.assertEqual(bundle.row_ids, ("a", "b"))
        rebuilt = DeleteGroups(("A", "B")).build(ledger.snapshot())
        refit = Calibrator.fit(rebuilt, ("population",), (3,), 0.5, 2)
        report = Reporter.ratios(refit, {"all": {"a", "b"}}).estimates[0]
        self.assertAlmostEqual(report.ratio, 1)
        self.assertTrue(all(abs(value - 1) < 1e-8
                            for value in report.replicate_ratios))
        self.assertAlmostEqual(report.variance, 0)
        # Reordering alone preserves weights keyed by stable row ID.
        ledger.apply(order=["a", "b"])
        ordered = Calibrator.fit(
            DeleteGroups(("A", "B")).build(ledger.snapshot()),
            ("population",), (3,), 0.5, 2)
        mapping = dict(zip(ordered.row_ids, ordered.weights))
        for row_id, weight in zip(bundle.row_ids, bundle.weights):
            self.assertAlmostEqual(mapping[row_id], weight)

    def test_input_validation_and_immutability(self):
        rows = [Row("a", 1, (1, 0), 0, "g"),
                Row("b", 1, (1, 1), 1, "g")]
        snapshot = ReplicateSnapshot(
            Ledger(("n", "x"), rows).snapshot(),
            ("r",), ((1, 1),), (0.5,))
        with self.assertRaises(ValueError):
            Calibrator.fit(snapshot, ("unknown",), (1,), 0.5, 2)
        with self.assertRaises(ValueError):
            Calibrator.fit(snapshot, ("n", "n"), (1, 1), 0.5, 2)
        with self.assertRaises(ValueError):
            Calibrator.fit(snapshot, ("n",), (1,), 0.0, 2)
        with self.assertRaises(ValueError):
            Calibrator.fit(snapshot, ("n",), (1,), 0.5, 11)
        with self.assertRaises(ValueError):
            Calibrator.fit(snapshot, ("n",), (1,), 1.5, 2)
        with self.assertRaises(ValueError):
            Calibrator.fit(snapshot, ("n",), (True,), 0.5, 2)
        with self.assertRaises(ValueError):
            Calibrator.fit(snapshot, (), (), 0.5, 2)
        with self.assertRaises(ValueError):
            Calibrator.fit(snapshot, ("n",), (1, 2), 0.5, 2)
        with self.assertRaises(ValueError):
            Reporter.ratios(snapshot, {"d": ["a"]})
        bundle = Calibrator.fit(snapshot, ("n", "x"), (2, 1), 0.5, 2)
        with self.assertRaises((TypeError, AttributeError)):
            bundle.weights = (1, 1)
        domains = {"d": ["a"]}
        report = Reporter.ratios(bundle, domains)
        domains["d"].append("b")
        self.assertEqual(report.estimates[0].row_ids, ("a",))
        with self.assertRaises(ValueError):
            Reporter.ratios(bundle, {"d": ["missing"]})
        with self.assertRaises(ValueError):
            Reporter.ratios(bundle, {"": []})

    def test_sixty_rows_four_controls_eight_replicates(self):
        rows = [
            Row(str(i), 1 + i % 5,
                (1, i % 7 - 3, i % 3 - 1, 2 * (i % 3) - 2),
                i % 2, str(i % 8))
            for i in range(60)
        ]
        names = ("n", "x", "z", "twice_z")
        snapshot = DeleteGroups(tuple(map(str, range(8)))).build(
            Ledger(names, rows).snapshot())
        targets = (
            sum(row.controls[0] * row.weight for row in rows),
            sum(row.controls[1] * row.weight for row in rows),
            sum(row.controls[2] * row.weight for row in rows),
            sum(row.controls[3] * row.weight for row in rows),
        )
        bundle = Calibrator.fit(snapshot, names, targets, 0.5, 2)
        self.assertEqual(len(bundle.replicate_weights), 8)
        self.assertTrue(all(len(column) == 60
                            for column in bundle.replicate_weights))
        report = Reporter.ratios(
            bundle, {"all": snapshot.row_ids, "empty": []})
        overall, empty = report.estimates
        self.assertAlmostEqual(overall.ratio, 0.5)
        self.assertIsNone(empty.ratio)
        self.assertIsNone(empty.variance)
        self.assertEqual(len(empty.missing_replicate_ids), 8)
        for r, column in enumerate(bundle.replicate_weights):
            base = snapshot.replicate_weights[r]
            for j, target in enumerate(targets):
                residual = sum(rows[i].controls[j] * column[i]
                               for i in range(60)) - target
                denominator = max(
                    1.0, abs(target),
                    sum(abs(rows[i].controls[j]) * base[i]
                        for i in range(60)))
                self.assertLessEqual(abs(residual) / denominator, 1e-8)
            for i, weight in enumerate(column):
                if base[i] == 0:
                    self.assertEqual(weight, 0)
                else:
                    self.assertTrue(0.5 * base[i] - 1e-8 <= weight <=
                                    2 * base[i] + 1e-8)

    def test_large_scale_null_space_optimality(self):
        """No feasible equality-preserving perturbation lowers the objective."""
        from survey import qp
        rng = random.Random(1618)
        n = 48
        group_count = 6
        rows = [
            Row(str(i), rng.randrange(1, 9),
                tuple(rng.randrange(-5, 6) for _ in range(3)),
                rng.randrange(-8, 9), str(i % group_count))
            for i in range(n)
        ]
        names = ("n", "x", "z")
        snapshot = DeleteGroups(
            tuple(str(g) for g in range(group_count))).build(
            Ledger(names, rows).snapshot())
        targets = tuple(
            round(sum(row.controls[j] * row.weight for row in rows) * 1.2)
            for j in range(3))
        bundle = Calibrator.fit(snapshot, names, targets, 0.5, 2)
        columns = [(bundle.weights, snapshot.weights)] + list(zip(
            bundle.replicate_weights, snapshot.replicate_weights))
        for fitted, base in columns:
            support = [i for i in range(n) if base[i] > 0]
            k = len(support)
            d = [base[i] for i in support]
            w = [fitted[i] for i in support]
            A = [[float(rows[i].controls[j]) for i in support]
                 for j in range(3)]
            lo = [0.5 * value for value in d]
            hi = [2.0 * value for value in d]
            kept = qp.independent_rows(A, k)
            Ar = [[A[r][p] for p in range(k)] for r in kept]
            rank = len(kept)
            objective = sum((w[p] - d[p]) ** 2 / d[p] for p in range(k))

            def null_vector(seed_vector):
                kkt = [[0.0] * (k + rank) for _ in range(k + rank)]
                for p in range(k):
                    kkt[p][p] = 1.0
                for t in range(rank):
                    for p in range(k):
                        kkt[p][k + t] = Ar[t][p]
                        kkt[k + t][p] = Ar[t][p]
                return qp.solve_linear(kkt, seed_vector + [0.0] * rank)[:k]

            for _ in range(20):
                direction = null_vector(
                    [rng.uniform(-1, 1) for _ in range(k)])
                step = rng.uniform(-0.05, 0.05) * min(
                    high - low for low, high in zip(lo, hi))
                candidate = [
                    min(hi[p], max(lo[p], w[p] + step * direction[p]))
                    for p in range(k)]
                residual = max(
                    abs(sum(Ar[t][p] * candidate[p] for p in range(k))
                        - sum(Ar[t][p] * w[p] for p in range(k)))
                    for t in range(rank))
                if residual > 1e-6:
                    continue
                candidate_objective = sum(
                    (candidate[p] - d[p]) ** 2 / d[p] for p in range(k))
                self.assertGreaterEqual(
                    candidate_objective,
                    objective - 1e-7 * max(1.0, abs(objective)))

    def test_calibrated_variance_matches_independent_formula(self):
        rng = random.Random(44201)
        group_count = 4
        rows = [
            Row("r" + str(i), rng.randrange(1, 8) / 2,
                (1, rng.randrange(-3, 4), rng.randrange(-2, 3)),
                rng.randrange(-5, 6), "g" + str(i % group_count))
            for i in range(12)
        ]
        names = ("n", "x", "z")
        snapshot = DeleteGroups(tuple("g" + str(g)
                                      for g in range(group_count))).build(
            Ledger(names, rows).snapshot())
        targets = tuple(
            round(sum(row.controls[j] * row.weight for row in rows))
            for j in range(3))
        bundle = Calibrator.fit(snapshot, names, targets, 0.5, 2)
        domains = {"d" + str(k): {row.row_id for row in rows
                                  if rng.randrange(2) == 0}
                   for k in range(5)}
        report = Reporter.ratios(bundle, domains)
        for estimate, members in zip(report.estimates, domains.values()):
            indices = [i for i, row in enumerate(rows)
                       if row.row_id in members]
            denominator = sum(bundle.weights[i] for i in indices)
            numerator = sum(bundle.weights[i] * rows[i].value
                            for i in indices)
            self.assertAlmostEqual(estimate.denominator, denominator)
            self.assertAlmostEqual(estimate.numerator, numerator)
            replicate_ratios = []
            missing = []
            for r, column in enumerate(bundle.replicate_weights):
                den = sum(column[i] for i in indices)
                if den == 0:
                    replicate_ratios.append(None)
                    missing.append(bundle.replicate_ids[r])
                    continue
                replicate_ratios.append(
                    sum(column[i] * rows[i].value for i in indices) / den)
            self.assertEqual(estimate.missing_replicate_ids, tuple(missing))
            for actual, expected in zip(estimate.replicate_ratios,
                                        replicate_ratios):
                if expected is None:
                    self.assertIsNone(actual)
                else:
                    self.assertAlmostEqual(actual, expected)
            if denominator == 0 or missing:
                self.assertIsNone(estimate.variance)
            else:
                theta = numerator / denominator
                expected_variance = sum(
                    bundle.coefficients[r] *
                    (replicate_ratios[r] - theta) ** 2
                    for r in range(group_count))
                self.assertAlmostEqual(
                    estimate.variance, expected_variance)


if __name__ == "__main__":
    unittest.main()
