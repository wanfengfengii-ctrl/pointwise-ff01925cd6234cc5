"""Bounded equality calibration of base and replicate weights."""
from dataclasses import dataclass

from . import qp
from .design import ReplicateSnapshot
from .ledger import _name, _number

_RESIDUAL_TOLERANCE = 1e-9
_BOUND_TOLERANCE = 1e-9


class CalibrationInfeasible(ValueError):
    """The calibration equations cannot be met within the weight bounds."""


@dataclass(frozen=True)
class CalibrationBundle:
    """Immutable package pinning one revision to its calibrated weights."""

    source: ReplicateSnapshot
    weights: tuple[float, ...]
    replicate_weights: tuple[tuple[float, ...], ...]

    @property
    def revision(self):
        return self.source.revision

    @property
    def rows(self):
        return self.source.rows

    @property
    def row_ids(self):
        return self.source.row_ids

    @property
    def replicate_ids(self):
        return self.source.replicate_ids

    @property
    def coefficients(self):
        return self.source.coefficients


def _finite_number(value):
    return _number(value)


def _calibrate(base, matrix, targets, lower, upper):
    """Unique bounded QP solution on one support; None when infeasible."""
    n = len(base)
    support = [i for i, weight in enumerate(base) if weight > 0.0]
    if not support:
        if all(target == 0.0 for target in targets):
            return [0.0] * n
        return None
    sub = [[matrix[j][i] for i in support] for j in range(len(targets))]
    lo = [lower * base[i] for i in support]
    hi = [upper * base[i] for i in support]
    denominators = [
        max(1.0, abs(targets[j]),
            sum(abs(sub[j][k]) * base[support[k]] for k in range(len(support))))
        for j in range(len(targets))
    ]
    diagonal = [1.0 / base[i] for i in support]
    # Objective 0.5 sum (w-d)^2/d has gradient D w - 1, so q is all ones.
    linear = [1.0] * len(support)
    try:
        point = qp.bounded_minimum(
            diagonal, linear, sub, list(targets), lo, hi,
            denominators,
        )
    except qp.InfeasibleError:
        return None
    weights = [0.0] * n
    for k, i in enumerate(support):
        value = point[k]
        if abs(value - lo[k]) <= _BOUND_TOLERANCE:
            value = lo[k]
        elif abs(value - hi[k]) <= _BOUND_TOLERANCE:
            value = hi[k]
        weights[i] = value
    return weights


def _verify(weights, base, matrix, targets, lower, upper):
    for i, weight in enumerate(weights):
        if base[i] == 0.0:
            if weight != 0.0:
                raise CalibrationInfeasible("zero support must remain zero")
            continue
        if weight < lower * base[i] - 1e-8 or weight > upper * base[i] + 1e-8:
            raise CalibrationInfeasible("calibrated weight outside its bounds")
    for j, target in enumerate(targets):
        residual = sum(matrix[j][i] * weights[i] for i in range(len(base))) - target
        denom = max(1.0, abs(target),
                    sum(abs(matrix[j][i]) * base[i] for i in range(len(base))))
        if abs(residual) > 1e-8 * denom:
            raise CalibrationInfeasible("calibration equation is not satisfied")


class Calibrator:
    @staticmethod
    def fit(snapshot, controls, targets, lower, upper):
        """Calibrate the population and every replicate independently."""
        if not isinstance(snapshot, ReplicateSnapshot):
            raise ValueError("expected a ReplicateSnapshot")
        rows = snapshot.rows
        if len(rows) > 60:
            raise ValueError("at most sixty rows are accepted")
        if isinstance(controls, (str, bytes)):
            raise ValueError("controls must be a name sequence")
        names = tuple(_name(name) for name in controls)
        if not 1 <= len(names) <= 4 or len(set(names)) != len(names):
            raise ValueError("expected one to four distinct control names")
        unknown = [name for name in names if name not in snapshot.control_names]
        if unknown:
            raise ValueError("unknown control name")
        if isinstance(targets, (str, bytes)):
            raise ValueError("targets must be a numeric sequence")
        target_values = tuple(_finite_number(value) for value in targets)
        if len(target_values) != len(names):
            raise ValueError("targets must match the controls")
        bound_low = _finite_number(lower)
        bound_high = _finite_number(upper)
        if not 0.0 < bound_low <= 1.0 <= bound_high <= 10.0:
            raise ValueError("require 0 < lower <= 1 <= upper <= 10")
        columns = [snapshot.control_names.index(name) for name in names]
        matrix = [
            [float(rows[i].controls[column]) for i in range(len(rows))]
            for column in columns
        ]

        def accepted(weights):
            for i, weight in enumerate(weights):
                if weight < 0.0:
                    raise ValueError("snapshot weights must be nonnegative")
                if weight > 0.0 and not 0.1 <= weight <= 100.0:
                    raise ValueError("positive fit weights must lie in [0.1, 100]")

        population = tuple(snapshot.weights)
        accepted(population)
        for column in snapshot.replicate_weights:
            accepted(column)
        if not 1 <= len(snapshot.replicate_ids) <= 8:
            raise ValueError("expected one to eight replicates")

        fitted_population = _calibrate(
            population, matrix, target_values, bound_low, bound_high)
        if fitted_population is None:
            raise CalibrationInfeasible("population calibration is infeasible")
        _verify(fitted_population, population, matrix, target_values,
                bound_low, bound_high)
        fitted_replicates = []
        for replicate_id, base in zip(snapshot.replicate_ids,
                                      snapshot.replicate_weights):
            fitted = _calibrate(base, matrix, target_values,
                                bound_low, bound_high)
            if fitted is None:
                raise CalibrationInfeasible(
                    "replicate calibration is infeasible: " + replicate_id)
            _verify(fitted, base, matrix, target_values, bound_low, bound_high)
            fitted_replicates.append(tuple(fitted))
        return CalibrationBundle(
            snapshot, tuple(fitted_population), tuple(fitted_replicates))
