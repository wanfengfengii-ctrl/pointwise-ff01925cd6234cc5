"""Point estimates using base weights from one fixed replication snapshot."""
from collections.abc import Mapping
from dataclasses import dataclass
from math import fsum

from .calibration import CalibrationBundle
from .design import ReplicateSnapshot
from .ledger import _name


@dataclass(frozen=True)
class DomainRatio:
    name: str
    row_ids: tuple[str, ...]
    count: int
    numerator: float
    denominator: float
    ratio: float | None


@dataclass(frozen=True)
class CalibratedDomainRatio:
    name: str
    row_ids: tuple[str, ...]
    count: int
    numerator: float
    denominator: float
    ratio: float | None
    replicate_ratios: tuple[float | None, ...]
    variance: float | None
    missing_replicate_ids: tuple[str, ...]


@dataclass(frozen=True)
class RawReport:
    revision: int
    estimates: tuple[DomainRatio, ...]


@dataclass(frozen=True)
class CalibratedReport:
    revision: int
    estimates: tuple[CalibratedDomainRatio, ...]


def _prepare_domains(domains, known_ids):
    if not isinstance(domains, Mapping):
        raise ValueError("domains must be a mapping")
    prepared = []
    for name, row_ids in domains.items():
        _name(name)
        if isinstance(row_ids, (str, bytes)):
            raise ValueError("domain membership must be an ID iterable")
        members = frozenset(_name(row_id) for row_id in row_ids)
        if not members.issubset(known_ids):
            raise ValueError("unknown domain row ID")
        prepared.append((name, members))
    return prepared


class Reporter:
    @staticmethod
    def raw_ratios(snapshot, domains):
        """Preserve domain mapping order; repeated member IDs denote a set.

        Unknown IDs invalidate the whole call. An empty domain has numerator
        and denominator zero and ratio None. Results contain no live aliases.
        """
        if not isinstance(snapshot, ReplicateSnapshot):
            raise ValueError("expected a ReplicateSnapshot")
        prepared = _prepare_domains(domains, set(snapshot.row_ids))
        estimates = []
        for name, members in prepared:
            denominator = fsum(
                row.weight for row in snapshot.rows if row.row_id in members
            )
            numerator = fsum(
                row.weight * row.value
                for row in snapshot.rows if row.row_id in members
            )
            estimates.append(DomainRatio(
                name, tuple(sorted(members)), len(members),
                numerator, denominator,
                None if denominator == 0 else numerator / denominator,
            ))
        return RawReport(snapshot.revision, tuple(estimates))

    @staticmethod
    def ratios(bundle, domains):
        """Calibrated domain ratios with delete-group replicate variance.

        Replicates whose domain denominator is zero keep their slot with a
        None ratio; any missing replicate makes the domain variance None.
        A zero population denominator makes both ratio and variance None.
        """
        if not isinstance(bundle, CalibrationBundle):
            raise ValueError("expected a CalibrationBundle")
        prepared = _prepare_domains(domains, set(bundle.row_ids))
        rows = bundle.rows
        estimates = []
        for name, members in prepared:
            ordered = tuple(sorted(members))
            population_den = fsum(
                weight for row, weight in zip(rows, bundle.weights)
                if row.row_id in members
            )
            population_num = fsum(
                weight * row.value
                for row, weight in zip(rows, bundle.weights)
                if row.row_id in members
            )
            theta = None if population_den == 0 else population_num / population_den
            replicate_ratios = []
            missing = []
            for replicate_id, column in zip(bundle.replicate_ids,
                                            bundle.replicate_weights):
                denominator = fsum(
                    weight for row, weight in zip(rows, column)
                    if row.row_id in members
                )
                if denominator == 0:
                    replicate_ratios.append(None)
                    missing.append(replicate_id)
                    continue
                numerator = fsum(
                    weight * row.value
                    for row, weight in zip(rows, column)
                    if row.row_id in members
                )
                replicate_ratios.append(numerator / denominator)
            if theta is None or missing:
                variance = None
            else:
                variance = fsum(
                    coefficient * (replicate_ratio - theta) ** 2
                    for coefficient, replicate_ratio
                    in zip(bundle.coefficients, replicate_ratios)
                )
            estimates.append(CalibratedDomainRatio(
                name, ordered, len(members), population_num, population_den,
                theta, tuple(replicate_ratios), variance, tuple(missing),
            ))
        return CalibratedReport(bundle.revision, tuple(estimates))
