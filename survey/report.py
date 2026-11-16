"""Point estimates using base weights from one fixed replication snapshot."""
from collections.abc import Mapping
from dataclasses import dataclass
from math import fsum

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
class RawReport:
    revision: int
    estimates: tuple[DomainRatio, ...]


class Reporter:
    @staticmethod
    def raw_ratios(snapshot, domains):
        """Preserve domain mapping order; repeated member IDs denote a set.

        Unknown IDs invalidate the whole call. An empty domain has numerator
        and denominator zero and ratio None. Results contain no live aliases.
        """
        if not isinstance(snapshot, ReplicateSnapshot):
            raise ValueError("expected a ReplicateSnapshot")
        if not isinstance(domains, Mapping):
            raise ValueError("domains must be a mapping")
        known = set(snapshot.row_ids)
        prepared = []
        for name, row_ids in domains.items():
            _name(name)
            if isinstance(row_ids, (str, bytes)):
                raise ValueError("domain membership must be an ID iterable")
            members = frozenset(_name(row_id) for row_id in row_ids)
            if not members.issubset(known):
                raise ValueError("unknown domain row ID")
            prepared.append((name, members))
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
