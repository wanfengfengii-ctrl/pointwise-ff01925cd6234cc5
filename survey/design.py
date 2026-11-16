"""Public delete-one-group replication of an immutable ledger snapshot."""
from dataclasses import dataclass

from .ledger import LedgerSnapshot, _names, _number


@dataclass(frozen=True)
class ReplicateSnapshot:
    source: LedgerSnapshot
    replicate_ids: tuple[str, ...]
    replicate_weights: tuple[tuple[float, ...], ...]
    coefficients: tuple[float, ...]

    def __post_init__(self):
        if not isinstance(self.source, LedgerSnapshot):
            raise ValueError("source must be a LedgerSnapshot")
        ids = _names(self.replicate_ids)
        weights = tuple(
            tuple(_number(w) for w in column)
            for column in self.replicate_weights
        )
        coefficients = tuple(_number(c) for c in self.coefficients)
        if not 1 <= len(ids) <= 8:
            raise ValueError("expected one to eight replicates")
        if len(weights) != len(ids) or len(coefficients) != len(ids):
            raise ValueError("replicate metadata length mismatch")
        if any(len(column) != len(self.rows) for column in weights):
            raise ValueError("replicate row length mismatch")
        if any(w < 0 for column in weights for w in column):
            raise ValueError("replicate weights must be nonnegative")
        if any(c < 0 for c in coefficients):
            raise ValueError("coefficients must be nonnegative")
        object.__setattr__(self, "replicate_ids", ids)
        object.__setattr__(self, "replicate_weights", weights)
        object.__setattr__(self, "coefficients", coefficients)

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
    def control_names(self):
        return self.source.control_names

    @property
    def weights(self):
        return tuple(row.weight for row in self.rows)


@dataclass(frozen=True)
class DeleteGroups:
    """Unstratified equal-group jackknife, with configured group order."""
    groups: tuple[str, ...]

    def __post_init__(self):
        groups = _names(self.groups)
        if not 2 <= len(groups) <= 8:
            raise ValueError("expected two to eight deletion groups")
        object.__setattr__(self, "groups", groups)

    def build(self, snapshot):
        if not isinstance(snapshot, LedgerSnapshot):
            raise ValueError("expected a LedgerSnapshot")
        present = {row.group for row in snapshot.rows}
        if present != set(self.groups):
            raise ValueError("every configured group must occur, and no others")
        count = len(self.groups)
        factor = count / (count - 1)
        columns = tuple(
            tuple(
                0.0 if row.group == group else row.weight * factor
                for row in snapshot.rows
            )
            for group in self.groups
        )
        return ReplicateSnapshot(
            snapshot, self.groups, columns,
            ((count - 1) / count,) * count,
        )
