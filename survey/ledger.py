"""Atomic row batches and immutable, ordered survey revisions."""
from dataclasses import dataclass
from math import isfinite
from typing import Iterable


def _name(value):
    if not isinstance(value, str) or not value:
        raise ValueError("identifiers must be nonempty strings")
    return value


def _number(value):
    if type(value) not in (int, float):
        raise ValueError("expected a finite int or float, excluding bool")
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError("number is not finite") from exc
    if not isfinite(result):
        raise ValueError("number is not finite")
    return result


def _names(values):
    if isinstance(values, (str, bytes)):
        raise ValueError("expected an identifier sequence")
    result = tuple(_name(v) for v in values)
    if len(set(result)) != len(result):
        raise ValueError("duplicate identifiers")
    return result


@dataclass(frozen=True)
class Row:
    row_id: str
    weight: float
    controls: tuple[float, ...]
    value: float
    group: str

    def __post_init__(self):
        _name(self.row_id)
        _name(self.group)
        weight = _number(self.weight)
        value = _number(self.value)
        controls = tuple(_number(x) for x in self.controls)
        if not 0.1 <= weight <= 100:
            raise ValueError("base weight must be in [0.1, 100]")
        if not 1 <= len(controls) <= 4:
            raise ValueError("expected one to four controls")
        if abs(value) > 100 or any(abs(x) > 100 for x in controls):
            raise ValueError("outcomes and controls must be in [-100, 100]")
        object.__setattr__(self, "weight", weight)
        object.__setattr__(self, "controls", controls)
        object.__setattr__(self, "value", value)


@dataclass(frozen=True)
class LedgerSnapshot:
    revision: int
    control_names: tuple[str, ...]
    rows: tuple[Row, ...]

    def __post_init__(self):
        if type(self.revision) is not int or self.revision < 0:
            raise ValueError("revision must be a nonnegative integer")
        names = _names(self.control_names)
        rows = tuple(self.rows)
        if not 1 <= len(names) <= 4:
            raise ValueError("expected one to four named controls")
        if any(not isinstance(row, Row) for row in rows):
            raise ValueError("rows must contain Row objects")
        if len({r.row_id for r in rows}) != len(rows):
            raise ValueError("duplicate row IDs")
        if any(len(r.controls) != len(names) for r in rows):
            raise ValueError("control dimension mismatch")
        object.__setattr__(self, "control_names", names)
        object.__setattr__(self, "rows", rows)

    @property
    def row_ids(self):
        return tuple(row.row_id for row in self.rows)


class Ledger:
    """Each accepted apply call advances revision, including a no-op batch."""

    def __init__(self, control_names, rows: Iterable[Row] = ()):
        self._current = LedgerSnapshot(0, control_names, rows)

    @property
    def revision(self):
        return self._current.revision

    def snapshot(self):
        return self._current

    def apply(self, upserts=(), deletes=(), order=None):
        """Validate a complete candidate before publishing one new revision.

        Existing IDs keep their slots; new IDs append in upsert order.
        order, when supplied, must name every final row exactly once.
        Source iterator exceptions propagate without publishing the batch.
        """
        incoming = tuple(upserts)
        removed = _names(deletes)
        if any(not isinstance(row, Row) for row in incoming):
            raise ValueError("upserts must contain Row objects")
        ids = tuple(row.row_id for row in incoming)
        if len(set(ids)) != len(ids):
            raise ValueError("duplicate upserts")
        if set(ids).intersection(removed):
            raise ValueError("an ID cannot be deleted and upserted together")
        candidate = {row.row_id: row for row in self._current.rows}
        if not set(removed).issubset(candidate):
            raise ValueError("cannot delete an unknown ID")
        for row_id in removed:
            del candidate[row_id]
        for row in incoming:
            candidate[row.row_id] = row
        if order is not None:
            requested = _names(order)
            if set(requested) != set(candidate):
                raise ValueError("order must be a permutation of final IDs")
            candidate = {row_id: candidate[row_id] for row_id in requested}
        next_snapshot = LedgerSnapshot(
            self.revision + 1, self._current.control_names,
            tuple(candidate.values()),
        )
        self._current = next_snapshot
        return next_snapshot
