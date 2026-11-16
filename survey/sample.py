"""A finite household survey used by the CLI and examples."""
from . import DeleteGroups, Ledger, Row


def household():
    ledger = Ledger(("population", "adults"), [
        Row("h1", 1, (1, 1), 0, "north"),
        Row("h2", 2, (1, 2), 1, "south"),
        Row("h3", 1.5, (1, 1), 1, "north"),
        Row("h4", 0.5, (1, 3), 0, "south"),
    ])
    design = DeleteGroups(("north", "south"))
    domains = {
        "all": {"h1", "h2", "h3", "h4"},
        "one-adult": {"h1", "h3"},
        "empty": set(),
    }
    return ledger, design, domains
