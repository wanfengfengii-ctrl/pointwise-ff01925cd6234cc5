"""Deterministic demonstration; no files are written."""
import argparse
import json
from dataclasses import asdict

from . import Reporter, Row
from .sample import household


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("demo",))
    parser.parse_args(argv)
    ledger, design, domains = household()
    original = design.build(ledger.snapshot())
    before = Reporter.raw_ratios(original, domains)
    ledger.apply(
        upserts=[Row("h1", 1, (1, 1), 1, "north")],
        order=["h4", "h3", "h2", "h1"],
    )
    current = design.build(ledger.snapshot())
    payload = {
        "before": asdict(before),
        "old_snapshot_after_edit": asdict(
            Reporter.raw_ratios(original, domains)
        ),
        "current": asdict(Reporter.raw_ratios(current, domains)),
        "current_row_ids": current.row_ids,
        "replicate_ids": current.replicate_ids,
        "replicate_weights": current.replicate_weights,
        "coefficients": current.coefficients,
    }
    print(json.dumps(payload, ensure_ascii=False, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
