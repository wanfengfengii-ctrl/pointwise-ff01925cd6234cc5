"""Deterministic demonstrations; no files are written."""
import argparse
import json
from dataclasses import asdict

from . import Calibrator, DeleteGroups, Ledger, Reporter, Row
from .sample import household


def _demo():
    """Original uncalibrated demonstration, unchanged in semantics."""
    ledger, design, domains = household()
    original = design.build(ledger.snapshot())
    before = Reporter.raw_ratios(original, domains)
    ledger.apply(
        upserts=[Row("h1", 1, (1, 1), 1, "north")],
        order=["h4", "h3", "h2", "h1"],
    )
    current = design.build(ledger.snapshot())
    return {
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


def _calibration_demo():
    """Calibrated point estimates, replicate variance and missing slots."""
    # Six rows in two deletion groups keep every replicate feasible for
    # the population-style targets used below.
    ledger = Ledger(("population", "adults"), [
        Row("n1", 1.0, (1, 1), 0, "north"),
        Row("n2", 2.0, (1, 2), 1, "north"),
        Row("n3", 1.5, (1, 1), 1, "north"),
        Row("s1", 1.0, (1, 2), 0, "south"),
        Row("s2", 2.5, (1, 3), 1, "south"),
        Row("s3", 0.5, (1, 1), 0, "south"),
    ])
    design = DeleteGroups(("north", "south"))
    snapshot = design.build(ledger.snapshot())
    # A single population-total control leaves every replicate feasible
    # while still exercising bounds, redundant information and variance.
    targets = (
        round(1.2 * sum(row.weight for row in snapshot.rows), 4),
    )
    bundle = Calibrator.fit(
        snapshot, ("population",), targets, 0.5, 2.0)
    domains = {
        "all": snapshot.row_ids,
        "north-only": {"n1", "n2", "n3"},
        "empty": set(),
    }
    calibrated = Reporter.ratios(bundle, domains)
    # A correction must not alter the already-published calibrated package.
    ledger.apply(
        upserts=[Row("n1", 1.0, (1, 1), 1, "north")],
        order=["s3", "s2", "s1", "n3", "n2", "n1"],
    )
    after_edit = Reporter.ratios(bundle, domains)
    rebuilt = design.build(ledger.snapshot())
    refit = Calibrator.fit(
        rebuilt, ("population",), targets, 0.5, 2.0)
    return {
        "targets": targets,
        "calibrated_population_weights": bundle.weights,
        "calibrated_replicate_weights": bundle.replicate_weights,
        "calibrated_report": asdict(calibrated),
        "old_bundle_after_edit": asdict(after_edit),
        "refit_revision": refit.revision,
        "refit_report": asdict(Reporter.ratios(refit, domains)),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command",
                        choices=("demo", "calibration-demo", "all"))
    args = parser.parse_args(argv)
    if args.command == "demo":
        print(json.dumps(_demo(), ensure_ascii=False, sort_keys=True,
                         indent=2))
    elif args.command == "calibration-demo":
        print(json.dumps(_calibration_demo(), ensure_ascii=False,
                         sort_keys=True, indent=2))
    else:
        # One valid JSON document combining both demonstrations.
        print(json.dumps(
            {"uncalibrated_demo": _demo(),
             "calibration_demo": _calibration_demo()},
            ensure_ascii=False, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
