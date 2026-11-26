#!/usr/bin/env python3
"""Thin CLI entrypoint for analyzing MFEM NURBS Kronecker benchmark results."""

from __future__ import annotations

import argparse
import dataclasses
import sys
from pathlib import Path
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11 fallback
    import tomli as tomllib  # type: ignore[no-redef]

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from analysis.ingest import apply_filters, load_dataset
from analysis.metrics import (
    apply_speedups,
    compute_derived_metrics,
    parse_baseline_spec,
    save_summary_csv,
)
from analysis.plotting import generate_all_plots
from analysis.report import print_console_report

_VALID_FORMATS = frozenset({"png", "pdf", "svg"})


def _parse_formats(raw_format: str) -> list[str]:
    """Parse and validate a comma-separated list of output plot formats."""
    formats: list[str] = []
    for token in raw_format.split(","):
        fmt = token.strip().lower().lstrip(".")
        if not fmt:
            continue
        if fmt not in _VALID_FORMATS:
            valid_list = ", ".join(sorted(_VALID_FORMATS))
            raise ValueError(
                f"Unsupported plot format '{fmt}'. Valid formats: {valid_list}"
            )
        if fmt not in formats:
            formats.append(fmt)
    if not formats:
        raise ValueError("At least one plot format must be specified (e.g. 'png').")
    return formats


def _load_plots_config(config_path: Path) -> list[dict[str, Any]]:
    """Load ``[[plots]]`` specifications from an external TOML file."""
    if not config_path.exists():
        raise FileNotFoundError(f"Plots config file not found: {config_path}")
    if not config_path.is_file():
        raise ValueError(f"Plots config path is not a file: {config_path}")

    with config_path.open("rb") as fh:
        data = tomllib.load(fh)

    raw_plots = data.get("plots", [])
    if not isinstance(raw_plots, list) or not all(
        isinstance(entry, dict) for entry in raw_plots
    ):
        raise ValueError(
            f"Invalid [[plots]] section in '{config_path}': expected an array of tables."
        )
    return list(raw_plots)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="analyze_results.py",
        description=(
            "Aggregate MFEM NURBS Kronecker benchmark JSON results, compute "
            "speedups, export summary CSV, generate plots/tables, and print a console report."
        ),
    )
    parser.add_argument(
        "result_dir",
        type=Path,
        help="Path to the benchmark results directory containing run JSON logs.",
    )
    parser.add_argument(
        "--plots-config",
        type=Path,
        default=None,
        metavar="PATH",
        help=(
            "Optional TOML file whose [[plots]] entries override config_snapshot.toml "
            "so plots and tables can be iterated without re-running benchmarks."
        ),
    )
    parser.add_argument(
        "--include-defaults",
        action="store_true",
        help="Also generate the minimal built-in academic default plot set.",
    )
    parser.add_argument(
        "--baseline",
        action="append",
        default=None,
        metavar="SPEC",
        help=(
            "Repeatable speedup baseline specification, e.g. "
            "'--baseline prec=jacobi' or compound '--baseline assembly=legacy,prec=none'."
        ),
    )
    parser.add_argument(
        "--baseline-prec",
        default="none",
        metavar="PREC",
        help="Default baseline preconditioner when --baseline is not given (default: none).",
    )
    parser.add_argument(
        "--baseline-assembly",
        default="legacy",
        metavar="ASM",
        help="Default baseline assembly mode when --baseline is not given (default: legacy).",
    )
    parser.add_argument(
        "--filter",
        action="append",
        default=None,
        metavar="EXPR",
        help=(
            "Repeatable subset filter, e.g. '--filter mesh_name=square_dist' "
            "or '--filter order=5,7'."
        ),
    )
    parser.add_argument(
        "--compare-axis",
        default=None,
        metavar="AXIS",
        help="Pivot default plot hue and stacked-bar series onto any parameter column.",
    )
    parser.add_argument(
        "--format",
        default="png",
        metavar="FORMATS",
        help="Comma-separated plot output formats from {png, pdf, svg} (default: png).",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip plot/table generation and only write summary CSV and console report.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        formats = _parse_formats(args.format)
        baseline_specs = (
            [parse_baseline_spec(spec) for spec in args.baseline]
            if args.baseline
            else None
        )
        override_plot_specs = (
            _load_plots_config(args.plots_config)
            if args.plots_config is not None
            else None
        )

        dataset = load_dataset(args.result_dir)
        if override_plot_specs is not None:
            dataset = dataclasses.replace(
                dataset, custom_plot_specs=override_plot_specs
            )

        if args.filter:
            dataset = apply_filters(dataset, args.filter)

        dataset = compute_derived_metrics(dataset)
        dataset = apply_speedups(
            dataset,
            baseline_specs=baseline_specs,
            default_prec=args.baseline_prec,
            default_assembly=args.baseline_assembly,
        )

        generated_paths: list[Path] = [save_summary_csv(dataset)]

        if not args.no_plots:
            plots_dir = dataset.result_dir / "plots"
            plot_paths = generate_all_plots(
                dataset,
                plots_dir=plots_dir,
                formats=formats,
                compare_axis=args.compare_axis,
                include_defaults=args.include_defaults,
            )
            generated_paths.extend(plot_paths)

        print_console_report(dataset, generated_paths)
        return 0
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130
    except (
        FileNotFoundError,
        NotADirectoryError,
        ValueError,
        KeyError,
        RuntimeError,
        tomllib.TOMLDecodeError,
    ) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
