"""Core analysis package for the MFEM NURBS Kronecker-product benchmark suite."""

from __future__ import annotations

from .ingest import BenchmarkDataset, apply_filters, load_dataset
from .metrics import (
    apply_speedups,
    compute_derived_metrics,
    parse_baseline_spec,
    save_summary_csv,
)
from .report import print_console_report

__all__ = [
    "BenchmarkDataset",
    "load_dataset",
    "apply_filters",
    "compute_derived_metrics",
    "apply_speedups",
    "parse_baseline_spec",
    "save_summary_csv",
    "print_console_report",
]
