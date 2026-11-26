"""Adaptive terminal reporting for MFEM NURBS benchmark results."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .ingest import (
    BenchmarkDataset,
    fmt_int_like,
    fmt_speedup,
    fmt_time,
    is_valid_number,
)

_VERIF_THRESHOLD = 1e-8
_REASON_MAX_LEN = 120


def _is_explicit_false(val: Any) -> bool:
    return isinstance(val, (bool, np.bool_)) and (not bool(val))


def _fmt_param_val(val: Any) -> str:
    if val is None:
        return "--"
    if isinstance(val, (bool, np.bool_)):
        return "true" if bool(val) else "false"
    if is_valid_number(val):
        return fmt_int_like(val, thousands=False)
    try:
        if pd.isna(val):
            return "--"
    except (TypeError, ValueError):
        pass
    return str(val)


def _extract_concise_reason(stderr_tail: Any) -> str:
    if stderr_tail is None:
        return "(no stderr)"
    try:
        if pd.isna(stderr_tail):
            return "(no stderr)"
    except (TypeError, ValueError):
        pass

    text = str(stderr_tail).strip()
    if not text:
        return "(no stderr)"

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return "(no stderr)"

    chosen = lines[-1]
    for line in lines:
        low = line.lower()
        if (
            "verification failed:" in low
            or "mfem abort" in low
            or "error" in low
        ):
            chosen = line
            break

    if len(chosen) > _REASON_MAX_LEN:
        return chosen[: _REASON_MAX_LEN - 3] + "..."
    return chosen


def _print_header(dataset: BenchmarkDataset) -> None:
    df = dataset.df
    exp_name = dataset.result_dir.name or str(dataset.result_dir)
    print(f"=== Benchmark Report: {exp_name} ===")

    total = len(df)
    status_series = (
        df["status"].astype(str)
        if "status" in df.columns
        else pd.Series(dtype=str)
    )
    ok_mask = status_series == "ok"
    n_ok = int(ok_mask.sum())

    n_nc = 0
    if "converged" in df.columns and total > 0:
        n_nc = int(
            sum(
                1
                for is_ok, conv in zip(ok_mask.tolist(), df["converged"].tolist())
                if is_ok and _is_explicit_false(conv)
            )
        )

    n_timeout = int((status_series == "timeout").sum())
    n_oom = int((status_series == "oom").sum())
    n_crash = int((status_series == "crash").sum())

    known = {"ok", "timeout", "oom", "crash"}
    other_counts: dict[str, int] = {}
    for st, cnt in status_series.value_counts().items():
        if st not in known:
            other_counts[str(st)] = int(cnt)

    status_parts = [
        f"total={total}",
        f"ok={n_ok}",
        f"NC={n_nc}",
        f"timeout={n_timeout}",
        f"oom={n_oom}",
        f"crash={n_crash}",
    ]
    for st, cnt in sorted(other_counts.items()):
        status_parts.append(f"{st}={cnt}")
    print("Runs:   " + "  ".join(status_parts))

    meta = dataset.env_metadata
    env_parts: list[str] = []
    git_commit = meta.get("git_commit")
    if git_commit:
        short_sha = str(git_commit)[:8]
        if meta.get("git_dirty"):
            short_sha += "*"
        env_parts.append(f"git={short_sha}")
    cpu_model = meta.get("cpu_model")
    if cpu_model:
        env_parts.append(f"cpu={cpu_model}")
    if env_parts:
        print("Env:    " + "  |  ".join(env_parts))

    if dataset.constant_params:
        fixed_str = ", ".join(
            f"{k}={_fmt_param_val(v)}"
            for k, v in dataset.constant_params.items()
        )
        print(f"Fixed:  {fixed_str}")


def _print_summary_table(dataset: BenchmarkDataset) -> None:
    df = dataset.df
    print()
    if df.empty:
        print("(No benchmark runs to display.)")
        return

    active_speedups = [
        (col, label)
        for col, label in dataset.speedup_cols
        if col in df.columns and df[col].notna().any()
    ]

    axis_cols = [c for c in dataset.varying_axes if c in df.columns]

    headers: list[str] = [
        *axis_cols,
        "DoFs",
        "status",
        "iters",
        "asm_s",
        "setup_s",
        "solve_s",
        *[label for _, label in active_speedups],
    ]

    align_right: list[bool] = []
    for col in axis_cols:
        is_num = pd.api.types.is_numeric_dtype(df[col]) and not pd.api.types.is_bool_dtype(df[col])
        align_right.append(bool(is_num))
    align_right.extend([True, False, True, True, True, True])
    align_right.extend([True] * len(active_speedups))

    has_dofs = "num_free_dofs" in df.columns
    has_conv = "converged" in df.columns
    has_iters = "iterations" in df.columns
    has_asm = "operator_assembly" in df.columns
    has_setup = "prec_setup_total" in df.columns
    has_solve = "solve" in df.columns

    formatted_rows: list[list[str]] = []
    for idx in df.index:
        status = str(df.at[idx, "status"]) if "status" in df.columns else "unknown"
        is_ok = status == "ok"
        is_nc = is_ok and has_conv and _is_explicit_false(df.at[idx, "converged"])

        row_cells: list[str] = [_fmt_param_val(df.at[idx, col]) for col in axis_cols]

        dofs_val = df.at[idx, "num_free_dofs"] if (is_ok and has_dofs) else None
        row_cells.append(fmt_int_like(dofs_val, thousands=True))

        row_cells.append("ok(NC)" if is_nc else status)

        if not is_ok:
            row_cells.append("--")
        elif is_nc:
            row_cells.append("NC")
        else:
            iters_val = df.at[idx, "iterations"] if has_iters else None
            row_cells.append(fmt_int_like(iters_val, thousands=False))

        asm_val = df.at[idx, "operator_assembly"] if (is_ok and has_asm) else None
        setup_val = df.at[idx, "prec_setup_total"] if (is_ok and has_setup) else None
        solve_val = df.at[idx, "solve"] if (is_ok and has_solve) else None
        row_cells.append(fmt_time(asm_val))
        row_cells.append(fmt_time(setup_val))
        row_cells.append(fmt_time(solve_val))

        for sp_col, _ in active_speedups:
            sp_val = df.at[idx, sp_col] if is_ok else None
            row_cells.append(fmt_speedup(sp_val))

        formatted_rows.append(row_cells)

    widths = [len(h) for h in headers]
    for r in formatted_rows:
        for i, cell in enumerate(r):
            if len(cell) > widths[i]:
                widths[i] = len(cell)

    def _render_row(cells: list[str]) -> str:
        parts = [
            cell.rjust(widths[i]) if align_right[i] else cell.ljust(widths[i])
            for i, cell in enumerate(cells)
        ]
        return "  ".join(parts)

    print(_render_row(headers))
    print("  ".join("-" * w for w in widths))
    for r in formatted_rows:
        print(_render_row(r))


def _print_verification_digest(dataset: BenchmarkDataset) -> None:
    df = dataset.df
    if df.empty or "status" not in df.columns:
        return

    check_cols = [
        c for c in ("prec_sym_err", "prec_stale_y_err") if c in df.columns
    ]
    if not check_cols:
        return

    anomalies: list[tuple[str, list[str]]] = []
    for idx in df.index:
        if str(df.at[idx, "status"]) != "ok":
            continue
        issues: list[str] = []
        for col in check_cols:
            val = df.at[idx, col]
            if is_valid_number(val) and float(val) > _VERIF_THRESHOLD:
                issues.append(f"{col}={float(val):.2e}")
        if issues:
            rid = str(df.at[idx, "run_id"]) if "run_id" in df.columns else str(idx)
            anomalies.append((rid, issues))

    if not anomalies:
        return

    print()
    print(
        f"Verification Warnings ({len(anomalies)} ok run(s) exceeded {_VERIF_THRESHOLD:.0e}):"
    )
    for rid, issues in anomalies[:10]:
        print(f"  - {rid}: {', '.join(issues)}")
    if len(anomalies) > 10:
        print(f"  ... and {len(anomalies) - 10} more")


def _print_failure_digest(dataset: BenchmarkDataset) -> None:
    df = dataset.df
    if df.empty or "status" not in df.columns:
        return

    failed_groups: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    nc_runs: list[str] = []

    has_rc = "error_return_code" in df.columns
    has_stderr = "error_stderr_tail" in df.columns
    has_conv = "converged" in df.columns
    has_rid = "run_id" in df.columns

    for idx in df.index:
        status = str(df.at[idx, "status"])
        rid = str(df.at[idx, "run_id"]) if has_rid else str(idx)
        if status == "ok":
            if has_conv and _is_explicit_false(df.at[idx, "converged"]):
                nc_runs.append(rid)
            continue

        rc_val = df.at[idx, "error_return_code"] if has_rc else None
        rc_str = fmt_int_like(rc_val, thousands=False)
        reason = _extract_concise_reason(
            df.at[idx, "error_stderr_tail"] if has_stderr else None
        )
        failed_groups[(status, rc_str, reason)].append(rid)

    if not failed_groups and not nc_runs:
        return

    print()
    print("Failure & Non-Convergence Digest:")

    sorted_groups = sorted(
        failed_groups.items(),
        key=lambda item: (-len(item[1]), item[0][0], item[0][1], item[0][2]),
    )
    for (status, rc_str, reason), rids in sorted_groups:
        rc_tag = f" rc={rc_str}" if rc_str != "--" else ""
        examples = ", ".join(rids[:3])
        extra = f" (+{len(rids) - 3} more)" if len(rids) > 3 else ""
        print(f"  - [{status}{rc_tag}] x{len(rids)}: {reason}")
        print(f"      runs: {examples}{extra}")

    if nc_runs:
        examples = ", ".join(nc_runs[:3])
        extra = f" (+{len(nc_runs) - 3} more)" if len(nc_runs) > 3 else ""
        print(f"  - [ok(NC)] x{len(nc_runs)}: solver did not converge")
        print(f"      runs: {examples}{extra}")


def _print_generated_outputs(generated_paths: list[Path]) -> None:
    print()
    print("Generated Outputs:")
    if not generated_paths:
        print("  (none)")
        return
    for path in generated_paths:
        print(f"  - {path}")


def print_console_report(
    dataset: BenchmarkDataset,
    generated_paths: list[Path],
) -> None:
    """Print an adaptive summary report, digests, and output manifest."""
    _print_header(dataset)
    _print_summary_table(dataset)
    _print_verification_digest(dataset)
    _print_failure_digest(dataset)
    _print_generated_outputs(generated_paths)
