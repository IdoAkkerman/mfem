#!/usr/bin/env python3
"""Post-processing, aggregation, plotting, and LaTeX table generation for the
MFEM NURBS Kronecker-product preconditioner benchmark suite.

Reads JSON v1.0 run records from <result_dir>/runs/*.json and produces:
  - <result_dir>/summary.csv (and summary.parquet if pyarrow is installed)
  - <result_dir>/tables/summary_table.tex
  - <result_dir>/tables/iterations_pivot.tex
  - <result_dir>/plots/{iterations_vs_dofs,timing_breakdown,convergence}_*.{png,pdf}
"""

from __future__ import annotations

import argparse
import json
import math
import re
import statistics
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants & Styling
# ---------------------------------------------------------------------------

GROUP_KEYS = [
    "mesh_name",
    "problem",
    "ref_strategy",
    "order",
    "ref_levels",
    "assembly",
    "solver",
]

SORT_KEYS = [
    "mesh_name",
    "ref_strategy",
    "order",
    "ref_levels",
    "solver",
    "prec",
    "run_id",
]

# Colorblind-safe Okabe-Ito inspired palette + distinct markers for print legibility
PALETTE = [
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#332288",  # indigo
    "#882255",  # wine
    "#44AA99",  # teal
    "#999933",  # olive
]

MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">"]
LINESTYLES = ["-", "--", "-.", ":"]

TIMING_STAGES: list[tuple[str, str, str]] = [
    ("operator_assembly", "Operator assembly", "#4C78A8"),
    ("prec_univariate_assembly", "Prec: 1D assembly", "#F58518"),
    ("prec_eigendecomp", "Prec: eigendecomp", "#E45756"),
    ("prec_fd_construct", "Prec: FD construct", "#72B7B2"),
    ("prec_other_setup", "Prec: other setup", "#B279A2"),
    ("solve", "Krylov solve", "#54A24B"),
]


def configure_matplotlib() -> None:
    """Configure clean, publication-ready Matplotlib defaults."""
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8.5,
            "legend.frameon": True,
            "legend.framealpha": 0.9,
            "legend.edgecolor": "#CCCCCC",
            "axes.grid": True,
            "grid.linestyle": "--",
            "grid.linewidth": 0.6,
            "grid.alpha": 0.45,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


# ---------------------------------------------------------------------------
# Generic Helpers
# ---------------------------------------------------------------------------


def slugify(value: Any) -> str:
    """Convert arbitrary identifiers into safe filename slugs."""
    text = str(value).strip()
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "unknown"


def latex_escape(value: Any) -> str:
    """Escape special LaTeX characters in plain text cells."""
    text = str(value)
    replacements = [
        ("\\", r"\textbackslash{}"),
        ("&", r"\&"),
        ("%", r"\%"),
        ("$", r"\$"),
        ("#", r"\#"),
        ("_", r"\_"),
        ("{", r"\{"),
        ("}", r"\}"),
        ("~", r"\textasciitilde{}"),
        ("^", r"\textasciicircum{}"),
    ]
    for old, new in replacements:
        text = text.replace(old, new)
    return text


def format_sizes_1d(sizes: Any) -> str:
    """Format a 1D size array such as [19, 19] into a compact string '19x19'."""
    if isinstance(sizes, list) and len(sizes) > 0:
        return "x".join(str(int(x)) for x in sizes)
    return ""


def is_valid_number(val: Any) -> bool:
    """Return True if val is a finite numeric value (excluding booleans)."""
    if val is None or isinstance(val, bool):
        return False
    try:
        return math.isfinite(float(val))
    except (TypeError, ValueError):
        return False


def fmt_int_like(val: Any, thousands: bool = False) -> str:
    """Format integer-like values without trailing '.0'."""
    if not is_valid_number(val):
        return "--"
    ival = int(round(float(val)))
    return f"{ival:,}" if thousands else str(ival)


def fmt_time(val: Any) -> str:
    """Format wall-clock time in seconds cleanly."""
    if not is_valid_number(val):
        return "--"
    fval = float(val)
    if fval < 10.0:
        return f"{fval:.4f}"
    return f"{fval:.2f}"


def fmt_speedup(val: Any, latex: bool = False) -> str:
    """Format speedup factor cleanly."""
    if not is_valid_number(val):
        return "--"
    fval = float(val)
    if latex:
        return f"{fval:.2f}$\\times$"
    return f"{fval:.2f}x"


def fmt_rss(val: Any) -> str:
    """Format memory in MB."""
    if not is_valid_number(val):
        return "--"
    return f"{float(val):.1f}"


def safe_div(num: Any, den: Any) -> float:
    """Return num / den as float, or NaN on invalid/zero denominator."""
    if not is_valid_number(num) or not is_valid_number(den):
        return float("nan")
    d = float(den)
    if d == 0.0:
        return float("nan")
    return float(num) / d


def compute_sample_std(values: Any) -> float:
    """Sample standard deviation of solve_repeats (0.0 if len <= 1)."""
    if not isinstance(values, list):
        return 0.0
    nums = [float(v) for v in values if is_valid_number(v)]
    if len(nums) <= 1:
        return 0.0
    return float(statistics.stdev(nums))


def order_preconditioners(precs: list[str], baseline_prec: str) -> list[str]:
    """Order preconditioners with baseline first, then alphabetically."""
    unique = sorted({str(p) for p in precs if pd.notna(p)})
    if baseline_prec in unique:
        unique.remove(baseline_prec)
        return [baseline_prec] + unique
    return unique


def build_prec_style_map(
    precs: list[str], baseline_prec: str
) -> dict[str, tuple[str, str]]:
    """Assign a consistent (color, marker) pair to each preconditioner."""
    ordered = order_preconditioners(precs, baseline_prec)
    style_map: dict[str, tuple[str, str]] = {}
    for idx, prec in enumerate(ordered):
        color = PALETTE[idx % len(PALETTE)]
        marker = MARKERS[idx % len(MARKERS)]
        style_map[prec] = (color, marker)
    return style_map


def status_token(status: Any, converged: Any, iterations: Any) -> str:
    """Return a concise status/iteration token for tables."""
    st = str(status).lower() if pd.notna(status) else "unknown"
    if st == "ok":
        if converged is False or (pd.notna(converged) and not bool(converged)):
            return "NC"
        return fmt_int_like(iterations)
    mapping = {
        "timeout": "TIMEOUT",
        "oom": "OOM",
        "crash": "CRASH",
    }
    return mapping.get(st, st.upper())


# ---------------------------------------------------------------------------
# Ingestion & Aggregation
# ---------------------------------------------------------------------------


def load_run_records(
    runs_dir: Path,
) -> tuple[list[dict[str, Any]], dict[str, list[float]]]:
    """Load and flatten all JSON v1.0 files under runs_dir."""
    json_files = sorted(runs_dir.glob("*.json"))
    records: list[dict[str, Any]] = []
    residual_histories: dict[str, list[float]] = {}

    for path in json_files:
        try:
            with path.open("r", encoding="utf-8") as fh:
                raw = json.load(fh)
        except (OSError, json.JSONDecodeError) as exc:
            print(f"[WARN] Skipping unreadable JSON {path}: {exc}", file=sys.stderr)
            continue

        if not isinstance(raw, dict):
            print(f"[WARN] Skipping non-object JSON {path}", file=sys.stderr)
            continue

        run_id = str(raw.get("run_id") or path.stem)
        status = str(raw.get("status", "unknown"))
        cfg = raw.get("config") or {}
        stats = raw.get("problem_stats") or {}
        timings = raw.get("timings_sec") or {}
        conv = raw.get("convergence") or {}
        verif = raw.get("verification") or {}
        sys_info = raw.get("system") or {}
        err = raw.get("error") or {}

        is_ok = status == "ok"

        solve_repeats = timings.get("solve_repeats") if is_ok else None
        res_hist = conv.get("residual_history") if is_ok else None
        if isinstance(res_hist, list) and len(res_hist) > 0:
            clean_hist = [float(v) for v in res_hist if is_valid_number(v)]
            if clean_hist:
                residual_histories[run_id] = clean_hist

        effective_order = (
            stats.get("order")
            if is_ok and stats.get("order") is not None
            else cfg.get("order")
        )

        record: dict[str, Any] = {
            # Top-level metadata
            "schema_version": raw.get("schema_version"),
            "run_id": run_id,
            "status": status,
            "source_file": path.name,
            # Config
            "sweep_name": cfg.get("sweep_name"),
            "mesh_name": cfg.get("mesh_name"),
            "mesh": cfg.get("mesh"),
            "problem": cfg.get("problem"),
            "ref_strategy": cfg.get("ref_strategy"),
            "requested_order": cfg.get("requested_order"),
            "config_order": cfg.get("order"),
            "order": effective_order,
            "ref_levels": cfg.get("ref_levels"),
            "assembly": cfg.get("assembly"),
            "solver": cfg.get("solver"),
            "prec": cfg.get("prec"),
            "rtol": cfg.get("rtol"),
            "atol": cfg.get("atol"),
            "max_iter": cfg.get("max_iter"),
            "gmres_restart": cfg.get("gmres_restart"),
            "repeats": cfg.get("repeats"),
            "verify": cfg.get("verify"),
            "record_residuals": cfg.get("record_residuals"),
            # Problem stats
            "dim": stats.get("dim") if is_ok else np.nan,
            "base_order": stats.get("base_order") if is_ok else np.nan,
            "num_elements": stats.get("num_elements") if is_ok else np.nan,
            "num_tdofs": stats.get("num_tdofs") if is_ok else np.nan,
            "num_free_dofs": stats.get("num_free_dofs") if is_ok else np.nan,
            "num_ess_dofs": stats.get("num_ess_dofs") if is_ok else np.nan,
            "full_sizes_1d": format_sizes_1d(stats.get("full_sizes_1d"))
            if is_ok
            else "",
            "free_sizes_1d": format_sizes_1d(stats.get("free_sizes_1d"))
            if is_ok
            else "",
            "matrix_nnz": stats.get("matrix_nnz") if is_ok else np.nan,
            # Timings
            "mesh_load": timings.get("mesh_load") if is_ok else np.nan,
            "refinement": timings.get("refinement") if is_ok else np.nan,
            "fespace_setup": timings.get("fespace_setup") if is_ok else np.nan,
            "operator_assembly": timings.get("operator_assembly") if is_ok else np.nan,
            "rhs_assembly": timings.get("rhs_assembly") if is_ok else np.nan,
            "prec_setup_total": timings.get("prec_setup_total") if is_ok else np.nan,
            "prec_univariate_assembly": (
                timings.get("prec_univariate_assembly") if is_ok else np.nan
            ),
            "prec_eigendecomp": timings.get("prec_eigendecomp") if is_ok else np.nan,
            "prec_fd_construct": timings.get("prec_fd_construct") if is_ok else np.nan,
            "solve": timings.get("solve") if is_ok else np.nan,
            "solve_repeats_count": (
                len(solve_repeats) if isinstance(solve_repeats, list) else 0
            )
            if is_ok
            else 0,
            "solve_std_sec": compute_sample_std(solve_repeats) if is_ok else np.nan,
            "solve_per_iter": timings.get("solve_per_iter") if is_ok else np.nan,
            "total": timings.get("total") if is_ok else np.nan,
            # Convergence
            "converged": conv.get("converged") if is_ok else None,
            "iterations": conv.get("iterations") if is_ok else np.nan,
            "initial_norm": conv.get("initial_norm") if is_ok else np.nan,
            "final_norm": conv.get("final_norm") if is_ok else np.nan,
            "avg_reduction_factor": (
                conv.get("avg_reduction_factor") if is_ok else np.nan
            ),
            "residual_history_len": (
                len(res_hist) if isinstance(res_hist, list) else 0
            )
            if is_ok
            else 0,
            # Verification
            "verification_enabled": verif.get("enabled") if is_ok else None,
            "prec_sym_err": verif.get("prec_sym_err") if is_ok else np.nan,
            "prec_stale_y_err": verif.get("prec_stale_y_err") if is_ok else np.nan,
            # System
            "peak_rss_kb": sys_info.get("peak_rss_kb") if is_ok else np.nan,
            # Error details
            "error_return_code": err.get("return_code") if not is_ok else np.nan,
            "error_wall_time_sec": err.get("wall_time_sec") if not is_ok else np.nan,
            "error_stderr_tail": err.get("stderr_tail") if not is_ok else "",
        }
        records.append(record)

    return records, residual_histories


def build_summary_dataframe(
    records: list[dict[str, Any]], baseline_prec: str
) -> pd.DataFrame:
    """Construct the tidy summary DataFrame with derived metrics and deterministic sort."""
    df = pd.DataFrame.from_records(records)

    numeric_cols = [
        "requested_order",
        "config_order",
        "order",
        "ref_levels",
        "rtol",
        "atol",
        "max_iter",
        "gmres_restart",
        "repeats",
        "dim",
        "base_order",
        "num_elements",
        "num_tdofs",
        "num_free_dofs",
        "num_ess_dofs",
        "matrix_nnz",
        "mesh_load",
        "refinement",
        "fespace_setup",
        "operator_assembly",
        "rhs_assembly",
        "prec_setup_total",
        "prec_univariate_assembly",
        "prec_eigendecomp",
        "prec_fd_construct",
        "solve",
        "solve_std_sec",
        "solve_per_iter",
        "total",
        "iterations",
        "initial_norm",
        "final_norm",
        "avg_reduction_factor",
        "prec_sym_err",
        "prec_stale_y_err",
        "peak_rss_kb",
        "error_return_code",
        "error_wall_time_sec",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    ok_mask = df["status"] == "ok"

    # Derived metrics for status == "ok" rows
    df["peak_rss_mb"] = np.where(ok_mask, df["peak_rss_kb"] / 1024.0, np.nan)

    df["time_per_free_dof"] = [
        safe_div(s, d) if ok and is_valid_number(d) and float(d) > 0 else float("nan")
        for ok, s, d in zip(ok_mask, df["solve"], df["num_free_dofs"], strict=True)
    ]

    df["setup_to_solve_ratio"] = [
        safe_div(setup, s) if ok else float("nan")
        for ok, setup, s in zip(
            ok_mask, df["prec_setup_total"], df["solve"], strict=True
        )
    ]

    df["speedup_vs_baseline"] = np.nan
    df["total_speedup_vs_baseline"] = np.nan

    for _, group in df.groupby(GROUP_KEYS, dropna=False, sort=False):
        base_rows = group[
            (group["prec"] == baseline_prec)
            & (group["status"] == "ok")
            & group["solve"].notna()
        ]
        if base_rows.empty:
            continue

        base_solve = float(base_rows["solve"].median())
        base_setup_vals = base_rows["prec_setup_total"].fillna(0.0)
        base_total = float((base_setup_vals + base_rows["solve"]).median())

        for idx, row in group.iterrows():
            if row["status"] != "ok":
                continue
            row_solve = row["solve"]
            if is_valid_number(base_solve) and base_solve > 0 and is_valid_number(row_solve):
                df.at[idx, "speedup_vs_baseline"] = safe_div(base_solve, row_solve)

            row_setup = (
                float(row["prec_setup_total"])
                if is_valid_number(row["prec_setup_total"])
                else 0.0
            )
            if is_valid_number(row_solve):
                row_total = row_setup + float(row_solve)
                if is_valid_number(base_total) and base_total > 0 and row_total > 0:
                    df.at[idx, "total_speedup_vs_baseline"] = safe_div(
                        base_total, row_total
                    )

    df = df.sort_values(by=SORT_KEYS, kind="mergesort", na_position="last").reset_index(
        drop=True
    )
    return df


def save_tabular_outputs(df: pd.DataFrame, result_dir: Path) -> list[Path]:
    """Write summary.csv and (if pyarrow is available) summary.parquet."""
    generated: list[Path] = []
    csv_path = result_dir / "summary.csv"
    df.to_csv(csv_path, index=False)
    generated.append(csv_path)

    try:
        import pyarrow  # noqa: F401

        parquet_path = result_dir / "summary.parquet"
        df.to_parquet(parquet_path, index=False, engine="pyarrow")
        generated.append(parquet_path)
    except ImportError:
        print(
            "[INFO] pyarrow not installed; skipping summary.parquet (summary.csv written).",
            file=sys.stderr,
        )
    except Exception as exc:
        print(
            f"[WARN] Failed to write summary.parquet ({exc}); summary.csv written.",
            file=sys.stderr,
        )

    return generated


# ---------------------------------------------------------------------------
# LaTeX Table Generation
# ---------------------------------------------------------------------------


def generate_summary_table_tex(
    df: pd.DataFrame, tables_dir: Path, baseline_prec: str
) -> Path:
    """Create booktabs summary table at <result_dir>/tables/summary_table.tex."""
    out_path = tables_dir / "summary_table.tex"

    lines: list[str] = [
        "% Auto-generated by benchmarks/analyze_results.py",
        "% Requires \\usepackage{booktabs}",
        "\\begin{table}[htbp]",
        "\\centering",
        "\\small",
        (
            "\\caption{Benchmark summary across meshes, refinement strategies, "
            f"orders, and preconditioners (speedup relative to \\texttt{{{latex_escape(baseline_prec)}}}).}}"
        ),
        "\\label{tab:benchmark_summary}",
        "\\begin{tabular}{llrrrllrrrrr}",
        "\\toprule",
        "Mesh & Strat & $p$ & Ref & DoFs & Solver & Prec & Iters & Setup (s) & Solve (s) & Speedup & RSS (MB) \\\\",
        "\\midrule",
    ]

    prev_group: tuple[Any, ...] | None = None
    for _, row in df.iterrows():
        cur_group = (
            row["mesh_name"],
            row["ref_strategy"],
            row["order"],
            row["ref_levels"],
        )
        if prev_group is not None and cur_group != prev_group:
            lines.append("\\addlinespace[0.4em]")
        prev_group = cur_group

        mesh_s = latex_escape(row["mesh_name"] if pd.notna(row["mesh_name"]) else "--")
        strat_s = latex_escape(
            row["ref_strategy"] if pd.notna(row["ref_strategy"]) else "--"
        )
        order_s = fmt_int_like(row["order"])
        ref_s = fmt_int_like(row["ref_levels"])
        dofs_s = fmt_int_like(row["num_free_dofs"], thousands=True)
        solver_s = latex_escape(row["solver"] if pd.notna(row["solver"]) else "--")
        prec_s = latex_escape(row["prec"] if pd.notna(row["prec"]) else "--")

        iters_tok = status_token(row["status"], row["converged"], row["iterations"])
        if iters_tok in {"TIMEOUT", "OOM", "CRASH", "NC"}:
            iters_s = f"\\texttt{{{iters_tok}}}"
        else:
            iters_s = iters_tok

        if row["status"] == "ok":
            setup_s = fmt_time(row["prec_setup_total"])
            solve_s = fmt_time(row["solve"])
            speedup_s = fmt_speedup(row["speedup_vs_baseline"], latex=True)
            rss_s = fmt_rss(row["peak_rss_mb"])
        else:
            setup_s = "--"
            solve_s = "--"
            speedup_s = "--"
            rss_s = "--"

        lines.append(
            f"{mesh_s} & {strat_s} & {order_s} & {ref_s} & {dofs_s} & "
            f"{solver_s} & {prec_s} & {iters_s} & {setup_s} & {solve_s} & "
            f"{speedup_s} & {rss_s} \\\\"
        )

    lines.extend(
        [
            "\\bottomrule",
            "\\end{tabular}",
            "\\end{table}",
            "",
        ]
    )

    out_path.write_text("\n".join(lines), encoding="utf-8")
    return out_path


def generate_iterations_pivot_tex(
    df: pd.DataFrame, tables_dir: Path, baseline_prec: str
) -> Path:
    """Create booktabs iteration-count pivot table at <result_dir>/tables/iterations_pivot.tex."""
    out_path = tables_dir / "iterations_pivot.tex"

    work = df.copy()
    # Backfill num_free_dofs within each problem configuration so failed runs
    # stay aligned on the same pivot row as successful sibling preconditioners.
    dof_lookup_keys = ["mesh_name", "ref_strategy", "order", "ref_levels"]
    work["_group_dofs"] = work.groupby(dof_lookup_keys, dropna=False)[
        "num_free_dofs"
    ].transform("max")

    precs = order_preconditioners(work["prec"].dropna().tolist(), baseline_prec)
    multi_solver = work["solver"].nunique(dropna=True) > 1

    row_keys = ["mesh_name", "ref_strategy", "order", "ref_levels"]
    if multi_solver:
        row_keys.append("solver")

    # Build lookup (row_key -> {prec: formatted_token}) and DOF map
    pivot_cells: dict[tuple[Any, ...], dict[str, str]] = {}
    pivot_dofs: dict[tuple[Any, ...], Any] = {}
    ordered_row_keys: list[tuple[Any, ...]] = []

    for _, row in work.iterrows():
        key = tuple(row[k] for k in row_keys)
        if key not in pivot_cells:
            pivot_cells[key] = {}
            pivot_dofs[key] = row["_group_dofs"]
            ordered_row_keys.append(key)
        prec = str(row["prec"])
        tok = status_token(row["status"], row["converged"], row["iterations"])
        if tok in {"TIMEOUT", "OOM", "CRASH", "NC"}:
            tok = f"\\texttt{{{tok}}}"
        pivot_cells[key][prec] = tok

    header_cols = ["Mesh", "Strat", "$p$", "Ref", "DoFs"]
    col_spec = "llrrr"
    if multi_solver:
        header_cols.append("Solver")
        col_spec += "l"
    header_cols.extend(latex_escape(p) for p in precs)
    col_spec += "r" * len(precs)

    lines: list[str] = [
        "% Auto-generated by benchmarks/analyze_results.py",
        "% Requires \\usepackage{booktabs}",
        "\\begin{table}[htbp]",
        "\\centering",
        "\\small",
        "\\caption{Krylov iteration counts by preconditioner across mesh, refinement strategy, order, and refinement level.}",
        "\\label{tab:iterations_pivot}",
        f"\\begin{{tabular}}{{{col_spec}}}",
        "\\toprule",
        " & ".join(header_cols) + " \\\\",
        "\\midrule",
    ]

    prev_mesh_strat: tuple[Any, Any] | None = None
    for key in ordered_row_keys:
        mesh_name, ref_strategy, order, ref_levels = key[:4]
        if prev_mesh_strat is not None and (mesh_name, ref_strategy) != prev_mesh_strat:
            lines.append("\\addlinespace[0.4em]")
        prev_mesh_strat = (mesh_name, ref_strategy)

        cells = [
            latex_escape(mesh_name if pd.notna(mesh_name) else "--"),
            latex_escape(ref_strategy if pd.notna(ref_strategy) else "--"),
            fmt_int_like(order),
            fmt_int_like(ref_levels),
            fmt_int_like(pivot_dofs[key], thousands=True),
        ]
        if multi_solver:
            solver_val = key[4]
            cells.append(latex_escape(solver_val if pd.notna(solver_val) else "--"))

        row_map = pivot_cells[key]
        for prec in precs:
            cells.append(row_map.get(prec, "--"))

        lines.append(" & ".join(cells) + " \\\\")

    lines.extend(
        [
            "\\bottomrule",
            "\\end{tabular}",
            "\\end{table}",
            "",
        ]
    )

    out_path.write_text("\n".join(lines), encoding="utf-8")
    return out_path


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def save_figure(fig: plt.Figure, base_path: Path) -> list[Path]:
    """Save a Matplotlib figure as both 300 DPI PNG and vector PDF."""
    png_path = base_path.with_suffix(".png")
    pdf_path = base_path.with_suffix(".pdf")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return [png_path, pdf_path]


def plot_iterations_vs_dofs(
    df: pd.DataFrame,
    plots_dir: Path,
    baseline_prec: str,
    style_map: dict[str, tuple[str, str]],
) -> list[Path]:
    """Plot iterations vs. free DoFs (log-x) per (mesh_name, ref_strategy)."""
    generated: list[Path] = []
    valid = df[
        (df["status"] == "ok")
        & df["num_free_dofs"].notna()
        & (df["num_free_dofs"] > 0)
        & df["iterations"].notna()
    ]
    if valid.empty:
        return generated

    for (mesh_name, ref_strategy), sub in valid.groupby(
        ["mesh_name", "ref_strategy"], dropna=False
    ):
        orders = sorted(o for o in sub["order"].dropna().unique())
        if not orders:
            continue

        solvers = sorted(s for s in sub["solver"].dropna().unique())
        multi_solver = len(solvers) > 1
        precs = order_preconditioners(sub["prec"].dropna().tolist(), baseline_prec)

        # 1-2 orders: single axes with hue=prec, linestyle=order
        # 3+ orders: one subplot per order p
        if len(orders) <= 2:
            fig, ax = plt.subplots(figsize=(6.8, 4.4))
            plotted_any = False

            for o_idx, order in enumerate(orders):
                o_sub = sub[sub["order"] == order]
                ls = LINESTYLES[o_idx % len(LINESTYLES)]
                for prec in precs:
                    p_sub = o_sub[o_sub["prec"] == prec]
                    for s_idx, solver in enumerate(solvers or [None]):
                        s_sub = (
                            p_sub[p_sub["solver"] == solver]
                            if solver is not None
                            else p_sub
                        )
                        s_sub = s_sub.sort_values("num_free_dofs")
                        if s_sub.empty:
                            continue

                        color, marker = style_map.get(prec, ("#333333", "o"))
                        label_parts = [prec]
                        if len(orders) == 2:
                            label_parts.append(f"p={fmt_int_like(order)}")
                        if multi_solver and solver is not None:
                            label_parts.append(str(solver))
                        label = ", ".join(label_parts)

                        cur_ls = (
                            LINESTYLES[(o_idx + s_idx) % len(LINESTYLES)]
                            if multi_solver
                            else ls
                        )
                        ax.plot(
                            s_sub["num_free_dofs"],
                            s_sub["iterations"],
                            color=color,
                            marker=marker,
                            linestyle=cur_ls,
                            linewidth=1.8,
                            markersize=6,
                            label=label,
                        )
                        # Highlight non-converged points with hollow markers
                        nc = s_sub[s_sub["converged"] == False]  # noqa: E712
                        if not nc.empty:
                            ax.scatter(
                                nc["num_free_dofs"],
                                nc["iterations"],
                                facecolors="white",
                                edgecolors=color,
                                marker=marker,
                                s=48,
                                linewidths=1.5,
                                zorder=5,
                            )
                        plotted_any = True

            if not plotted_any:
                plt.close(fig)
                continue

            ax.set_xscale("log")
            ax.set_xlabel("Free DoFs")
            ax.set_ylabel("Iterations")
            title_suffix = (
                f" (p={fmt_int_like(orders[0])})" if len(orders) == 1 else ""
            )
            ax.set_title(
                f"Iterations vs. Free DoFs — {mesh_name} ({ref_strategy}){title_suffix}"
            )
            ax.set_ylim(bottom=0)
            ax.legend(loc="best", ncol=1 if len(precs) <= 5 else 2)
        else:
            ncols = min(len(orders), 3)
            nrows = math.ceil(len(orders) / ncols)
            fig, axes = plt.subplots(
                nrows, ncols, figsize=(4.6 * ncols, 3.8 * nrows), squeeze=False
            )
            plotted_any = False

            for idx, order in enumerate(orders):
                ax = axes[idx // ncols][idx % ncols]
                o_sub = sub[sub["order"] == order]
                for prec in precs:
                    p_sub = o_sub[o_sub["prec"] == prec]
                    for s_idx, solver in enumerate(solvers or [None]):
                        s_sub = (
                            p_sub[p_sub["solver"] == solver]
                            if solver is not None
                            else p_sub
                        )
                        s_sub = s_sub.sort_values("num_free_dofs")
                        if s_sub.empty:
                            continue

                        color, marker = style_map.get(prec, ("#333333", "o"))
                        label = (
                            f"{prec} ({solver})"
                            if multi_solver and solver is not None
                            else prec
                        )
                        ax.plot(
                            s_sub["num_free_dofs"],
                            s_sub["iterations"],
                            color=color,
                            marker=marker,
                            linestyle=LINESTYLES[s_idx % len(LINESTYLES)],
                            linewidth=1.8,
                            markersize=5.5,
                            label=label,
                        )
                        nc = s_sub[s_sub["converged"] == False]  # noqa: E712
                        if not nc.empty:
                            ax.scatter(
                                nc["num_free_dofs"],
                                nc["iterations"],
                                facecolors="white",
                                edgecolors=color,
                                marker=marker,
                                s=42,
                                linewidths=1.5,
                                zorder=5,
                            )
                        plotted_any = True

                ax.set_xscale("log")
                ax.set_title(f"$p = {fmt_int_like(order)}$")
                ax.set_xlabel("Free DoFs")
                if idx % ncols == 0:
                    ax.set_ylabel("Iterations")
                ax.set_ylim(bottom=0)

            for empty_idx in range(len(orders), nrows * ncols):
                axes[empty_idx // ncols][empty_idx % ncols].set_visible(False)

            if not plotted_any:
                plt.close(fig)
                continue

            # Deduplicated shared legend at top
            handles, labels = [], []
            seen: set[str] = set()
            for ax_row in axes:
                for ax in ax_row:
                    h_list, l_list = ax.get_legend_handles_labels()
                    for h, lbl in zip(h_list, l_list, strict=True):
                        if lbl not in seen:
                            seen.add(lbl)
                            handles.append(h)
                            labels.append(lbl)

            fig.suptitle(
                f"Iterations vs. Free DoFs — {mesh_name} ({ref_strategy})",
                y=1.02,
                fontweight="bold",
            )
            if handles:
                fig.legend(
                    handles,
                    labels,
                    loc="upper center",
                    bbox_to_anchor=(0.5, 0.99),
                    ncol=min(len(labels), 5),
                )
            fig.tight_layout(rect=(0, 0, 1, 0.93))

        stem = f"iterations_vs_dofs_{slugify(mesh_name)}_{slugify(ref_strategy)}"
        generated.extend(save_figure(fig, plots_dir / stem))

    return generated


def plot_timing_breakdowns(
    df: pd.DataFrame, plots_dir: Path, baseline_prec: str
) -> list[Path]:
    """Grouped stacked bar charts across ref_levels and prec per (mesh, strat, order)."""
    generated: list[Path] = []
    valid = df[df["status"] == "ok"].copy()
    if valid.empty:
        return generated

    stage_cols = [
        "operator_assembly",
        "prec_univariate_assembly",
        "prec_eigendecomp",
        "prec_fd_construct",
        "solve",
    ]
    for col in stage_cols + ["prec_setup_total"]:
        valid[col] = valid[col].fillna(0.0).clip(lower=0.0)

    # Capture non-Kronecker preconditioner setup time if present (e.g., AMG, ILU)
    kronecker_sub = (
        valid["prec_univariate_assembly"]
        + valid["prec_eigendecomp"]
        + valid["prec_fd_construct"]
    )
    valid["prec_other_setup"] = (valid["prec_setup_total"] - kronecker_sub).clip(
        lower=0.0
    )

    for (mesh_name, ref_strategy, order), sub in valid.groupby(
        ["mesh_name", "ref_strategy", "order"], dropna=False
    ):
        if sub.empty or not is_valid_number(order):
            continue

        refs = sorted(r for r in sub["ref_levels"].dropna().unique())
        if not refs:
            continue

        solvers = sorted(s for s in sub["solver"].dropna().unique())
        multi_solver = len(solvers) > 1
        precs = order_preconditioners(sub["prec"].dropna().tolist(), baseline_prec)

        # Build ordered (prec, solver) bar series present in this subset
        series_keys: list[tuple[str, str | None]] = []
        for prec in precs:
            for solver in solvers or [None]:
                mask = sub["prec"] == prec
                if solver is not None:
                    mask &= sub["solver"] == solver
                if mask.any():
                    series_keys.append((prec, solver))

        if not series_keys:
            continue

        # Only include "prec_other_setup" in the stack/legend if non-negligible
        active_stages = [
            (col, label, color)
            for col, label, color in TIMING_STAGES
            if col != "prec_other_setup" or (sub["prec_other_setup"] > 1e-9).any()
        ]

        n_refs = len(refs)
        n_series = len(series_keys)
        group_width = 0.8
        bar_width = group_width / max(n_series, 1)

        fig_w = max(6.8, 1.5 * n_refs * max(1.0, n_series / 3.0))
        fig, ax = plt.subplots(figsize=(min(fig_w, 16.0), 4.8))

        x_centers = np.arange(n_refs, dtype=float)
        plotted_labels: set[str] = set()
        bar_positions: list[float] = []
        bar_totals: list[float] = []
        bar_names: list[str] = []

        for s_idx, (prec, solver) in enumerate(series_keys):
            offset = (s_idx - (n_series - 1) / 2.0) * bar_width
            bar_label_short = (
                f"{prec}\n({solver})" if multi_solver and solver is not None else prec
            )

            for r_idx, ref_lvl in enumerate(refs):
                mask = (sub["ref_levels"] == ref_lvl) & (sub["prec"] == prec)
                if solver is not None:
                    mask &= sub["solver"] == solver
                rows = sub[mask]
                if rows.empty:
                    continue

                # Use median if duplicate entries exist for the same key
                row_vals = rows[[col for col, _, _ in active_stages]].median()
                x_pos = x_centers[r_idx] + offset
                bottom = 0.0

                for col, stage_label, color in active_stages:
                    height = float(row_vals[col])
                    if height <= 0.0:
                        continue
                    lbl = stage_label if stage_label not in plotted_labels else "_nolegend_"
                    plotted_labels.add(stage_label)
                    ax.bar(
                        x_pos,
                        height,
                        width=bar_width * 0.88,
                        bottom=bottom,
                        color=color,
                        edgecolor="white",
                        linewidth=0.5,
                        label=lbl,
                    )
                    bottom += height

                bar_positions.append(x_pos)
                bar_totals.append(bottom)
                bar_names.append(bar_label_short)

        if not bar_positions:
            plt.close(fig)
            continue

        max_height = max(bar_totals) if bar_totals else 1.0
        for x_pos, total_h, bname in zip(
            bar_positions, bar_totals, bar_names, strict=True
        ):
            ax.text(
                x_pos,
                total_h + 0.015 * max_height,
                bname,
                ha="center",
                va="bottom",
                fontsize=7.5,
                rotation=45 if n_series > 3 else 0,
            )

        # Annotate ref level and DOF count on x-axis
        x_tick_labels: list[str] = []
        for ref_lvl in refs:
            dof_vals = sub[sub["ref_levels"] == ref_lvl]["num_free_dofs"].dropna()
            dof_txt = (
                f"\n({fmt_int_like(dof_vals.iloc[0], thousands=True)} DoFs)"
                if not dof_vals.empty
                else ""
            )
            x_tick_labels.append(f"ref={fmt_int_like(ref_lvl)}{dof_txt}")

        ax.set_xticks(x_centers)
        ax.set_xticklabels(x_tick_labels)
        ax.set_ylabel("Wall time (s)")
        ax.set_ylim(0, max_height * 1.22 if max_height > 0 else 1.0)
        ax.set_title(
            f"Timing Breakdown — {mesh_name} ({ref_strategy}), $p={fmt_int_like(order)}$"
        )
        ax.legend(loc="upper left", ncol=min(len(plotted_labels), 3))

        stem = (
            f"timing_breakdown_{slugify(mesh_name)}_{slugify(ref_strategy)}"
            f"_p{fmt_int_like(order)}"
        )
        generated.extend(save_figure(fig, plots_dir / stem))

    return generated


def plot_convergence_histories(
    df: pd.DataFrame,
    residual_histories: dict[str, list[float]],
    plots_dir: Path,
    baseline_prec: str,
    style_map: dict[str, tuple[str, str]],
) -> list[Path]:
    """Plot relative residual r_k / r_0 vs. iteration for the finest refinement level
    with non-empty residual_history per (mesh_name, ref_strategy, order).
    """
    generated: list[Path] = []
    if not residual_histories:
        return generated

    has_hist = df[
        (df["status"] == "ok") & df["run_id"].isin(residual_histories.keys())
    ].copy()
    if has_hist.empty:
        return generated

    for (mesh_name, ref_strategy, order), sub in has_hist.groupby(
        ["mesh_name", "ref_strategy", "order"], dropna=False
    ):
        valid_refs = sub["ref_levels"].dropna()
        if valid_refs.empty or not is_valid_number(order):
            continue

        finest_ref = valid_refs.max()
        finest_sub = sub[sub["ref_levels"] == finest_ref]
        if finest_sub.empty:
            continue

        solvers = sorted(s for s in finest_sub["solver"].dropna().unique())
        multi_solver = len(solvers) > 1
        precs = order_preconditioners(
            finest_sub["prec"].dropna().tolist(), baseline_prec
        )

        fig, ax = plt.subplots(figsize=(6.8, 4.4))
        plotted_any = False

        for prec in precs:
            p_sub = finest_sub[finest_sub["prec"] == prec]
            for s_idx, solver in enumerate(solvers or [None]):
                s_sub = (
                    p_sub[p_sub["solver"] == solver]
                    if solver is not None
                    else p_sub
                )
                if s_sub.empty:
                    continue

                row = s_sub.iloc[0]
                hist = residual_histories.get(str(row["run_id"]), [])
                if not hist:
                    continue

                r0 = hist[0]
                if not is_valid_number(r0) or r0 <= 0.0:
                    init_norm = row.get("initial_norm")
                    if is_valid_number(init_norm) and float(init_norm) > 0.0:
                        r0 = float(init_norm)
                    else:
                        continue

                rel = np.array([val / r0 for val in hist], dtype=float)
                iters = np.arange(len(rel), dtype=int)
                pos_mask = np.isfinite(rel) & (rel > 0.0)
                if not pos_mask.any():
                    continue

                color, marker = style_map.get(prec, ("#333333", "o"))
                label = (
                    f"{prec} ({solver})"
                    if multi_solver and solver is not None
                    else prec
                )
                markevery = max(1, int(pos_mask.sum()) // 12)
                ax.semilogy(
                    iters[pos_mask],
                    rel[pos_mask],
                    color=color,
                    marker=marker,
                    markevery=markevery,
                    linestyle=LINESTYLES[s_idx % len(LINESTYLES)],
                    linewidth=1.7,
                    markersize=5,
                    label=label,
                )
                plotted_any = True

        if not plotted_any:
            plt.close(fig)
            continue

        rtols = finest_sub["rtol"].dropna().unique()
        if len(rtols) == 1 and is_valid_number(rtols[0]) and float(rtols[0]) > 0:
            rtol_val = float(rtols[0])
            ax.axhline(
                rtol_val,
                color="#666666",
                linestyle=":",
                linewidth=1.1,
                label=f"rtol = {rtol_val:.0e}",
            )

        dof_vals = finest_sub["num_free_dofs"].dropna()
        dof_str = (
            f", {fmt_int_like(dof_vals.iloc[0], thousands=True)} DoFs"
            if not dof_vals.empty
            else ""
        )
        ax.set_xlabel("Iteration $k$")
        ax.set_ylabel("Relative residual $\\|r_k\\| / \\|r_0\\|$")
        ax.set_title(
            f"Convergence — {mesh_name} ({ref_strategy}), "
            f"$p={fmt_int_like(order)}$, ref={fmt_int_like(finest_ref)}{dof_str}"
        )
        ax.legend(loc="best")

        stem = (
            f"convergence_{slugify(mesh_name)}_{slugify(ref_strategy)}"
            f"_p{fmt_int_like(order)}_r{fmt_int_like(finest_ref)}"
        )
        generated.extend(save_figure(fig, plots_dir / stem))

    return generated


def generate_all_plots(
    df: pd.DataFrame,
    residual_histories: dict[str, list[float]],
    plots_dir: Path,
    baseline_prec: str,
) -> list[Path]:
    """Generate all applicable plots and return their file paths."""
    configure_matplotlib()
    plots_dir.mkdir(parents=True, exist_ok=True)
    style_map = build_prec_style_map(df["prec"].dropna().tolist(), baseline_prec)

    generated: list[Path] = []
    generated.extend(plot_iterations_vs_dofs(df, plots_dir, baseline_prec, style_map))
    generated.extend(plot_timing_breakdowns(df, plots_dir, baseline_prec))
    generated.extend(
        plot_convergence_histories(
            df, residual_histories, plots_dir, baseline_prec, style_map
        )
    )
    return generated


# ---------------------------------------------------------------------------
# Console Summary
# ---------------------------------------------------------------------------


def print_console_summary(
    df: pd.DataFrame, generated_paths: list[Path], result_dir: Path
) -> None:
    """Print a clean, readable summary table and all generated output paths."""
    headers = [
        "mesh",
        "p",
        "ref",
        "DoFs",
        "prec",
        "status",
        "iters",
        "setup_s",
        "solve_s",
        "speedup",
    ]
    right_aligned = {"p", "ref", "DoFs", "iters", "setup_s", "solve_s", "speedup"}

    rows: list[list[str]] = []
    for _, row in df.iterrows():
        status_str = str(row["status"])
        if status_str == "ok" and row["converged"] is False:
            status_str = "ok(NC)"

        iters_str = status_token(row["status"], row["converged"], row["iterations"])
        setup_str = fmt_time(row["prec_setup_total"]) if row["status"] == "ok" else "--"
        solve_str = fmt_time(row["solve"]) if row["status"] == "ok" else "--"
        speedup_str = (
            fmt_speedup(row["speedup_vs_baseline"], latex=False)
            if row["status"] == "ok"
            else "--"
        )

        rows.append(
            [
                str(row["mesh_name"]) if pd.notna(row["mesh_name"]) else "--",
                fmt_int_like(row["order"]),
                fmt_int_like(row["ref_levels"]),
                fmt_int_like(row["num_free_dofs"], thousands=True),
                str(row["prec"]) if pd.notna(row["prec"]) else "--",
                status_str,
                iters_str,
                setup_str,
                solve_str,
                speedup_str,
            ]
        )

    widths = [
        max(len(h), *(len(r[i]) for r in rows)) if rows else len(h)
        for i, h in enumerate(headers)
    ]

    def fmt_row(vals: list[str]) -> str:
        parts = []
        for h, val, w in zip(headers, vals, widths, strict=True):
            parts.append(val.rjust(w) if h in right_aligned else val.ljust(w))
        return "  ".join(parts)

    sep = "  ".join("-" * w for w in widths)
    print(f"\nBenchmark Summary ({len(df)} runs in {result_dir}):")
    print(fmt_row(headers))
    print(sep)
    for r in rows:
        print(fmt_row(r))

    print(f"\nGenerated outputs ({len(generated_paths)} files):")
    for p in generated_paths:
        print(f"  - {p}")


# ---------------------------------------------------------------------------
# CLI Entrypoint
# ---------------------------------------------------------------------------


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate, analyze, plot, and generate LaTeX tables from MFEM NURBS "
            "Kronecker-product preconditioner benchmark results."
        )
    )
    parser.add_argument(
        "result_dir",
        type=Path,
        help="Path to a benchmark result directory containing runs/*.json.",
    )
    parser.add_argument(
        "--baseline-prec",
        type=str,
        default="none",
        help='Reference preconditioner for speedup calculations (default: "none").',
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip plot generation and only produce CSV/Parquet and LaTeX tables.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    result_dir: Path = args.result_dir.resolve()
    runs_dir = result_dir / "runs"

    if not runs_dir.is_dir():
        print(f"[ERROR] Runs directory not found: {runs_dir}", file=sys.stderr)
        return 1

    records, residual_histories = load_run_records(runs_dir)
    if not records:
        print(f"[ERROR] No valid run JSON files found in: {runs_dir}", file=sys.stderr)
        return 1

    df = build_summary_dataframe(records, baseline_prec=args.baseline_prec)

    generated_paths: list[Path] = []
    generated_paths.extend(save_tabular_outputs(df, result_dir))

    tables_dir = result_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    generated_paths.append(
        generate_summary_table_tex(df, tables_dir, baseline_prec=args.baseline_prec)
    )
    generated_paths.append(
        generate_iterations_pivot_tex(df, tables_dir, baseline_prec=args.baseline_prec)
    )

    if not args.no_plots:
        plots_dir = result_dir / "plots"
        generated_paths.extend(
            generate_all_plots(
                df,
                residual_histories,
                plots_dir=plots_dir,
                baseline_prec=args.baseline_prec,
            )
        )

    print_console_summary(df, generated_paths, result_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
