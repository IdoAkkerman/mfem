"""Derived metric calculations, generic speedup analysis, and CSV export."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .ingest import (
    BenchmarkDataset,
    coerce_match_value,
    is_valid_number,
    series_matches_value,
    slugify,
)

_KNOWN_NUMERIC_COLS: tuple[str, ...] = (
    "dim",
    "base_order",
    "order",
    "num_elements",
    "num_tdofs",
    "num_free_dofs",
    "num_ess_dofs",
    "matrix_nnz",
    "solve_repeats_count",
    "solve_std_sec",
    "iterations",
    "initial_norm",
    "final_norm",
    "avg_reduction_factor",
    "residual_history_len",
    "prec_sym_err",
    "prec_stale_y_err",
    "peak_rss_kb",
    "error_return_code",
    "error_wall_time_sec",
)

_NAN_SENTINEL = object()


def parse_baseline_spec(spec: str) -> dict[str, str]:
    """Parse a 'key=val' or 'k1=v1,k2=v2' baseline spec string into a dict."""
    cleaned = spec.strip()
    if not cleaned:
        raise ValueError("Empty baseline specification.")

    result: dict[str, str] = {}
    for part in cleaned.split(","):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            raise ValueError(
                f"Invalid baseline specification {spec!r}: expected 'key=value' pairs."
            )
        k_raw, v_raw = part.split("=", 1)
        key = k_raw.strip()
        val = v_raw.strip()
        if not key or not val:
            raise ValueError(
                f"Invalid baseline specification {spec!r}: empty key or value in {part!r}."
            )
        result[key] = val

    if not result:
        raise ValueError(f"Invalid baseline specification {spec!r}: no pairs found.")
    return result


def _num_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col in df.columns:
        return pd.to_numeric(df[col], errors="coerce").astype(float)
    return pd.Series(np.nan, index=df.index, dtype=float)


def _coerce_numeric_columns(dataset: BenchmarkDataset) -> None:
    df = dataset.df
    if df.empty:
        return

    candidates: set[str] = set(_KNOWN_NUMERIC_COLS) | set(dataset.timing_cols)
    for col in df.columns:
        if col in candidates:
            if not pd.api.types.is_bool_dtype(df[col]):
                df[col] = pd.to_numeric(df[col], errors="coerce")
            continue

        if pd.api.types.is_object_dtype(df[col]):
            non_null = df[col].dropna()
            if len(non_null) > 0 and all(is_valid_number(v) for v in non_null.tolist()):
                df[col] = pd.to_numeric(df[col], errors="coerce")


def compute_derived_metrics(dataset: BenchmarkDataset) -> BenchmarkDataset:
    """Coerce numeric columns and compute derived metrics on status=='ok' rows."""
    _coerce_numeric_columns(dataset)
    df = dataset.df

    if "status" in df.columns:
        ok_mask = (df["status"] == "ok").fillna(False)
    else:
        ok_mask = pd.Series(False, index=df.index, dtype=bool)

    peak_rss_kb = _num_series(df, "peak_rss_kb")
    num_free_dofs = _num_series(df, "num_free_dofs")
    op_asm = _num_series(df, "operator_assembly")
    prec_setup = _num_series(df, "prec_setup_total")
    prec_univ = _num_series(df, "prec_univariate_assembly")
    prec_eig = _num_series(df, "prec_eigendecomp")
    prec_fd = _num_series(df, "prec_fd_construct")
    solve = _num_series(df, "solve")

    df["peak_rss_mb"] = np.where(
        ok_mask & peak_rss_kb.notna(),
        peak_rss_kb / 1024.0,
        np.nan,
    )

    valid_dofs = ok_mask & num_free_dofs.notna() & (num_free_dofs > 0.0)
    df["time_per_free_dof"] = np.where(
        valid_dofs & solve.notna(),
        solve / num_free_dofs,
        np.nan,
    )
    df["assembly_time_per_free_dof"] = np.where(
        valid_dofs & op_asm.notna(),
        op_asm / num_free_dofs,
        np.nan,
    )

    df["setup_to_solve_ratio"] = np.where(
        ok_mask & prec_setup.notna() & solve.notna() & (solve > 0.0),
        prec_setup / solve,
        np.nan,
    )

    df["pipeline_time_sec"] = np.where(
        ok_mask & op_asm.notna() & solve.notna(),
        op_asm + prec_setup.fillna(0.0) + solve,
        np.nan,
    )

    sub_setup = prec_univ.fillna(0.0) + prec_eig.fillna(0.0) + prec_fd.fillna(0.0)
    df["prec_other_setup"] = np.where(
        ok_mask & prec_setup.notna(),
        np.clip(prec_setup - sub_setup, a_min=0.0, a_max=None),
        np.nan,
    )

    df["prec_and_solve_sec"] = np.where(
        ok_mask & solve.notna(),
        prec_setup.fillna(0.0) + solve,
        np.nan,
    )

    return dataset


def _normalize_key_elem(val: Any) -> Any:
    if val is None:
        return _NAN_SENTINEL
    if isinstance(val, (bool, np.bool_)):
        return ("bool", bool(val))
    if is_valid_number(val):
        return ("num", round(float(val), 12))
    try:
        if pd.isna(val):
            return _NAN_SENTINEL
    except (TypeError, ValueError):
        pass
    return ("str", str(val))


def _row_group_key(df: pd.DataFrame, idx: Any, keys: list[str]) -> tuple[Any, ...]:
    return tuple(_normalize_key_elem(df.at[idx, k]) for k in keys)


def _evaluate_numerator(
    df: pd.DataFrame,
    metric: str,
    numerator_expr: str | None,
) -> pd.Series:
    if numerator_expr is None:
        return _num_series(df, metric)
    expr = numerator_expr.strip()
    if expr in df.columns:
        return _num_series(df, expr)
    evaluated = df.eval(expr)
    return pd.to_numeric(evaluated, errors="coerce").astype(float)


def compute_speedup(
    dataset: BenchmarkDataset,
    baseline: dict[str, Any],
    metric: str,
    out_col: str,
    *,
    fallback_drop_keys: tuple[str, ...] = (),
    numerator_expr: str | None = None,
) -> BenchmarkDataset:
    """Compute relative speedup (baseline_median / row_val) into *out_col*."""
    df = dataset.df
    if df.empty or metric not in df.columns:
        df[out_col] = pd.Series(np.nan, index=df.index, dtype=float)
        return dataset

    if any(k not in df.columns for k in baseline):
        df[out_col] = pd.Series(np.nan, index=df.index, dtype=float)
        return dataset

    ok_mask = (
        (df["status"] == "ok").fillna(False)
        if "status" in df.columns
        else pd.Series(False, index=df.index, dtype=bool)
    )

    metric_vals = _num_series(df, metric)
    num_vals = _evaluate_numerator(df, metric, numerator_expr)

    base_mask = ok_mask.copy()
    for col, raw_target in baseline.items():
        coerced = (
            coerce_match_value(raw_target, df[col])
            if isinstance(raw_target, str)
            else raw_target
        )
        base_mask &= series_matches_value(df[col], coerced)

    valid_base_mask = (
        base_mask
        & num_vals.notna()
        & pd.Series(np.isfinite(num_vals), index=df.index)
        & (num_vals > 0.0)
    )

    group_keys = [c for c in dataset.varying_axes if c not in baseline and c in df.columns]
    fallback_keys = [c for c in group_keys if c not in fallback_drop_keys]
    use_fallback = bool(fallback_drop_keys) and (fallback_keys != group_keys)

    primary_buckets: dict[tuple[Any, ...], list[float]] = {}
    fallback_buckets: dict[tuple[Any, ...], list[float]] = {}

    for idx in df.index[valid_base_mask]:
        val = float(num_vals.at[idx])
        pkey = _row_group_key(df, idx, group_keys)
        primary_buckets.setdefault(pkey, []).append(val)
        if use_fallback:
            fkey = _row_group_key(df, idx, fallback_keys)
            fallback_buckets.setdefault(fkey, []).append(val)

    primary_medians = {k: float(np.median(v)) for k, v in primary_buckets.items()}
    fallback_medians = {k: float(np.median(v)) for k, v in fallback_buckets.items()}

    out_vals = np.full(len(df), np.nan, dtype=float)
    for pos, idx in enumerate(df.index):
        if not bool(ok_mask.at[idx]):
            continue
        row_val = metric_vals.at[idx]
        if not is_valid_number(row_val) or float(row_val) <= 0.0:
            continue

        pkey = _row_group_key(df, idx, group_keys)
        base_val = primary_medians.get(pkey)
        if base_val is None and use_fallback:
            fkey = _row_group_key(df, idx, fallback_keys)
            base_val = fallback_medians.get(fkey)

        if base_val is not None and math.isfinite(base_val) and base_val > 0.0:
            out_vals[pos] = base_val / float(row_val)

    df[out_col] = out_vals
    return dataset


def _column_has_value(df: pd.DataFrame, col: str, value: Any) -> bool:
    if col not in df.columns or df.empty:
        return False
    coerced = coerce_match_value(str(value), df[col]) if isinstance(value, str) else value
    return bool(series_matches_value(df[col], coerced).any())


def _compact_key_name(key: str) -> str:
    replacements = (
        ("wq_interior_pts", "wq_pts"),
        ("wq_symmetrize", "wq_sym"),
        ("gmres_restart", "gmres_k"),
        ("ref_strategy", "ref"),
        ("ref_levels", "lvl"),
    )
    for full, short in replacements:
        if key == full:
            return short
    return key


def _compact_baseline_label(spec: dict[str, str], kind: str) -> str:
    if len(spec) == 1:
        only_k, only_v = next(iter(spec.items()))
        if kind == "solve" and only_k == "prec":
            return f"solve_x({only_v})"
        if kind == "asm" and only_k == "assembly":
            return f"asm_x({only_v})"
        short_k = _compact_key_name(only_k)
        prefix = "vs_" if kind == "solve" else "asm_vs_"
        return f"{prefix}{short_k}={only_v}"

    joined = ",".join(f"{_compact_key_name(k)}={v}" for k, v in spec.items())
    prefix = "vs_" if kind == "solve" else "asm_vs_"
    return f"{prefix}{joined}"


def _append_speedup_col(
    dataset: BenchmarkDataset,
    col_name: str,
    label: str,
) -> None:
    if all(existing_col != col_name for existing_col, _ in dataset.speedup_cols):
        dataset.speedup_cols.append((col_name, label))


def _sort_dataset(dataset: BenchmarkDataset) -> None:
    df = dataset.df
    if df.empty:
        return
    sort_cols = [c for c in [*dataset.varying_axes, "run_id"] if c in df.columns]
    if not sort_cols:
        return
    try:
        dataset.df = df.sort_values(
            by=sort_cols,
            kind="mergesort",
            na_position="last",
        ).reset_index(drop=True)
    except TypeError:
        dataset.df = df.sort_values(
            by=sort_cols,
            kind="mergesort",
            na_position="last",
            key=lambda s: s.map(
                lambda x: (1, "") if pd.isna(x) else (0, str(x))
            ),
        ).reset_index(drop=True)


def apply_speedups(
    dataset: BenchmarkDataset,
    baseline_specs: list[dict[str, str]] | None = None,
    default_prec: str = "none",
    default_assembly: str = "legacy",
) -> BenchmarkDataset:
    """Compute default and user-specified speedup columns and sort *dataset*."""
    df = dataset.df
    dataset.speedup_cols = []

    has_default_prec = _column_has_value(df, "prec", default_prec)
    if has_default_prec:
        compute_speedup(
            dataset,
            baseline={"prec": default_prec},
            metric="solve",
            out_col="speedup_vs_baseline",
        )
        _append_speedup_col(dataset, "speedup_vs_baseline", "solve_x")

        compute_speedup(
            dataset,
            baseline={"prec": default_prec},
            metric="prec_and_solve_sec",
            out_col="total_speedup_vs_baseline",
        )
        _append_speedup_col(dataset, "total_speedup_vs_baseline", "tot_x")

    has_varying_asm = "assembly" in dataset.varying_axes
    has_default_asm = has_varying_asm and _column_has_value(
        df, "assembly", default_assembly
    )
    if has_default_asm:
        compute_speedup(
            dataset,
            baseline={"assembly": default_assembly},
            metric="operator_assembly",
            out_col="assembly_speedup_vs_baseline",
            fallback_drop_keys=("prec", "solver"),
        )
        _append_speedup_col(dataset, "assembly_speedup_vs_baseline", "asm_x")

        if has_default_prec:
            compute_speedup(
                dataset,
                baseline={"assembly": default_assembly, "prec": default_prec},
                metric="pipeline_time_sec",
                out_col="pipeline_speedup_vs_baseline",
            )
            _append_speedup_col(dataset, "pipeline_speedup_vs_baseline", "pipe_x")

    if baseline_specs:
        for spec in baseline_specs:
            if not spec:
                raise ValueError("Baseline specification dictionary must not be empty.")
            missing = [k for k in spec if k not in df.columns]
            if missing:
                avail = ", ".join(sorted(str(c) for c in df.columns))
                raise ValueError(
                    f"Baseline column(s) {missing} not found in dataset. "
                    f"Available: [{avail}]"
                )

            slug = slugify("_".join(f"{k}_{v}" for k, v in spec.items()))
            solve_col = f"speedup_solve_vs_{slug}"
            solve_label = _compact_baseline_label(spec, kind="solve")
            compute_speedup(
                dataset,
                baseline=spec,
                metric="solve",
                out_col=solve_col,
            )
            _append_speedup_col(dataset, solve_col, solve_label)

            is_asm_spec = all(
                k == "assembly" or k.startswith("wq_") or k.startswith("op_")
                for k in spec
            )
            if is_asm_spec:
                asm_col = f"speedup_asm_vs_{slug}"
                asm_label = _compact_baseline_label(spec, kind="asm")
                compute_speedup(
                    dataset,
                    baseline=spec,
                    metric="operator_assembly",
                    out_col=asm_col,
                    fallback_drop_keys=("prec", "solver"),
                )
                _append_speedup_col(dataset, asm_col, asm_label)

    _sort_dataset(dataset)
    return dataset


def save_summary_csv(dataset: BenchmarkDataset) -> Path:
    """Write dataset.df to '<result_dir>/summary.csv' and return its Path."""
    out_path = dataset.result_dir / "summary.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    dataset.df.to_csv(out_path, index=False)
    return out_path
