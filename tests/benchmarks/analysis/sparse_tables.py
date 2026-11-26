"""Publication-grade sparse aggregate table renderer (Matplotlib booktabs + Markdown)."""

from __future__ import annotations

import math
import textwrap
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

from .ingest import (
    BenchmarkDataset,
    fmt_int_like,
    fmt_speedup,
    is_valid_number,
    slugify,
)

VALID_VALUE_FMTS: frozenset[str] = frozenset(
    {
        "auto",
        "int",
        "time",
        "speedup",
        "float",
        "iters_and_time",
        "time_and_speedup",
    }
)
VALID_AGGS: frozenset[str] = frozenset({"median", "mean", "min", "max"})
VALID_HIGHLIGHTS: frozenset[str] = frozenset({"min", "max"})

_PREC_ORDER: dict[str, int] = {
    "none": 0,
    "jacobi": 1,
    "gs": 2,
    "kronecker-fd": 3,
}
_ASSEMBLY_ORDER: dict[str, int] = {
    "legacy": 0,
    "wq": 1,
}
_AXIS_DISPLAY_NAMES: dict[str, str] = {
    "mesh_name": "Mesh",
    "mesh": "Mesh",
    "problem": "Problem",
    "ref_strategy": "Strat",
    "order": "p",
    "ref_levels": "Ref",
    "assembly": "Asm",
    "solver": "Solver",
    "prec": "Prec",
    "wq_interior_pts": "WQ Pts",
    "wq_symmetrize": "WQ Sym",
    "gmres_restart": "Restart",
    "rtol": "rtol",
    "atol": "atol",
    "max_iter": "MaxIter",
}
_INT_AUTO_COLS: frozenset[str] = frozenset(
    {"iterations", "num_free_dofs", "ref_levels"}
)
_STATUS_LABELS: dict[str, str] = {
    "timeout": "TIMEOUT",
    "oom": "OOM",
    "crash": "CRASH",
}


@dataclass
class SparseTableSpec:
    name: str
    row_axes: list[str]
    col_axis: str
    value: str
    secondary_value: str | None = None
    value_fmt: str = "auto"  # "auto" | "int" | "time" | "speedup" | "float" | "iters_and_time" | "time_and_speedup"
    title: str = ""
    missing: str = "—"
    agg: str = "median"  # "median" | "mean" | "min" | "max"
    row_group_axes: list[str] = field(default_factory=list)
    highlight_best: str | None = None  # "min" | "max" | None
    col_header_prefix: str = ""  # e.g. "p=" -> "p=2", "p=3"
    split_by: list[str] = field(default_factory=list)


@dataclass
class _CellData:
    text: str
    kind: str  # "ok" | "missing" | "fail" | "nc"
    num_val: float | None = None
    is_best: bool = False


def _coerce_str_list(raw_val: Any, field_name: str, *, required: bool = False) -> list[str]:
    if raw_val is None:
        if required:
            raise ValueError(f"SparseTableSpec missing required field {field_name!r}.")
        return []
    if isinstance(raw_val, str):
        cleaned = raw_val.strip()
        if not cleaned:
            if required:
                raise ValueError(f"SparseTableSpec field {field_name!r} must not be empty.")
            return []
        return [cleaned]
    if isinstance(raw_val, (list, tuple)):
        out: list[str] = []
        for item in raw_val:
            if not isinstance(item, str) or not item.strip():
                raise ValueError(
                    f"SparseTableSpec field {field_name!r} must contain non-empty strings."
                )
            out.append(item.strip())
        if required and not out:
            raise ValueError(f"SparseTableSpec field {field_name!r} must not be empty.")
        return out
    raise ValueError(
        f"SparseTableSpec field {field_name!r} must be a string or list of strings."
    )


def _require_nonempty_str(raw: dict[str, Any], field_name: str) -> str:
    val = raw.get(field_name)
    if not isinstance(val, str) or not val.strip():
        raise ValueError(
            f"SparseTableSpec missing or invalid required string field {field_name!r}."
        )
    return val.strip()


def parse_sparse_table_spec(raw: dict[str, Any]) -> SparseTableSpec:
    """Validate and construct a SparseTableSpec from a TOML/dict specification."""
    if not isinstance(raw, dict):
        raise ValueError("SparseTableSpec raw specification must be a dict.")

    name = _require_nonempty_str(raw, "name")
    row_axes = _coerce_str_list(raw.get("row_axes"), "row_axes", required=True)
    col_axis = _require_nonempty_str(raw, "col_axis")
    value = _require_nonempty_str(raw, "value")

    raw_fmt = raw.get("value_fmt", "auto")
    if not isinstance(raw_fmt, str) or raw_fmt.strip() not in VALID_VALUE_FMTS:
        raise ValueError(
            f"Invalid value_fmt {raw_fmt!r}; expected one of {sorted(VALID_VALUE_FMTS)}."
        )
    value_fmt = raw_fmt.strip()

    raw_secondary = raw.get("secondary_value")
    if raw_secondary is None or (isinstance(raw_secondary, str) and not raw_secondary.strip()):
        secondary_value: str | None = None
    elif isinstance(raw_secondary, str):
        secondary_value = raw_secondary.strip()
    else:
        raise ValueError("SparseTableSpec field 'secondary_value' must be a string or None.")

    if value_fmt == "iters_and_time" and secondary_value is None:
        secondary_value = "solve"
    elif value_fmt == "time_and_speedup" and secondary_value is None:
        secondary_value = "speedup_vs_baseline"

    raw_agg = raw.get("agg", "median")
    if not isinstance(raw_agg, str) or raw_agg.strip() not in VALID_AGGS:
        raise ValueError(
            f"Invalid agg {raw_agg!r}; expected one of {sorted(VALID_AGGS)}."
        )
    agg = raw_agg.strip()

    raw_hl = raw.get("highlight_best")
    if raw_hl is None:
        highlight_best: str | None = None
    elif isinstance(raw_hl, str) and raw_hl.strip() in VALID_HIGHLIGHTS:
        highlight_best = raw_hl.strip()
    else:
        raise ValueError(
            f"Invalid highlight_best {raw_hl!r}; expected 'min', 'max', or None."
        )

    title = str(raw.get("title", "") or "")
    missing = str(raw.get("missing", "—"))
    col_header_prefix = str(raw.get("col_header_prefix", "") or "")
    row_group_axes = _coerce_str_list(raw.get("row_group_axes"), "row_group_axes")
    split_by = _coerce_str_list(raw.get("split_by"), "split_by")

    return SparseTableSpec(
        name=name,
        row_axes=row_axes,
        col_axis=col_axis,
        value=value,
        secondary_value=secondary_value,
        value_fmt=value_fmt,
        title=title,
        missing=missing,
        agg=agg,
        row_group_axes=row_group_axes,
        highlight_best=highlight_best,
        col_header_prefix=col_header_prefix,
        split_by=split_by,
    )


def _is_na(val: Any) -> bool:
    if val is None:
        return True
    if isinstance(val, (bool, np.bool_)):
        return False
    if isinstance(val, (float, np.floating)):
        return not math.isfinite(float(val))
    try:
        return bool(pd.isna(val))
    except (TypeError, ValueError):
        return False


def _norm_key(val: Any) -> tuple[str, Any]:
    if _is_na(val):
        return ("na", None)
    if isinstance(val, (bool, np.bool_)):
        return ("bool", bool(val))
    if is_valid_number(val):
        return ("num", round(float(val), 12))
    return ("str", str(val))


def _axis_sort_key(axis: str, val: Any) -> tuple[int, float, str]:
    if _is_na(val):
        return (9, 0.0, "")
    if axis == "prec":
        low = str(val).strip().lower()
        if low in _PREC_ORDER:
            return (0, float(_PREC_ORDER[low]), low)
        return (1, 0.0, low)
    if axis == "assembly":
        low = str(val).strip().lower()
        if low in _ASSEMBLY_ORDER:
            return (0, float(_ASSEMBLY_ORDER[low]), low)
        return (1, 0.0, low)
    if isinstance(val, (bool, np.bool_)):
        return (0, float(int(bool(val))), "")
    if is_valid_number(val):
        return (0, float(val), "")
    return (1, 0.0, str(val).strip().lower())


def _axis_display_name(axis: str) -> str:
    if axis in _AXIS_DISPLAY_NAMES:
        return _AXIS_DISPLAY_NAMES[axis]
    return axis.replace("_", " ").title()


def _fmt_axis_val(val: Any, missing: str = "—") -> str:
    if _is_na(val):
        return missing
    if isinstance(val, (bool, np.bool_)):
        return "true" if bool(val) else "false"
    if is_valid_number(val):
        return fmt_int_like(val, thousands=False)
    return str(val)


def _fmt_compact_time(val: float) -> str:
    aval = abs(val)
    if aval < 0.01:
        return f"{val:.4f}s"
    if aval < 10.0:
        return f"{val:.3f}s"
    return f"{val:.1f}s"


def _resolve_value_fmt(
    value_fmt: str,
    value_col: str,
    timing_cols: list[str],
) -> str:
    if value_fmt != "auto":
        return value_fmt
    if value_col in _INT_AUTO_COLS:
        return "int"
    if "speedup" in value_col:
        return "speedup"
    if value_col in timing_cols or value_col.endswith("_sec"):
        return "time"
    return "float"


def _format_cell_value(
    primary: float,
    secondary: float | None,
    resolved_fmt: str,
    value_col: str,
) -> str:
    if resolved_fmt == "int":
        return fmt_int_like(primary, thousands=(value_col == "num_free_dofs"))
    if resolved_fmt == "time":
        return _fmt_compact_time(primary)
    if resolved_fmt == "speedup":
        return fmt_speedup(primary)
    if resolved_fmt == "iters_and_time":
        p_str = fmt_int_like(primary, thousands=False)
        if secondary is not None and is_valid_number(secondary):
            return f"{p_str} ({_fmt_compact_time(secondary)})"
        return p_str
    if resolved_fmt == "time_and_speedup":
        p_str = _fmt_compact_time(primary)
        if secondary is not None and is_valid_number(secondary):
            return f"{p_str} ({fmt_speedup(secondary)})"
        return p_str
    return f"{primary:.3g}"


def _aggregate(values: list[float], agg: str) -> float:
    arr = np.asarray(values, dtype=float)
    if agg == "median":
        return float(np.median(arr))
    if agg == "mean":
        return float(np.mean(arr))
    if agg == "min":
        return float(np.min(arr))
    if agg == "max":
        return float(np.max(arr))
    raise ValueError(f"Unsupported aggregation: {agg!r}")


def _is_explicit_false(val: Any) -> bool:
    return isinstance(val, (bool, np.bool_)) and (not bool(val))


def _format_fail_status(statuses: list[str]) -> str:
    counts = Counter(str(s).strip().lower() for s in statuses if str(s).strip())
    if not counts:
        return "FAIL"
    top_status = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
    return _STATUS_LABELS.get(top_status, top_status.upper())


def _extract_finite_floats(df: pd.DataFrame, col: str) -> list[float]:
    if col not in df.columns or df.empty:
        return []
    out: list[float] = []
    for v in df[col].tolist():
        if is_valid_number(v):
            out.append(float(v))
    return out


def _build_cell(
    cell_df: pd.DataFrame,
    spec: SparseTableSpec,
    resolved_fmt: str,
) -> _CellData:
    if cell_df.empty:
        return _CellData(text=spec.missing, kind="missing")

    if "status" in cell_df.columns:
        ok_df = cell_df[cell_df["status"].astype(str) == "ok"]
    else:
        ok_df = cell_df

    if ok_df.empty:
        fail_statuses = (
            cell_df["status"].astype(str).tolist()
            if "status" in cell_df.columns
            else ["fail"]
        )
        return _CellData(text=_format_fail_status(fail_statuses), kind="fail")

    if (
        spec.value == "iterations"
        and "converged" in ok_df.columns
        and all(_is_explicit_false(v) for v in ok_df["converged"].tolist())
    ):
        return _CellData(text="NC", kind="nc")

    primary_vals = _extract_finite_floats(ok_df, spec.value)
    if not primary_vals:
        return _CellData(text=spec.missing, kind="missing")

    primary_agg = _aggregate(primary_vals, spec.agg)
    secondary_agg: float | None = None
    if spec.secondary_value is not None and spec.secondary_value in ok_df.columns:
        sec_vals = _extract_finite_floats(ok_df, spec.secondary_value)
        if sec_vals:
            secondary_agg = _aggregate(sec_vals, spec.agg)

    text = _format_cell_value(primary_agg, secondary_agg, resolved_fmt, spec.value)
    return _CellData(text=text, kind="ok", num_val=primary_agg)


def _iter_splits(
    df: pd.DataFrame,
    split_cols: list[str],
) -> list[tuple[list[tuple[str, Any]], pd.DataFrame]]:
    if not split_cols:
        return [([], df)]

    buckets: dict[tuple[tuple[str, Any], ...], list[int]] = {}
    representatives: dict[tuple[tuple[str, Any], ...], list[Any]] = {}

    for pos, idx in enumerate(df.index):
        raw_vals = [df.at[idx, c] for c in split_cols]
        if any(_is_na(v) for v in raw_vals):
            continue
        key = tuple(_norm_key(v) for v in raw_vals)
        if key not in buckets:
            buckets[key] = []
            representatives[key] = raw_vals
        buckets[key].append(pos)

    ordered_keys = sorted(
        buckets.keys(),
        key=lambda k: tuple(
            _axis_sort_key(c, v) for c, v in zip(split_cols, representatives[k])
        ),
    )

    result: list[tuple[list[tuple[str, Any]], pd.DataFrame]] = []
    for k in ordered_keys:
        pairs = list(zip(split_cols, representatives[k]))
        sub_df = df.iloc[buckets[k]].reset_index(drop=True)
        result.append((pairs, sub_df))
    return result


def _assemble_table_grid(
    sub_df: pd.DataFrame,
    spec: SparseTableSpec,
    timing_cols: list[str],
) -> tuple[
    list[str],
    list[str],
    list[list[str]],
    list[list[str]],
    list[list[_CellData]],
    list[int],
] | None:
    """Build headers, row labels (full and group-suppressed), cells, and group ids."""
    col_rep: dict[tuple[str, Any], Any] = {}
    row_rep: dict[tuple[tuple[str, Any], ...], list[Any]] = {}
    cell_indices: dict[
        tuple[tuple[tuple[str, Any], ...], tuple[str, Any]], list[int]
    ] = {}

    for pos, idx in enumerate(sub_df.index):
        c_raw = sub_df.at[idx, spec.col_axis]
        if _is_na(c_raw):
            continue
        c_key = _norm_key(c_raw)
        if c_key not in col_rep:
            col_rep[c_key] = c_raw

        r_raw = [sub_df.at[idx, a] for a in spec.row_axes]
        r_key = tuple(_norm_key(v) for v in r_raw)
        if r_key not in row_rep:
            row_rep[r_key] = r_raw

        cell_indices.setdefault((r_key, c_key), []).append(pos)

    if not col_rep or not row_rep:
        return None

    sorted_col_keys = sorted(
        col_rep.keys(),
        key=lambda ck: _axis_sort_key(spec.col_axis, col_rep[ck]),
    )
    sorted_row_keys = sorted(
        row_rep.keys(),
        key=lambda rk: tuple(
            _axis_sort_key(a, v) for a, v in zip(spec.row_axes, row_rep[rk])
        ),
    )

    row_headers = [_axis_display_name(a) for a in spec.row_axes]
    col_headers = [
        f"{spec.col_header_prefix}{_fmt_axis_val(col_rep[ck], spec.missing)}"
        for ck in sorted_col_keys
    ]

    resolved_fmt = _resolve_value_fmt(spec.value_fmt, spec.value, timing_cols)

    cells: list[list[_CellData]] = []
    for rk in sorted_row_keys:
        row_cells: list[_CellData] = []
        for ck in sorted_col_keys:
            idxs = cell_indices.get((rk, ck))
            if not idxs:
                row_cells.append(_CellData(text=spec.missing, kind="missing"))
            else:
                cell_df = sub_df.iloc[idxs]
                row_cells.append(_build_cell(cell_df, spec, resolved_fmt))
        cells.append(row_cells)

    # Determine row grouping axes
    if spec.row_group_axes:
        group_axis_indices = [
            i for i, a in enumerate(spec.row_axes) if a in spec.row_group_axes
        ]
    elif len(spec.row_axes) >= 2:
        group_axis_indices = list(range(len(spec.row_axes) - 1))
    else:
        group_axis_indices = []

    group_ids: list[int] = []
    group_key_to_rows: dict[int, list[int]] = {}
    prev_gkey: tuple[tuple[str, Any], ...] | None = None
    cur_gid = -1
    for r_idx, rk in enumerate(sorted_row_keys):
        gkey = tuple(rk[i] for i in group_axis_indices)
        if gkey != prev_gkey:
            cur_gid += 1
            prev_gkey = gkey
        group_ids.append(cur_gid)
        group_key_to_rows.setdefault(cur_gid, []).append(r_idx)

    # Highlight best cell per (row_group, col) when >= 2 numeric values exist
    if spec.highlight_best in VALID_HIGHLIGHTS:
        n_cols = len(sorted_col_keys)
        for gid_rows in group_key_to_rows.values():
            for c_idx in range(n_cols):
                candidates = [
                    (r_idx, cells[r_idx][c_idx].num_val)
                    for r_idx in gid_rows
                    if cells[r_idx][c_idx].kind == "ok"
                    and cells[r_idx][c_idx].num_val is not None
                    and math.isfinite(cells[r_idx][c_idx].num_val)
                ]
                if len(candidates) < 2:
                    continue
                vals = [v for _, v in candidates if v is not None]
                best_val = min(vals) if spec.highlight_best == "min" else max(vals)
                for r_idx, val in candidates:
                    if val is not None and math.isclose(
                        val, best_val, rel_tol=1e-9, abs_tol=1e-12
                    ):
                        cells[r_idx][c_idx].is_best = True

    full_row_labels: list[list[str]] = [
        [_fmt_axis_val(v, spec.missing) for v in row_rep[rk]]
        for rk in sorted_row_keys
    ]

    display_row_labels: list[list[str]] = []
    group_idx_set = set(group_axis_indices)
    for r_idx, labels in enumerate(full_row_labels):
        if r_idx > 0 and group_ids[r_idx] == group_ids[r_idx - 1]:
            suppressed = [
                "" if col_i in group_idx_set else lbl
                for col_i, lbl in enumerate(labels)
            ]
            display_row_labels.append(suppressed)
        else:
            display_row_labels.append(list(labels))

    return (
        row_headers,
        col_headers,
        full_row_labels,
        display_row_labels,
        cells,
        group_ids,
    )


def _estimate_text_width_in(text: str, fontsize_pt: float, *, bold: bool = False) -> float:
    if not text:
        return 0.0
    em_in = fontsize_pt / 72.0
    units = 0.0
    for ch in text:
        if ch in "ilI1.,:;|!()[]'` ":
            units += 0.44
        elif ch in "MWmw@%—–":
            units += 0.94
        elif ch.isupper():
            units += 0.76
        else:
            units += 0.66
    weight_factor = 1.07 if bold else 1.0
    return units * em_in * weight_factor * 1.12


def _render_matplotlib_table(
    title: str,
    row_headers: list[str],
    col_headers: list[str],
    display_row_labels: list[list[str]],
    cells: list[list[_CellData]],
    group_ids: list[int],
    out_dir: Path,
    stem: str,
    formats: list[str],
) -> list[Path]:
    n_row_cols = len(row_headers)
    n_data_cols = len(col_headers)
    n_rows = len(display_row_labels)
    n_groups = (group_ids[-1] + 1) if group_ids else 1

    title_fs = 11.5
    header_fs = 9.5
    body_fs = 9.2
    fail_fs = 8.0

    pad_left = 0.10
    pad_right = 0.10
    gutter = 0.16 if (n_row_cols > 0 and n_data_cols > 0) else 0.0

    col_widths: list[float] = []
    for j in range(n_row_cols):
        w = _estimate_text_width_in(row_headers[j], header_fs, bold=True)
        for r_idx in range(n_rows):
            w = max(w, _estimate_text_width_in(display_row_labels[r_idx][j], body_fs))
        col_widths.append(max(w + pad_left + pad_right, 0.52))

    for j in range(n_data_cols):
        w = _estimate_text_width_in(col_headers[j], header_fs, bold=True)
        for r_idx in range(n_rows):
            cell = cells[r_idx][j]
            if cell.kind in ("fail", "nc"):
                cw = _estimate_text_width_in(cell.text, fail_fs, bold=True) * 1.05
            else:
                cw = _estimate_text_width_in(cell.text, body_fs, bold=cell.is_best)
            w = max(w, cw)
        col_widths.append(max(w + pad_left + pad_right, 0.56))

    table_w = sum(col_widths) + gutter

    # Wrap long titles gently so narrow tables stay well-proportioned
    title_lines: list[str] = []
    if title.strip():
        max_title_w = max(table_w * 1.25, 4.8)
        est_title_w = _estimate_text_width_in(title.strip(), title_fs, bold=True)
        if est_title_w > max_title_w:
            avg_char_w = max(est_title_w / max(len(title.strip()), 1), 0.06)
            wrap_chars = max(int(max_title_w / avg_char_w), 28)
            title_lines = textwrap.wrap(title.strip(), width=wrap_chars)
        else:
            title_lines = [title.strip()]

    max_rendered_title_w = (
        max(_estimate_text_width_in(line, title_fs, bold=True) for line in title_lines)
        if title_lines
        else 0.0
    )

    if max_rendered_title_w > table_w and n_data_cols > 0:
        extra = max_rendered_title_w - table_w
        per_data_col = extra / n_data_cols
        for j in range(n_row_cols, n_row_cols + n_data_cols):
            col_widths[j] += per_data_col
        table_w = sum(col_widths) + gutter

    margin_x = 0.24
    margin_top = 0.20
    margin_bottom = 0.20
    title_line_h = 0.25
    title_gap = 0.12
    header_h = 0.34
    top_body_pad = 0.04
    bottom_body_pad = 0.05
    row_h = 0.27
    group_gap = 0.07

    title_block_h = (
        len(title_lines) * title_line_h + title_gap if title_lines else 0.0
    )
    body_h = (
        top_body_pad
        + n_rows * row_h
        + max(0, n_groups - 1) * group_gap
        + bottom_body_pad
    )

    fig_w = max(table_w, max_rendered_title_w) + 2.0 * margin_x
    fig_h = margin_top + title_block_h + header_h + body_h + margin_bottom

    x_left = margin_x
    x_right = x_left + table_w

    col_starts: list[float] = []
    col_ends: list[float] = []
    x_cur = x_left
    for j, w in enumerate(col_widths):
        if j == n_row_cols and gutter > 0.0:
            x_cur += gutter
        col_starts.append(x_cur)
        col_ends.append(x_cur + w)
        x_cur += w

    rc_overrides = {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Liberation Sans", "Arial", "Helvetica"],
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    }

    written: list[Path] = []
    with plt.rc_context(rc_overrides):
        fig = plt.figure(figsize=(fig_w, fig_h), dpi=300)
        try:
            ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
            ax.set_xlim(0.0, fig_w)
            ax.set_ylim(0.0, fig_h)
            ax.axis("off")

            y_cur = fig_h - margin_top

            # Title
            if title_lines:
                for line in title_lines:
                    y_cur -= title_line_h
                    ax.text(
                        x_left,
                        y_cur + title_line_h * 0.25,
                        line,
                        fontsize=title_fs,
                        fontweight="bold",
                        color="#111111",
                        ha="left",
                        va="baseline",
                    )
                y_cur -= title_gap

            # Header band + top/mid rules
            y_top_rule = y_cur
            y_mid_rule = y_top_rule - header_h
            ax.add_patch(
                Rectangle(
                    (x_left, y_mid_rule),
                    x_right - x_left,
                    header_h,
                    facecolor="#F4F6F9",
                    edgecolor="none",
                    zorder=0,
                )
            )
            ax.hlines(
                y_top_rule,
                x_left,
                x_right,
                colors="#1A1A1A",
                linewidth=1.4,
                zorder=3,
            )
            ax.hlines(
                y_mid_rule,
                x_left,
                x_right,
                colors="#333333",
                linewidth=0.9,
                zorder=3,
            )

            y_hdr_mid = y_mid_rule + header_h * 0.5
            for j, h_text in enumerate(row_headers):
                ax.text(
                    col_starts[j] + pad_left,
                    y_hdr_mid,
                    h_text,
                    fontsize=header_fs,
                    fontweight="bold",
                    color="#1A1A1A",
                    ha="left",
                    va="center",
                    zorder=2,
                )
            for j, h_text in enumerate(col_headers):
                col_i = n_row_cols + j
                ax.text(
                    col_ends[col_i] - pad_right,
                    y_hdr_mid,
                    h_text,
                    fontsize=header_fs,
                    fontweight="bold",
                    color="#1A1A1A",
                    ha="right",
                    va="center",
                    zorder=2,
                )

            # Table body rows
            y_cur = y_mid_rule - top_body_pad
            use_group_shading = n_groups > 1
            use_row_shading = (not use_group_shading) and (n_rows >= 4)

            for r_idx in range(n_rows):
                if r_idx > 0 and group_ids[r_idx] != group_ids[r_idx - 1]:
                    rule_y = y_cur - group_gap * 0.5
                    ax.hlines(
                        rule_y,
                        x_left,
                        x_right,
                        colors="#D0D5DD",
                        linewidth=0.5,
                        zorder=2,
                    )
                    y_cur -= group_gap

                row_bottom = y_cur - row_h
                row_mid = row_bottom + row_h * 0.5

                shade = (
                    (group_ids[r_idx] % 2 == 1)
                    if use_group_shading
                    else (use_row_shading and (r_idx % 2 == 1))
                )
                if shade:
                    ax.add_patch(
                        Rectangle(
                            (x_left, row_bottom),
                            x_right - x_left,
                            row_h,
                            facecolor="#F9FAFB",
                            edgecolor="none",
                            zorder=0,
                        )
                    )

                for j, lbl in enumerate(display_row_labels[r_idx]):
                    if not lbl:
                        continue
                    ax.text(
                        col_starts[j] + pad_left,
                        row_mid,
                        lbl,
                        fontsize=body_fs,
                        fontweight="normal",
                        color="#222222",
                        ha="left",
                        va="center",
                        zorder=2,
                    )

                for j, cell in enumerate(cells[r_idx]):
                    col_i = n_row_cols + j
                    x_pos = col_ends[col_i] - pad_right
                    if cell.kind == "missing":
                        ax.text(
                            x_pos,
                            row_mid,
                            cell.text,
                            fontsize=body_fs,
                            fontweight="normal",
                            color="#999999",
                            ha="right",
                            va="center",
                            zorder=2,
                        )
                    elif cell.kind in ("fail", "nc"):
                        ax.text(
                            x_pos,
                            row_mid,
                            cell.text,
                            fontsize=fail_fs,
                            fontfamily="monospace",
                            fontweight="medium",
                            color="#B03030",
                            ha="right",
                            va="center",
                            zorder=2,
                        )
                    elif cell.is_best:
                        ax.text(
                            x_pos,
                            row_mid,
                            cell.text,
                            fontsize=body_fs,
                            fontweight="bold",
                            color="#111111",
                            ha="right",
                            va="center",
                            zorder=2,
                        )
                    else:
                        ax.text(
                            x_pos,
                            row_mid,
                            cell.text,
                            fontsize=body_fs,
                            fontweight="normal",
                            color="#222222",
                            ha="right",
                            va="center",
                            zorder=2,
                        )

                y_cur = row_bottom

            # Bottom rule
            y_bottom_rule = y_cur - bottom_body_pad
            ax.hlines(
                y_bottom_rule,
                x_left,
                x_right,
                colors="#1A1A1A",
                linewidth=1.4,
                zorder=3,
            )

            out_dir.mkdir(parents=True, exist_ok=True)
            for fmt in formats:
                ext = fmt.strip().lower().lstrip(".")
                if not ext:
                    continue
                out_path = out_dir / f"{stem}.{ext}"
                save_kwargs: dict[str, Any] = {
                    "bbox_inches": "tight",
                    "pad_inches": 0.06,
                    "facecolor": "white",
                }
                if ext == "png":
                    save_kwargs["dpi"] = 300
                fig.savefig(out_path, **save_kwargs)
                written.append(out_path)
        finally:
            plt.close(fig)

    return written


def _md_escape(text: str) -> str:
    return text.replace("|", "\\|")


def _write_markdown_table(
    title: str,
    row_headers: list[str],
    col_headers: list[str],
    full_row_labels: list[list[str]],
    cells: list[list[_CellData]],
    out_dir: Path,
    stem: str,
) -> Path:
    n_row_cols = len(row_headers)
    headers = [_md_escape(h) for h in [*row_headers, *col_headers]]

    rows_str: list[list[str]] = []
    for r_labels, r_cells in zip(full_row_labels, cells):
        row_out = [_md_escape(lbl) for lbl in r_labels]
        for cell in r_cells:
            escaped = _md_escape(cell.text)
            row_out.append(f"**{escaped}**" if cell.is_best else escaped)
        rows_str.append(row_out)

    widths = [max(len(h), 3) for h in headers]
    for row in rows_str:
        for j, val in enumerate(row):
            if len(val) > widths[j]:
                widths[j] = len(val)

    def _fmt_md_row(vals: list[str]) -> str:
        padded: list[str] = []
        for j, val in enumerate(vals):
            if j < n_row_cols:
                padded.append(val.ljust(widths[j]))
            else:
                padded.append(val.rjust(widths[j]))
        return "| " + " | ".join(padded) + " |"

    sep_cells: list[str] = []
    for j, w in enumerate(widths):
        if j < n_row_cols:
            sep_cells.append(":" + "-" * max(w - 1, 2))
        else:
            sep_cells.append("-" * max(w - 1, 2) + ":")
    sep_line = "| " + " | ".join(sep_cells) + " |"

    lines: list[str] = []
    if title.strip():
        lines.append(f"## {title.strip()}")
        lines.append("")
    lines.append(_fmt_md_row(headers))
    lines.append(sep_line)
    for row in rows_str:
        lines.append(_fmt_md_row(row))
    lines.append("")

    out_dir.mkdir(parents=True, exist_ok=True)
    md_path = out_dir / f"{stem}.md"
    md_path.write_text("\n".join(lines), encoding="utf-8")
    return md_path


def render_sparse_table(
    dataset: BenchmarkDataset,
    spec: SparseTableSpec,
    out_dir: Path,
    formats: list[str] | None = None,
) -> list[Path]:
    """Render a sparse aggregate table to Matplotlib image(s) and companion Markdown."""
    df = dataset.df
    if df.empty:
        return []

    required_cols = [*spec.row_axes, spec.col_axis, spec.value]
    if any(col not in df.columns for col in required_cols):
        return []

    out_dir = Path(out_dir)
    img_formats = list(formats) if formats is not None else ["png"]
    split_cols = [c for c in spec.split_by if c in df.columns]
    splits = _iter_splits(df, split_cols)

    base_title = spec.title.strip()
    if not base_title:
        base_title = (
            f"{spec.value} by {', '.join(_axis_display_name(a) for a in spec.row_axes)} "
            f"vs {_axis_display_name(spec.col_axis)}"
        )

    generated: list[Path] = []
    for split_pairs, sub_df in splits:
        if sub_df.empty:
            continue

        grid = _assemble_table_grid(sub_df, spec, dataset.timing_cols)
        if grid is None:
            continue

        (
            row_headers,
            col_headers,
            full_row_labels,
            display_row_labels,
            cells,
            group_ids,
        ) = grid

        if split_pairs:
            split_label = ", ".join(
                f"{k}={_fmt_axis_val(v, spec.missing)}" for k, v in split_pairs
            )
            table_title = f"{base_title} ({split_label})"
            split_slug = "__".join(
                f"{slugify(k)}_{slugify(_fmt_axis_val(v, spec.missing))}"
                for k, v in split_pairs
            )
            stem = f"{slugify(spec.name)}__{split_slug}"
        else:
            table_title = base_title
            stem = slugify(spec.name)

        if img_formats:
            generated.extend(
                _render_matplotlib_table(
                    title=table_title,
                    row_headers=row_headers,
                    col_headers=col_headers,
                    display_row_labels=display_row_labels,
                    cells=cells,
                    group_ids=group_ids,
                    out_dir=out_dir,
                    stem=stem,
                    formats=img_formats,
                )
            )

        md_path = _write_markdown_table(
            title=table_title,
            row_headers=row_headers,
            col_headers=col_headers,
            full_row_labels=full_row_labels,
            cells=cells,
            out_dir=out_dir,
            stem=stem,
        )
        generated.append(md_path)

    return generated
