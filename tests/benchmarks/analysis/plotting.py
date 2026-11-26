"""Declarative plotting engine for the MFEM NURBS Kronecker benchmark suite.

Provides three reusable plot primitives (faceted line plots, stacked stage
breakdowns, and residual convergence plots), sparse-table dispatch via
``analysis.sparse_tables``, color-first series encoding, strict per-subplot
legends, an Anti-Collision Rule for arbitrary sweep parameters, per-plot
filtering, and manual-first orchestration of ``[[plots]]`` specifications.
"""

from __future__ import annotations

import math
import re
import sys
from collections.abc import Collection, Iterable, Sequence
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any, TypeVar

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from .ingest import (
    BenchmarkDataset,
    apply_filters,
    fmt_int_like,
    is_valid_number,
    slugify,
)

try:
    from .sparse_tables import parse_sparse_table_spec, render_sparse_table
except ImportError:  # pragma: no cover - allows standalone import before module exists
    parse_sparse_table_spec = None  # type: ignore[assignment]
    render_sparse_table = None  # type: ignore[assignment]

T = TypeVar("T")

# ---------------------------------------------------------------------------
# Publication Styling & Visual Encoding Constants
# ---------------------------------------------------------------------------

# 10-color colorblind-safe palette (Okabe-Ito + Paul Tol high-contrast tones)
PALETTE: tuple[str, ...] = (
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#332288",  # indigo
    "#882255",  # wine
    "#44AA99",  # teal
    "#000000",  # black
)

MARKERS: tuple[str, ...] = (
    "o",
    "s",
    "^",
    "D",
    "v",
    "P",
    "X",
    "*",
    "<",
    ">",
    "h",
    "p",
)

LINESTYLES: tuple[str | tuple[int, tuple[int, ...]], ...] = (
    "-",
    "--",
    "-.",
    ":",
    (0, (5, 1, 1, 1)),
    (0, (3, 1, 1, 1, 1, 1)),
)

DEFAULT_TIMING_STAGES: tuple[tuple[str, str, str], ...] = (
    ("operator_assembly", "Operator Assembly", "#0072B2"),
    ("prec_univariate_assembly", "Univariate Assembly", "#56B4E9"),
    ("prec_eigendecomp", "Eigendecomposition", "#E69F00"),
    ("prec_fd_construct", "FD Construct", "#CC79A7"),
    ("prec_other_setup", "Other Prec Setup", "#882255"),
    ("solve", "Solve", "#009E73"),
)

_PREC_ORDER: dict[str, int] = {
    "none": 0,
    "jacobi": 1,
    "gs": 2,
    "gauss-seidel": 2,
    "ilu": 3,
    "amg": 4,
    "boomeramg": 4,
    "kronecker-fd": 5,
    "kronecker": 6,
}

_ASSEMBLY_ORDER: dict[str, int] = {
    "legacy": 0,
    "wq": 1,
}

_REF_STRATEGY_ORDER: dict[str, int] = {
    "h": 0,
    "p": 1,
    "k": 2,
}

_NAKED_PRIMARY_AXES = frozenset({"prec", "assembly", "solver", "mesh_name"})
_NAKED_SECONDARY_AXES = frozenset({"prec", "assembly", "solver", "ref_strategy"})

_STAGE_EPS = 1e-12
_DWARF_RATIO_THRESHOLD = 20.0
_NATURAL_SPLIT_RE = re.compile(r"(\d+)")


def _configure_matplotlib() -> None:
    """Apply publication-ready Matplotlib defaults."""
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": [
                "DejaVu Serif",
                "Times New Roman",
                "Computer Modern Roman",
                "serif",
            ],
            "mathtext.fontset": "dejavuserif",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8.5,
            "figure.titlesize": 12,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.35,
            "grid.linestyle": "--",
            "grid.linewidth": 0.6,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


_configure_matplotlib()


# ---------------------------------------------------------------------------
# Value Normalization, Sorting, Formatting & Export Helpers
# ---------------------------------------------------------------------------


def _is_missing(val: Any) -> bool:
    """Return True if ``val`` is None or NaN/NA."""
    if val is None:
        return True
    try:
        return bool(pd.isna(val))
    except (TypeError, ValueError):
        return False


def _canonical_key(val: Any) -> tuple[str, Any]:
    """Produce a hashable, type-safe key so bools and ints never collide."""
    if _is_missing(val):
        return ("na", None)
    if isinstance(val, (bool, np.bool_)):
        return ("b", bool(val))
    if isinstance(val, (int, np.integer)):
        return ("n", int(val))
    if isinstance(val, (float, np.floating)):
        fval = float(val)
        if math.isfinite(fval) and fval.is_integer():
            return ("n", int(fval))
        return ("n", fval)
    return ("s", str(val))


def _natural_sort_key(text: str) -> list[tuple[int, Any]]:
    """Split text into numeric and alphabetic runs for natural sorting."""
    parts = _NATURAL_SPLIT_RE.split(text.lower())
    key: list[tuple[int, Any]] = []
    for part in parts:
        if not part:
            continue
        if part.isdigit():
            key.append((0, int(part)))
        else:
            key.append((1, part))
    return key


def _value_sort_key(axis: str, val: Any) -> tuple[Any, ...]:
    """Return a canonical sort key for known categorical and numeric axes."""
    if _is_missing(val):
        return (9, 0, "")

    if isinstance(val, (bool, np.bool_)):
        return (1, int(bool(val)), str(bool(val)).lower())

    sval = str(val).strip().lower()
    if axis == "prec":
        if sval in _PREC_ORDER:
            return (0, _PREC_ORDER[sval], _natural_sort_key(str(val)))
        return (2, 0, _natural_sort_key(str(val)))

    if axis == "assembly":
        if sval in _ASSEMBLY_ORDER:
            return (0, _ASSEMBLY_ORDER[sval], _natural_sort_key(str(val)))
        return (2, 0, _natural_sort_key(str(val)))

    if axis == "ref_strategy":
        if sval in _REF_STRATEGY_ORDER:
            return (0, _REF_STRATEGY_ORDER[sval], _natural_sort_key(str(val)))
        return (2, 0, _natural_sort_key(str(val)))

    if is_valid_number(val):
        return (1, float(val), "")

    try:
        fval = float(str(val))
        if math.isfinite(fval):
            return (1, fval, "")
    except ValueError:
        pass

    return (2, 0, _natural_sort_key(str(val)))


def sort_axis_values(axis: str, values: Iterable[Any]) -> list[Any]:
    """Deduplicate and canonically sort unique values for an axis."""
    seen: dict[tuple[str, Any], Any] = {}
    for val in values:
        key = _canonical_key(val)
        if key not in seen:
            seen[key] = val
    return sorted(seen.values(), key=lambda v: _value_sort_key(axis, v))


def build_style_map(
    df: pd.DataFrame,
    axis: str | None,
    styles: Sequence[T],
) -> dict[tuple[str, Any], T]:
    """Build a consistent canonical-value -> visual-style map from a DataFrame."""
    if not axis or axis not in df.columns or not styles:
        return {}
    ordered = sort_axis_values(axis, df[axis].tolist())
    return {
        _canonical_key(val): styles[idx % len(styles)]
        for idx, val in enumerate(ordered)
    }


def _fmt_val(val: Any) -> str:
    """Format a scalar parameter value cleanly for titles, legends, and filenames."""
    if _is_missing(val):
        return "na"
    if isinstance(val, (bool, np.bool_)):
        return "true" if bool(val) else "false"
    if isinstance(val, (int, np.integer)):
        return str(int(val))
    if isinstance(val, (float, np.floating)):
        fval = float(val)
        if not math.isfinite(fval):
            return "na"
        if fval.is_integer():
            return fmt_int_like(fval)
        return f"{fval:g}"
    return str(val)


def _fmt_bar_total(total_sec: float) -> str:
    """Format a positive wall-time total compactly for bar-top annotations."""
    if not math.isfinite(total_sec) or total_sec < 0.0:
        return "na"
    if total_sec < 0.001:
        return f"{total_sec:.4f}s"
    if total_sec < 0.1:
        return f"{total_sec:.3f}s"
    if total_sec < 10.0:
        return f"{total_sec:.2f}s"
    if total_sec < 100.0:
        return f"{total_sec:.1f}s"
    return f"{total_sec:.0f}s"


def _humanize_col_name(col: str) -> str:
    """Convert a snake_case column name into a clean display label."""
    if col == "num_free_dofs":
        return "Free DoFs"
    if col == "ref_levels":
        return "Refinement Level"
    if col == "order":
        return "Order"
    if col == "pipeline_time_sec":
        return "Pipeline Time (s)"
    if col == "operator_assembly":
        return "Operator Assembly Time (s)"
    if col == "solve":
        return "Solve Time (s)"
    if col == "iterations":
        return "Iterations"
    return col.replace("_", " ").title()


def _format_split_context(
    split_items: Sequence[tuple[str, Any]],
) -> tuple[str, str]:
    """Return ``(title_suffix, filename_slug_suffix)`` for a split combination."""
    if not split_items:
        return "", ""
    title_parts: list[str] = []
    slug_parts: list[str] = []
    for axis, val in split_items:
        val_str = _fmt_val(val)
        title_parts.append(f"{axis}={val_str}")
        slug_parts.append(f"{slugify(axis)}={slugify(val_str)}")
    return ", ".join(title_parts), "__".join(slug_parts)


def _compose_title(
    base_title: str,
    context_parts: Sequence[str],
) -> str:
    """Join a base title with optional parenthesized context strings."""
    active = [part for part in context_parts if part]
    if not active:
        return base_title
    return f"{base_title} ({', '.join(active)})"


def _is_explicitly_false(val: Any) -> bool:
    """Return True if a value explicitly represents boolean False."""
    if _is_missing(val):
        return False
    if isinstance(val, (bool, np.bool_)):
        return not bool(val)
    if isinstance(val, (int, np.integer, float, np.floating)):
        return float(val) == 0.0
    if isinstance(val, str):
        return val.strip().lower() in {"false", "0", "no"}
    return False


def _filter_ok_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Filter DataFrame to rows with ``status == 'ok'`` if a status column exists."""
    if df.empty or "status" not in df.columns:
        return df.copy()
    return df.loc[df["status"].astype(str).str.lower() == "ok"].copy()


def _iter_split_groups(
    df: pd.DataFrame,
    split_axes: Sequence[str],
) -> list[tuple[list[tuple[str, Any]], pd.DataFrame]]:
    """Partition ``df`` into canonically sorted groups for ``split_axes``."""
    active_splits = [a for a in split_axes if a in df.columns]
    if not active_splits:
        return [([], df)]

    key_cols: list[str] = []
    work = df.copy()
    for idx, axis in enumerate(active_splits):
        kcol = f"__split_key_{idx}"
        work[kcol] = work[axis].map(_canonical_key)
        key_cols.append(kcol)

    grouped: list[tuple[tuple[Any, ...], list[tuple[str, Any]], pd.DataFrame]] = []
    for _, grp in work.groupby(key_cols, sort=False, dropna=False):
        split_items: list[tuple[str, Any]] = []
        sort_tuple: list[Any] = []
        for axis in active_splits:
            rep_val = grp[axis].iloc[0]
            split_items.append((axis, rep_val))
            sort_tuple.append(_value_sort_key(axis, rep_val))
        clean_grp = grp.drop(columns=key_cols)
        grouped.append((tuple(sort_tuple), split_items, clean_grp))

    grouped.sort(key=lambda item: item[0])
    return [(items, grp) for _, items, grp in grouped]


def save_figure(
    fig: Figure,
    base_path: Path,
    formats: Sequence[str] = ("png",),
) -> list[Path]:
    """Save ``fig`` to ``base_path.<ext>`` for each requested format and close it."""
    saved: list[Path] = []
    try:
        base_path.parent.mkdir(parents=True, exist_ok=True)
        norm_formats: list[str] = []
        for fmt in formats:
            ext = fmt.strip().lower().lstrip(".")
            if ext and ext not in norm_formats:
                norm_formats.append(ext)
        if not norm_formats:
            norm_formats = ["png"]

        for ext in norm_formats:
            out_path = base_path.with_suffix(f".{ext}")
            fig.savefig(out_path, format=ext, dpi=300, bbox_inches="tight")
            saved.append(out_path)
    finally:
        plt.close(fig)
    return saved


# ---------------------------------------------------------------------------
# Color-First Series Visual Mapping
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _SeriesStyle:
    color: str
    marker: str
    linestyle: str | tuple[int, tuple[int, ...]]
    label: str


def _normalize_hue_cols(
    hue: str | Sequence[str] | None,
    df: pd.DataFrame,
    exclude: Collection[str] = (),
) -> list[str]:
    """Resolve ``hue`` (single column or list of columns) to active DataFrame columns."""
    if hue is None:
        return []
    candidates = [hue] if isinstance(hue, str) else [str(c) for c in hue]
    excluded = set(exclude)
    resolved: list[str] = []
    for col in candidates:
        col_str = col.strip()
        if (
            col_str
            and col_str not in excluded
            and col_str not in resolved
            and col_str in df.columns
            and df[col_str].notna().any()
        ):
            resolved.append(col_str)
    return resolved


def _format_series_label(
    cols: Sequence[str],
    vals: Sequence[Any],
) -> str:
    """Build a clean, compact legend label from one or more series columns.

    Examples:
        - ``["prec"]`` -> ``"kronecker-fd"``
        - ``["prec", "ref_strategy"]`` -> ``"kronecker-fd (k)"``
        - ``["prec", "order"]`` -> ``"kronecker-fd (order=3)"``
    """
    if not cols:
        return "value"
    first_col, first_val = cols[0], _fmt_val(vals[0])
    first_str = (
        first_val
        if first_col in _NAKED_PRIMARY_AXES
        else f"{first_col}={first_val}"
    )
    if len(cols) == 1:
        return first_str

    rest_parts: list[str] = []
    for col, val in zip(cols[1:], vals[1:]):
        v_str = _fmt_val(val)
        if col in _NAKED_SECONDARY_AXES:
            rest_parts.append(v_str)
        else:
            rest_parts.append(f"{col}={v_str}")
    return f"{first_str} ({', '.join(rest_parts)})"


def _extract_sorted_combos(
    df: pd.DataFrame,
    cols: Sequence[str],
) -> list[tuple[tuple[tuple[str, Any], ...], tuple[Any, ...]]]:
    """Return canonically sorted ``(canonical_key_tuple, raw_value_tuple)`` pairs."""
    if not cols:
        return [((), ())]

    seen: dict[tuple[tuple[str, Any], ...], tuple[Any, ...]] = {}
    for row_vals in df[list(cols)].itertuples(index=False, name=None):
        key_tup = tuple(_canonical_key(v) for v in row_vals)
        if key_tup not in seen:
            seen[key_tup] = row_vals

    items = list(seen.items())
    items.sort(
        key=lambda item: tuple(
            _value_sort_key(col, val) for col, val in zip(cols, item[1])
        )
    )
    return items


def _build_series_visual_map(
    df: pd.DataFrame,
    hue_cols: Sequence[str],
    style_col: str | None,
    use_linestyles: bool = False,
) -> tuple[list[str], dict[tuple[tuple[str, Any], ...], _SeriesStyle]]:
    """Assign color-first visual styles to every unique series combination.

    While the number of unique series is ``<= len(PALETTE)`` and
    ``use_linestyles`` is False, every unique series receives its own distinct
    color from ``PALETTE``, a distinct marker, and a solid line (``"-"``).
    Non-solid linestyles are only cycled when total series exceed
    ``len(PALETTE)`` or ``use_linestyles`` is explicitly True.
    """
    series_cols: list[str] = list(hue_cols)
    has_style = style_col is not None and style_col in df.columns and style_col not in series_cols
    if has_style:
        assert style_col is not None
        series_cols.append(style_col)

    combos = _extract_sorted_combos(df, series_cols)
    n_series = len(combos)
    visual_map: dict[tuple[tuple[str, Any], ...], _SeriesStyle] = {}

    if not use_linestyles and n_series <= len(PALETTE):
        for idx, (key_tup, val_tup) in enumerate(combos):
            visual_map[key_tup] = _SeriesStyle(
                color=PALETTE[idx],
                marker=MARKERS[idx % len(MARKERS)],
                linestyle="-",
                label=_format_series_label(series_cols, val_tup),
            )
        return series_cols, visual_map

    # Cycle linestyles only when total series > len(PALETTE) or use_linestyles=True
    if has_style and hue_cols:
        hue_combos = _extract_sorted_combos(df, hue_cols)
        hue_idx_map = {k: i for i, (k, _) in enumerate(hue_combos)}
        assert style_col is not None
        style_vals = sort_axis_values(style_col, df[style_col].tolist())
        style_idx_map = {_canonical_key(v): i for i, v in enumerate(style_vals)}

        for key_tup, val_tup in combos:
            h_key = key_tup[:-1]
            s_key = key_tup[-1]
            h_idx = hue_idx_map.get(h_key, 0)
            s_idx = style_idx_map.get(s_key, 0)
            ls_idx = (s_idx + (h_idx // len(PALETTE))) % len(LINESTYLES)
            visual_map[key_tup] = _SeriesStyle(
                color=PALETTE[h_idx % len(PALETTE)],
                marker=MARKERS[h_idx % len(MARKERS)],
                linestyle=LINESTYLES[ls_idx],
                label=_format_series_label(series_cols, val_tup),
            )
        return series_cols, visual_map

    for idx, (key_tup, val_tup) in enumerate(combos):
        ls_idx = (
            idx % len(LINESTYLES)
            if use_linestyles
            else (idx // len(PALETTE)) % len(LINESTYLES)
        )
        visual_map[key_tup] = _SeriesStyle(
            color=PALETTE[idx % len(PALETTE)],
            marker=MARKERS[idx % len(MARKERS)],
            linestyle=LINESTYLES[ls_idx],
            label=_format_series_label(series_cols, val_tup),
        )
    return series_cols, visual_map


# ---------------------------------------------------------------------------
# Anti-Collision Rule
# ---------------------------------------------------------------------------


def resolve_split_axes(
    varying_axes: list[str],
    mapped_axes: set[str],
    agg_over: set[str],
    explicit_split: list[str],
    df_columns: Collection[str] | None = None,
) -> list[str]:
    """Determine figure split axes so no unmapped varying axis ever collides.

    Starts with axes from ``explicit_split`` that exist in ``varying_axes`` (or
    ``df_columns`` when provided) and are not already consumed by
    ``mapped_axes`` or ``agg_over``. Any remaining axis in ``varying_axes`` not
    covered by ``mapped_axes``, ``agg_over``, or the explicit split list is
    automatically appended.
    """
    consumed = set(mapped_axes) | set(agg_over)
    varying_set = set(varying_axes)
    col_set = set(df_columns) if df_columns is not None else set()

    result: list[str] = []
    seen: set[str] = set()

    for axis in explicit_split:
        if axis in seen or axis in consumed:
            continue
        if axis in varying_set or (df_columns is not None and axis in col_set):
            result.append(axis)
            seen.add(axis)

    for axis in varying_axes:
        if axis not in consumed and axis not in seen:
            result.append(axis)
            seen.add(axis)

    return result


# ---------------------------------------------------------------------------
# Plot Spec Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class LinePlotSpec:
    """Declarative specification for a faceted or single-axis line plot."""

    name: str
    x: str
    y: str
    xlabel: str = ""
    ylabel: str = ""
    title_prefix: str = ""
    title: str = ""
    xscale: str = "linear"
    yscale: str = "linear"
    hue: str | list[str] | None = None
    style: str | None = None
    facet_col: str | None = None
    split_by: list[str] = field(default_factory=list)
    agg_over: set[str] = field(default_factory=set)
    hollow_non_converged: bool = False
    ylim_bottom_zero: bool = False
    xticks_from_data: bool = False
    use_linestyles: bool = False
    min_points_per_panel: int = 1
    drop_empty_facets: bool = True
    reference_slopes: list[float] = field(default_factory=list)


@dataclass
class StackedBarSpec:
    """Declarative specification for a grouped stacked-bar breakdown plot."""

    name: str
    x: str = "ref_levels"
    series: list[str] = field(default_factory=list)
    stages: list[tuple[str, str, str]] = field(default_factory=list)
    ylabel: str | None = None
    title_prefix: str = "Timing Breakdown"
    title: str = ""
    split_by: list[str] = field(default_factory=list)
    agg_over: set[str] = field(default_factory=set)
    annotate_dofs: bool = True
    normalize: bool = True
    facet_col: str | None = None
    drop_empty_facets: bool = True


@dataclass
class ConvergencePlotSpec:
    """Declarative specification for a semilogy residual convergence plot."""

    name: str
    hue: str | list[str] | None = None
    style: str | None = None
    facet_col: str | None = None
    split_by: list[str] = field(default_factory=list)
    title_prefix: str = "Residual Convergence"
    title: str = ""
    use_linestyles: bool = False


# ---------------------------------------------------------------------------
# Primitive 1: Faceted Line Plot
# ---------------------------------------------------------------------------


def _format_slope_label(slope: float) -> str:
    """Format a log-log reference slope as clean LaTeX complexity notation."""
    if abs(slope - round(slope)) < 1e-2:
        k = int(round(slope))
        if k == 1:
            return r"$\mathcal{O}(N)$"
        return rf"$\mathcal{{O}}(N^{{{k}}})$"

    frac = Fraction(slope).limit_denominator(12)
    if abs(slope - float(frac)) < 0.015:
        if frac.denominator == 1:
            if frac.numerator == 1:
                return r"$\mathcal{O}(N)$"
            return rf"$\mathcal{{O}}(N^{{{frac.numerator}}})$"
        return rf"$\mathcal{{O}}(N^{{{frac.numerator}/{frac.denominator}}})$"

    return rf"$\mathcal{{O}}(N^{{{slope:.2g}}})$"


def _draw_reference_slopes(
    ax: Axes,
    panel_df: pd.DataFrame,
    spec: LinePlotSpec,
) -> list[tuple[str, Any]]:
    """Overlay subtle gray dotted reference slope lines on a log-log axes."""
    if not spec.reference_slopes or spec.xscale != "log" or spec.yscale != "log":
        return []

    xs = panel_df[spec.x].to_numpy(dtype=float)
    ys = panel_df[spec.y].to_numpy(dtype=float)
    pos = np.isfinite(xs) & np.isfinite(ys) & (xs > 0.0) & (ys > 0.0)
    xs, ys = xs[pos], ys[pos]
    if xs.size < 2:
        return []

    x_min, x_max = float(np.min(xs)), float(np.max(xs))
    if x_max <= x_min * (1.0 + 1e-9):
        return []

    ymin_cur, ymax_cur = ax.get_ylim()
    if not (math.isfinite(ymin_cur) and math.isfinite(ymax_cur) and ymin_cur > 0.0 and ymax_cur > ymin_cur):
        ymin_cur = float(np.min(ys)) * 0.8
        ymax_cur = float(np.max(ys)) * 1.2

    lx0, lx1 = math.log10(x_min), math.log10(x_max)
    lx_span = lx1 - lx0
    seg_lx_a = lx0 + 0.12 * lx_span if lx_span >= 0.35 else lx0
    seg_lx_b = lx1 - 0.04 * lx_span if lx_span >= 0.35 else lx1
    lx_anchor = 0.5 * (seg_lx_a + seg_lx_b)

    ly_min_data = math.log10(float(np.min(ys)))
    ly_max_data = math.log10(float(np.max(ys)))
    ly_lo_bound = math.log10(ymin_cur) + 0.04 * (math.log10(ymax_cur) - math.log10(ymin_cur))
    ly_hi_bound = math.log10(ymax_cur) - 0.04 * (math.log10(ymax_cur) - math.log10(ymin_cur))

    ref_styles = (":", "--", "-.")
    ref_grays = ("#666666", "#888888", "#4D4D4D")
    legend_items: list[tuple[str, Any]] = []

    n_slopes = len(spec.reference_slopes)
    for idx, slope in enumerate(spec.reference_slopes):
        if not is_valid_number(slope):
            continue
        s_val = float(slope)

        # Distribute anchor heights neatly across the lower-to-mid data band
        frac_pos = 0.25 if n_slopes == 1 else (0.18 + 0.55 * idx / max(n_slopes - 1, 1))
        ly_anchor = ly_min_data + frac_pos * max(ly_max_data - ly_min_data, 0.4)
        ly_anchor = min(max(ly_anchor, ly_lo_bound), ly_hi_bound)

        lx_a, lx_b = seg_lx_a, seg_lx_b
        if abs(s_val) > 1e-12:
            t1 = lx_anchor + (ly_lo_bound - ly_anchor) / s_val
            t2 = lx_anchor + (ly_hi_bound - ly_anchor) / s_val
            lx_a = max(lx_a, min(t1, t2))
            lx_b = min(lx_b, max(t1, t2))

        if lx_b <= lx_a + 1e-6:
            continue

        x_ref = np.array([10.0**lx_a, 10.0**lx_b], dtype=float)
        y_ref = np.array(
            [
                10.0 ** (ly_anchor + s_val * (lx_a - lx_anchor)),
                10.0 ** (ly_anchor + s_val * (lx_b - lx_anchor)),
            ],
            dtype=float,
        )

        label = _format_slope_label(s_val)
        (ref_line,) = ax.plot(
            x_ref,
            y_ref,
            color=ref_grays[idx % len(ref_grays)],
            linestyle=ref_styles[idx % len(ref_styles)],
            linewidth=1.35,
            alpha=0.85,
            label=label,
            zorder=2,
        )
        legend_items.append((label, ref_line))

    ax.set_ylim(ymin_cur, ymax_cur)
    return legend_items


def _aggregate_line_df(
    df: pd.DataFrame,
    spec: LinePlotSpec,
    varying_axes: list[str],
    split_axes: Sequence[str],
    hue_cols: Sequence[str],
    style_col: str | None,
    facet_col: str | None,
    agg_over: set[str],
) -> pd.DataFrame:
    """Aggregate ``x`` and ``y`` via median across ``agg_over`` while preserving other axes."""
    if not agg_over:
        return df

    group_cols: list[str] = []

    def _add_group_col(col: str | None) -> None:
        if col and col in df.columns and col not in agg_over and col not in group_cols:
            group_cols.append(col)

    for col in split_axes:
        _add_group_col(col)
    for col in hue_cols:
        _add_group_col(col)
    _add_group_col(style_col)
    _add_group_col(facet_col)

    for col in varying_axes:
        _add_group_col(col)

    if spec.x == "num_free_dofs":
        if "ref_levels" in df.columns and "ref_levels" not in agg_over:
            _add_group_col("ref_levels")
        else:
            _add_group_col(spec.x)
    else:
        _add_group_col(spec.x)

    if not group_cols:
        agg_row: dict[str, Any] = {
            spec.x: float(df[spec.x].median()),
            spec.y: float(df[spec.y].median()),
        }
        if spec.hollow_non_converged and "converged" in df.columns:
            agg_row["converged"] = not any(
                _is_explicitly_false(v) for v in df["converged"].tolist()
            )
        return pd.DataFrame([agg_row])

    work = df.copy()
    key_cols: list[str] = []
    for idx, col in enumerate(group_cols):
        kcol = f"__agg_key_{idx}"
        work[kcol] = work[col].map(_canonical_key)
        key_cols.append(kcol)

    records: list[dict[str, Any]] = []
    for _, grp in work.groupby(key_cols, sort=False, dropna=False):
        rec: dict[str, Any] = {col: grp[col].iloc[0] for col in group_cols}
        if spec.x not in group_cols:
            rec[spec.x] = float(grp[spec.x].median())
        rec[spec.y] = float(grp[spec.y].median())
        if (
            spec.hollow_non_converged
            and "converged" in grp.columns
            and "converged" not in group_cols
        ):
            rec["converged"] = not any(
                _is_explicitly_false(v) for v in grp["converged"].tolist()
            )
        records.append(rec)

    return pd.DataFrame(records)


def _draw_line_panel(
    ax: Axes,
    panel_df: pd.DataFrame,
    spec: LinePlotSpec,
    series_cols: Sequence[str],
    visual_map: dict[tuple[tuple[str, Any], ...], _SeriesStyle],
) -> int:
    """Draw curves, optional reference slopes, and a per-subplot legend onto ``ax``.

    Returns the number of data curves drawn on ``ax``.
    """
    if panel_df.empty:
        return 0

    panel_combos = _extract_sorted_combos(panel_df, series_cols)
    key_series_list = [panel_df[col].map(_canonical_key) for col in series_cols]

    legend_handles: list[Any] = []
    legend_labels: list[str] = []
    seen_labels: set[str] = set()
    has_non_converged = False
    curves_drawn = 0

    for key_tup, val_tup in panel_combos:
        mask = pd.Series(True, index=panel_df.index)
        for s_series, k_val in zip(key_series_list, key_tup):
            mask &= s_series == k_val

        curve_df = panel_df.loc[mask].sort_values(by=spec.x, ascending=True)
        if curve_df.empty:
            continue

        style = visual_map.get(
            key_tup,
            _SeriesStyle(
                color=PALETTE[0],
                marker=MARKERS[0],
                linestyle="-",
                label=_format_series_label(series_cols, val_tup),
            ),
        )

        xs = curve_df[spec.x].to_numpy(dtype=float)
        ys = curve_df[spec.y].to_numpy(dtype=float)

        (line,) = ax.plot(
            xs,
            ys,
            color=style.color,
            marker=style.marker,
            linestyle=style.linestyle,
            linewidth=1.8,
            markersize=5.5,
            label=style.label,
            zorder=3,
        )
        curves_drawn += 1

        if style.label not in seen_labels:
            seen_labels.add(style.label)
            legend_handles.append(line)
            legend_labels.append(style.label)

        if spec.hollow_non_converged and "converged" in curve_df.columns:
            nc_mask = np.array(
                [_is_explicitly_false(v) for v in curve_df["converged"].tolist()],
                dtype=bool,
            )
            if np.any(nc_mask):
                has_non_converged = True
                ax.scatter(
                    xs[nc_mask],
                    ys[nc_mask],
                    facecolors="white",
                    edgecolors=style.color,
                    marker=style.marker,
                    s=48,
                    linewidths=1.5,
                    zorder=5,
                )

    if curves_drawn == 0:
        return 0

    ax.set_xscale(spec.xscale)
    ax.set_yscale(spec.yscale)

    if spec.ylim_bottom_zero and spec.yscale == "linear":
        _, ymax = ax.get_ylim()
        if ymax <= 0.0:
            ax.set_ylim(0.0, 1.0)
        else:
            ax.set_ylim(bottom=0.0)

    if spec.xticks_from_data:
        unique_x = sorted(
            {
                float(v)
                for v in panel_df[spec.x].dropna().tolist()
                if math.isfinite(float(v))
            }
        )
        if unique_x:
            ax.set_xticks(unique_x)
            ax.set_xticklabels([_fmt_val(v) for v in unique_x])

    # Overlay reference slopes on log-log plots after limits are established
    for ref_label, ref_handle in _draw_reference_slopes(ax, panel_df, spec):
        if ref_label not in seen_labels:
            seen_labels.add(ref_label)
            legend_handles.append(ref_handle)
            legend_labels.append(ref_label)

    if has_non_converged:
        legend_handles.append(
            Line2D(
                [],
                [],
                color="#444444",
                marker="o",
                linestyle="None",
                markerfacecolor="white",
                markeredgewidth=1.5,
                markersize=5.5,
            )
        )
        legend_labels.append("Non-converged (hollow)")

    ax.grid(True, which="major", alpha=0.35, linestyle="--")
    if spec.xscale == "log" or spec.yscale == "log":
        ax.grid(True, which="minor", alpha=0.15, linestyle=":")

    if legend_handles:
        ax.legend(
            legend_handles,
            legend_labels,
            loc="best",
            frameon=True,
            framealpha=0.88,
        )

    return curves_drawn


def render_faceted_line_plot(
    dataset: BenchmarkDataset,
    spec: LinePlotSpec,
    plots_dir: Path,
    formats: Sequence[str] = ("png",),
) -> list[Path]:
    """Render one or more single-axis or multi-panel line figures for ``spec``."""
    df = _filter_ok_rows(dataset.df)
    if df.empty or spec.x not in df.columns or spec.y not in df.columns:
        return []

    df[spec.x] = pd.to_numeric(df[spec.x], errors="coerce")
    df[spec.y] = pd.to_numeric(df[spec.y], errors="coerce")
    valid_mask = np.isfinite(df[spec.x]) & np.isfinite(df[spec.y])
    if spec.xscale == "log":
        valid_mask &= df[spec.x] > 0.0
    if spec.yscale == "log":
        valid_mask &= df[spec.y] > 0.0
    df = df.loc[valid_mask].copy()
    if df.empty:
        return []

    xlabel = spec.xlabel or _humanize_col_name(spec.x)
    ylabel = spec.ylabel or _humanize_col_name(spec.y)
    base_title = spec.title or spec.title_prefix or f"{ylabel} vs. {xlabel}"

    hue_cols = _normalize_hue_cols(spec.hue, df, exclude={spec.x})
    style_col = (
        spec.style
        if (
            spec.style
            and spec.style in df.columns
            and spec.style not in {spec.x, *hue_cols}
            and df[spec.style].notna().any()
            and df[spec.style].dropna().nunique() > 1
        )
        else None
    )
    facet_col = (
        spec.facet_col
        if (
            spec.facet_col
            and spec.facet_col in df.columns
            and spec.facet_col not in {spec.x, *hue_cols, style_col}
            and df[spec.facet_col].notna().any()
        )
        else None
    )

    mapped_axes: set[str] = {spec.x, *hue_cols}
    if spec.x == "num_free_dofs":
        mapped_axes.add("ref_levels")
    for col in (style_col, facet_col):
        if col is not None:
            mapped_axes.add(col)

    effective_agg = {a for a in spec.agg_over if a not in mapped_axes}
    split_axes = resolve_split_axes(
        dataset.varying_axes,
        mapped_axes,
        effective_agg,
        spec.split_by,
        df_columns=df.columns,
    )

    df = _aggregate_line_df(
        df=df,
        spec=spec,
        varying_axes=dataset.varying_axes,
        split_axes=split_axes,
        hue_cols=hue_cols,
        style_col=style_col,
        facet_col=facet_col,
        agg_over=effective_agg,
    )
    if df.empty:
        return []

    min_pts = max(1, int(spec.min_points_per_panel))
    saved_paths: list[Path] = []

    for split_items, sub_df in _iter_split_groups(df, split_axes):
        if sub_df.empty:
            continue

        title_ctx, slug_ctx = _format_split_context(split_items)
        base_name = f"{spec.name}__{slug_ctx}" if slug_ctx else spec.name
        base_path = plots_dir / base_name

        effective_style_col = (
            style_col
            if (style_col is not None and sub_df[style_col].dropna().nunique() > 1)
            else None
        )

        # Build color-first visual map consistently across all panels of this figure
        series_cols, visual_map = _build_series_visual_map(
            sub_df,
            hue_cols=hue_cols,
            style_col=effective_style_col,
            use_linestyles=spec.use_linestyles,
        )

        if facet_col is None:
            if len(sub_df) < min_pts:
                continue

            fig, ax = plt.subplots(figsize=(7.0, 4.8))
            drawn = _draw_line_panel(
                ax=ax,
                panel_df=sub_df,
                spec=spec,
                series_cols=series_cols,
                visual_map=visual_map,
            )
            if drawn == 0:
                plt.close(fig)
                continue

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(_compose_title(base_title, [title_ctx]), pad=10)
            fig.tight_layout()
            saved_paths.extend(save_figure(fig, base_path, formats))
            continue

        # Faceted rendering: NEVER collapse k >= 2 facets onto 1 axes
        facet_vals = sort_axis_values(facet_col, sub_df[facet_col].dropna().tolist())
        facet_keys_series = sub_df[facet_col].map(_canonical_key)

        active_panels: list[tuple[Any, pd.DataFrame]] = []
        for f_val in facet_vals:
            f_key = _canonical_key(f_val)
            panel_df = sub_df.loc[facet_keys_series == f_key]
            if spec.drop_empty_facets:
                if len(panel_df) < min_pts:
                    continue
            elif panel_df.empty:
                continue
            active_panels.append((f_val, panel_df))

        if not active_panels:
            continue

        if len(active_panels) == 1:
            f_val, panel_df = active_panels[0]
            fig, ax = plt.subplots(figsize=(7.0, 4.8))
            drawn = _draw_line_panel(
                ax=ax,
                panel_df=panel_df,
                spec=spec,
                series_cols=series_cols,
                visual_map=visual_map,
            )
            if drawn == 0:
                plt.close(fig)
                continue

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            full_title = _compose_title(
                base_title,
                [title_ctx, f"{facet_col}={_fmt_val(f_val)}"],
            )
            ax.set_title(full_title, pad=10)
            fig.tight_layout()
            saved_paths.extend(save_figure(fig, base_path, formats))
            continue

        n_panels = len(active_panels)
        ncols = min(n_panels, 3)
        nrows = math.ceil(n_panels / ncols)
        fig_w = 5.0 * ncols
        fig_h = 4.1 * nrows + 0.55
        fig, axes_grid = plt.subplots(
            nrows,
            ncols,
            figsize=(fig_w, fig_h),
            squeeze=False,
        )

        any_plotted = False
        for panel_idx, (f_val, panel_df) in enumerate(active_panels):
            row_i, col_i = divmod(panel_idx, ncols)
            ax = axes_grid[row_i][col_i]

            drawn = _draw_line_panel(
                ax=ax,
                panel_df=panel_df,
                spec=spec,
                series_cols=series_cols,
                visual_map=visual_map,
            )
            if drawn > 0:
                any_plotted = True

            ax.set_title(f"{facet_col} = {_fmt_val(f_val)}", pad=8)
            if col_i == 0:
                ax.set_ylabel(ylabel)
            if (panel_idx + ncols) >= n_panels:
                ax.set_xlabel(xlabel)

        for empty_idx in range(n_panels, nrows * ncols):
            row_i, col_i = divmod(empty_idx, ncols)
            axes_grid[row_i][col_i].set_visible(False)

        if not any_plotted:
            plt.close(fig)
            continue

        fig.suptitle(_compose_title(base_title, [title_ctx]), y=0.985)
        top_rect = 0.93 if nrows == 1 else 0.955
        fig.tight_layout(rect=(0.0, 0.0, 1.0, top_rect))
        saved_paths.extend(save_figure(fig, base_path, formats))

    return saved_paths


# ---------------------------------------------------------------------------
# Primitive 2: Stacked Stage Breakdown
# ---------------------------------------------------------------------------


def _format_series_bar_label(
    active_series: Sequence[str],
    series_vals: Sequence[Any],
) -> str:
    """Build a compact label displayed above a stacked bar."""
    if not active_series:
        return ""
    if len(active_series) == 1:
        axis = active_series[0]
        val_str = _fmt_val(series_vals[0])
        if axis in _NAKED_PRIMARY_AXES or axis in _NAKED_SECONDARY_AXES:
            return val_str
        return f"{axis}={val_str}"
    parts: list[str] = []
    for axis, val in zip(active_series, series_vals):
        val_str = _fmt_val(val)
        if axis in _NAKED_PRIMARY_AXES or axis in _NAKED_SECONDARY_AXES:
            parts.append(val_str)
        else:
            parts.append(f"{axis}={val_str}")
    return "/".join(parts)


def _draw_stacked_panel(
    ax: Axes,
    panel_df: pd.DataFrame,
    spec: StackedBarSpec,
    active_stages: Sequence[tuple[str, str, str]],
    active_series: Sequence[str],
    ylabel: str,
) -> tuple[bool, list[float], int]:
    """Draw a grouped stacked-bar panel onto ``ax`` with its own legend.

    Returns ``(any_bar_drawn, positive_bar_totals, n_series)``.
    """
    if panel_df.empty:
        return False, [], 0

    x_vals = sort_axis_values(spec.x, panel_df[spec.x].dropna().tolist())
    if not x_vals:
        return False, [], 0

    if active_series:
        series_combos = [
            val_tup for _, val_tup in _extract_sorted_combos(panel_df, active_series)
        ]
    else:
        series_combos = [()]

    n_x = len(x_vals)
    n_series = len(series_combos)

    x_key_series = panel_df[spec.x].map(_canonical_key)
    series_key_cols = [panel_df[c].map(_canonical_key) for c in active_series]

    dofs_by_x: dict[tuple[str, Any], float] = {}
    if (
        spec.annotate_dofs
        and spec.x == "ref_levels"
        and "num_free_dofs" in panel_df.columns
    ):
        dofs_num = pd.to_numeric(panel_df["num_free_dofs"], errors="coerce")
        for xv in x_vals:
            xk = _canonical_key(xv)
            vals = dofs_num.loc[x_key_series == xk].dropna()
            if not vals.empty:
                dofs_by_x[xk] = float(vals.median())

    group_width = 0.78
    bar_width = min(0.45, group_width / max(n_series, 1))

    labeled_stages: set[str] = set()
    max_bar_height = 0.0
    positive_totals: list[float] = []
    any_bar_drawn = False

    for xi, xv in enumerate(x_vals):
        xk = _canonical_key(xv)
        x_mask = x_key_series == xk

        for si, s_vals in enumerate(series_combos):
            mask = x_mask
            for s_series, sv in zip(series_key_cols, s_vals):
                mask = mask & (s_series == _canonical_key(sv))

            bar_df = panel_df.loc[mask]
            if bar_df.empty:
                continue

            stage_values: list[tuple[str, str, str, float]] = []
            total = 0.0
            for col, stage_label, stage_color in active_stages:
                med = bar_df[col].median()
                val = float(med) if is_valid_number(med) and float(med) > 0.0 else 0.0
                if val > 0.0:
                    stage_values.append((col, stage_label, stage_color, val))
                    total += val

            if total <= 0.0:
                continue

            xpos = xi + (si - (n_series - 1) / 2.0) * bar_width
            bottom = 0.0
            scale = (100.0 / total) if spec.normalize else 1.0

            for col, stage_label, stage_color, raw_val in stage_values:
                height = raw_val * scale
                lbl = stage_label if col not in labeled_stages else "_nolegend_"
                ax.bar(
                    xpos,
                    height,
                    width=bar_width * 0.9,
                    bottom=bottom,
                    color=stage_color,
                    edgecolor="white",
                    linewidth=0.5,
                    label=lbl,
                    zorder=3,
                )
                labeled_stages.add(col)
                bottom += height

            any_bar_drawn = True
            positive_totals.append(total)
            max_bar_height = max(max_bar_height, bottom)

            bar_lbl = _format_series_bar_label(active_series, s_vals)
            total_str = _fmt_bar_total(total)
            rotate = n_series > 3
            y_text = (101.2 if spec.normalize else bottom * 1.015)

            if bar_lbl:
                ann_text = f" {bar_lbl} ({total_str})" if rotate else f"{bar_lbl}\n{total_str}"
            else:
                ann_text = f" {total_str}" if rotate else total_str

            ax.text(
                xpos,
                y_text,
                ann_text,
                ha="left" if rotate else "center",
                va="bottom",
                rotation=40 if rotate else 0,
                fontsize=7.0 if rotate else 7.6,
                color="#222222",
            )

    if not any_bar_drawn:
        return False, [], n_series

    xtick_labels: list[str] = []
    for xv in x_vals:
        xk = _canonical_key(xv)
        if spec.x == "ref_levels":
            base_lbl = f"ref {_fmt_val(xv)}"
            if xk in dofs_by_x:
                base_lbl += f"\n({fmt_int_like(dofs_by_x[xk], thousands=True)} DoFs)"
            xtick_labels.append(base_lbl)
        else:
            xtick_labels.append(_fmt_val(xv))

    ax.set_xticks(range(n_x))
    ax.set_xticklabels(xtick_labels)
    ax.set_xlabel(_humanize_col_name(spec.x))
    ax.set_ylabel(ylabel)

    if spec.normalize:
        ax.set_ylim(0.0, 122.0)
        ax.set_yticks([0, 20, 40, 60, 80, 100])
    else:
        ax.set_ylim(0.0, max_bar_height * (1.34 if n_series > 3 else 1.22))

    ax.grid(True, axis="y", alpha=0.35, linestyle="--")
    ax.grid(False, axis="x")

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(
            handles,
            labels,
            loc="best",
            frameon=True,
            framealpha=0.88,
            fontsize=8.0,
        )

    return True, positive_totals, n_series


def render_stacked_breakdown(
    dataset: BenchmarkDataset,
    spec: StackedBarSpec,
    plots_dir: Path,
    formats: Sequence[str] = ("png",),
) -> list[Path]:
    """Render grouped stacked-bar timing breakdowns for ``spec``."""
    df = _filter_ok_rows(dataset.df)
    stages_to_use = list(spec.stages) if spec.stages else list(DEFAULT_TIMING_STAGES)
    if df.empty or spec.x not in df.columns or not stages_to_use:
        return []

    present_stages = [
        (col, label, color)
        for col, label, color in stages_to_use
        if col in df.columns
    ]
    if not present_stages:
        return []

    for col, _, _ in present_stages:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    facet_col = (
        spec.facet_col
        if (
            spec.facet_col
            and spec.facet_col in df.columns
            and df[spec.facet_col].notna().any()
        )
        else None
    )

    active_series: list[str] = []
    for s_col in spec.series:
        if (
            s_col
            and s_col in df.columns
            and s_col != spec.x
            and s_col != facet_col
            and s_col not in active_series
            and df[s_col].notna().any()
        ):
            active_series.append(s_col)

    mapped_axes: set[str] = {spec.x, *active_series}
    if facet_col is not None:
        mapped_axes.add(facet_col)
    if spec.x == "num_free_dofs" or facet_col == "num_free_dofs":
        mapped_axes.add("ref_levels")

    effective_agg = {a for a in spec.agg_over if a not in mapped_axes}

    split_axes = resolve_split_axes(
        dataset.varying_axes,
        mapped_axes,
        effective_agg,
        spec.split_by,
        df_columns=df.columns,
    )

    ylabel = (
        spec.ylabel
        if spec.ylabel
        else ("Share of wall time (%)" if spec.normalize else "Time (s)")
    )
    base_title = spec.title or spec.title_prefix or "Timing Breakdown"

    saved_paths: list[Path] = []

    for split_items, sub_df in _iter_split_groups(df, split_axes):
        if sub_df.empty:
            continue

        # Dynamically drop stages that are all missing or <= _STAGE_EPS
        active_stages: list[tuple[str, str, str]] = []
        for col, label, color in present_stages:
            col_vals = sub_df[col].dropna()
            if not col_vals.empty and (col_vals > _STAGE_EPS).any():
                active_stages.append((col, label, color))
        if not active_stages:
            continue

        title_ctx, slug_ctx = _format_split_context(split_items)
        base_name = f"{spec.name}__{slug_ctx}" if slug_ctx else spec.name
        base_path = plots_dir / base_name

        # Build panel list if facet_col is set
        active_panels: list[tuple[Any, pd.DataFrame]] = []
        if facet_col is not None:
            facet_vals = sort_axis_values(facet_col, sub_df[facet_col].dropna().tolist())
            facet_keys_series = sub_df[facet_col].map(_canonical_key)
            for f_val in facet_vals:
                panel_df = sub_df.loc[facet_keys_series == _canonical_key(f_val)]
                if panel_df.empty and spec.drop_empty_facets:
                    continue
                if not panel_df.empty:
                    active_panels.append((f_val, panel_df))

        if facet_col is None or len(active_panels) <= 1:
            panel_df = active_panels[0][1] if len(active_panels) == 1 else sub_df
            extra_ctx = (
                [title_ctx, f"{facet_col}={_fmt_val(active_panels[0][0])}"]
                if len(active_panels) == 1
                else [title_ctx]
            )

            n_x = len(sort_axis_values(spec.x, panel_df[spec.x].dropna().tolist()))
            n_series_est = (
                len(_extract_sorted_combos(panel_df, active_series))
                if active_series
                else 1
            )
            fig_w = max(7.0, min(15.0, 2.2 + n_x * max(1.25, 0.65 * n_series_est)))
            fig_h = 5.2 if n_series_est > 3 else 4.8
            fig, ax = plt.subplots(figsize=(fig_w, fig_h))

            drawn, pos_totals, _ = _draw_stacked_panel(
                ax=ax,
                panel_df=panel_df,
                spec=spec,
                active_stages=active_stages,
                active_series=active_series,
                ylabel=ylabel,
            )
            if not drawn:
                plt.close(fig)
                continue

            if not spec.normalize and len(pos_totals) >= 2:
                ratio = max(pos_totals) / min(pos_totals)
                if ratio > _DWARF_RATIO_THRESHOLD:
                    print(
                        f"[warn] Stacked bar plot '{spec.name}' has a {ratio:.1f}x "
                        f"bar-height ratio on a single linear y-axis; smaller bars may "
                        f"be dwarfed. Consider setting 'normalize = true' or "
                        f"'facet_col = \"{spec.x}\"'.",
                        file=sys.stderr,
                    )

            ax.set_title(_compose_title(base_title, extra_ctx), pad=10)
            fig.tight_layout()
            saved_paths.extend(save_figure(fig, base_path, formats))
            continue

        # Multi-panel faceted stacked bars with independent y-axes (sharey=False)
        n_panels = len(active_panels)
        ncols = min(n_panels, 3)
        nrows = math.ceil(n_panels / ncols)
        fig_w = 5.2 * ncols
        fig_h = 4.3 * nrows + 0.55
        fig, axes_grid = plt.subplots(
            nrows,
            ncols,
            figsize=(fig_w, fig_h),
            squeeze=False,
            sharey=False,
        )

        any_plotted = False
        for panel_idx, (f_val, panel_df) in enumerate(active_panels):
            row_i, col_i = divmod(panel_idx, ncols)
            ax = axes_grid[row_i][col_i]

            drawn, _, _ = _draw_stacked_panel(
                ax=ax,
                panel_df=panel_df,
                spec=spec,
                active_stages=active_stages,
                active_series=active_series,
                ylabel=ylabel if (col_i == 0 or not spec.normalize) else "",
            )
            if drawn:
                any_plotted = True
            ax.set_title(f"{facet_col} = {_fmt_val(f_val)}", pad=8)

        for empty_idx in range(n_panels, nrows * ncols):
            row_i, col_i = divmod(empty_idx, ncols)
            axes_grid[row_i][col_i].set_visible(False)

        if not any_plotted:
            plt.close(fig)
            continue

        fig.suptitle(_compose_title(base_title, [title_ctx]), y=0.985)
        top_rect = 0.93 if nrows == 1 else 0.955
        fig.tight_layout(rect=(0.0, 0.0, 1.0, top_rect))
        saved_paths.extend(save_figure(fig, base_path, formats))

    return saved_paths


# ---------------------------------------------------------------------------
# Primitive 3: Residual Convergence Plot
# ---------------------------------------------------------------------------


def _extract_relative_history(raw_hist: Sequence[float]) -> np.ndarray | None:
    """Convert a residual history into positive relative residuals ``||r_k|| / ||r_0||``."""
    if not raw_hist:
        return None
    arr = np.asarray(raw_hist, dtype=float)
    if arr.size == 0 or not np.isfinite(arr[0]) or arr[0] <= 0.0:
        pos_idx = np.where(np.isfinite(arr) & (arr > 0.0))[0]
        if pos_idx.size == 0:
            return None
        r0 = float(arr[pos_idx[0]])
    else:
        r0 = float(arr[0])

    rel = arr / r0
    rel = np.where(np.isfinite(rel) & (rel > 0.0), rel, np.nan)
    if not np.any(np.isfinite(rel)):
        return None
    return rel


def _draw_convergence_panel(
    ax: Axes,
    panel_df: pd.DataFrame,
    residual_histories: dict[str, list[float]],
    series_cols: Sequence[str],
    visual_map: dict[tuple[tuple[str, Any], ...], _SeriesStyle],
) -> bool:
    """Draw semilogy residual curves and a per-subplot legend onto ``ax``."""
    if panel_df.empty:
        return False

    combos = _extract_sorted_combos(panel_df, series_cols)
    key_series_list = [panel_df[col].map(_canonical_key) for col in series_cols]

    plotted_any = False
    for key_tup, val_tup in combos:
        mask = pd.Series(True, index=panel_df.index)
        for s_series, k_val in zip(key_series_list, key_tup):
            mask &= s_series == k_val

        matches = panel_df.loc[mask]
        if matches.empty:
            continue

        rid = str(matches["run_id"].iloc[0])
        rel = _extract_relative_history(residual_histories[rid])
        if rel is None:
            continue

        style = visual_map.get(
            key_tup,
            _SeriesStyle(
                color=PALETTE[0],
                marker=MARKERS[0],
                linestyle="-",
                label=_format_series_label(series_cols, val_tup),
            ),
        )

        iters = np.arange(rel.size, dtype=int)
        markevery = max(1, rel.size // 12)
        ax.semilogy(
            iters,
            rel,
            color=style.color,
            marker=style.marker,
            linestyle=style.linestyle,
            markevery=markevery,
            linewidth=1.7,
            markersize=5.0,
            label=style.label,
            zorder=3,
        )
        plotted_any = True

    if not plotted_any:
        return False

    if "rtol" in panel_df.columns:
        rtol_vals = pd.to_numeric(panel_df["rtol"], errors="coerce").dropna()
        rtol_pos = rtol_vals.loc[rtol_vals > 0.0]
        unique_rtols = sorted({float(v) for v in rtol_pos.tolist()})
        if len(unique_rtols) == 1:
            rtol_val = unique_rtols[0]
            ax.axhline(
                rtol_val,
                color="#555555",
                linestyle=":",
                linewidth=1.3,
                label=f"rtol = {rtol_val:.1e}",
                zorder=2,
            )

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(True, which="major", alpha=0.35, linestyle="--")
    ax.grid(True, which="minor", alpha=0.15, linestyle=":")
    ax.legend(loc="best", frameon=True, framealpha=0.88)
    return True


def render_convergence_plot(
    dataset: BenchmarkDataset,
    spec: ConvergencePlotSpec,
    plots_dir: Path,
    formats: Sequence[str] = ("png",),
) -> list[Path]:
    """Render semilogy residual convergence curves for the finest ``ref_levels`` per split."""
    if (
        dataset.df.empty
        or not dataset.residual_histories
        or "run_id" not in dataset.df.columns
    ):
        return []

    df = _filter_ok_rows(dataset.df)
    if df.empty:
        return []

    def _has_valid_hist(rid: Any) -> bool:
        if _is_missing(rid):
            return False
        hist = dataset.residual_histories.get(str(rid))
        return hist is not None and _extract_relative_history(hist) is not None

    df = df.loc[df["run_id"].map(_has_valid_hist)].copy()
    if df.empty:
        return []

    hue_cols = _normalize_hue_cols(spec.hue, df)
    style_col = (
        spec.style
        if (
            spec.style
            and spec.style in df.columns
            and spec.style not in hue_cols
            and df[spec.style].notna().any()
            and df[spec.style].dropna().nunique() > 1
        )
        else None
    )
    facet_col = (
        spec.facet_col
        if (
            spec.facet_col
            and spec.facet_col in df.columns
            and spec.facet_col not in {*hue_cols, style_col}
            and df[spec.facet_col].notna().any()
        )
        else None
    )

    mapped_axes: set[str] = set(hue_cols)
    for col in (style_col, facet_col):
        if col is not None:
            mapped_axes.add(col)

    auto_select_finest = (
        "ref_levels" in df.columns
        and "ref_levels" not in spec.split_by
        and "ref_levels" not in mapped_axes
    )
    agg_over: set[str] = {"ref_levels"} if auto_select_finest else set()

    split_axes = resolve_split_axes(
        dataset.varying_axes,
        mapped_axes,
        agg_over,
        spec.split_by,
        df_columns=df.columns,
    )

    base_title = spec.title or spec.title_prefix or "Residual Convergence"
    saved_paths: list[Path] = []

    for split_items, sub_df in _iter_split_groups(df, split_axes):
        if sub_df.empty:
            continue

        selected_ref_str: str | None = None
        selected_dofs_str: str | None = None

        if auto_select_finest:
            ref_num = pd.to_numeric(sub_df["ref_levels"], errors="coerce")
            if ref_num.notna().any():
                max_ref = ref_num.max()
                sub_df = sub_df.loc[ref_num == max_ref].copy()
                selected_ref_str = _fmt_val(max_ref)
                if "num_free_dofs" in sub_df.columns:
                    dofs_vals = pd.to_numeric(
                        sub_df["num_free_dofs"], errors="coerce"
                    ).dropna()
                    if not dofs_vals.empty:
                        selected_dofs_str = fmt_int_like(
                            float(dofs_vals.median()), thousands=True
                        )

        if sub_df.empty:
            continue

        sub_style_col = (
            style_col
            if (style_col is not None and sub_df[style_col].dropna().nunique() > 1)
            else None
        )
        series_cols, visual_map = _build_series_visual_map(
            sub_df,
            hue_cols=hue_cols,
            style_col=sub_style_col,
            use_linestyles=spec.use_linestyles,
        )

        title_ctx, slug_ctx = _format_split_context(split_items)
        ctx_parts: list[str] = []
        if title_ctx:
            ctx_parts.append(title_ctx)
        if selected_ref_str is not None:
            ref_info = f"ref={selected_ref_str}"
            if selected_dofs_str is not None:
                ref_info += f", {selected_dofs_str} DoFs"
            ctx_parts.append(ref_info)

        base_name = f"{spec.name}__{slug_ctx}" if slug_ctx else spec.name
        base_path = plots_dir / base_name

        active_panels: list[tuple[Any, pd.DataFrame]] = []
        if facet_col is not None:
            facet_vals = sort_axis_values(facet_col, sub_df[facet_col].dropna().tolist())
            facet_keys_series = sub_df[facet_col].map(_canonical_key)
            for f_val in facet_vals:
                panel_df = sub_df.loc[facet_keys_series == _canonical_key(f_val)]
                if not panel_df.empty:
                    active_panels.append((f_val, panel_df))

        if facet_col is None or len(active_panels) <= 1:
            panel_df = active_panels[0][1] if len(active_panels) == 1 else sub_df
            if len(active_panels) == 1:
                ctx_parts.append(f"{facet_col}={_fmt_val(active_panels[0][0])}")

            fig, ax = plt.subplots(figsize=(7.0, 4.8))
            plotted = _draw_convergence_panel(
                ax=ax,
                panel_df=panel_df,
                residual_histories=dataset.residual_histories,
                series_cols=series_cols,
                visual_map=visual_map,
            )
            if not plotted:
                plt.close(fig)
                continue

            ax.set_xlabel("Iteration $k$")
            ax.set_ylabel(r"Relative Residual $\|r_k\| / \|r_0\|$")
            ax.set_title(_compose_title(base_title, ctx_parts), pad=10)
            fig.tight_layout()
            saved_paths.extend(save_figure(fig, base_path, formats))
            continue

        n_panels = len(active_panels)
        ncols = min(n_panels, 3)
        nrows = math.ceil(n_panels / ncols)
        fig_w = 5.0 * ncols
        fig_h = 4.1 * nrows + 0.55
        fig, axes_grid = plt.subplots(
            nrows,
            ncols,
            figsize=(fig_w, fig_h),
            squeeze=False,
        )

        any_plotted = False
        for panel_idx, (f_val, panel_df) in enumerate(active_panels):
            row_i, col_i = divmod(panel_idx, ncols)
            ax = axes_grid[row_i][col_i]
            plotted = _draw_convergence_panel(
                ax=ax,
                panel_df=panel_df,
                residual_histories=dataset.residual_histories,
                series_cols=series_cols,
                visual_map=visual_map,
            )
            if plotted:
                any_plotted = True
            ax.set_title(f"{facet_col} = {_fmt_val(f_val)}", pad=8)
            if col_i == 0:
                ax.set_ylabel(r"Relative Residual $\|r_k\| / \|r_0\|$")
            if (panel_idx + ncols) >= n_panels:
                ax.set_xlabel("Iteration $k$")

        for empty_idx in range(n_panels, nrows * ncols):
            row_i, col_i = divmod(empty_idx, ncols)
            axes_grid[row_i][col_i].set_visible(False)

        if not any_plotted:
            plt.close(fig)
            continue

        fig.suptitle(_compose_title(base_title, ctx_parts), y=0.985)
        top_rect = 0.93 if nrows == 1 else 0.955
        fig.tight_layout(rect=(0.0, 0.0, 1.0, top_rect))
        saved_paths.extend(save_figure(fig, base_path, formats))

    return saved_paths


# ---------------------------------------------------------------------------
# Custom TOML [[plots]] Parsing & Per-Plot Filtering
# ---------------------------------------------------------------------------


def _as_str_list(val: Any) -> list[str]:
    """Coerce a TOML scalar or list into a list of strings."""
    if val is None:
        return []
    if isinstance(val, str):
        return [val]
    if isinstance(val, (list, tuple, set)):
        return [str(v) for v in val]
    raise ValueError(f"Expected string or list of strings, got {type(val).__name__}")


def _as_hue_spec(val: Any) -> str | list[str] | None:
    """Coerce a TOML ``hue`` field into ``str | list[str] | None``."""
    if val is None:
        return None
    if isinstance(val, str):
        return val if val.strip() else None
    if isinstance(val, (list, tuple)):
        items = [str(v).strip() for v in val if str(v).strip()]
        return items if items else None
    raise ValueError(
        f"Expected 'hue' to be a string or list of strings, got {type(val).__name__}"
    )


def _as_float_list(val: Any) -> list[float]:
    """Coerce a TOML scalar or list into a list of finite floats."""
    if val is None:
        return []
    raw_items = [val] if isinstance(val, (int, float)) else list(val)
    out: list[float] = []
    for item in raw_items:
        fval = float(item)
        if not math.isfinite(fval):
            raise ValueError(f"Non-finite slope value: {item!r}")
        out.append(fval)
    return out


def _fmt_filter_val(val: Any) -> str:
    """Serialize a TOML scalar filter value into the ``apply_filters`` grammar."""
    if isinstance(val, (bool, np.bool_)):
        return "true" if bool(val) else "false"
    if isinstance(val, (int, np.integer)):
        return str(int(val))
    if isinstance(val, (float, np.floating)):
        fval = float(val)
        if math.isfinite(fval) and fval.is_integer():
            return str(int(fval))
        return str(fval)
    return str(val)


def _filter_dict_to_specs(raw_filter: Any) -> list[str]:
    """Convert a TOML ``filter`` dict (or list of strings) into ``apply_filters`` specs."""
    if raw_filter is None:
        return []
    if isinstance(raw_filter, list):
        return [str(item) for item in raw_filter if str(item).strip()]
    if not isinstance(raw_filter, dict):
        raise ValueError(
            f"'filter' must be a table/dict or list of filter strings, "
            f"got {type(raw_filter).__name__}"
        )

    specs: list[str] = []
    for key, val in raw_filter.items():
        key_str = str(key).strip()
        if not key_str:
            raise ValueError("Empty key in plot 'filter' specification")
        if isinstance(val, (list, tuple, set)):
            vals = [_fmt_filter_val(v) for v in val]
            if not vals:
                raise ValueError(f"Empty filter list for key '{key_str}'")
            specs.append(f"{key_str}={','.join(vals)}")
        else:
            specs.append(f"{key_str}={_fmt_filter_val(val)}")
    return specs


def _apply_plot_filter(
    dataset: BenchmarkDataset,
    raw_filter: Any,
    spec_label: str,
) -> BenchmarkDataset | None:
    """Apply a per-plot filter dict, returning ``None`` cleanly when 0 rows remain."""
    filter_specs = _filter_dict_to_specs(raw_filter)
    if not filter_specs:
        return dataset

    try:
        filtered = apply_filters(dataset, filter_specs)
    except ValueError as exc:
        msg = str(exc).lower()
        if any(tok in msg for tok in ("0 rows", "no rows", "empty", "matched no")):
            print(
                f"[warn] Skipping '{spec_label}': per-plot filter {filter_specs} matched 0 rows.",
                file=sys.stderr,
            )
            return None
        raise

    if filtered.df.empty:
        print(
            f"[warn] Skipping '{spec_label}': per-plot filter {filter_specs} matched 0 rows.",
            file=sys.stderr,
        )
        return None

    return filtered


def _parse_custom_stages(raw_stages: Any) -> list[tuple[str, str, str]]:
    """Parse declarative ``stages`` definitions from a TOML custom stacked plot."""
    if raw_stages is None:
        return list(DEFAULT_TIMING_STAGES)
    if not isinstance(raw_stages, list) or not raw_stages:
        raise ValueError("'stages' must be a non-empty list")

    parsed: list[tuple[str, str, str]] = []
    for idx, entry in enumerate(raw_stages):
        default_color = PALETTE[idx % len(PALETTE)]
        if isinstance(entry, str):
            parsed.append((entry, _humanize_col_name(entry), default_color))
        elif isinstance(entry, (list, tuple)):
            if len(entry) == 2:
                col, label = str(entry[0]), str(entry[1])
                parsed.append((col, label, default_color))
            elif len(entry) == 3:
                col, label, color = str(entry[0]), str(entry[1]), str(entry[2])
                parsed.append((col, label, color))
            else:
                raise ValueError(
                    f"Stage entry #{idx + 1} must have 2 or 3 elements, got {len(entry)}"
                )
        elif isinstance(entry, dict):
            col_val = entry.get("col") or entry.get("column")
            if not col_val:
                raise ValueError(f"Stage dict #{idx + 1} is missing 'col'")
            col = str(col_val)
            label = str(entry.get("label") or _humanize_col_name(col))
            color = str(entry.get("color") or default_color)
            parsed.append((col, label, color))
        else:
            raise ValueError(
                f"Unsupported stage entry type at index {idx}: {type(entry).__name__}"
            )
    return parsed


def _resolve_tables_dir(plots_dir: Path) -> Path:
    """Choose the sibling ``tables/`` directory when ``plots_dir`` is named ``plots``."""
    if plots_dir.name == "plots":
        return plots_dir.parent / "tables"
    return plots_dir / "tables"


def _render_custom_plot_spec(
    dataset: BenchmarkDataset,
    raw_spec: dict[str, Any],
    idx: int,
    plots_dir: Path,
    formats: Sequence[str],
) -> list[Path]:
    """Validate, filter, and render a single custom TOML ``[[plots]]`` entry."""
    if not isinstance(raw_spec, dict):
        raise ValueError(
            f"Plot spec must be a table/dict, got {type(raw_spec).__name__}"
        )

    kind = str(raw_spec.get("kind", "")).strip().lower()
    if kind not in {"line", "stacked", "convergence", "table"}:
        raise ValueError(
            f"Unsupported or missing plot kind '{kind}' "
            "(expected 'line', 'stacked', 'convergence', or 'table')"
        )

    name = str(raw_spec.get("name") or f"custom_{kind}_{idx + 1}")
    plot_dataset = _apply_plot_filter(dataset, raw_spec.get("filter"), name)
    if plot_dataset is None:
        return []

    if kind == "table":
        if parse_sparse_table_spec is None or render_sparse_table is None:
            raise ImportError(
                "analysis.sparse_tables is required to render kind='table' specs."
            )
        table_raw = {
            k: v for k, v in raw_spec.items() if k not in {"filter", "formats"}
        }
        table_raw.setdefault("name", name)
        table_spec = parse_sparse_table_spec(table_raw)
        table_formats = (
            _as_str_list(raw_spec["formats"])
            if raw_spec.get("formats") is not None
            else None
        )
        tables_dir = _resolve_tables_dir(plots_dir)
        return render_sparse_table(
            plot_dataset,
            table_spec,
            tables_dir,
            formats=table_formats,
        )

    split_by = _as_str_list(raw_spec.get("split_by"))
    agg_over = set(_as_str_list(raw_spec.get("agg_over")))

    if kind == "line":
        x = raw_spec.get("x")
        y = raw_spec.get("y")
        if not x or not y:
            raise ValueError("Line plot spec requires both 'x' and 'y'")
        x_str, y_str = str(x), str(y)
        spec = LinePlotSpec(
            name=name,
            x=x_str,
            y=y_str,
            xlabel=str(raw_spec.get("xlabel") or _humanize_col_name(x_str)),
            ylabel=str(raw_spec.get("ylabel") or _humanize_col_name(y_str)),
            title=str(raw_spec.get("title") or ""),
            title_prefix=str(
                raw_spec.get("title_prefix")
                or f"{_humanize_col_name(y_str)} vs. {_humanize_col_name(x_str)}"
            ),
            xscale=str(raw_spec.get("xscale", "linear")),
            yscale=str(raw_spec.get("yscale", "linear")),
            hue=_as_hue_spec(raw_spec.get("hue")),
            style=str(raw_spec["style"]) if raw_spec.get("style") else None,
            facet_col=(
                str(raw_spec["facet_col"]) if raw_spec.get("facet_col") else None
            ),
            split_by=split_by,
            agg_over=agg_over,
            hollow_non_converged=bool(raw_spec.get("hollow_non_converged", False)),
            ylim_bottom_zero=bool(raw_spec.get("ylim_bottom_zero", False)),
            xticks_from_data=bool(raw_spec.get("xticks_from_data", False)),
            use_linestyles=bool(raw_spec.get("use_linestyles", False)),
            min_points_per_panel=int(raw_spec.get("min_points_per_panel", 1)),
            drop_empty_facets=bool(raw_spec.get("drop_empty_facets", True)),
            reference_slopes=_as_float_list(raw_spec.get("reference_slopes")),
        )
        return render_faceted_line_plot(plot_dataset, spec, plots_dir, formats)

    if kind == "stacked":
        stages = _parse_custom_stages(raw_spec.get("stages"))
        raw_ylabel = raw_spec.get("ylabel")
        spec_stacked = StackedBarSpec(
            name=name,
            x=str(raw_spec.get("x") or "ref_levels"),
            series=_as_str_list(raw_spec.get("series")),
            stages=stages,
            ylabel=str(raw_ylabel) if raw_ylabel else None,
            title=str(raw_spec.get("title") or ""),
            title_prefix=str(
                raw_spec.get("title_prefix") or _humanize_col_name(name)
            ),
            split_by=split_by,
            agg_over=agg_over,
            annotate_dofs=bool(raw_spec.get("annotate_dofs", True)),
            normalize=bool(raw_spec.get("normalize", True)),
            facet_col=(
                str(raw_spec["facet_col"]) if raw_spec.get("facet_col") else None
            ),
            drop_empty_facets=bool(raw_spec.get("drop_empty_facets", True)),
        )
        return render_stacked_breakdown(plot_dataset, spec_stacked, plots_dir, formats)

    # kind == "convergence"
    spec_conv = ConvergencePlotSpec(
        name=name,
        hue=_as_hue_spec(raw_spec.get("hue")),
        style=str(raw_spec["style"]) if raw_spec.get("style") else None,
        facet_col=(
            str(raw_spec["facet_col"]) if raw_spec.get("facet_col") else None
        ),
        split_by=split_by,
        title=str(raw_spec.get("title") or ""),
        title_prefix=str(raw_spec.get("title_prefix") or "Residual Convergence"),
        use_linestyles=bool(raw_spec.get("use_linestyles", False)),
    )
    return render_convergence_plot(plot_dataset, spec_conv, plots_dir, formats)


# ---------------------------------------------------------------------------
# Minimal Academic Defaults & Manual-First Orchestrator
# ---------------------------------------------------------------------------


def _render_minimal_defaults(
    dataset: BenchmarkDataset,
    plots_dir: Path,
    export_formats: Sequence[str],
    compare_axis: str | None,
) -> list[Path]:
    """Render a concise set of standard academic default plots."""
    ok_df = _filter_ok_rows(dataset.df)
    if ok_df.empty:
        return []

    varying_set = set(dataset.varying_axes)
    primary_hue: str | None = compare_axis or (
        "prec" if "prec" in ok_df.columns else None
    )
    if "assembly" in varying_set and primary_hue != "assembly":
        primary_style: str | None = "assembly"
    elif "solver" in varying_set and primary_hue != "solver":
        primary_style = "solver"
    else:
        primary_style = None

    primary_facet: str | None = (
        "order"
        if primary_hue != "order"
        else ("prec" if "prec" in varying_set else None)
    )

    generated: list[Path] = []

    # 1. Iterations vs. DoFs
    generated.extend(
        render_faceted_line_plot(
            dataset,
            LinePlotSpec(
                name="iterations_vs_dofs",
                x="num_free_dofs",
                y="iterations",
                xlabel="Free DoFs",
                ylabel="Iterations",
                title_prefix="Solver Iterations vs. DoFs",
                xscale="log",
                yscale="linear",
                hue=primary_hue,
                style=primary_style,
                facet_col=primary_facet,
                split_by=["mesh_name", "ref_strategy"],
                hollow_non_converged=True,
                ylim_bottom_zero=True,
            ),
            plots_dir,
            export_formats,
        )
    )

    # 2. Solve time vs. DoFs (with O(N) reference slope)
    generated.extend(
        render_faceted_line_plot(
            dataset,
            LinePlotSpec(
                name="solve_time_vs_dofs",
                x="num_free_dofs",
                y="solve",
                xlabel="Free DoFs",
                ylabel="Solve Time (s)",
                title_prefix="Solve Time vs. DoFs",
                xscale="log",
                yscale="log",
                hue=primary_hue,
                style=primary_style,
                facet_col=primary_facet,
                split_by=["mesh_name", "ref_strategy"],
                reference_slopes=[1.0],
            ),
            plots_dir,
            export_formats,
        )
    )

    # 3. Normalized stage timing breakdown
    breakdown_series: list[str] = []
    if compare_axis and compare_axis != "ref_levels":
        breakdown_series.append(compare_axis)
    for candidate in ("prec", "assembly", "solver"):
        if (
            candidate in varying_set
            and candidate in ok_df.columns
            and candidate not in breakdown_series
        ):
            breakdown_series.append(candidate)
    if not breakdown_series and "prec" in ok_df.columns:
        breakdown_series.append("prec")

    generated.extend(
        render_stacked_breakdown(
            dataset,
            StackedBarSpec(
                name="timing_breakdown",
                x="ref_levels",
                series=breakdown_series,
                stages=list(DEFAULT_TIMING_STAGES),
                title_prefix="Stage Timing Breakdown",
                split_by=["mesh_name", "ref_strategy", "order"],
                annotate_dofs=True,
                normalize=True,
            ),
            plots_dir,
            export_formats,
        )
    )

    # 4. Residual convergence
    if dataset.residual_histories:
        generated.extend(
            render_convergence_plot(
                dataset,
                ConvergencePlotSpec(
                    name="convergence",
                    hue=primary_hue,
                    style=primary_style,
                    split_by=["mesh_name", "ref_strategy", "order"],
                    title_prefix="Residual Convergence",
                ),
                plots_dir,
                export_formats,
            )
        )

    return generated


def generate_all_plots(
    dataset: BenchmarkDataset,
    plots_dir: Path,
    formats: list[str] | None = None,
    compare_axis: str | None = None,
    include_defaults: bool = False,
) -> list[Path]:
    """Generate plots and tables using a manual-first workflow.

    Behavior:
      - If ``dataset.custom_plot_specs`` is non-empty, renders those specs
        (and renders minimal built-in defaults only if ``include_defaults=True``).
      - If ``dataset.custom_plot_specs`` is empty and ``include_defaults=False``,
        prints an informational hint to stderr and returns ``[]``.
      - If ``include_defaults=True``, renders the minimal academic default suite.

    Args:
        dataset: Loaded and metric-enriched ``BenchmarkDataset``.
        plots_dir: Destination directory for rendered figures.
        formats: Image/vector formats to export (defaults to ``["png"]``).
        compare_axis: Optional DataFrame column name onto which default plot
            ``hue`` and stacked-bar series are pivoted.
        include_defaults: Whether to render built-in minimal default plots.

    Returns:
        List of all generated plot and table file paths.

    Raises:
        ValueError: If ``compare_axis`` is provided but not found in ``dataset.df``.
    """
    if compare_axis is not None and compare_axis not in dataset.df.columns:
        available = ", ".join(sorted(map(str, dataset.df.columns)))
        raise ValueError(
            f"Invalid --compare-axis '{compare_axis}': column not found in dataset. "
            f"Available columns: {available}"
        )

    export_formats = list(formats) if formats else ["png"]
    has_custom = bool(dataset.custom_plot_specs)

    if not has_custom and not include_defaults:
        print(
            "[INFO] No [[plots]] entries found; skipping plot generation. "
            "Pass --plots-config <file.toml> or --include-defaults to generate plots.",
            file=sys.stderr,
        )
        return []

    if dataset.df.empty:
        return []

    generated: list[Path] = []

    if include_defaults:
        generated.extend(
            _render_minimal_defaults(
                dataset=dataset,
                plots_dir=plots_dir,
                export_formats=export_formats,
                compare_axis=compare_axis,
            )
        )

    for idx, raw_spec in enumerate(dataset.custom_plot_specs):
        spec_label = (
            raw_spec.get("name", f"#{idx + 1}")
            if isinstance(raw_spec, dict)
            else f"#{idx + 1}"
        )
        try:
            generated.extend(
                _render_custom_plot_spec(
                    dataset, raw_spec, idx, plots_dir, export_formats
                )
            )
        except Exception as exc:  # noqa: BLE001
            print(
                f"[warn] Skipping malformed custom plot spec '{spec_label}': {exc}",
                file=sys.stderr,
            )

    return generated
