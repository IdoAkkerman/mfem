"""Dynamic ingestion, formatting helpers, and filtering for benchmark results."""

from __future__ import annotations

import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11 fallback
    try:
        import tomli as tomllib  # type: ignore[no-redef]
    except ModuleNotFoundError:
        tomllib = None  # type: ignore[assignment]

METADATA_CONFIG_KEYS: frozenset[str] = frozenset(
    {
        "sweep_name",
        "requested_order",
        "config_order",
        "repeats",
        "verify",
        "record_residuals",
    }
)

CANONICAL_AXIS_ORDER: tuple[str, ...] = (
    "mesh_name",
    "problem",
    "ref_strategy",
    "order",
    "ref_levels",
    "assembly",
    "solver",
    "prec",
)

_TRUE_STRINGS: frozenset[str] = frozenset({"true", "1", "yes", "y", "on"})
_FALSE_STRINGS: frozenset[str] = frozenset({"false", "0", "no", "n", "off"})


@dataclass
class BenchmarkDataset:
    df: pd.DataFrame
    result_dir: Path
    config_cols: list[str]
    varying_axes: list[str]
    constant_params: dict[str, Any]
    timing_cols: list[str]
    residual_histories: dict[str, list[float]]
    env_metadata: dict[str, Any]
    custom_plot_specs: list[dict[str, Any]]
    speedup_cols: list[tuple[str, str]] = field(default_factory=list)


def is_valid_number(val: Any) -> bool:
    """Return True if *val* is a finite non-boolean numeric value."""
    if val is None or isinstance(val, (bool, np.bool_)):
        return False
    if isinstance(val, (int, float, np.integer, np.floating)):
        return math.isfinite(float(val))
    return False


def fmt_int_like(val: Any, thousands: bool = False) -> str:
    """Format an integer-like numeric value cleanly, returning '--' if invalid."""
    if not is_valid_number(val):
        return "--"
    fval = float(val)
    if abs(fval - round(fval)) < 1e-9:
        ival = int(round(fval))
        return f"{ival:,}" if thousands else f"{ival}"
    return f"{fval:,.4g}" if thousands else f"{fval:.4g}"


def fmt_time(val: Any) -> str:
    """Format seconds with 4 decimals (<10s) or 2 decimals (>=10s)."""
    if not is_valid_number(val):
        return "--"
    fval = float(val)
    if abs(fval) < 10.0:
        return f"{fval:.4f}"
    return f"{fval:.2f}"


def fmt_speedup(val: Any) -> str:
    """Format a speedup ratio as e.g. '2.34x', returning '--' if invalid."""
    if not is_valid_number(val):
        return "--"
    return f"{float(val):.2f}x"


def slugify(value: Any) -> str:
    """Convert *value* into a filesystem- and identifier-safe lowercase slug."""
    text = str(value).strip().lower()
    slug = re.sub(r"[^a-z0-9]+", "_", text).strip("_")
    return slug or "value"


def safe_div(num: Any, den: Any) -> float:
    """Safely divide two finite numbers, returning NaN on invalid input or zero."""
    if not is_valid_number(num) or not is_valid_number(den):
        return float("nan")
    den_f = float(den)
    if den_f == 0.0:
        return float("nan")
    res = float(num) / den_f
    return res if math.isfinite(res) else float("nan")


def _parse_bool_token(token: str) -> bool | None:
    low = token.strip().lower()
    if low in _TRUE_STRINGS:
        return True
    if low in _FALSE_STRINGS:
        return False
    return None


def coerce_match_value(raw: str, series: pd.Series) -> Any:
    """Coerce a raw string token to match the runtime type of *series*."""
    token = str(raw).strip()

    if pd.api.types.is_bool_dtype(series):
        parsed_bool = _parse_bool_token(token)
        return parsed_bool if parsed_bool is not None else token

    if pd.api.types.is_integer_dtype(series):
        try:
            return int(token)
        except ValueError:
            try:
                fval = float(token)
                return int(fval) if fval.is_integer() else fval
            except ValueError:
                return token

    if pd.api.types.is_float_dtype(series):
        try:
            return float(token)
        except ValueError:
            return token

    non_null = series.dropna()
    if len(non_null) == 0:
        return token

    sample_vals = non_null.tolist()
    if all(isinstance(v, (bool, np.bool_)) for v in sample_vals):
        parsed_bool = _parse_bool_token(token)
        return parsed_bool if parsed_bool is not None else token

    if all(is_valid_number(v) for v in sample_vals):
        if all(isinstance(v, (int, np.integer)) for v in sample_vals):
            try:
                return int(token)
            except ValueError:
                try:
                    fval = float(token)
                    return int(fval) if fval.is_integer() else fval
                except ValueError:
                    return token
        try:
            return float(token)
        except ValueError:
            return token

    return token


def series_matches_value(series: pd.Series, target: Any) -> pd.Series:
    """Return a boolean Series indicating elements matching *target* robustly."""
    if len(series) == 0:
        return pd.Series(dtype=bool, index=series.index)

    if isinstance(target, (bool, np.bool_)):
        target_b = bool(target)
        if pd.api.types.is_bool_dtype(series):
            return series.fillna(False).astype(bool) == target_b
        return series.map(
            lambda x: isinstance(x, (bool, np.bool_)) and bool(x) == target_b
        ).astype(bool)

    if is_valid_number(target):
        target_f = float(target)
        if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series):
            vals = series.astype(float)
            return pd.Series(
                np.isclose(vals, target_f, rtol=1e-9, atol=1e-12, equal_nan=False),
                index=series.index,
                dtype=bool,
            )

        def _match_num(x: Any) -> bool:
            if is_valid_number(x):
                return math.isclose(float(x), target_f, rel_tol=1e-9, abs_tol=1e-12)
            return False

        return series.map(_match_num).astype(bool)

    target_s = str(target)
    return (series.notna() & (series.astype(str) == target_s)).astype(bool)


def _is_scalar(val: Any) -> bool:
    return val is None or isinstance(val, (str, bool, int, float, np.generic))


def _flatten_dict(
    data: dict[str, Any],
    prefix: str = "",
) -> list[tuple[str, Any]]:
    """Recursively flatten nested dicts with '_' separators."""
    items: list[tuple[str, Any]] = []
    for key, val in data.items():
        full_key = f"{prefix}_{key}" if prefix else str(key)
        if isinstance(val, dict):
            items.extend(_flatten_dict(val, prefix=full_key))
        else:
            items.append((full_key, val))
    return items


def _to_python_scalar(val: Any) -> Any:
    if isinstance(val, np.bool_):
        return bool(val)
    if isinstance(val, np.integer):
        return int(val)
    if isinstance(val, np.floating):
        fval = float(val)
        return int(round(fval)) if math.isfinite(fval) and abs(fval - round(fval)) < 1e-9 and fval.is_integer() else fval
    return val


def compute_axis_partition(
    df: pd.DataFrame,
    config_cols: list[str],
) -> tuple[list[str], dict[str, Any]]:
    """Partition config columns into ordered varying_axes and constant_params."""
    present_config = [c for c in config_cols if c in df.columns]
    has_mesh_name = "mesh_name" in present_config

    excluded = set(METADATA_CONFIG_KEYS)
    if has_mesh_name:
        excluded.add("mesh")

    varying_set: set[str] = set()
    const_raw: dict[str, Any] = {}

    for col in present_config:
        if col in excluded:
            continue
        non_null = df[col].dropna()
        n_unique = int(non_null.nunique())
        if n_unique > 1:
            varying_set.add(col)
        elif n_unique == 1:
            const_raw[col] = _to_python_scalar(non_null.iloc[0])

    canonical_set = set(CANONICAL_AXIS_ORDER)
    varying_axes = [c for c in CANONICAL_AXIS_ORDER if c in varying_set]
    varying_axes.extend(sorted(c for c in varying_set if c not in canonical_set))

    constant_params: dict[str, Any] = {}
    for c in CANONICAL_AXIS_ORDER:
        if c in const_raw:
            constant_params[c] = const_raw[c]
    for c in sorted(k for k in const_raw if k not in canonical_set):
        constant_params[c] = const_raw[c]

    return varying_axes, constant_params


def _load_env_metadata(result_dir: Path) -> dict[str, Any]:
    meta_path = result_dir / "env_metadata.json"
    if not meta_path.is_file():
        return {}
    try:
        loaded = json.loads(meta_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, UnicodeDecodeError) as exc:
        print(f"Warning: failed to read {meta_path}: {exc}", file=sys.stderr)
        return {}
    return loaded if isinstance(loaded, dict) else {}


def _load_custom_plot_specs(result_dir: Path) -> list[dict[str, Any]]:
    toml_path = result_dir / "config_snapshot.toml"
    if not toml_path.is_file():
        return []
    if tomllib is None:
        print(
            f"Warning: tomllib/tomli unavailable; skipping {toml_path}",
            file=sys.stderr,
        )
        return []
    try:
        raw = tomllib.loads(toml_path.read_text(encoding="utf-8"))
    except Exception as exc:
        print(f"Warning: failed to parse {toml_path}: {exc}", file=sys.stderr)
        return []
    plots = raw.get("plots", [])
    if not isinstance(plots, list):
        return []
    return [entry for entry in plots if isinstance(entry, dict)]


def _register_unique(target: list[str], seen: set[str], key: str) -> None:
    if key not in seen:
        seen.add(key)
        target.append(key)


def load_dataset(result_dir: Path) -> BenchmarkDataset:
    """Load and dynamically flatten all run JSON files from *result_dir*."""
    result_dir = Path(result_dir)
    runs_dir = result_dir / "runs"
    if not runs_dir.is_dir():
        raise FileNotFoundError(f"Runs directory not found: {runs_dir}")

    rows: list[dict[str, Any]] = []
    config_cols: list[str] = []
    config_seen: set[str] = set()
    timing_cols: list[str] = []
    timing_seen: set[str] = set()
    residual_histories: dict[str, list[float]] = {}

    for json_path in sorted(runs_dir.glob("*.json")):
        try:
            record = json.loads(json_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError, UnicodeDecodeError) as exc:
            print(f"Warning: skipping unreadable JSON {json_path}: {exc}", file=sys.stderr)
            continue

        if not isinstance(record, dict):
            print(
                f"Warning: skipping non-object JSON {json_path}",
                file=sys.stderr,
            )
            continue

        run_id = str(record.get("run_id") or json_path.stem)
        status = str(record.get("status", "unknown"))
        is_ok = status == "ok"

        row: dict[str, Any] = {
            "schema_version": record.get("schema_version"),
            "run_id": run_id,
            "status": status,
        }

        # 1. Config (always ingested)
        cfg = record.get("config")
        config_order_val: Any = None
        if isinstance(cfg, dict):
            config_order_val = cfg.get("order")
            for key, val in _flatten_dict(cfg):
                if _is_scalar(val):
                    row[key] = val
                    _register_unique(config_cols, config_seen, key)
                elif isinstance(val, list):
                    row[key] = ",".join(str(x) for x in val)
                    _register_unique(config_cols, config_seen, key)

        # 2. Measured sections (only for status == "ok")
        if is_ok:
            # problem_stats
            pstats = record.get("problem_stats")
            if isinstance(pstats, dict):
                ps_order = pstats.get("order")
                if ps_order is not None:
                    row["order"] = ps_order
                    _register_unique(config_cols, config_seen, "order")
                elif config_order_val is not None:
                    row["order"] = config_order_val
                    _register_unique(config_cols, config_seen, "order")

                for key, val in _flatten_dict(pstats):
                    if key == "order":
                        continue
                    if key in ("full_sizes_1d", "free_sizes_1d") and isinstance(val, list):
                        row[key] = "x".join(str(x) for x in val) if val else None
                    elif _is_scalar(val):
                        row[key] = val
                    elif isinstance(val, list):
                        row[key] = ",".join(str(x) for x in val) if val else None

            # timings_sec
            timings = record.get("timings_sec")
            if isinstance(timings, dict):
                for key, val in _flatten_dict(timings):
                    if key == "solve_repeats":
                        if isinstance(val, list):
                            reps = [float(x) for x in val if is_valid_number(x)]
                            row["solve_repeats_count"] = len(reps)
                            row["solve_std_sec"] = (
                                float(np.std(reps, ddof=1)) if len(reps) > 1 else 0.0
                            )
                        continue
                    if _is_scalar(val):
                        row[key] = val
                        _register_unique(timing_cols, timing_seen, key)

            # convergence
            conv = record.get("convergence")
            if isinstance(conv, dict):
                for key, val in _flatten_dict(conv):
                    if key == "residual_history":
                        if isinstance(val, list):
                            hist = [float(x) for x in val if is_valid_number(x)]
                            if hist:
                                residual_histories[run_id] = hist
                            row["residual_history_len"] = len(hist)
                        continue
                    if _is_scalar(val):
                        row[key] = val

            # verification
            verif = record.get("verification")
            if isinstance(verif, dict):
                for key, val in _flatten_dict(verif):
                    target_key = "verification_enabled" if key == "enabled" else key
                    if _is_scalar(val):
                        row[target_key] = val

            # system
            sys_info = record.get("system")
            if isinstance(sys_info, dict):
                for key, val in _flatten_dict(sys_info):
                    if _is_scalar(val):
                        row[key] = val
        else:
            # Ensure order falls back cleanly to config["order"] for non-ok runs
            if config_order_val is not None:
                row["order"] = config_order_val
                _register_unique(config_cols, config_seen, "order")

            # error
            err = record.get("error")
            if isinstance(err, dict):
                for key, val in _flatten_dict(err):
                    prefixed = key if key.startswith("error_") else f"error_{key}"
                    if _is_scalar(val):
                        row[prefixed] = val

        rows.append(row)

    if rows:
        df = pd.DataFrame(rows)
    else:
        df = pd.DataFrame(columns=["run_id", "status"])

    varying_axes, constant_params = compute_axis_partition(df, config_cols)
    env_metadata = _load_env_metadata(result_dir)
    custom_plot_specs = _load_custom_plot_specs(result_dir)

    return BenchmarkDataset(
        df=df,
        result_dir=result_dir,
        config_cols=config_cols,
        varying_axes=varying_axes,
        constant_params=constant_params,
        timing_cols=timing_cols,
        residual_histories=residual_histories,
        env_metadata=env_metadata,
        custom_plot_specs=custom_plot_specs,
    )


def apply_filters(
    dataset: BenchmarkDataset,
    filter_specs: list[str],
) -> BenchmarkDataset:
    """Filter *dataset* using 'key=v1,v2' specs and recompute axis partitions."""
    df = dataset.df
    if not filter_specs:
        filtered_df = df.copy().reset_index(drop=True)
        varying_axes, constant_params = compute_axis_partition(
            filtered_df, dataset.config_cols
        )
        return BenchmarkDataset(
            df=filtered_df,
            result_dir=dataset.result_dir,
            config_cols=list(dataset.config_cols),
            varying_axes=varying_axes,
            constant_params=constant_params,
            timing_cols=list(dataset.timing_cols),
            residual_histories=dict(dataset.residual_histories),
            env_metadata=dict(dataset.env_metadata),
            custom_plot_specs=list(dataset.custom_plot_specs),
            speedup_cols=list(dataset.speedup_cols),
        )

    mask = pd.Series(True, index=df.index, dtype=bool)

    for raw_spec in filter_specs:
        if "=" not in raw_spec:
            raise ValueError(
                f"Invalid filter syntax {raw_spec!r}: expected 'key=v1,v2'."
            )
        key_raw, vals_raw = raw_spec.split("=", 1)
        key = key_raw.strip()
        if not key:
            raise ValueError(
                f"Invalid filter syntax {raw_spec!r}: empty column key."
            )
        if key not in df.columns:
            avail = ", ".join(sorted(str(c) for c in df.columns))
            raise ValueError(
                f"Filter column {key!r} not found in dataset. Available: [{avail}]"
            )

        tokens = [tok.strip() for tok in vals_raw.split(",") if tok.strip() != ""]
        if not tokens:
            raise ValueError(
                f"Invalid filter syntax {raw_spec!r}: no filter values provided."
            )

        series = df[key]
        col_mask = pd.Series(False, index=df.index, dtype=bool)
        for tok in tokens:
            coerced = coerce_match_value(tok, series)
            col_mask |= series_matches_value(series, coerced)
        mask &= col_mask

    filtered_df = df.loc[mask].copy().reset_index(drop=True)
    kept_ids = (
        set(filtered_df["run_id"].astype(str).tolist())
        if "run_id" in filtered_df.columns
        else set()
    )
    filtered_histories = {
        rid: hist
        for rid, hist in dataset.residual_histories.items()
        if rid in kept_ids
    }
    varying_axes, constant_params = compute_axis_partition(
        filtered_df, dataset.config_cols
    )

    return BenchmarkDataset(
        df=filtered_df,
        result_dir=dataset.result_dir,
        config_cols=list(dataset.config_cols),
        varying_axes=varying_axes,
        constant_params=constant_params,
        timing_cols=list(dataset.timing_cols),
        residual_histories=filtered_histories,
        env_metadata=dict(dataset.env_metadata),
        custom_plot_specs=list(dataset.custom_plot_specs),
        speedup_cols=list(dataset.speedup_cols),
    )
