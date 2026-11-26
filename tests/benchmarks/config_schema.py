from __future__ import annotations

from dataclasses import dataclass
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11 fallback
    import tomli as tomllib  # type: ignore[no-redef]


ScalarValue = str | int | float | bool

ALLOWED_TOP_LEVEL_KEYS: frozenset[str] = frozenset(
    {"experiment", "meshes", "defaults", "sweeps", "cases", "plots"}
)

REQUIRED_EXPERIMENT_KEYS: frozenset[str] = frozenset(
    {
        "name",
        "executable",
        "working_dir",
        "output_root",
        "timeout_sec",
        "openblas_threads",
        "omp_threads",
    }
)

RESERVED_PARAM_KEYS: frozenset[str] = frozenset(
    {"out_file", "run_id", "quiet", "sweep_name"}
)

COMMON_KEY_ORDER: tuple[str, ...] = (
    "mesh_name",
    "problem",
    "ref_strategy",
    "order",
    "ref_levels",
    "assembly",
    "solver",
    "prec",
)

SWEEP_AXIS_PRIORITY: tuple[str, ...] = (
    "mesh",
    "mesh_name",
    "problem",
    "ref_strategy",
    "order",
    "ref_levels",
    "assembly",
    "solver",
    "prec",
)

RUN_ID_PREFIX_KEYS: tuple[str, ...] = (
    "mesh_name",
    "ref_strategy",
    "ref_levels",
    "order",
    "assembly",
    "solver",
    "prec",
)

_SAFE_TOKEN_RE = re.compile(r"[^A-Za-z0-9._-]+")
_VALID_PARAM_KEY_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
_MISSING = object()


@dataclass(frozen=True)
class ExperimentConfig:
    name: str
    executable: str
    working_dir: str
    output_root: str
    timeout_sec: float
    openblas_threads: int
    omp_threads: int


@dataclass(frozen=True)
class RunSpec:
    sweep_name: str
    params: dict[str, Any]
    run_id: str


@dataclass(frozen=True)
class BenchmarkPlan:
    experiment: ExperimentConfig
    meshes: dict[str, str]
    runs: list[RunSpec]
    excluded_count: int
    source_path: Path
    varying_keys: list[str]


def _validate_non_empty_str(value: Any, ctx: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{ctx} must be a non-empty string, got {value!r}")
    return value


def _validate_positive_int(value: Any, ctx: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{ctx} must be an integer >= 1, got {value!r}")
    return value


def _validate_positive_float(value: Any, ctx: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{ctx} must be a finite positive number, got {value!r}")
    fval = float(value)
    if not math.isfinite(fval) or fval <= 0.0:
        raise ValueError(f"{ctx} must be a finite positive number, got {value!r}")
    return fval


def _validate_param_key(key: Any, ctx: str) -> str:
    if not isinstance(key, str) or not _VALID_PARAM_KEY_RE.match(key):
        raise ValueError(
            f"{ctx}: invalid parameter key {key!r} "
            "(expected snake_case identifier starting with a letter)"
        )
    if key in RESERVED_PARAM_KEYS:
        raise ValueError(
            f"{ctx}: parameter key {key!r} is managed automatically by the runner "
            "and must not be set in the TOML config"
        )
    return key


def _validate_scalar(value: Any, ctx: str) -> ScalarValue:
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"{ctx} must be a finite float, got {value!r}")
        return value
    if isinstance(value, str):
        if not value.strip():
            raise ValueError(f"{ctx} must be a non-empty string, got {value!r}")
        return value
    raise ValueError(
        f"{ctx} must be a scalar TOML value (str, int, float, or bool), "
        f"got {type(value).__name__}: {value!r}"
    )


def _scalar_equal(lhs: Any, rhs: Any) -> bool:
    if isinstance(lhs, bool) or isinstance(rhs, bool):
        return type(lhs) is type(rhs) and lhs == rhs
    return lhs == rhs


def _parse_experiment(raw: Any) -> ExperimentConfig:
    if not isinstance(raw, dict):
        raise ValueError("[experiment] section is required and must be a table")

    keys = set(raw.keys())
    missing = sorted(REQUIRED_EXPERIMENT_KEYS - keys)
    if missing:
        raise ValueError(f"[experiment] missing required keys: {', '.join(missing)}")
    extra = sorted(keys - REQUIRED_EXPERIMENT_KEYS)
    if extra:
        raise ValueError(f"[experiment] contains unknown keys: {', '.join(extra)}")

    return ExperimentConfig(
        name=_validate_non_empty_str(raw["name"], "experiment.name"),
        executable=_validate_non_empty_str(raw["executable"], "experiment.executable"),
        working_dir=_validate_non_empty_str(raw["working_dir"], "experiment.working_dir"),
        output_root=_validate_non_empty_str(raw["output_root"], "experiment.output_root"),
        timeout_sec=_validate_positive_float(raw["timeout_sec"], "experiment.timeout_sec"),
        openblas_threads=_validate_positive_int(
            raw["openblas_threads"], "experiment.openblas_threads"
        ),
        omp_threads=_validate_positive_int(raw["omp_threads"], "experiment.omp_threads"),
    )


def _parse_meshes(raw: Any) -> dict[str, str]:
    if raw is None:
        return {}
    if not isinstance(raw, dict):
        raise ValueError("[meshes] must be a table mapping alias -> mesh path")
    meshes: dict[str, str] = {}
    for alias, path in raw.items():
        alias_str = _validate_non_empty_str(alias, "meshes key")
        path_str = _validate_non_empty_str(path, f"meshes.{alias_str}")
        meshes[alias_str] = path_str
    return meshes


def _parse_defaults(raw: Any) -> dict[str, ScalarValue]:
    if raw is None:
        return {}
    if not isinstance(raw, dict):
        raise ValueError("[defaults] must be a table of scalar parameter values")
    defaults: dict[str, ScalarValue] = {}
    for key, val in raw.items():
        valid_key = _validate_param_key(key, "defaults")
        defaults[valid_key] = _validate_scalar(val, f"defaults.{valid_key}")
    return defaults


def _validate_plots(raw: Any) -> None:
    if raw is None:
        return
    if not isinstance(raw, list) or not all(isinstance(item, dict) for item in raw):
        raise ValueError("[[plots]] must be a list of tables")


def _parse_exclude_rules(raw: Any, ctx: str) -> list[dict[str, ScalarValue]]:
    if raw is None:
        return []
    if not isinstance(raw, list):
        raise ValueError(f"{ctx}.exclude must be a list of tables")
    rules: list[dict[str, ScalarValue]] = []
    for idx, entry in enumerate(raw):
        entry_ctx = f"{ctx}.exclude[{idx}]"
        if not isinstance(entry, dict) or not entry:
            raise ValueError(f"{entry_ctx} must be a non-empty table")
        parsed_rule: dict[str, ScalarValue] = {}
        for key, val in entry.items():
            valid_key = _validate_param_key(key, entry_ctx)
            parsed_rule[valid_key] = _validate_scalar(val, f"{entry_ctx}.{valid_key}")
        rules.append(parsed_rule)
    return rules


def _resolve_mesh_alias(
    params: dict[str, Any], meshes: dict[str, str], ctx: str
) -> dict[str, Any]:
    resolved = dict(params)
    if "mesh" not in resolved:
        return resolved

    mesh_val = resolved["mesh"]
    if not isinstance(mesh_val, str) or not mesh_val.strip():
        raise ValueError(f"{ctx}.mesh must be a non-empty string, got {mesh_val!r}")

    if meshes:
        if mesh_val not in meshes:
            available = ", ".join(sorted(meshes.keys()))
            raise ValueError(
                f"{ctx}.mesh alias {mesh_val!r} is not defined in [meshes] "
                f"(available aliases: {available})"
            )
        resolved["mesh"] = meshes[mesh_val]
        if "mesh_name" not in resolved:
            resolved["mesh_name"] = mesh_val
    return resolved


def _matches_exclude_rule(
    rule: dict[str, ScalarValue],
    raw_params: dict[str, Any],
    resolved_params: dict[str, Any],
) -> bool:
    for key, target in rule.items():
        if key == "mesh":
            matched = (
                ("mesh" in raw_params and _scalar_equal(raw_params["mesh"], target))
                or (
                    "mesh" in resolved_params
                    and _scalar_equal(resolved_params["mesh"], target)
                )
                or (
                    "mesh_name" in resolved_params
                    and _scalar_equal(resolved_params["mesh_name"], target)
                )
            )
            if not matched:
                return False
        elif key == "mesh_name":
            matched = (
                "mesh_name" in resolved_params
                and _scalar_equal(resolved_params["mesh_name"], target)
            ) or (
                "mesh" in raw_params and _scalar_equal(raw_params["mesh"], target)
            )
            if not matched:
                return False
        else:
            if key not in resolved_params or not _scalar_equal(
                resolved_params[key], target
            ):
                return False
    return True


def _ordered_sweep_keys(keys: list[str]) -> list[str]:
    priority = {k: i for i, k in enumerate(SWEEP_AXIS_PRIORITY)}
    return sorted(keys, key=lambda k: (priority.get(k, len(priority)), k))


def _sanitize_token(value: str) -> str:
    cleaned = _SAFE_TOKEN_RE.sub("_", value.strip())
    return cleaned.strip("_")


def compute_run_id(params: dict[str, Any]) -> str:
    """Return a deterministic, filesystem-safe run identifier for a parameter dict."""
    payload = {k: v for k, v in params.items() if k != "sweep_name"}
    canonical = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    )
    hash8 = hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:8]

    parts: list[str] = []
    for key in RUN_ID_PREFIX_KEYS:
        if key not in payload:
            continue
        val = payload[key]
        if key == "ref_levels":
            raw_tok = f"r{val}"
        elif key == "order":
            raw_tok = f"p{val}"
        elif isinstance(val, bool):
            raw_tok = f"{key}_{'true' if val else 'false'}"
        else:
            raw_tok = str(val)
        token = _sanitize_token(raw_tok)
        if token:
            parts.append(token)

    if not parts:
        return f"run__{hash8}"

    prefix = "_".join(parts)
    if len(prefix) > 120:
        prefix = prefix[:120].rstrip("_")
    return f"{prefix}__{hash8}"


def _distinct_key_value(value: Any) -> Any:
    if value is _MISSING:
        return ("missing", None)
    if isinstance(value, bool):
        return ("bool", value)
    return ("val", value)


def _compute_varying_keys(runs: list[RunSpec]) -> list[str]:
    if len(runs) <= 1:
        return []

    all_keys: set[str] = set()
    for run in runs:
        all_keys.update(run.params.keys())

    if "mesh_name" in all_keys:
        all_keys.discard("mesh")

    varying: list[str] = []
    for key in all_keys:
        distinct = {
            _distinct_key_value(run.params.get(key, _MISSING)) for run in runs
        }
        if len(distinct) > 1:
            varying.append(key)

    priority = {k: i for i, k in enumerate(COMMON_KEY_ORDER)}
    return sorted(varying, key=lambda k: (priority.get(k, len(priority)), k))


def load_and_expand_config(config_path: str | Path) -> BenchmarkPlan:
    """Load, structurally validate, and expand a TOML benchmark config."""
    source_path = Path(config_path).resolve()
    if not source_path.is_file():
        raise FileNotFoundError(f"Benchmark config file not found: {source_path}")

    with source_path.open("rb") as f:
        raw = tomllib.load(f)

    if not isinstance(raw, dict):
        raise ValueError("TOML root document must be a table")

    unknown_top = sorted(set(raw.keys()) - ALLOWED_TOP_LEVEL_KEYS)
    if unknown_top:
        raise ValueError(
            f"Unknown top-level TOML section(s): {', '.join(unknown_top)} "
            f"(allowed: {', '.join(sorted(ALLOWED_TOP_LEVEL_KEYS))})"
        )

    experiment = _parse_experiment(raw.get("experiment"))
    meshes = _parse_meshes(raw.get("meshes"))
    defaults = _parse_defaults(raw.get("defaults"))
    _validate_plots(raw.get("plots"))

    raw_sweeps = raw.get("sweeps", [])
    raw_cases = raw.get("cases", [])

    if not isinstance(raw_sweeps, list) or not all(
        isinstance(s, dict) for s in raw_sweeps
    ):
        raise ValueError("[[sweeps]] must be a list of tables")
    if not isinstance(raw_cases, list) or not all(
        isinstance(c, dict) for c in raw_cases
    ):
        raise ValueError("[[cases]] must be a list of tables")
    if not raw_sweeps and not raw_cases:
        raise ValueError("Config must define at least one [[sweeps]] or [[cases]] entry")

    runs: list[RunSpec] = []
    seen_run_ids: set[str] = set()
    excluded_count = 0

    # Expand [[sweeps]]
    for sweep_idx, sweep in enumerate(raw_sweeps):
        ctx = f"sweeps[{sweep_idx}]"
        raw_name = sweep.get("name")
        if raw_name is not None:
            sweep_name = _validate_non_empty_str(raw_name, f"{ctx}.name")
        else:
            sweep_name = f"sweep_{sweep_idx}"

        exclude_rules = _parse_exclude_rules(sweep.get("exclude"), ctx)

        axes: dict[str, list[ScalarValue]] = {}
        for key, val in sweep.items():
            if key in ("name", "exclude"):
                continue
            valid_key = _validate_param_key(key, ctx)
            if isinstance(val, list):
                if not val:
                    raise ValueError(f"{ctx}.{valid_key} list must not be empty")
                axes[valid_key] = [
                    _validate_scalar(item, f"{ctx}.{valid_key}[{i}]")
                    for i, item in enumerate(val)
                ]
            else:
                axes[valid_key] = [_validate_scalar(val, f"{ctx}.{valid_key}")]

        ordered_keys = _ordered_sweep_keys(list(axes.keys()))
        value_lists = [axes[k] for k in ordered_keys]

        for combo in itertools.product(*value_lists):
            raw_params: dict[str, Any] = dict(defaults)
            for k, v in zip(ordered_keys, combo):
                raw_params[k] = v

            resolved_params = _resolve_mesh_alias(raw_params, meshes, ctx)

            if any(
                _matches_exclude_rule(rule, raw_params, resolved_params)
                for rule in exclude_rules
            ):
                excluded_count += 1
                continue

            run_id = compute_run_id(resolved_params)
            if run_id in seen_run_ids:
                continue
            seen_run_ids.add(run_id)
            runs.append(
                RunSpec(
                    sweep_name=sweep_name,
                    params=resolved_params,
                    run_id=run_id,
                )
            )

    # Expand [[cases]]
    for case_idx, case in enumerate(raw_cases):
        ctx = f"cases[{case_idx}]"
        raw_name = case.get("name")
        if raw_name is not None:
            case_name = _validate_non_empty_str(raw_name, f"{ctx}.name")
        else:
            case_name = f"case_{case_idx}"

        raw_params = dict(defaults)
        for key, val in case.items():
            if key == "name":
                continue
            valid_key = _validate_param_key(key, ctx)
            raw_params[valid_key] = _validate_scalar(val, f"{ctx}.{valid_key}")

        resolved_params = _resolve_mesh_alias(raw_params, meshes, ctx)
        run_id = compute_run_id(resolved_params)
        if run_id in seen_run_ids:
            continue
        seen_run_ids.add(run_id)
        runs.append(
            RunSpec(
                sweep_name=case_name,
                params=resolved_params,
                run_id=run_id,
            )
        )

    varying_keys = _compute_varying_keys(runs)

    return BenchmarkPlan(
        experiment=experiment,
        meshes=meshes,
        runs=runs,
        excluded_count=excluded_count,
        source_path=source_path,
        varying_keys=varying_keys,
    )


def build_cli_args(
    executable: Path,
    spec: RunSpec,
    out_file: Path,
    quiet: bool = True,
) -> list[str]:
    """Build C++ CLI arguments from an arbitrary RunSpec parameter dictionary."""
    args: list[str] = [str(executable)]

    for key, val in sorted(spec.params.items()):
        if not isinstance(key, str) or not key.strip():
            raise ValueError(f"Invalid CLI parameter key: {key!r}")
        flag_stem = key.replace("_", "-")
        flag = f"--{flag_stem}"

        if isinstance(val, bool):
            args.append(flag if val else f"--no-{flag_stem}")
        elif isinstance(val, float):
            if not math.isfinite(val):
                raise ValueError(f"Non-finite float for parameter {key!r}: {val!r}")
            args.extend([flag, repr(val)])
        elif isinstance(val, (int, str)):
            args.extend([flag, str(val)])
        else:
            raise ValueError(
                f"Unsupported non-scalar value for parameter {key!r}: "
                f"{type(val).__name__} ({val!r})"
            )

    args.extend(
        [
            "--out-file",
            str(out_file),
            "--run-id",
            spec.run_id,
            "--quiet" if quiet else "--no-quiet",
        ]
    )
    return args
