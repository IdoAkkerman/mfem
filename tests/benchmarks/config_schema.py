"""Configuration schema, validator, Cartesian sweep expander, deterministic run-ID
generator, and CLI command builder for MFEM NURBS Kronecker benchmarks.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11 fallback
    import tomli as tomllib  # type: ignore[no-redef]


ALLOWED_TOP_LEVEL_KEYS: frozenset[str] = frozenset(
    {"experiment", "meshes", "defaults", "sweeps", "cases"}
)

EXPERIMENT_REQUIRED_KEYS: frozenset[str] = frozenset(
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

ALLOWED_ENUMS: dict[str, frozenset[str]] = {
    "problem": frozenset({"const-1"}),
    "ref_strategy": frozenset({"k", "h", "p"}),
    "assembly": frozenset({"legacy"}),
    "solver": frozenset({"cg", "gmres", "minres"}),
    "prec": frozenset({"none", "jacobi", "gs", "kronecker-fd"}),
}

RUN_PARAM_ORDER: tuple[str, ...] = (
    "mesh_name",
    "problem",
    "ref_strategy",
    "order",
    "ref_levels",
    "assembly",
    "solver",
    "prec",
    "rtol",
    "atol",
    "max_iter",
    "gmres_restart",
    "repeats",
    "verify",
    "record_residuals",
)

REQUIRED_RUN_PARAM_KEYS: frozenset[str] = frozenset(RUN_PARAM_ORDER)
ALLOWED_TOML_PARAM_KEYS: frozenset[str] = REQUIRED_RUN_PARAM_KEYS | frozenset({"mesh"})

RUN_ID_KEYS: tuple[str, ...] = (
    "mesh_name",
    "mesh",
    "problem",
    "ref_strategy",
    "order",
    "ref_levels",
    "assembly",
    "solver",
    "prec",
    "rtol",
    "atol",
    "max_iter",
    "gmres_restart",
)


@dataclass(frozen=True)
class ExperimentConfig:
    name: str
    executable: str
    working_dir: str
    output_root: str
    timeout_sec: int
    openblas_threads: int
    omp_threads: int


@dataclass(frozen=True)
class RunSpec:
    sweep_name: str
    mesh_name: str
    mesh: str
    problem: str
    ref_strategy: str
    order: int
    ref_levels: int
    assembly: str
    solver: str
    prec: str
    rtol: float
    atol: float
    max_iter: int
    gmres_restart: int
    repeats: int
    verify: bool
    record_residuals: bool
    run_id: str


@dataclass(frozen=True)
class BenchmarkPlan:
    experiment: ExperimentConfig
    meshes: dict[str, str]
    runs: list[RunSpec]
    excluded_count: int
    source_path: Path


def _is_strict_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _is_strict_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _is_finite(value: float) -> bool:
    return value == value and value != float("inf") and value != float("-inf")


def _require_non_empty_str(value: Any, context: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{context} must be a non-empty string, got {value!r}")
    return value


def _require_int(value: Any, context: str, *, min_value: int | None = None) -> int:
    if not _is_strict_int(value):
        raise ValueError(
            f"{context} must be an integer, got {type(value).__name__} ({value!r})"
        )
    if min_value is not None and value < min_value:
        raise ValueError(f"{context} must be >= {min_value}, got {value}")
    return value


def _require_float(
    value: Any,
    context: str,
    *,
    min_value: float | None = None,
    strictly_greater: bool = False,
) -> float:
    if not _is_strict_number(value):
        raise ValueError(
            f"{context} must be a number, got {type(value).__name__} ({value!r})"
        )
    fval = float(value)
    if not _is_finite(fval):
        raise ValueError(f"{context} must be finite, got {value!r}")
    if min_value is not None:
        if strictly_greater and not (fval > min_value):
            raise ValueError(f"{context} must be > {min_value}, got {fval}")
        if not strictly_greater and not (fval >= min_value):
            raise ValueError(f"{context} must be >= {min_value}, got {fval}")
    return fval


def _require_bool(value: Any, context: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(
            f"{context} must be a boolean, got {type(value).__name__} ({value!r})"
        )
    return value


def _sanitize_mesh_name(mesh_name: str) -> str:
    sanitized = "".join(
        ch if (ch.isascii() and (ch.isalnum() or ch in "_-")) else "_"
        for ch in mesh_name
    )
    if not sanitized:
        raise ValueError("mesh_name cannot be empty")
    return sanitized


def _validate_run_param_value(
    key: str,
    value: Any,
    meshes: dict[str, str],
    context: str,
) -> Any:
    if key == "mesh_name":
        alias = _require_non_empty_str(value, context)
        if alias not in meshes:
            available = ", ".join(sorted(meshes)) or "<none>"
            raise ValueError(
                f"{context} references unknown mesh alias {alias!r} "
                f"(available: {available})"
            )
        return alias

    if key in ALLOWED_ENUMS:
        sval = _require_non_empty_str(value, context)
        allowed = ALLOWED_ENUMS[key]
        if sval not in allowed:
            raise ValueError(
                f"{context} has invalid value {sval!r}; allowed: {sorted(allowed)}"
            )
        return sval

    if key == "order":
        ival = _require_int(value, context)
        if ival != -1 and ival < 1:
            raise ValueError(f"{context} must be -1 or >= 1, got {ival}")
        return ival

    if key == "ref_levels":
        return _require_int(value, context, min_value=0)

    if key == "rtol":
        return _require_float(value, context, min_value=0.0, strictly_greater=True)

    if key == "atol":
        return _require_float(value, context, min_value=0.0, strictly_greater=False)

    if key in ("max_iter", "gmres_restart", "repeats"):
        return _require_int(value, context, min_value=1)

    if key in ("verify", "record_residuals"):
        return _require_bool(value, context)

    raise ValueError(f"Unknown parameter {key!r} in {context}")


def _normalize_param_keys(
    raw: dict[str, Any],
    allowed_extra_keys: frozenset[str],
    context: str,
) -> dict[str, Any]:
    allowed = ALLOWED_TOML_PARAM_KEYS | allowed_extra_keys
    unknown = sorted(set(raw.keys()) - allowed)
    if unknown:
        raise ValueError(f"{context} contains unknown key(s): {unknown}")

    if "mesh" in raw and "mesh_name" in raw:
        raise ValueError(
            f"{context} cannot specify both 'mesh' and 'mesh_name'; use one"
        )

    normalized: dict[str, Any] = {}
    for key, value in raw.items():
        target_key = "mesh_name" if key == "mesh" else key
        normalized[target_key] = value
    return normalized


def _parse_experiment(raw: Any) -> ExperimentConfig:
    if not isinstance(raw, dict):
        raise ValueError("[experiment] must be a TOML table")

    unknown = sorted(set(raw.keys()) - EXPERIMENT_REQUIRED_KEYS)
    if unknown:
        raise ValueError(f"[experiment] contains unknown key(s): {unknown}")

    missing = sorted(EXPERIMENT_REQUIRED_KEYS - set(raw.keys()))
    if missing:
        raise ValueError(f"[experiment] is missing required key(s): {missing}")

    return ExperimentConfig(
        name=_require_non_empty_str(raw["name"], "[experiment].name"),
        executable=_require_non_empty_str(raw["executable"], "[experiment].executable"),
        working_dir=_require_non_empty_str(
            raw["working_dir"], "[experiment].working_dir"
        ),
        output_root=_require_non_empty_str(
            raw["output_root"], "[experiment].output_root"
        ),
        timeout_sec=_require_int(
            raw["timeout_sec"], "[experiment].timeout_sec", min_value=1
        ),
        openblas_threads=_require_int(
            raw["openblas_threads"], "[experiment].openblas_threads", min_value=1
        ),
        omp_threads=_require_int(
            raw["omp_threads"], "[experiment].omp_threads", min_value=1
        ),
    )


def _parse_meshes(raw: Any) -> dict[str, str]:
    if not isinstance(raw, dict):
        raise ValueError("[meshes] must be a TOML table")
    if not raw:
        raise ValueError("[meshes] must define at least one mesh alias")

    meshes: dict[str, str] = {}
    for alias, path_str in raw.items():
        clean_alias = _require_non_empty_str(alias, "[meshes] alias")
        _sanitize_mesh_name(clean_alias)
        clean_path = _require_non_empty_str(path_str, f"[meshes].{clean_alias}")
        meshes[clean_alias] = clean_path
    return meshes


def _parse_defaults(raw: Any, meshes: dict[str, str]) -> dict[str, Any]:
    if raw is None:
        return {}
    if not isinstance(raw, dict):
        raise ValueError("[defaults] must be a TOML table")

    normalized = _normalize_param_keys(raw, frozenset(), "[defaults]")
    defaults: dict[str, Any] = {}
    for key, value in normalized.items():
        if isinstance(value, list):
            raise ValueError(
                f"[defaults].{key} must be a scalar value, not a list"
            )
        defaults[key] = _validate_run_param_value(
            key, value, meshes, f"[defaults].{key}"
        )
    return defaults


def _parse_exclude_rules(
    raw: Any,
    meshes: dict[str, str],
    sweep_context: str,
) -> list[dict[str, Any]]:
    if raw is None:
        return []
    if not isinstance(raw, list):
        raise ValueError(f"{sweep_context}.exclude must be a list of tables")

    rules: list[dict[str, Any]] = []
    for idx, entry in enumerate(raw):
        rule_ctx = f"{sweep_context}.exclude[{idx}]"
        if not isinstance(entry, dict):
            raise ValueError(f"{rule_ctx} must be a table/dict")
        if not entry:
            raise ValueError(f"{rule_ctx} must not be empty")

        normalized = _normalize_param_keys(entry, frozenset(), rule_ctx)
        validated_rule: dict[str, Any] = {}
        for key, value in normalized.items():
            if isinstance(value, list):
                raise ValueError(
                    f"{rule_ctx}.{key} must be a scalar value, not a list"
                )
            validated_rule[key] = _validate_run_param_value(
                key, value, meshes, f"{rule_ctx}.{key}"
            )
        rules.append(validated_rule)
    return rules


def _matches_exclude_rule(perm: dict[str, Any], rule: dict[str, Any]) -> bool:
    return all(perm[key] == expected for key, expected in rule.items())


def _build_run_spec(
    sweep_name: str,
    params: dict[str, Any],
    meshes: dict[str, str],
) -> RunSpec:
    mesh_name = params["mesh_name"]
    mesh_path = meshes[mesh_name]
    full_params = {**params, "mesh": mesh_path}
    run_id = compute_run_id(full_params)
    return RunSpec(
        sweep_name=sweep_name,
        mesh_name=mesh_name,
        mesh=mesh_path,
        problem=params["problem"],
        ref_strategy=params["ref_strategy"],
        order=params["order"],
        ref_levels=params["ref_levels"],
        assembly=params["assembly"],
        solver=params["solver"],
        prec=params["prec"],
        rtol=params["rtol"],
        atol=params["atol"],
        max_iter=params["max_iter"],
        gmres_restart=params["gmres_restart"],
        repeats=params["repeats"],
        verify=params["verify"],
        record_residuals=params["record_residuals"],
        run_id=run_id,
    )


def compute_run_id(params: dict) -> str:
    """Compute a deterministic, filesystem-safe run identifier from canonical parameters."""
    missing = [key for key in RUN_ID_KEYS if key not in params]
    if missing:
        raise ValueError(
            f"compute_run_id missing required parameter(s): {missing}"
        )

    mesh_name = _require_non_empty_str(params["mesh_name"], "params['mesh_name']")
    mesh = _require_non_empty_str(params["mesh"], "params['mesh']")
    problem = _require_non_empty_str(params["problem"], "params['problem']")
    ref_strategy = _require_non_empty_str(
        params["ref_strategy"], "params['ref_strategy']"
    )
    order = _require_int(params["order"], "params['order']")
    ref_levels = _require_int(params["ref_levels"], "params['ref_levels']", min_value=0)
    assembly = _require_non_empty_str(params["assembly"], "params['assembly']")
    solver = _require_non_empty_str(params["solver"], "params['solver']")
    prec = _require_non_empty_str(params["prec"], "params['prec']")
    rtol = _require_float(params["rtol"], "params['rtol']", min_value=0.0, strictly_greater=True)
    atol = _require_float(params["atol"], "params['atol']", min_value=0.0, strictly_greater=False)
    max_iter = _require_int(params["max_iter"], "params['max_iter']", min_value=1)
    gmres_restart = _require_int(
        params["gmres_restart"], "params['gmres_restart']", min_value=1
    )

    canonical_dict = {
        "mesh_name": mesh_name,
        "mesh": mesh,
        "problem": problem,
        "ref_strategy": ref_strategy,
        "order": order,
        "ref_levels": ref_levels,
        "assembly": assembly,
        "solver": solver,
        "prec": prec,
        "rtol": rtol,
        "atol": atol,
        "max_iter": max_iter,
        "gmres_restart": gmres_restart,
    }

    payload = json.dumps(canonical_dict, sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )
    hash6 = hashlib.sha256(payload).hexdigest()[:6]
    safe_mesh_name = _sanitize_mesh_name(mesh_name)
    order_tag = "base" if order == -1 else str(order)

    return (
        f"{safe_mesh_name}__{ref_strategy}_r{ref_levels}_p{order_tag}"
        f"__{assembly}__{solver}_{prec}__{hash6}"
    )


def load_and_expand_config(toml_path: Path) -> BenchmarkPlan:
    """Load, strictly validate, expand, exclude, and deduplicate a benchmark TOML config."""
    path = Path(toml_path)
    with path.open("rb") as fh:
        raw_config = tomllib.load(fh)

    if not isinstance(raw_config, dict):
        raise ValueError("Top-level TOML document must be a table")

    unknown_top = sorted(set(raw_config.keys()) - ALLOWED_TOP_LEVEL_KEYS)
    if unknown_top:
        raise ValueError(f"Unknown top-level section(s): {unknown_top}")

    if "experiment" not in raw_config:
        raise ValueError("Missing required [experiment] section")
    if "meshes" not in raw_config:
        raise ValueError("Missing required [meshes] section")
    if "sweeps" not in raw_config and "cases" not in raw_config:
        raise ValueError("Config must define at least one [[sweeps]] or [[cases]] entry")

    experiment = _parse_experiment(raw_config["experiment"])
    meshes = _parse_meshes(raw_config["meshes"])
    defaults = _parse_defaults(raw_config.get("defaults"), meshes)

    runs: list[RunSpec] = []
    seen_run_ids: set[str] = set()
    excluded_count = 0

    raw_sweeps = raw_config.get("sweeps", [])
    if not isinstance(raw_sweeps, list):
        raise ValueError("[[sweeps]] must be an array of tables")

    for idx, raw_sweep in enumerate(raw_sweeps):
        sweep_ctx = f"sweeps[{idx}]"
        if not isinstance(raw_sweep, dict):
            raise ValueError(f"{sweep_ctx} must be a table")

        normalized = _normalize_param_keys(
            raw_sweep, frozenset({"name", "exclude"}), sweep_ctx
        )
        if "name" not in normalized:
            raise ValueError(f"{sweep_ctx} is missing required 'name' field")
        sweep_name = _require_non_empty_str(normalized["name"], f"{sweep_ctx}.name")
        named_ctx = f"sweeps[{idx}] ({sweep_name!r})"

        exclude_rules = _parse_exclude_rules(
            normalized.get("exclude"), meshes, named_ctx
        )

        axes: dict[str, list[Any]] = {k: [v] for k, v in defaults.items()}
        for key, raw_val in normalized.items():
            if key in ("name", "exclude"):
                continue
            if isinstance(raw_val, list):
                if not raw_val:
                    raise ValueError(f"{named_ctx}.{key} list must not be empty")
                axes[key] = [
                    _validate_run_param_value(
                        key, elem, meshes, f"{named_ctx}.{key}[{i}]"
                    )
                    for i, elem in enumerate(raw_val)
                ]
            else:
                axes[key] = [
                    _validate_run_param_value(
                        key, raw_val, meshes, f"{named_ctx}.{key}"
                    )
                ]

        missing = sorted(REQUIRED_RUN_PARAM_KEYS - set(axes.keys()))
        if missing:
            raise ValueError(
                f"{named_ctx} is missing required parameter(s): {missing}"
            )

        value_lists = [axes[k] for k in RUN_PARAM_ORDER]
        for combo in itertools.product(*value_lists):
            perm = dict(zip(RUN_PARAM_ORDER, combo))
            if any(_matches_exclude_rule(perm, rule) for rule in exclude_rules):
                excluded_count += 1
                continue

            spec = _build_run_spec(sweep_name, perm, meshes)
            if spec.run_id in seen_run_ids:
                continue
            seen_run_ids.add(spec.run_id)
            runs.append(spec)

    raw_cases = raw_config.get("cases", [])
    if not isinstance(raw_cases, list):
        raise ValueError("[[cases]] must be an array of tables")

    for idx, raw_case in enumerate(raw_cases):
        case_ctx = f"cases[{idx}]"
        if not isinstance(raw_case, dict):
            raise ValueError(f"{case_ctx} must be a table")

        normalized = _normalize_param_keys(
            raw_case, frozenset({"sweep_name"}), case_ctx
        )
        sweep_name = _require_non_empty_str(
            normalized.get("sweep_name", "explicit"), f"{case_ctx}.sweep_name"
        )

        merged = dict(defaults)
        for key, raw_val in normalized.items():
            if key == "sweep_name":
                continue
            if isinstance(raw_val, list):
                raise ValueError(
                    f"{case_ctx}.{key} must be a scalar value, not a list"
                )
            merged[key] = _validate_run_param_value(
                key, raw_val, meshes, f"{case_ctx}.{key}"
            )

        missing = sorted(REQUIRED_RUN_PARAM_KEYS - set(merged.keys()))
        if missing:
            raise ValueError(
                f"{case_ctx} is missing required parameter(s): {missing}"
            )

        spec = _build_run_spec(sweep_name, merged, meshes)
        if spec.run_id in seen_run_ids:
            continue
        seen_run_ids.add(spec.run_id)
        runs.append(spec)

    return BenchmarkPlan(
        experiment=experiment,
        meshes=meshes,
        runs=runs,
        excluded_count=excluded_count,
        source_path=path,
    )


def build_cli_args(
    executable: Path,
    spec: RunSpec,
    out_file: Path,
    quiet: bool = True,
) -> list[str]:
    """Build exact mfem::OptionsParser CLI arguments for a single benchmark run."""
    return [
        str(executable),
        "-m",
        spec.mesh,
        "-mn",
        spec.mesh_name,
        "-prob",
        spec.problem,
        "-rs",
        spec.ref_strategy,
        "-o",
        str(spec.order),
        "-r",
        str(spec.ref_levels),
        "-asm",
        spec.assembly,
        "-sol",
        spec.solver,
        "-prec",
        spec.prec,
        "-rtol",
        repr(spec.rtol),
        "-atol",
        repr(spec.atol),
        "-mi",
        str(spec.max_iter),
        "-kdim",
        str(spec.gmres_restart),
        "-rep",
        str(spec.repeats),
        "--verify" if spec.verify else "--no-verify",
        "--record-residuals" if spec.record_residuals else "--no-record-residuals",
        "-out",
        str(out_file),
        "-id",
        spec.run_id,
        "--quiet" if quiet else "--no-quiet",
    ]
