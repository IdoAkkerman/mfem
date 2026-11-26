#!/usr/bin/env python3
"""Sequential CLI benchmark runner for MFEM NURBS Kronecker experiments.

Parses a TOML sweep configuration via sibling ``config_schema.py``, validates
the execution environment, runs the compiled ``nurbs_kronecker`` binary
sequentially with strict threading/timing hygiene, captures failures as v1.0
JSON stubs, and writes structured per-run results.
"""

from __future__ import annotations

import argparse
import dataclasses
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import socket
import subprocess
import sys
import time
from typing import Any

# Ensure sibling config_schema.py resolves cleanly regardless of CWD.
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from config_schema import (  # noqa: E402
    BenchmarkPlan,
    ExperimentConfig,
    RunSpec,
    build_cli_args,
    load_and_expand_config,
)

SCHEMA_VERSION = "1.0"
STDERR_TAIL_LIMIT = 2000
ACCEPTABLE_CMAKE_BUILD_TYPES = frozenset({"Release", "RelWithDebInfo"})
RUN_SPEC_CONFIG_FIELDS = (
    "sweep_name",
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
    "repeats",
    "verify",
    "record_residuals",
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run sequential MFEM NURBS Kronecker benchmarks from a TOML "
            "sweep configuration."
        )
    )
    parser.add_argument(
        "config",
        type=Path,
        help="Path to the TOML benchmark sweep configuration file.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Validate config and pre-flight checks, print run summary and "
            "first 5 CLI commands, then exit without creating directories or "
            "executing runs."
        ),
    )
    parser.add_argument(
        "--resume",
        type=Path,
        default=None,
        metavar="DIR",
        help=(
            "Reuse an existing result directory; skip runs whose "
            "runs/<run_id>.json already exists with status == 'ok'."
        ),
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=None,
        metavar="DIR",
        help=(
            "Base directory for resolving relative paths in [experiment] "
            "(defaults to current working directory)."
        ),
    )
    parser.add_argument(
        "--verbose-cpp",
        action="store_true",
        help="Pass quiet=False to build_cli_args and echo C++ stdout/stderr on failure.",
    )
    return parser.parse_args(argv)


def resolve_under(base: Path, raw_path: str | Path) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path.resolve()
    return (base / path).resolve()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    tmp_path = path.with_name(f"{path.name}.{os.getpid()}.tmp")
    try:
        with tmp_path.open("w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, sort_keys=False)
            fh.write("\n")
            fh.flush()
            os.fsync(fh.fileno())
        tmp_path.replace(path)
    except Exception:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass
        raise


def discover_cmake_cache(executable: Path, repo_root: Path) -> Path | None:
    candidates: list[Path] = []
    seen: set[Path] = set()

    for parent in executable.parents:
        candidate = parent / "CMakeCache.txt"
        if candidate not in seen:
            seen.add(candidate)
            candidates.append(candidate)

    default_build_cache = (repo_root / "build" / "CMakeCache.txt").resolve()
    if default_build_cache not in seen:
        candidates.append(default_build_cache)

    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def inspect_cmake_build_type(executable: Path, repo_root: Path) -> str | None:
    cache_path = discover_cmake_cache(executable, repo_root)
    if cache_path is None:
        return None

    pattern = re.compile(r"^CMAKE_BUILD_TYPE(?::[^=]+)?=(.*)$")
    try:
        with cache_path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                match = pattern.match(line.strip())
                if match:
                    build_type = match.group(1).strip()
                    if build_type not in ACCEPTABLE_CMAKE_BUILD_TYPES:
                        shown = build_type if build_type else "<unset>"
                        print(
                            "\n"
                            "========================================================================\n"
                            f"WARNING: Non-optimized CMake build detected: CMAKE_BUILD_TYPE={shown}\n"
                            f"         Cache file: {cache_path}\n"
                            "         Benchmark timings are unreliable unless built with\n"
                            "         Release or RelWithDebInfo.\n"
                            "========================================================================\n",
                            file=sys.stderr,
                            flush=True,
                        )
                    return build_type
    except OSError:
        return None
    return None


def run_preflight_checks(
    plan: BenchmarkPlan,
    repo_root: Path,
) -> tuple[Path, Path, Path]:
    exp: ExperimentConfig = plan.experiment
    executable = resolve_under(repo_root, exp.executable)
    working_dir = resolve_under(repo_root, exp.working_dir)
    output_root = resolve_under(repo_root, exp.output_root)

    errors: list[str] = []

    if not executable.exists():
        errors.append(f"Executable does not exist: {executable}")
    elif not executable.is_file():
        errors.append(f"Executable path is not a regular file: {executable}")
    elif not os.access(executable, os.X_OK):
        errors.append(f"Executable is not executable (missing +x): {executable}")

    if not working_dir.exists():
        errors.append(f"Working directory does not exist: {working_dir}")
    elif not working_dir.is_dir():
        errors.append(f"Working directory is not a directory: {working_dir}")

    if working_dir.is_dir():
        unique_meshes = sorted({spec.mesh for spec in plan.runs})
        missing_meshes: list[str] = []
        for mesh_rel in unique_meshes:
            mesh_path = working_dir / mesh_rel
            if not mesh_path.is_file():
                missing_meshes.append(f"  - {mesh_rel} (resolved: {mesh_path})")
        if missing_meshes:
            errors.append(
                "Missing mesh file(s) required by benchmark plan:\n"
                + "\n".join(missing_meshes)
            )

    if errors:
        print("PRE-FLIGHT CHECK FAILED:", file=sys.stderr)
        for err in errors:
            print(f"  * {err}", file=sys.stderr)
        sys.exit(1)

    inspect_cmake_build_type(executable, repo_root)
    return executable, working_dir, output_root


def detect_cpu_model() -> str:
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        try:
            with cpuinfo.open("r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if ":" not in line:
                        continue
                    key, val = line.split(":", 1)
                    key_norm = key.strip().lower()
                    if key_norm in {"model name", "hardware", "model"}:
                        model = val.strip()
                        if model:
                            return model
        except OSError:
            pass

    proc = platform.processor().strip()
    if proc:
        return proc
    return platform.machine() or "unknown"


def _run_git_cmd(repo_root: Path, git_args: list[str]) -> str | None:
    try:
        completed = subprocess.run(
            ["git", *git_args],
            cwd=repo_root,
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode != 0:
        return None
    return completed.stdout.strip()


def collect_env_metadata(
    exp: ExperimentConfig,
    executable: Path,
    working_dir: Path,
    repo_root: Path,
) -> dict[str, Any]:
    git_commit = _run_git_cmd(repo_root, ["rev-parse", "HEAD"])
    git_branch = _run_git_cmd(repo_root, ["rev-parse", "--abbrev-ref", "HEAD"])
    git_status = _run_git_cmd(repo_root, ["status", "--porcelain"])
    git_dirty: bool | None = None if git_status is None else bool(git_status)

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "cpu_model": detect_cpu_model(),
        "cpu_count": os.cpu_count(),
        "git_commit": git_commit,
        "git_dirty": git_dirty,
        "git_branch": git_branch,
        "openblas_threads": exp.openblas_threads,
        "omp_threads": exp.omp_threads,
        "executable": str(executable),
        "working_dir": str(working_dir),
    }


def make_result_dir_name(experiment_name: str) -> str:
    stamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", experiment_name).strip("_")
    if not safe_name:
        safe_name = "benchmark"
    return f"{stamp}_{safe_name}"


def spec_to_config_dict(spec: RunSpec) -> dict[str, Any]:
    if dataclasses.is_dataclass(spec):
        raw = dataclasses.asdict(spec)
        raw.pop("run_id", None)
        return raw
    return {field: getattr(spec, field) for field in RUN_SPEC_CONFIG_FIELDS}


def is_cached_ok(out_file: Path) -> bool:
    if not out_file.is_file():
        return False
    try:
        with out_file.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return False
    return isinstance(data, dict) and data.get("status") == "ok"


def decode_output_chunk(chunk: str | bytes | None) -> str:
    if chunk is None:
        return ""
    if isinstance(chunk, bytes):
        return chunk.decode("utf-8", errors="replace")
    return chunk


def build_failure_stub(
    spec: RunSpec,
    status: str,
    return_code: int | None,
    wall_time_sec: float,
    stderr_text: str,
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "run_id": spec.run_id,
        "status": status,
        "config": spec_to_config_dict(spec),
        "error": {
            "return_code": return_code,
            "wall_time_sec": round(wall_time_sec, 6),
            "stderr_tail": stderr_text[-STDERR_TAIL_LIMIT:],
        },
    }


def build_clean_env(exp: ExperimentConfig) -> dict[str, str]:
    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = str(exp.openblas_threads)
    env["OMP_NUM_THREADS"] = str(exp.omp_threads)
    env["MKL_NUM_THREADS"] = str(exp.openblas_threads)
    env["NUMEXPR_NUM_THREADS"] = str(exp.openblas_threads)
    return env


def format_run_prefix(index: int, total: int, spec: RunSpec) -> str:
    width = max(1, len(str(total)))
    return (
        f"[{index:>{width}}/{total}] {spec.mesh_name} | "
        f"{spec.ref_strategy} r{spec.ref_levels} p{spec.order} | "
        f"{spec.solver} + {spec.prec}"
    )


def format_ok_details(data: dict[str, Any], wall_time_sec: float) -> str:
    convergence = data.get("convergence")
    if not isinstance(convergence, dict):
        convergence = {}
    timings = data.get("timings_sec")
    if not isinstance(timings, dict):
        timings = {}

    iterations = data.get("iterations")
    if iterations is None:
        iterations = convergence.get("iterations")
    converged = convergence.get("converged")
    solve_sec = timings.get("solve")
    total_sec = timings.get("total")
    if not isinstance(total_sec, (int, float)):
        total_sec = wall_time_sec

    parts: list[str] = []
    if iterations is not None:
        parts.append(f"{iterations} iters")
    if converged is False:
        parts.append("NOT CONVERGED")
    if isinstance(solve_sec, (int, float)):
        parts.append(f"solve={solve_sec:.3f}s")
    parts.append(f"total={total_sec:.3f}s")
    return ", ".join(parts)


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(Path.cwd().resolve()))
    except ValueError:
        return str(path)


def echo_verbose_failure(stdout_text: str, stderr_text: str) -> None:
    if stdout_text.strip():
        print("  --- C++ stdout ---", file=sys.stderr)
        print(stdout_text.rstrip(), file=sys.stderr)
    if stderr_text.strip():
        print("  --- C++ stderr ---", file=sys.stderr)
        print(stderr_text.rstrip(), file=sys.stderr)
    sys.stderr.flush()


def execute_single_run(
    spec: RunSpec,
    executable: Path,
    working_dir: Path,
    out_file: Path,
    env: dict[str, str],
    timeout_sec: float | None,
    verbose_cpp: bool,
) -> tuple[str, str]:
    """Run one benchmark spec sequentially and persist its JSON result.

    Returns ``(outcome, log_message)`` where ``outcome`` is ``"ok"`` or
    ``"failed"``.
    """
    # Clear any stale output from an earlier failed attempt before launching.
    if out_file.exists():
        try:
            out_file.unlink()
        except OSError:
            pass

    cmd = build_cli_args(
        executable=executable,
        spec=spec,
        out_file=out_file,
        quiet=not verbose_cpp,
    )

    timeout_arg = timeout_sec if (timeout_sec is not None and timeout_sec > 0) else None
    t_start = time.perf_counter()

    try:
        completed = subprocess.run(
            cmd,
            cwd=working_dir,
            env=env,
            timeout=timeout_arg,
            capture_output=True,
            text=True,
            check=False,
        )
        wall_time = time.perf_counter() - t_start
    except subprocess.TimeoutExpired as exc:
        wall_time = time.perf_counter() - t_start
        stdout_text = decode_output_chunk(exc.stdout)
        stderr_text = decode_output_chunk(exc.stderr)
        stub = build_failure_stub(
            spec=spec,
            status="timeout",
            return_code=None,
            wall_time_sec=wall_time,
            stderr_text=stderr_text,
        )
        atomic_write_json(out_file, stub)
        if verbose_cpp:
            echo_verbose_failure(stdout_text, stderr_text)
        return "failed", f"FAIL (timeout, {wall_time:.1f}s)"
    except OSError as exc:
        wall_time = time.perf_counter() - t_start
        err_msg = f"Failed to spawn executable: {exc}"
        stub = build_failure_stub(
            spec=spec,
            status="crash",
            return_code=None,
            wall_time_sec=wall_time,
            stderr_text=err_msg,
        )
        atomic_write_json(out_file, stub)
        if verbose_cpp:
            echo_verbose_failure("", err_msg)
        return "failed", f"FAIL (crash, {wall_time:.1f}s)"

    stdout_text = completed.stdout or ""
    stderr_text = completed.stderr or ""
    rc = completed.returncode

    if rc in (-9, 137):
        stub = build_failure_stub(
            spec=spec,
            status="oom",
            return_code=rc,
            wall_time_sec=wall_time,
            stderr_text=stderr_text,
        )
        atomic_write_json(out_file, stub)
        if verbose_cpp:
            echo_verbose_failure(stdout_text, stderr_text)
        return "failed", f"FAIL (oom, {wall_time:.1f}s)"

    if rc != 0:
        stub = build_failure_stub(
            spec=spec,
            status="crash",
            return_code=rc,
            wall_time_sec=wall_time,
            stderr_text=stderr_text,
        )
        atomic_write_json(out_file, stub)
        if verbose_cpp:
            echo_verbose_failure(stdout_text, stderr_text)
        return "failed", f"FAIL (crash, rc={rc}, {wall_time:.1f}s)"

    # Return code == 0: verify C++ wrote valid JSON.
    if not out_file.is_file():
        diag = (
            f"{stderr_text}\n[runner] Expected JSON output missing: {out_file}"
        ).strip()
        stub = build_failure_stub(
            spec=spec,
            status="crash",
            return_code=rc,
            wall_time_sec=wall_time,
            stderr_text=diag,
        )
        atomic_write_json(out_file, stub)
        if verbose_cpp:
            echo_verbose_failure(stdout_text, diag)
        return "failed", f"FAIL (crash, missing json, {wall_time:.1f}s)"

    try:
        with out_file.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
        if not isinstance(data, dict):
            raise ValueError("Top-level JSON value is not an object.")
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        diag = (
            f"{stderr_text}\n[runner] Invalid JSON in {out_file}: {exc}"
        ).strip()
        stub = build_failure_stub(
            spec=spec,
            status="crash",
            return_code=rc,
            wall_time_sec=wall_time,
            stderr_text=diag,
        )
        atomic_write_json(out_file, stub)
        if verbose_cpp:
            echo_verbose_failure(stdout_text, diag)
        return "failed", f"FAIL (crash, invalid json, {wall_time:.1f}s)"

    # Ensure sweep_name is recorded under config and write back atomically.
    config_obj = data.get("config")
    if not isinstance(config_obj, dict):
        config_obj = spec_to_config_dict(spec)
        data["config"] = config_obj
        atomic_write_json(out_file, data)
    elif "sweep_name" not in config_obj:
        config_obj["sweep_name"] = spec.sweep_name
        atomic_write_json(out_file, data)

    if data.get("status") not in (None, "ok"):
        status_val = str(data.get("status"))
        if verbose_cpp:
            echo_verbose_failure(stdout_text, stderr_text)
        return "failed", f"FAIL ({status_val}, {wall_time:.1f}s)"

    details = format_ok_details(data, wall_time)
    return "ok", f"OK ({details})"


def handle_dry_run(
    plan: BenchmarkPlan,
    executable: Path,
    output_root: Path,
    resume_dir: Path | None,
    verbose_cpp: bool,
) -> int:
    preview_dir = (
        resume_dir.resolve()
        if resume_dir is not None
        else output_root / make_result_dir_name(plan.experiment.name)
    )
    preview_runs_dir = preview_dir / "runs"

    print(f"Dry run validated: {plan.source_path}")
    print(f"Experiment:        {plan.experiment.name}")
    print(f"Total runs:        {len(plan.runs)}")
    print(f"Excluded runs:     {plan.excluded_count}")
    print(f"Target dir:        {display_path(preview_dir)} (not created)")
    print()

    sample_count = min(5, len(plan.runs))
    print(f"First {sample_count} CLI command(s):")
    for idx, spec in enumerate(plan.runs[:sample_count], start=1):
        out_file = (preview_runs_dir / f"{spec.run_id}.json").resolve()
        cmd = build_cli_args(
            executable=executable,
            spec=spec,
            out_file=out_file,
            quiet=not verbose_cpp,
        )
        print(f"  {idx}. {shlex.join(cmd)}")

    return 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    config_path = args.config.resolve()
    repo_root = (
        args.repo_root.resolve()
        if args.repo_root is not None
        else Path.cwd().resolve()
    )

    if not config_path.is_file():
        print(f"ERROR: Config file not found: {config_path}", file=sys.stderr)
        return 1

    try:
        plan: BenchmarkPlan = load_and_expand_config(config_path)
    except Exception as exc:
        print(f"ERROR: Failed to load config {config_path}: {exc}", file=sys.stderr)
        return 1

    executable, working_dir, output_root = run_preflight_checks(plan, repo_root)

    if args.dry_run:
        return handle_dry_run(
            plan=plan,
            executable=executable,
            output_root=output_root,
            resume_dir=args.resume,
            verbose_cpp=args.verbose_cpp,
        )

    # Set up result directory layout.
    if args.resume is not None:
        result_dir = args.resume.resolve()
        if not result_dir.is_dir():
            print(
                f"ERROR: --resume directory does not exist: {result_dir}",
                file=sys.stderr,
            )
            return 1
        runs_dir = result_dir / "runs"
        runs_dir.mkdir(parents=True, exist_ok=True)
    else:
        result_dir = output_root / make_result_dir_name(plan.experiment.name)
        runs_dir = result_dir / "runs"
        runs_dir.mkdir(parents=True, exist_ok=True)

    snapshot_path = result_dir / "config_snapshot.toml"
    if not snapshot_path.exists():
        shutil.copy2(config_path, snapshot_path)

    metadata = collect_env_metadata(
        exp=plan.experiment,
        executable=executable,
        working_dir=working_dir,
        repo_root=repo_root,
    )
    atomic_write_json(result_dir / "env_metadata.json", metadata)

    env = build_clean_env(plan.experiment)
    total_runs = len(plan.runs)
    ok_count = 0
    skipped_count = 0
    failed_count = 0

    print(
        f"Starting sweep '{plan.experiment.name}': "
        f"{total_runs} runs ({plan.excluded_count} excluded)",
        flush=True,
    )
    print(f"Results directory: {display_path(result_dir)}", flush=True)

    sweep_start = time.perf_counter()
    interrupted = False

    try:
        for idx, spec in enumerate(plan.runs, start=1):
            prefix = format_run_prefix(idx, total_runs, spec)
            out_file = (runs_dir / f"{spec.run_id}.json").resolve()

            if args.resume is not None and is_cached_ok(out_file):
                skipped_count += 1
                print(f"{prefix} -> SKIP (cached)", flush=True)
                continue

            outcome, message = execute_single_run(
                spec=spec,
                executable=executable,
                working_dir=working_dir,
                out_file=out_file,
                env=env,
                timeout_sec=plan.experiment.timeout_sec,
                verbose_cpp=args.verbose_cpp,
            )
            if outcome == "ok":
                ok_count += 1
            else:
                failed_count += 1
            print(f"{prefix} -> {message}", flush=True)
    except KeyboardInterrupt:
        interrupted = True
        print("\nInterrupted by user; partial results preserved.", file=sys.stderr)

    total_wall = time.perf_counter() - sweep_start
    rel_result_dir = display_path(result_dir)

    print()
    print(
        f"Summary: ok={ok_count}, skipped={skipped_count}, failed={failed_count} "
        f"(total={total_runs}) in {total_wall:.2f}s"
    )
    print(f"Analyze results with:\n  python3 tests/benchmarks/analyze_results.py {rel_result_dir}")

    return 130 if interrupted else 0


if __name__ == "__main__":
    raise SystemExit(main())
