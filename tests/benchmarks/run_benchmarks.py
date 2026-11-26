#!/usr/bin/env python3
from __future__ import annotations

import argparse
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
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
import threading
import time
from typing import Any

# Ensure sibling config_schema.py is importable regardless of current working directory.
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from config_schema import (  # noqa: E402
    COMMON_KEY_ORDER,
    BenchmarkPlan,
    ExperimentConfig,
    RunSpec,
    build_cli_args,
    load_and_expand_config,
)


_SAFE_DIR_RE = re.compile(r"[^A-Za-z0-9._-]+")
_STDERR_TAIL_LIMIT = 2000


def _positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid int value: {value!r}") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError(f"must be >= 1 (got {parsed})")
    return parsed


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Execute MFEM NURBS Kronecker preconditioner benchmark sweeps "
            "and persist per-run JSON results."
        )
    )
    parser.add_argument(
        "config",
        type=Path,
        help="Path to the benchmark TOML configuration file.",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=_positive_int,
        default=1,
        metavar="N",
        help="Number of benchmark runs to execute concurrently (default: 1).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate config, run pre-flight checks, and print the plan without executing.",
    )
    parser.add_argument(
        "--resume",
        type=Path,
        default=None,
        metavar="DIR",
        help="Resume into an existing results directory, skipping runs with status == 'ok'.",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=None,
        metavar="DIR",
        help="Override repository root directory (default: auto-detected from script path).",
    )
    parser.add_argument(
        "--verbose-cpp",
        action="store_true",
        help="Pass --no-quiet to the C++ binary and echo its stdout/stderr.",
    )
    return parser.parse_args(argv)


def detect_repo_root(override: Path | None) -> Path:
    if override is not None:
        return override.resolve()
    # Expected layout: <repo_root>/tests/benchmarks/run_benchmarks.py
    candidate = _SCRIPT_DIR.parent.parent
    if (candidate / ".git").exists() or (candidate / "CMakeLists.txt").exists():
        return candidate
    return Path.cwd().resolve()


def resolve_under(base: Path, raw_path: str) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path
    return (base / path).resolve()


def _unlink_quiet(path: Path) -> None:
    try:
        path.unlink(missing_ok=True)
    except OSError:
        pass


def atomic_write_json(target: Path, payload: dict[str, Any]) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = target.with_name(
        f".{target.name}.{os.getpid()}.{threading.get_ident()}.tmp"
    )
    try:
        serialized = json.dumps(payload, indent=2, sort_keys=False) + "\n"
        with tmp_path.open("w", encoding="utf-8") as f:
            f.write(serialized)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp_path, target)
    except Exception:
        _unlink_quiet(tmp_path)
        raise


def detect_cmake_build_type(executable: Path, repo_root: Path) -> str | None:
    current = executable.resolve().parent
    stop = repo_root.resolve()
    for _ in range(8):
        cache_file = current / "CMakeCache.txt"
        if cache_file.is_file():
            try:
                for line in cache_file.read_text(
                    encoding="utf-8", errors="replace"
                ).splitlines():
                    if line.startswith("CMAKE_BUILD_TYPE:"):
                        _, _, value = line.partition("=")
                        return value.strip() or None
            except OSError:
                return None
            return None
        if current == stop or current.parent == current:
            break
        current = current.parent
    return None


def detect_cpu_model() -> str:
    proc_cpuinfo = Path("/proc/cpuinfo")
    if proc_cpuinfo.is_file():
        try:
            for line in proc_cpuinfo.read_text(
                encoding="utf-8", errors="replace"
            ).splitlines():
                if line.lower().startswith("model name"):
                    _, _, value = line.partition(":")
                    model = value.strip()
                    if model:
                        return model
        except OSError:
            pass
    fallback = platform.processor().strip()
    return fallback if fallback else "unknown"


def _run_git(repo_root: Path, *args: str) -> str | None:
    try:
        proc = subprocess.run(
            ["git", *args],
            cwd=repo_root,
            capture_output=True,
            text=True,
            timeout=5.0,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if proc.returncode != 0:
        return None
    return proc.stdout.strip()


def collect_env_metadata(
    plan: BenchmarkPlan,
    executable: Path,
    working_dir: Path,
    repo_root: Path,
    cmake_build_type: str | None,
    jobs: int,
) -> dict[str, Any]:
    git_commit = _run_git(repo_root, "rev-parse", "HEAD")
    git_branch = _run_git(repo_root, "rev-parse", "--abbrev-ref", "HEAD")
    git_status = _run_git(repo_root, "status", "--porcelain")
    git_dirty = bool(git_status) if git_status is not None else None

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "experiment_name": plan.experiment.name,
        "config_source": str(plan.source_path),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "cpu_model": detect_cpu_model(),
        "cpu_count": os.cpu_count(),
        "git_commit": git_commit,
        "git_branch": git_branch,
        "git_dirty": git_dirty,
        "jobs": jobs,
        "openblas_threads": plan.experiment.openblas_threads,
        "omp_threads": plan.experiment.omp_threads,
        "executable": str(executable),
        "working_dir": str(working_dir),
        "cmake_build_type": cmake_build_type,
    }


def run_preflight_checks(
    plan: BenchmarkPlan,
    executable: Path,
    working_dir: Path,
    repo_root: Path,
) -> tuple[list[str], str | None]:
    errors: list[str] = []

    if not executable.is_file():
        errors.append(f"Executable not found: {executable}")
    elif not os.access(executable, os.X_OK):
        errors.append(f"Executable is not marked executable: {executable}")

    if not working_dir.is_dir():
        errors.append(f"Working directory does not exist: {working_dir}")

    seen_meshes: set[str] = set()
    for spec in plan.runs:
        mesh = spec.params.get("mesh")
        if not isinstance(mesh, str) or mesh in seen_meshes:
            continue
        seen_meshes.add(mesh)
        mesh_path = Path(mesh) if Path(mesh).is_absolute() else (working_dir / mesh)
        if not mesh_path.is_file():
            errors.append(f"Mesh file not found: {mesh_path} (from mesh={mesh!r})")

    cmake_build_type = detect_cmake_build_type(executable, repo_root)
    if cmake_build_type is not None and cmake_build_type not in {
        "Release",
        "RelWithDebInfo",
    }:
        print(
            f"[WARNING] CMake build type is {cmake_build_type!r} "
            "(expected 'Release' or 'RelWithDebInfo'); timings may be distorted.",
            file=sys.stderr,
        )

    return errors, cmake_build_type


def build_pinned_env(exp: ExperimentConfig) -> dict[str, str]:
    env = os.environ.copy()
    blas_threads = str(exp.openblas_threads)
    omp_threads = str(exp.omp_threads)
    env["OPENBLAS_NUM_THREADS"] = blas_threads
    env["OMP_NUM_THREADS"] = omp_threads
    env["MKL_NUM_THREADS"] = blas_threads
    env["NUMEXPR_NUM_THREADS"] = blas_threads
    return env


def _format_scalar_compact(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:g}"
    return str(value)


def format_run_prefix(
    idx: int,
    total: int,
    spec: RunSpec,
    varying_keys: list[str],
) -> str:
    """Format a compact progress prefix based on varying keys (or core fallback)."""
    width = len(str(max(total, 1)))
    counter = f"[{idx:>{width}}/{total}]"

    if varying_keys:
        keys = list(varying_keys)
    else:
        keys = [k for k in COMMON_KEY_ORDER if k in spec.params]
        if not keys and "mesh" in spec.params:
            keys.append("mesh")
        if not keys:
            keys = sorted(spec.params.keys())[:6]

    key_set = set(keys)
    segments: list[str] = []
    consumed: set[str] = set()

    # Segment 1: mesh / problem / ref_strategy
    geo_parts: list[str] = []
    if "mesh_name" in key_set and "mesh_name" in spec.params:
        geo_parts.append(str(spec.params["mesh_name"]))
        consumed.add("mesh_name")
    elif "mesh" in key_set and "mesh" in spec.params:
        geo_parts.append(Path(str(spec.params["mesh"])).stem)
        consumed.add("mesh")
    if "problem" in key_set and "problem" in spec.params:
        geo_parts.append(str(spec.params["problem"]))
        consumed.add("problem")
    if "ref_strategy" in key_set and "ref_strategy" in spec.params:
        geo_parts.append(str(spec.params["ref_strategy"]))
        consumed.add("ref_strategy")
    if geo_parts:
        segments.append(" ".join(geo_parts))

    # Segment 2: discretization (order + ref_levels)
    disc_parts: list[str] = []
    if "order" in key_set and "order" in spec.params:
        disc_parts.append(f"p={spec.params['order']}")
        consumed.add("order")
    if "ref_levels" in key_set and "ref_levels" in spec.params:
        disc_parts.append(f"r={spec.params['ref_levels']}")
        consumed.add("ref_levels")
    if disc_parts:
        segments.append(" ".join(disc_parts))

    # Segment 3: assembly
    if "assembly" in key_set and "assembly" in spec.params:
        segments.append(str(spec.params["assembly"]))
        consumed.add("assembly")

    # Segment 4: solver + prec
    has_solver = "solver" in key_set and "solver" in spec.params
    has_prec = "prec" in key_set and "prec" in spec.params
    if has_solver and has_prec:
        segments.append(f"{spec.params['solver']}+{spec.params['prec']}")
        consumed.update({"solver", "prec"})
    elif has_solver:
        segments.append(str(spec.params["solver"]))
        consumed.add("solver")
    elif has_prec:
        segments.append(str(spec.params["prec"]))
        consumed.add("prec")

    # Remaining arbitrary varying keys in deterministic order
    for key in keys:
        if key in consumed or key not in spec.params:
            continue
        segments.append(f"{key}={_format_scalar_compact(spec.params[key])}")

    if not segments:
        return f"{counter} {spec.run_id}"
    return f"{counter} {' | '.join(segments)}"


def _coerce_text(data: str | bytes | None) -> str:
    if data is None:
        return ""
    if isinstance(data, bytes):
        return data.decode("utf-8", errors="replace")
    return data


def _build_failure_stub(
    spec: RunSpec,
    status: str,
    returncode: int | None,
    message: str,
    stderr_text: str,
    wall_sec: float,
    cmd: list[str],
) -> dict[str, Any]:
    return {
        "schema_version": "1.0",
        "run_id": spec.run_id,
        "status": status,
        "config": {**spec.params, "sweep_name": spec.sweep_name},
        "problem_stats": None,
        "timings_sec": {
            "wall_runner_sec": round(wall_sec, 6),
        },
        "convergence": None,
        "verification": None,
        "error": {
            "returncode": returncode,
            "message": message,
            "stderr_tail": stderr_text[-_STDERR_TAIL_LIMIT:],
            "command": cmd,
        },
    }


def _is_cached_ok(out_file: Path) -> bool:
    if not out_file.is_file():
        return False
    try:
        data = json.loads(out_file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return isinstance(data, dict) and data.get("status") == "ok"


def _format_ok_metrics(data: dict[str, Any], wall_sec: float) -> str:
    parts: list[str] = []
    timings = data.get("timings_sec")
    if isinstance(timings, dict):
        total = timings.get("total_sec")
        if isinstance(total, (int, float)):
            parts.append(f"total={total:.3f}s")
        else:
            parts.append(f"wall={wall_sec:.3f}s")
        asm = timings.get("assemble_sec")
        if isinstance(asm, (int, float)):
            parts.append(f"asm={asm:.3f}s")
        solve = timings.get("solve_median_sec")
        if not isinstance(solve, (int, float)):
            solve = timings.get("solve_sec")
        if isinstance(solve, (int, float)):
            parts.append(f"solve={solve:.3f}s")
    else:
        parts.append(f"wall={wall_sec:.3f}s")

    conv = data.get("convergence")
    if isinstance(conv, dict):
        iters = conv.get("iterations")
        if isinstance(iters, int):
            parts.append(f"iters={iters}")

    return " ".join(parts)


def _create_result_dir(output_root: Path, experiment_name: str) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_name = _SAFE_DIR_RE.sub("_", experiment_name.strip()).strip("_") or "exp"
    base_name = f"{ts}_{safe_name}"
    candidate = output_root / base_name
    suffix = 1
    while candidate.exists():
        candidate = output_root / f"{base_name}_{suffix}"
        suffix += 1
    candidate.mkdir(parents=True, exist_ok=False)
    (candidate / "runs").mkdir(parents=True, exist_ok=True)
    return candidate


def _relative_display_path(target: Path, repo_root: Path) -> str:
    try:
        return str(target.relative_to(Path.cwd().resolve()))
    except ValueError:
        pass
    try:
        return str(target.relative_to(repo_root))
    except ValueError:
        return str(target)


def _echo_verbose_output(
    stdout_text: str,
    stderr_text: str,
    print_lock: threading.Lock | None,
) -> None:
    if not stdout_text and not stderr_text:
        return
    if print_lock is not None:
        with print_lock:
            if stdout_text:
                print(stdout_text, end="", flush=True)
            if stderr_text:
                print(stderr_text, end="", file=sys.stderr, flush=True)
    else:
        if stdout_text:
            print(stdout_text, end="", flush=True)
        if stderr_text:
            print(stderr_text, end="", file=sys.stderr, flush=True)


def _terminate_active_procs(
    active_procs: set[subprocess.Popen[str]],
    procs_lock: threading.Lock,
    grace_sec: float = 2.0,
) -> None:
    with procs_lock:
        procs = list(active_procs)

    for proc in procs:
        try:
            proc.terminate()
        except OSError:
            pass

    deadline = time.monotonic() + grace_sec
    for proc in procs:
        remaining = max(0.0, deadline - time.monotonic())
        try:
            proc.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            try:
                proc.kill()
            except OSError:
                pass
        except OSError:
            pass


def _execute_single_run(
    spec: RunSpec,
    executable: Path,
    working_dir: Path,
    out_file: Path,
    env: dict[str, str],
    timeout_sec: float,
    verbose_cpp: bool,
    stop_event: threading.Event | None = None,
    active_procs: set[subprocess.Popen[str]] | None = None,
    procs_lock: threading.Lock | None = None,
    print_lock: threading.Lock | None = None,
) -> tuple[str, str]:
    """Execute one C++ benchmark run and persist its result or failure stub.

    Returns (status, short_detail_string).
    """
    if stop_event is not None and stop_event.is_set():
        return "cancelled", "CANCELLED"

    tmp_out = out_file.with_name(
        f".{spec.run_id}.{os.getpid()}.{threading.get_ident()}.cpp.tmp.json"
    )
    _unlink_quiet(tmp_out)

    cmd = build_cli_args(
        executable=executable,
        spec=spec,
        out_file=tmp_out,
        quiet=not verbose_cpp,
    )

    t0 = time.perf_counter()
    proc: subprocess.Popen[str] | None = None
    timed_out = False
    stdout_text = ""
    stderr_text = ""

    try:
        if procs_lock is not None and active_procs is not None:
            with procs_lock:
                if stop_event is not None and stop_event.is_set():
                    return "cancelled", "CANCELLED"
                proc = subprocess.Popen(
                    cmd,
                    cwd=working_dir,
                    env=env,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                active_procs.add(proc)
        else:
            proc = subprocess.Popen(
                cmd,
                cwd=working_dir,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

        try:
            raw_out, raw_err = proc.communicate(timeout=timeout_sec)
            stdout_text = _coerce_text(raw_out)
            stderr_text = _coerce_text(raw_err)
        except subprocess.TimeoutExpired:
            timed_out = True
            try:
                proc.kill()
            except OSError:
                pass
            try:
                raw_out, raw_err = proc.communicate(timeout=5.0)
                stdout_text = _coerce_text(raw_out)
                stderr_text = _coerce_text(raw_err)
            except subprocess.TimeoutExpired:
                pass
        wall_sec = time.perf_counter() - t0
    except KeyboardInterrupt:
        if proc is not None:
            try:
                proc.kill()
                proc.communicate(timeout=5.0)
            except (OSError, subprocess.TimeoutExpired):
                pass
        _unlink_quiet(tmp_out)
        raise
    finally:
        if proc is not None and procs_lock is not None and active_procs is not None:
            with procs_lock:
                active_procs.discard(proc)

    # If aborted via Ctrl+C while running, discard temp output instead of writing a failure stub.
    if stop_event is not None and stop_event.is_set():
        _unlink_quiet(tmp_out)
        return "cancelled", "CANCELLED"

    if verbose_cpp:
        _echo_verbose_output(stdout_text, stderr_text, print_lock)

    if timed_out:
        _unlink_quiet(tmp_out)
        stub = _build_failure_stub(
            spec=spec,
            status="timeout",
            returncode=None,
            message=f"Process exceeded timeout of {timeout_sec:.1f}s",
            stderr_text=stderr_text,
            wall_sec=wall_sec,
            cmd=cmd,
        )
        atomic_write_json(out_file, stub)
        return "timeout", f"TIMEOUT ({wall_sec:.1f}s)"

    # Non-zero exit code -> OOM or crash stub
    if proc.returncode != 0:
        _unlink_quiet(tmp_out)
        if proc.returncode in (-9, 137):
            status = "oom"
            msg = f"Process killed (likely OOM, rc={proc.returncode})"
            detail = f"OOM (rc={proc.returncode}, {wall_sec:.2f}s)"
        else:
            status = "crash"
            msg = f"Process exited with non-zero return code {proc.returncode}"
            detail = f"CRASH (rc={proc.returncode}, {wall_sec:.2f}s)"
        stub = _build_failure_stub(
            spec=spec,
            status=status,
            returncode=proc.returncode,
            message=msg,
            stderr_text=stderr_text,
            wall_sec=wall_sec,
            cmd=cmd,
        )
        atomic_write_json(out_file, stub)
        return status, detail

    # Exit code 0 -> verify JSON exists and parses cleanly
    if not tmp_out.is_file():
        stub = _build_failure_stub(
            spec=spec,
            status="crash",
            returncode=0,
            message="Process exited 0 but did not produce JSON output file",
            stderr_text=stderr_text,
            wall_sec=wall_sec,
            cmd=cmd,
        )
        atomic_write_json(out_file, stub)
        return "crash", f"CRASH (missing JSON, {wall_sec:.2f}s)"

    try:
        raw_text = tmp_out.read_text(encoding="utf-8")
        data = json.loads(raw_text)
        if not isinstance(data, dict):
            raise ValueError("JSON root is not an object")
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        _unlink_quiet(tmp_out)
        stub = _build_failure_stub(
            spec=spec,
            status="crash",
            returncode=0,
            message=f"Process exited 0 but produced invalid JSON: {exc}",
            stderr_text=stderr_text,
            wall_sec=wall_sec,
            cmd=cmd,
        )
        atomic_write_json(out_file, stub)
        return "crash", f"CRASH (invalid JSON, {wall_sec:.2f}s)"

    # Merge TOML params with C++ resolved config (C++ wins on overlapping keys)
    existing_cfg = data.get("config")
    cpp_config: dict[str, Any] = (
        existing_cfg if isinstance(existing_cfg, dict) else {}
    )
    merged_config: dict[str, Any] = {
        **spec.params,
        **cpp_config,
        "sweep_name": spec.sweep_name,
    }

    if "config" not in data or merged_config != existing_cfg:
        data["config"] = merged_config
        atomic_write_json(out_file, data)
        _unlink_quiet(tmp_out)
    else:
        os.replace(tmp_out, out_file)

    status_str = str(data.get("status", "ok"))
    metrics = _format_ok_metrics(data, wall_sec)
    return status_str, f"{status_str.upper()}  {metrics}"


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    repo_root = detect_repo_root(args.repo_root)

    try:
        plan = load_and_expand_config(args.config)
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"[ERROR] Failed to load config: {exc}", file=sys.stderr)
        return 1

    exp = plan.experiment
    executable = resolve_under(repo_root, exp.executable)
    working_dir = resolve_under(repo_root, exp.working_dir)
    output_root = resolve_under(repo_root, exp.output_root)

    total = len(plan.runs)
    print(f"Experiment : {exp.name}")
    print(f"Config     : {plan.source_path}")
    print(
        f"Runs       : {total} "
        f"(excluded: {plan.excluded_count})"
    )
    print(
        f"Varying    : {', '.join(plan.varying_keys) if plan.varying_keys else '(none)'}"
    )
    print(
        f"Threads    : jobs={args.jobs}, OPENBLAS={exp.openblas_threads}, "
        f"OMP={exp.omp_threads}, timeout={exp.timeout_sec:g}s"
    )

    if args.jobs > 1:
        threads_per_job = max(exp.openblas_threads, exp.omp_threads, 1)
        total_threads = args.jobs * threads_per_job
        cpu_count = os.cpu_count()
        if cpu_count is not None and total_threads > cpu_count:
            print(
                f"[WARNING] Running {args.jobs} parallel jobs x {threads_per_job} "
                f"threads/job = {total_threads} threads on {cpu_count} logical CPUs; "
                "oversubscription may severely distort timings.",
                file=sys.stderr,
            )
        else:
            print(
                f"[WARNING] Running {args.jobs} jobs in parallel; shared cache and "
                "memory bandwidth contention may introduce timing noise.",
                file=sys.stderr,
            )

    preflight_errors, cmake_build_type = run_preflight_checks(
        plan=plan,
        executable=executable,
        working_dir=working_dir,
        repo_root=repo_root,
    )

    if args.dry_run:
        print("\n[DRY-RUN] Planned runs:")
        sample_runs_dir = output_root / "<timestamp>_dry_run" / "runs"
        for idx, spec in enumerate(plan.runs, start=1):
            prefix = format_run_prefix(idx, total, spec, plan.varying_keys)
            print(f"  {prefix}  ({spec.run_id})")
        if plan.runs:
            sample_spec = plan.runs[0]
            sample_cmd = build_cli_args(
                executable=executable,
                spec=sample_spec,
                out_file=sample_runs_dir / f"{sample_spec.run_id}.json",
                quiet=not args.verbose_cpp,
            )
            print("\n[DRY-RUN] Sample command (run 1):")
            print(f"  {shlex.join(sample_cmd)}")
        if preflight_errors:
            print("\n[PRE-FLIGHT FAILED]", file=sys.stderr)
            for err in preflight_errors:
                print(f"  - {err}", file=sys.stderr)
            return 1
        print("\n[DRY-RUN] Pre-flight checks passed.")
        return 0

    if preflight_errors:
        print("[PRE-FLIGHT FAILED]", file=sys.stderr)
        for err in preflight_errors:
            print(f"  - {err}", file=sys.stderr)
        return 1

    if total == 0:
        print("[WARNING] Expanded benchmark plan contains 0 runs. Nothing to do.")
        return 0

    if args.resume is not None:
        result_dir = args.resume.resolve()
        if not result_dir.is_dir():
            print(
                f"[ERROR] Resume directory does not exist: {result_dir}",
                file=sys.stderr,
            )
            return 1
        runs_dir = result_dir / "runs"
        runs_dir.mkdir(parents=True, exist_ok=True)
        print(f"Resuming   : {result_dir}")
    else:
        result_dir = _create_result_dir(output_root, exp.name)
        runs_dir = result_dir / "runs"
        print(f"Output dir : {result_dir}")

    snapshot_path = result_dir / "config_snapshot.toml"
    if not snapshot_path.exists():
        shutil.copy2(plan.source_path, snapshot_path)

    env_meta_path = result_dir / "env_metadata.json"
    if not env_meta_path.exists():
        env_metadata = collect_env_metadata(
            plan=plan,
            executable=executable,
            working_dir=working_dir,
            repo_root=repo_root,
            cmake_build_type=cmake_build_type,
            jobs=args.jobs,
        )
        atomic_write_json(env_meta_path, env_metadata)

    env = build_pinned_env(exp)
    counts: dict[str, int] = {
        "ok": 0,
        "cached": 0,
        "timeout": 0,
        "oom": 0,
        "crash": 0,
        "other": 0,
    }
    suite_t0 = time.perf_counter()
    print()

    jobs_flag = f" -j {args.jobs}" if args.jobs > 1 else ""

    try:
        if args.jobs == 1:
            for idx, spec in enumerate(plan.runs, start=1):
                prefix = format_run_prefix(idx, total, spec, plan.varying_keys)
                out_file = runs_dir / f"{spec.run_id}.json"

                if args.resume is not None and _is_cached_ok(out_file):
                    counts["cached"] += 1
                    print(f"{prefix} ... CACHED", flush=True)
                    continue

                print(f"{prefix} ... ", end="", flush=True)
                status, detail = _execute_single_run(
                    spec=spec,
                    executable=executable,
                    working_dir=working_dir,
                    out_file=out_file,
                    env=env,
                    timeout_sec=exp.timeout_sec,
                    verbose_cpp=args.verbose_cpp,
                )
                if status in counts:
                    counts[status] += 1
                else:
                    counts["other"] += 1
                print(detail, flush=True)
        else:
            pending: list[tuple[str, RunSpec, Path]] = []
            for idx, spec in enumerate(plan.runs, start=1):
                prefix = format_run_prefix(idx, total, spec, plan.varying_keys)
                out_file = runs_dir / f"{spec.run_id}.json"
                if args.resume is not None and _is_cached_ok(out_file):
                    counts["cached"] += 1
                    print(f"{prefix} ... CACHED", flush=True)
                else:
                    pending.append((prefix, spec, out_file))

            if pending:
                stop_event = threading.Event()
                active_procs: set[subprocess.Popen[str]] = set()
                procs_lock = threading.Lock()
                print_lock = threading.Lock()
                max_workers = min(args.jobs, len(pending))

                executor = ThreadPoolExecutor(max_workers=max_workers)
                try:
                    future_to_prefix: dict[Future[tuple[str, str]], str] = {
                        executor.submit(
                            _execute_single_run,
                            spec=spec,
                            executable=executable,
                            working_dir=working_dir,
                            out_file=out_file,
                            env=env,
                            timeout_sec=exp.timeout_sec,
                            verbose_cpp=args.verbose_cpp,
                            stop_event=stop_event,
                            active_procs=active_procs,
                            procs_lock=procs_lock,
                            print_lock=print_lock,
                        ): prefix
                        for prefix, spec, out_file in pending
                    }

                    for future in as_completed(future_to_prefix):
                        prefix = future_to_prefix[future]
                        status, detail = future.result()
                        if status == "cancelled":
                            continue
                        if status in counts:
                            counts[status] += 1
                        else:
                            counts["other"] += 1
                        with print_lock:
                            print(f"{prefix} ... {detail}", flush=True)
                except KeyboardInterrupt:
                    stop_event.set()
                    _terminate_active_procs(active_procs, procs_lock)
                    executor.shutdown(wait=True, cancel_futures=True)
                    raise
                else:
                    executor.shutdown(wait=True)

    except KeyboardInterrupt:
        elapsed = time.perf_counter() - suite_t0
        rel_dir = _relative_display_path(result_dir, repo_root)
        print(
            f"\n[INTERRUPTED] Stopped after {elapsed:.1f}s "
            f"(ok={counts['ok']}, cached={counts['cached']}, "
            f"failed={counts['timeout'] + counts['oom'] + counts['crash'] + counts['other']}).",
            file=sys.stderr,
        )
        print(
            f"Resume with: python3 tests/benchmarks/run_benchmarks.py "
            f"{args.config}{jobs_flag} --resume {rel_dir}",
            file=sys.stderr,
        )
        return 130

    elapsed = time.perf_counter() - suite_t0
    failed = counts["timeout"] + counts["oom"] + counts["crash"] + counts["other"]
    rel_dir = _relative_display_path(result_dir, repo_root)

    print(
        f"\nCompleted in {elapsed:.1f}s | "
        f"ok={counts['ok']}, cached={counts['cached']}, "
        f"timeout={counts['timeout']}, oom={counts['oom']}, "
        f"crash={counts['crash']}"
        + (f", other={counts['other']}" if counts["other"] else "")
    )
    if failed > 0:
        print(
            f"Re-run failed cases: python3 tests/benchmarks/run_benchmarks.py "
            f"{args.config}{jobs_flag} --resume {rel_dir}"
        )
    print(f"python3 tests/benchmarks/analyze_results.py {rel_dir}")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
