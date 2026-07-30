from __future__ import annotations

import os
import shlex
import subprocess
import shutil
import sys
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from enum import Enum
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import yaml
from rich.console import Console

_RESERVED_CONFIG_KEY = "_snakemake"
_PROFILE_FILENAMES = ("profile.v9+.yaml", "profile.yaml", "config.yaml")
_PROTECTED_LONG_OPTIONS = {
    "--config",
    "--configfile",
    "--configfiles",
    "--cores",
    "--directory",
    "--dry-run",
    "--executor",
    "--jobs",
    "--printshellcmds",
    "--sdm",
    "--snakefile",
    "--skip-script-cleanup",
    "--software-deployment-method",
    "--use-apptainer",
    "--use-conda",
    "--use-singularity",
    "--workflow-profile",
}
_PROTECTED_SHORT_OPTIONS = {"-c", "-d", "-e", "-j", "-n", "-p", "-s"}


class LauncherError(RuntimeError):
    """An expected launcher failure suitable for CLI display."""


@dataclass(frozen=True)
class WorkflowProfile:
    name: str
    path: Path
    config_file: Path
    config: dict[str, Any]


@dataclass(frozen=True)
class ProfileSummary:
    name: str
    jobs: str
    cores: str
    partition: str
    qos: str
    runtime: str
    memory_per_cpu: str
    path: str


def workflow_directory() -> Path:
    return Path(__file__).resolve().parent / "workflow"


def resolve_workflow_profile(
    selection: str,
    *,
    invocation_directory: Path,
) -> WorkflowProfile:
    packaged_root = workflow_directory() / "profiles"
    if Path(selection).name == selection:
        packaged = packaged_root / selection
        if packaged.exists():
            return _load_profile(selection, packaged)

    external = Path(selection).expanduser()
    if not external.is_absolute():
        external = invocation_directory / external
    external = external.resolve()
    if not external.exists():
        available = ", ".join(profile.name for profile in packaged_profiles())
        suffix = f" Packaged profiles: {available}." if available else ""
        raise LauncherError(f"Workflow profile not found: {selection}.{suffix}")
    name = external.parent.name if external.is_file() else external.name
    return _load_profile(name, external)


def packaged_profiles() -> list[WorkflowProfile]:
    root = workflow_directory() / "profiles"
    if not root.exists():
        return []
    return [
        _load_profile(path.name, path)
        for path in sorted(root.iterdir(), key=lambda item: item.name)
        if path.is_dir()
    ]


def profile_summaries() -> list[ProfileSummary]:
    return [_summarize_profile(profile) for profile in packaged_profiles()]


def reject_passthrough_conflicts(arguments: Sequence[str]) -> None:
    for argument in arguments:
        if argument == "--":
            continue
        if argument.startswith("--"):
            option = argument.split("=", 1)[0]
            if option in _PROTECTED_LONG_OPTIONS:
                raise LauncherError(
                    f"Passthrough option {option!r} conflicts with a wrapper-owned "
                    "setting. Use the corresponding 'run --help' option."
                )
            continue
        for option in _PROTECTED_SHORT_OPTIONS:
            if argument == option or argument.startswith(option):
                raise LauncherError(
                    f"Passthrough option {argument!r} conflicts with wrapper-owned "
                    f"option {option!r}. Use the corresponding 'run --help' option."
                )


def normalize_config(value: Any, *, invocation_directory: Path) -> Any:
    return _normalize(value, path="$", invocation_directory=invocation_directory)


def launch_workflow(
    *,
    distribution_name: str,
    output: Path,
    workflow_config: Mapping[str, Any],
    profile: str | None,
    executor: str | None,
    jobs: int | None,
    cores: int | None,
    dry_run: bool,
    use_conda: bool,
    use_apptainer: bool,
    passthrough: Sequence[str],
    console: Console | None = None,
) -> int:
    invocation_directory = Path.cwd().resolve()
    output_directory = output.expanduser()
    if not output_directory.is_absolute():
        output_directory = invocation_directory / output_directory
    output_directory = output_directory.resolve()

    reject_passthrough_conflicts(passthrough)
    if _RESERVED_CONFIG_KEY in workflow_config:
        raise LauncherError(
            f"Workflow config key {_RESERVED_CONFIG_KEY!r} is reserved by the launcher."
        )

    selected_profile = (
        resolve_workflow_profile(profile, invocation_directory=invocation_directory)
        if profile is not None
        else None
    )
    effective_cores = (
        1
        if executor == "local" and cores is None and selected_profile is None
        else cores
    )
    normalized_workflow_config = normalize_config(
        workflow_config,
        invocation_directory=invocation_directory,
    )
    if not isinstance(normalized_workflow_config, dict):
        raise LauncherError("Workflow config must be a mapping with string keys.")

    output_directory.mkdir(parents=True, exist_ok=True)
    config_path = output_directory / "workflow-config.yaml"
    snakefile = workflow_directory() / "Snakefile"
    if not snakefile.is_file():
        raise LauncherError(f"Packaged Snakefile not found: {snakefile}")

    try:
        workflow_version = version(distribution_name)
    except PackageNotFoundError:
        workflow_version = "0+unknown"

    persisted_config = dict(normalized_workflow_config)
    persisted_config[_RESERVED_CONFIG_KEY] = {
        "generated_at": datetime.now(UTC).isoformat(timespec="seconds").replace(
            "+00:00", "Z"
        ),
        "workflow_version": workflow_version,
        "output": str(output_directory),
        "workflow_profile": selected_profile.name if selected_profile else None,
        "workflow_profile_path": (
            str(selected_profile.path) if selected_profile else None
        ),
        "workflow_profile_config": (
            normalize_config(
                selected_profile.config,
                invocation_directory=invocation_directory,
            )
            if selected_profile
            else None
        ),
        "executor": executor,
        "jobs": jobs,
        "cores": effective_cores,
        "dry_run": dry_run,
        "use_conda": use_conda,
        "use_apptainer": use_apptainer,
        "printshellcmds": True,
        "skip_script_cleanup": True,
    }
    _write_config_atomically(config_path, persisted_config)

    wrapper_arguments = _wrapper_arguments(
        snakefile=snakefile,
        config_path=config_path,
        profile=selected_profile,
        executor=executor,
        jobs=jobs,
        cores=effective_cores,
        dry_run=dry_run,
        use_conda=use_conda,
        use_apptainer=use_apptainer,
    )
    output_console = console or Console()
    output_console.print(
        f"Snakemake command: {shlex.join(wrapper_arguments)}",
        markup=False,
    )
    if passthrough:
        output_console.print(
            f"{len(passthrough)} passthrough argument(s) omitted from display.",
            style="dim",
        )

    try:
        process = subprocess.run(
            [*wrapper_arguments, *passthrough],
            cwd=output_directory,
            check=False,
        )
    except FileNotFoundError as error:
        raise LauncherError(
            "Could not launch Snakemake: the 'snakemake' executable was not found."
        ) from error
    except OSError as error:
        raise LauncherError(f"Could not launch Snakemake: {error}") from error

    if process.returncode < 0:
        return 128 - process.returncode
    return process.returncode


def _load_profile(name: str, path: Path) -> WorkflowProfile:
    config_file = _find_profile_file(path)
    try:
        parsed = yaml.safe_load(config_file.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as error:
        raise LauncherError(f"Could not read workflow profile {name!r}: {error}") from error
    if parsed is None:
        parsed = {}
    if not isinstance(parsed, dict) or not all(
        isinstance(key, str) for key in parsed
    ):
        raise LauncherError(
            f"Workflow profile {name!r} must contain a YAML mapping with string keys."
        )
    resolved_path = path.resolve()
    return WorkflowProfile(
        name=name,
        path=resolved_path,
        config_file=config_file.resolve(),
        config=parsed,
    )


def _find_profile_file(path: Path) -> Path:
    if path.is_file():
        if path.suffix.lower() not in {".yaml", ".yml"}:
            raise LauncherError(f"Workflow profile is not YAML: {path}")
        return path
    for filename in _PROFILE_FILENAMES:
        candidate = path / filename
        if candidate.is_file():
            return candidate
    expected = ", ".join(_PROFILE_FILENAMES)
    raise LauncherError(
        f"Workflow profile directory {path} contains none of: {expected}."
    )


def _summarize_profile(profile: WorkflowProfile) -> ProfileSummary:
    default_resources = _resource_mapping(profile.config.get("default-resources"))
    return ProfileSummary(
        name=profile.name,
        jobs=_display_value(profile.config.get("jobs")),
        cores=_display_value(profile.config.get("cores")),
        partition=_display_value(
            _first_value(default_resources, "slurm_partition", "partition")
        ),
        qos=_display_value(_first_value(default_resources, "slurm_qos", "qos")),
        runtime=_display_value(default_resources.get("runtime")),
        memory_per_cpu=_display_value(
            _first_value(default_resources, "mem_mb_per_cpu", "mem-per-cpu")
        ),
        path=str(profile.path),
    )


def _resource_mapping(value: Any) -> dict[str, Any]:
    if isinstance(value, dict):
        return value
    if isinstance(value, str):
        values = [value]
    elif isinstance(value, list):
        values = value
    else:
        return {}

    resources: dict[str, Any] = {}
    for item in values:
        if not isinstance(item, str):
            continue
        delimiter = "=" if "=" in item else ":" if ":" in item else None
        if delimiter:
            key, resource_value = item.split(delimiter, 1)
            resources[key.strip()] = resource_value.strip()
    return resources


def _first_value(mapping: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in mapping:
            return mapping[key]
    return None


def _display_value(value: Any) -> str:
    return "-" if value is None else str(value)


def _normalize(value: Any, *, path: str, invocation_directory: Path) -> Any:
    if isinstance(value, Enum):
        return _normalize(
            value.value,
            path=path,
            invocation_directory=invocation_directory,
        )
    if isinstance(value, Path):
        expanded = value.expanduser()
        if not expanded.is_absolute():
            expanded = invocation_directory / expanded
        return str(expanded.resolve())
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, Mapping):
        normalized: dict[str, Any] = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise LauncherError(f"Unsupported non-string config key at {path}: {key!r}")
            normalized[key] = _normalize(
                item,
                path=f"{path}.{key}",
                invocation_directory=invocation_directory,
            )
        return normalized
    if isinstance(value, (list, tuple)):
        return [
            _normalize(
                item,
                path=f"{path}[{index}]",
                invocation_directory=invocation_directory,
            )
            for index, item in enumerate(value)
        ]
    if isinstance(value, (set, frozenset)):
        normalized_items = [
            _normalize(
                item,
                path=f"{path}[*]",
                invocation_directory=invocation_directory,
            )
            for item in value
        ]
        return sorted(
            normalized_items,
            key=lambda item: yaml.safe_dump(item, sort_keys=True),
        )
    raise LauncherError(
        f"Unsupported config value at {path}: {type(value).__name__}."
    )


def _write_config_atomically(path: Path, config: Mapping[str, Any]) -> None:
    temporary_path: Path | None = None
    try:
        descriptor, temporary_name = tempfile.mkstemp(
            dir=path.parent,
            prefix=".workflow-config.",
            suffix=".tmp",
            text=True,
        )
        temporary_path = Path(temporary_name)
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            yaml.safe_dump(
                dict(config),
                handle,
                allow_unicode=True,
                default_flow_style=False,
                sort_keys=False,
            )
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    except (OSError, yaml.YAMLError) as error:
        raise LauncherError(f"Could not persist workflow config at {path}: {error}") from error
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _snakemake_executable() -> str:
    executable = shutil.which("snakemake")
    if executable is not None:
        return executable
    sibling = Path(sys.executable).parent / "snakemake"
    return str(sibling) if sibling.is_file() else "snakemake"


def _wrapper_arguments(
    *,
    snakefile: Path,
    config_path: Path,
    profile: WorkflowProfile | None,
    executor: str | None,
    jobs: int | None,
    cores: int | None,
    dry_run: bool,
    use_conda: bool,
    use_apptainer: bool,
) -> list[str]:
    arguments = [
        _snakemake_executable(),
        "--snakefile",
        str(snakefile.resolve()),
        "--configfile",
        str(config_path.resolve()),
        "-p",
        "--skip-script-cleanup",
    ]
    if profile is not None:
        arguments.extend(["--workflow-profile", str(profile.path)])
    else:
        arguments.extend(["--workflow-profile", "none"])
    if executor is not None:
        arguments.extend(["--executor", executor])
    if jobs is not None:
        arguments.extend(["--jobs", str(jobs)])
    if cores is not None:
        arguments.extend(["--cores", str(cores)])
    if dry_run:
        arguments.append("--dry-run")
    deployment_methods = []
    if use_conda:
        deployment_methods.append("conda")
    if use_apptainer:
        deployment_methods.append("apptainer")
    if deployment_methods:
        arguments.extend(["--software-deployment-method", *deployment_methods])
    return arguments
