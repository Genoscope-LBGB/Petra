# petra

PETRA - Predicted protein-enabled Taxonomic Rank Assignment

## Installation

```console
pip install .
```

For development, install the package in editable mode:

```console
pip install -e .
```

## Usage

The packaged `default` workflow profile uses Slurm with at most 100 concurrent jobs:

```console
petra run --output run-1
```

Run the starter workflow locally with:

```console
petra run \
  --output run-1 \
  --no-profile \
  --executor local
```

The starter writes `run-1/results/message.txt`. Change its message with `--message`.

Use `petra run --help` for workflow inputs and common execution controls. Targets and advanced Snakemake options follow `--`:

```console
petra run --output run-1 -- --targets hello
```

Raw passthrough values are neither displayed nor persisted because plugin arguments can contain credentials.

## Project structure

- `src/petra/cli.py`: workflow-specific Typer options and the explicit Snakemake config mapping.
- `src/petra/launcher.py`: profile, persistence, passthrough, and subprocess mechanics.
- `src/petra/workflow/Snakefile`: workflow entry point.
- `src/petra/workflow/rules/`: modular Snakemake rules.
- `src/petra/workflow/scripts/`: workflow scripts.
- `src/petra/workflow/envs/`: software environment definitions.
- `src/petra/workflow/profiles/`: packaged workflow-specific profiles.

To add a workflow option, declare a typed Typer parameter in `run()` and add its transformed value to `workflow_config`. The wrapper is the supported public interface; direct `snakemake` invocation does not receive the wrapper's friendly validation.

## Author

Benjamin Istace <benjamin.istace@genoscope.cns.fr>
