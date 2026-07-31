# PETRA - Predicted protein-enabled Taxonomic Rank Assignment

## Installation

```console
pip install .
```

## Usage

### Running the workflow

Supply every long-read file by repeating `--long-reads`. The sequencing technology
defaults to `pacbio` (`map-hifi`), select `ont` to use `map-ont`:

```console
petra run \
  --assembly assembly.fasta \
  --diamond-db database/petra.dmnd \
  --long-reads reads-1.fastq.gz \
  --long-reads reads-2.fastq.gz \
  --technology ont \
  --output run-1
```

### Building the local DIAMOND database

Build a taxonomy-aware DIAMOND database from local inputs:

```console
petra build-db \
  --fasta proteins.fasta \
  --taxon-map prot.accession2taxid.FULL.gz \
  --output database
```

The command requires `diamond` on `PATH` and creates `database/petra.dmnd`. If
an input is omitted, it downloads NCBI NR or the full NCBI protein taxon map
into the output directory, logs the source and destination, and removes the
download after the build. It fails rather than replacing an existing download
or `petra.dmnd`.
