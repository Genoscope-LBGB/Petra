# PETRA - Predicted protein-enabled Taxonomic Rank Assignment

## Installation

```console
pip install .
```

## Usage

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
