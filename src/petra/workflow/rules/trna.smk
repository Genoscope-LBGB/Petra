rule trnascan:
    """Annotate tRNA genes across the whole assembly.

    -G uses the general model, the only appropriate one for an assembly holding
    sequences of several domains. The scan runs on the assembly rather than on
    the cmscan chunks because tRNAscan-SE reports its statistics per run.
    """
    input:
        assembly=config["assembly"],
    output:
        predictions=temp("rna/predictions.trna.raw.tsv"),
        stats="rna/predictions.trna.stats",
    log:
        "logs/trnascan.log",
    benchmark:
        "benchmarks/trnascan.tsv"
    shadow:
        "minimal"
    shell:
        "tRNAscan-SE -G --thread {threads} -o {output.predictions:q} "
        "--stats {output.stats:q} {input.assembly:q} > {log:q} 2>&1"


rule filter_trna_predictions:
    """Drop pseudogene, undetermined and suppressor predictions."""
    input:
        predictions="rna/predictions.trna.raw.tsv",
    output:
        predictions="rna/predictions.trna.tsv",
    log:
        "logs/filter_trna_predictions.log",
    benchmark:
        "benchmarks/filter_trna_predictions.tsv"
    script:
        "../scripts/filter_trna_predictions.py"
