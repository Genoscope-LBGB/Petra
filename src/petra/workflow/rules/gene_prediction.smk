rule pyrodigal:
    """
    -m masks runs of N so genes cannot span them
    -c prevents genes from running off sequence edges
    -p meta uses pre-trained models from different organisms instead of training on the input
    """
    input:
        assembly=config["assembly"]
    output:
        gff=temp("predictions.gff"),
        proteins=temp("predictions.faa")
    benchmark:
        "benchmarks/pyrodigal.tsv"
    shell:
        "pyrodigal -m -c -f gff -p meta -j {threads} "
        "-i {input.assembly:q} -o {output.gff:q} -a {output.proteins:q}"

rule segmasker:
    """Mask low-complexity regions in the predicted protein sequences."""
    input:
        proteins="predictions.faa"
    output:
        proteins=temp("predictions.masked.faa")
    benchmark:
        "benchmarks/segmasker.tsv"
    shell:
        "segmasker -in {input.proteins:q} -out {output.proteins:q} "
        "-outfmt fasta"

rule filter_low_complexity:
    """Discard proteins exceeding the configured low-complexity threshold."""
    input:
        proteins="predictions.masked.faa"
    output:
        proteins=temp("predictions.filtered.faa")
    params:
        maximum=config["max_low_complexity"]
    benchmark:
        "benchmarks/filter_low_complexity.tsv"
    script:
        "../scripts/filter_low_complexity.py"

rule cdhit:
    """Cluster proteins independently within each originating assembly sequence and output at most ``max_clusters`` per assembly sequence."""
    input:
        proteins="predictions.filtered.faa"
    output:
        proteins=temp("predictions.clustered.faa")
    params:
        maximum=config["max_clusters_per_sequence"]
    benchmark:
        "benchmarks/cdhit.tsv"
    script:
        "../scripts/cluster_proteins.py"
