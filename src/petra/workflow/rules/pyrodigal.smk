rule pyrodigal:
    """
    -m masks runs of N so genes cannot span them
    -c prevents genes from running off sequence edges
    -p meta uses pre-trained models from different organisms instead of training on the input
    """
    input:
        assembly=config["assembly"]
    output:
        "predictions.gff"
    shell:
        "pyrodigal -m -c -f gff -p meta -j {threads} "
        "-i {input.assembly:q} -o {output:q}"
