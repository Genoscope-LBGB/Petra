rule whokaryote:
    """Classify assembly sequences as eukaryotic or prokaryotic."""
    input:
        assembly=config["assembly"],
        gff="predictions.gff",
    output:
        predictions="predictions.whokaryote.tsv",
    log:
        "logs/whokaryote.log",
    benchmark:
        "benchmarks/whokaryote.tsv"
    shadow:
        "minimal"
    params:
        output_directory="whokaryote",
    shell:
        "whokaryote.py --contigs {input.assembly:q} --gff {input.gff:q} "
        "--outdir {params.output_directory:q} --model T --minsize 5000 "
        "--threads {threads} > {log:q} 2>&1 && "
        "cp {params.output_directory:q}/whokaryote_predictions_T.tsv "
        "{output.predictions:q}"


rule control_taxonomy_with_whokaryote:
    """Remove PETRA assignments that disagree with Whokaryote."""
    input:
        raw_taxonomy="predictions.taxonomy.raw.tsv",
        whokaryote_predictions="predictions.whokaryote.tsv",
    output:
        taxonomy="predictions.final.tsv",
    log:
        "logs/control_taxonomy_with_whokaryote.log",
    benchmark:
        "benchmarks/control_taxonomy_with_whokaryote.tsv"
    script:
        "../scripts/control_taxonomy_with_whokaryote.py"


rule split_assembly_by_domain:
    """Split assembly sequences into FASTA files by final PETRA domain."""
    input:
        assembly=config["assembly"],
        taxonomy="predictions.final.tsv",
    output:
        eukaryota="predictions.eukaryota.fasta",
        bacteria="predictions.bacteria.fasta",
        viruses="predictions.viruses.fasta",
        archaea="predictions.archaea.fasta",
        unassigned="predictions.unassigned.fasta",
    log:
        "logs/split_assembly_by_domain.log",
    benchmark:
        "benchmarks/split_assembly_by_domain.tsv"
    script:
        "../scripts/split_assembly_by_domain.py"
