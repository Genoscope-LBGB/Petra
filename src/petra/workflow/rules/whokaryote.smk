rule whokaryote:
    """Classify assembly sequences as eukaryotic or prokaryotic."""
    input:
        assembly=config["assembly"],
        gff="predictions.gff",
    output:
        predictions=temp("predictions.whokaryote.tsv"),
        tiara_predictions=temp("predictions.tiara.tsv"),
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
        "{output.predictions:q} && "
        "cp {params.output_directory:q}/tiara_pred.txt "
        "{output.tiara_predictions:q}"


rule extract_organelles:
    """List sequences that Tiara classified as organelles."""
    input:
        predictions="predictions.tiara.tsv",
    output:
        organelles="predictions.organelles.tsv",
    log:
        "logs/extract_organelles.log",
    benchmark:
        "benchmarks/extract_organelles.tsv"
    script:
        "../scripts/extract_organelles.py"


rule control_taxonomy_with_whokaryote:
    """Report PETRA assignments that disagree with Whokaryote."""
    input:
        raw_taxonomy="predictions.taxonomy.raw.tsv",
        whokaryote_predictions="predictions.whokaryote.tsv",
    output:
        reasons=temp("predictions.unassigned_reasons.whokaryote.tsv"),
    log:
        "logs/control_taxonomy_with_whokaryote.log",
    benchmark:
        "benchmarks/control_taxonomy_with_whokaryote.tsv"
    script:
        "../scripts/control_taxonomy_with_whokaryote.py"


rule aggregate_taxonomy:
    """Combine raw taxonomy and tool-specific reasons into the final table."""
    input:
        raw_taxonomy="predictions.taxonomy.raw.tsv",
        coverage="predictions.coverage.tsv",
        assignment_reasons="predictions.unassigned_reasons.assign_taxonomy.tsv",
        whokaryote_reasons="predictions.unassigned_reasons.whokaryote.tsv",
    output:
        taxonomy="predictions.final.tsv",
    log:
        "logs/aggregate_taxonomy.log",
    benchmark:
        "benchmarks/aggregate_taxonomy.tsv"
    script:
        "../scripts/aggregate_taxonomy.py"


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
