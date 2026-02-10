
rule prepare_unify_references:
    input:
        genome=config["genome_path"],
        transcriptome=config["transcriptome_path"],
    output:
        genome="results/prepare/unify_references/genome.fna",
        transcriptome="results/prepare/unify_references/transcriptome.gtf",
    log:
        "logs/prepare/unify_references.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/unify_references.R"


rule prepare_install_bambu:
    output:
        "results/prepare/install_bambu/done.txt",
    params:
        version=config["bambu_version"],
    log:
        "logs/prepare/install_bambu/out.log",
    conda:
        "../envs/r/bambu.yaml"
    script:
        "../scripts/r/install_bambu.R"


rule prepare_install_transmogr:
    output:
        "results/prepare/install_transmogr/done.txt",
    params:
        version=config["transmogr_version"],
    log:
        "logs/prepare/install_transmogr/out.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/install_transmogr.R"


rule prepare_install_rcppgreedysetcover:
    output:
        "results/prepare/install_rcppgreedysetcover/done.txt",
    params:
        version=config["rcppgreedysetcover_version"],
    log:
        "logs/prepare/install_rcppgreedysetcover/out.log",
    conda:
        "../envs/r/protein_inference.yaml"
    script:
        "../scripts/r/install_rcppgreedysetcover.R"
