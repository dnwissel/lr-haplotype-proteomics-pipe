

rule discover_isoforms_run_bambu:
    input:
        "results/prepare/install_bambu/done.txt",
        reads=f"results/align/run_minimap2_genome/{config['sample_id']}.sorted.bam",
        annotations="results/prepare/unify_references/transcriptome.gtf",
        genome="results/prepare/unify_references/genome.fna",
    output:
        output_path=f"results/discover_isoforms/run_bambu/{config['sample_id']}/extended_annotations.gtf",
    threads: config["parallel_threads"]
    log:
        f"logs/discover_isoforms/run_bambu/{config['sample_id']}.log",
    conda:
        "../envs/r/bambu.yaml"
    script:
        "../scripts/r/run_bambu.R"


rule discover_isoforms_postprocess_transcriptome:
    input:
        genome="results/prepare/unify_references/genome.fna",
        transcriptome="results/prepare/unify_references/transcriptome.gtf",
    output:
        f"results/discover_isoforms/postprocess_transcriptome/transcriptome.gtf",
    log:
        f"logs/discover_isoforms/postprocess_transcriptome/{config['sample_id']}.log",
    conda:
        "../envs/standalone/gffread.yaml"
    shell:
        """
        conda list &> {log};
        gffread -g {input.genome} --adj-stop -J -T -o {output} \
            {input.transcriptome} &>> {log}
        """


rule discover_isoforms_postprocess_transcriptome_length:
    input:
        transcriptome=f"results/discover_isoforms/postprocess_transcriptome/transcriptome.gtf",
        unfiltered="results/prepare/unify_references/transcriptome.gtf",
    output:
        output_path=f"results/discover_isoforms/discover_isoforms_postprocess_transcriptome_length/transcriptome.gtf",
    params:
        minimum_orf_length=config["call_orfs_minimum_orf_length_nt"],
    log:
        f"logs/discover_isoforms/postprocess_transcriptome_length/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/postprocess_transcriptome_length.R"


rule discover_isoforms_postprocess_bambu:
    input:
        novel_transcriptome=f"results/discover_isoforms/run_bambu/{config['sample_id']}/extended_annotations.gtf",
        filtered_transcriptome=f"results/discover_isoforms/discover_isoforms_postprocess_transcriptome_length/transcriptome.gtf",
    output:
        output_path=f"results/discover_isoforms/postprocess_bambu/{config['sample_id']}.gtf",
    log:
        f"logs/discover_isoforms/postprocess_bambu/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/postprocess_bambu.R"
