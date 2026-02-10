

rule call_orfs_run_orfanage:
    input:
        query_transcriptome=f"results/discover_isoforms/postprocess_bambu/{config['sample_id']}.gtf",
        genome="results/prepare/unify_references/genome.fna",
        reference_transcriptome=f"results/discover_isoforms/discover_isoforms_postprocess_transcriptome_length/transcriptome.gtf",
    output:
        f"results/call_orfs/run_orfanage/{config['sample_id']}.gtf",
    params:
        mode=config["call_orfs_orfanage_mode"],
    threads: config["parallel_threads"]
    log:
        f"logs/call_orfs/run_orfanage/{config['sample_id']}.log",
    conda:
        "../envs/standalone/orfanage.yaml"
    shell:
        """
        conda list &> {log};
        orfanage \
            --mode {params.mode} \
            --reference {input.genome} \
            --query {input.query_transcriptome} \
            --output {output} \
            --threads {threads} \
            --keep_cds \
            {input.reference_transcriptome} \
            &>> {log}
        """


rule call_orfs_filter_orfanage:
    input:
        orfanage_gtf_file_path=f"results/call_orfs/run_orfanage/{config['sample_id']}.gtf",
    output:
        output_path_to_be_predicted=f"results/call_orfs/filter_orfanage/orfanage_no_cds_{config['sample_id']}.gtf",
        output_path_filtered=f"results/call_orfs/filter_orfanage/orfanage_cds_{config['sample_id']}.gtf",
    params:
        minimum_orf_length=config["call_orfs_minimum_orf_length_nt"],
    log:
        f"logs/call_orfs/filter_orfanage/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/filter_orfanage.R"
