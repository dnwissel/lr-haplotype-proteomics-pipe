

rule search_mass_spec_run_sage:
    input:
        spectra=config["ms_data_path"],
        config=config["search_mass_spec_sage_config_path"],
        proteome=f"results/create_mass_spec_dbs/concatenate_contaminants/{config['sample_id']}.fasta",
    output:
        f"results/search_mass_spec/run_sage/{config['sample_id']}/results.sage.tsv",
    params:
        outdir=f"results/search_mass_spec/run_sage/{config['sample_id']}",
    threads: config["parallel_threads"]
    log:
        f"logs/search_mass_spec/run_sage/{config['sample_id']}.log",
    conda:
        "../envs/standalone/sage.yaml"
    shell:
        """
        conda list &> {log};
        export SAGE_LOG=trace;
        export RAYON_NUM_THREADS={threads};
        sage {input.config} -f {input.proteome} \
            --write-pin --output_directory {params.outdir} {input.spectra} \
            &>> {log}
        """
