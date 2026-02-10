

rule create_mass_spec_dbs_convert_gtf_to_gff:
    input:
        f"results/call_orfs/filter_orfanage/orfanage_cds_{config['sample_id']}.gtf",
    output:
        standardized=f"results/create_mass_spec_dbs/convert_gtf_to_gff/{config['sample_id']}_standardized.gff",
        _sorted=f"results/create_mass_spec_dbs/convert_gtf_to_gff/{config['sample_id']}_sorted.gff",
    log:
        f"logs/create_mass_spec_dbs/convert_gtf_to_gff/{config['sample_id']}.log",
    conda:
        "../envs/standalone/gtfsort.yaml"
    shell:
        """
        conda list &> {log};
        gffread --keep-genes -E {input} -o {output.standardized} &>> {log};
        gff3sort.pl --precise --chr_order natural {output.standardized} > {output._sorted} 2>> {log}
        """


rule create_mass_spec_dbs_add_transcript_biotype:
    input:
        f"results/create_mass_spec_dbs/convert_gtf_to_gff/{config['sample_id']}_sorted.gff",
    output:
        f"results/create_mass_spec_dbs/add_transcript_biotype/{config['sample_id']}.gff",
    params:
        file_name=f"results/create_mass_spec_dbs/add_transcript_biotype/{config['sample_id']}",
    log:
        f"logs/create_mass_spec_dbs/add_transcript_biotype/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/add_transcript_biotype.R"


rule create_mass_spec_dbs_zip_gff:
    input:
        f"results/create_mass_spec_dbs/add_transcript_biotype/{config['sample_id']}.gff",
    output:
        f"results/create_mass_spec_dbs/zip_gff/{config['sample_id']}.gff.gz",
    log:
        f"logs/create_mass_spec_dbs/zip_gff/{config['sample_id']}.log",
    conda:
        "../envs/standalone/gtfsort.yaml"
    shell:
        """
        conda list &> {log};
        bgzip -c {input} > {output} 2>> {log};
        tabix {output} &>> {log}
        """


rule create_mass_spec_dbs_run_haplosaurus:
    input:
        variants=f"results/phase_variants/compress_postprocessed_phasing/{config['sample_id']}.vcf.gz",
        transcriptome=f"results/create_mass_spec_dbs/zip_gff/{config['sample_id']}.gff.gz",
        genome="results/prepare/unify_references/genome.fna",
    output:
        f"results/create_mass_spec_dbs/run_haplosaurus/{config['sample_id']}.txt",
    log:
        f"logs/create_mass_spec_dbs/create_protein_db_gencode/{config['sample_id']}.log",
    threads: config["parallel_threads"]
    container:
        "docker://ensemblorg/ensembl-vep:release_114.1"
    #conda:
    #    "../envs/standalone/vep.yaml"
    shell:
        """
        touch {log};
        haplo -i  {input.variants} --gff {input.transcriptome} \
            --fasta {input.genome} -o {output} --fork {threads} \
            --force_overwrite --json &>> {log}
        """


rule create_mass_spec_dbs_extract_splice_protein_db:
    input:
        genome="results/prepare/unify_references/genome.fna",
        transcriptome=f"results/call_orfs/filter_orfanage/orfanage_cds_{config['sample_id']}.gtf",
    output:
        protein="results/create_mass_spec_dbs/extract_splice_protein_db/splice.fasta",
        cds="results/create_mass_spec_dbs/extract_splice_protein_db/splice_cds.fasta",
    log:
        f"logs/create_mass_spec_dbs/create_protein_db_gencode/{config['sample_id']}.log",
    conda:
        "../envs/standalone/gffread.yaml"
    shell:
        """
        conda list &> {log};
        gffread -y {output.protein} -g {input.genome} \
            {input.transcriptome} &>> {log};
        gffread -x {output.cds} -g {input.genome} \
            {input.transcriptome} &>> {log}
        """


rule create_mass_spec_dbs_create_protein_db:
    input:
        splice_proteome="results/create_mass_spec_dbs/extract_splice_protein_db/splice.fasta",
        haplosaurus_json=f"results/create_mass_spec_dbs/run_haplosaurus/{config['sample_id']}.txt",
    output:
        f"results/create_mass_spec_dbs/create_protein_db/{config['sample_id']}.fasta",
    log:
        f"logs/create_mass_spec_dbs/create_protein_db/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/create_protein_db.R"


rule create_mass_spec_db_disambiguate_proteome:
    input:
        input_path_fasta=f"results/create_mass_spec_dbs/create_protein_db/{config['sample_id']}.fasta",
    output:
        output_path_tsv=f"results/create_mass_spec_dbs/disambiguate_proteome/{config['sample_id']}.tsv",
        output_path_fasta=f"results/create_mass_spec_dbs/disambiguate_proteome/{config['sample_id']}.fasta",
    log:
        f"logs/create_mass_spec_dbs/disambiguate_proteome/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/disambiguate_proteome.R"


rule create_mass_spec_db_create_decoys:
    input:
        f"results/create_mass_spec_dbs/disambiguate_proteome/{config['sample_id']}.fasta",
    output:
        f"results/create_mass_spec_dbs/create_decoys/{config['sample_id']}.fasta",
    params:
        enzyme=config["search_mass_spec_ms_enzyme"],
        cleavage_position=config["search_mass_spec_ms_cleavage_position"],
        max_missed_cleavages=config["search_mass_spec_ms_max_missed_cleavages"],
        min_peptide_length=config["search_mass_spec_ms_decoy_min_peptide_length"],
        max_peptide_length=config["search_mass_spec_ms_decoy_max_peptide_length"],
        max_iterations=config["search_mass_spec_ms_decoy_max_iterations"],
    threads: config["parallel_threads"]
    log:
        f"logs/create_mass_spec_dbs/create_decoys_gencode/{config['sample_id']}.log",
    conda:
        "../envs/standalone/pypgatk.yaml"
    shell:
        """
        conda list &> {log};
        pypgatk generate-decoy --input_database {input} \
            --output_database {output} --method decoypyrat \
            --decoy_prefix decoy_ \
            --enzyme {params.enzyme} \
            --cleavage_position {params.cleavage_position} \
            --max_missed_cleavages {params.max_missed_cleavages} \
            --min_peptide_length {params.min_peptide_length} \
            --max_peptide_length {params.max_peptide_length} \
            --max_iterations {params.max_iterations} &>> {log}
        """


rule create_mass_spec_db_concatenate_contaminants:
    input:
        protein_database=f"results/create_mass_spec_dbs/create_decoys/{config['sample_id']}.fasta",
        contaminants=config["contaminants_path"],
    output:
        f"results/create_mass_spec_dbs/concatenate_contaminants/{config['sample_id']}.fasta",
    log:
        f"logs/create_mass_spec_dbs/concatenate_contaminants_gencode/{config['sample_id']}.log",
    shell:
        """
        cat {input.protein_database} {input.contaminants} > {output} 2> {log}
        """
