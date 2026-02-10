

rule call_variants_run_clair3_rna:
    input:
        bam=f"results/align/run_minimap2_genome/{config['sample_id']}.sorted.bam",
        genome="results/prepare/unify_references/genome.fna",
    output:
        f"results/call_variants/run_clair3_rna/{config['sample_id']}/{config['sample_id']}.vcf.gz",
        results_dir=directory(
            f"results/call_variants/run_clair3_rna/{config['sample_id']}"
        ),
    params:
        sample_id=config["sample_id"],
        platform=config["call_variants_clair3_rna_mapping"][config["technology"]],
        qual=config["call_variants_min_qual"],
        snp_min_af=config["call_variants_snp_min_af"],
        indel_min_af=config["call_variants_indel_min_af"],
        min_coverage=config["call_variants_min_coverage"],
        min_mq=config["call_variants_min_mq"],
    threads: config["parallel_threads"]
    log:
        f"logs/call_variants/run_clair3_rna/{config['sample_id']}.log",
    container:
        f"docker://hkubal/clair3-rna:{config["clair3_rna_version"]}"
    shell:
        """
        /opt/bin/run_clair3_rna \
            --bam_fn {input.bam} \
            --ref_fn {input.genome} \
            --threads {threads} \
            --platform {params.platform} \
            --output_dir {output.results_dir} \
            --tag_variant_using_readiportal \
            --output_prefix {params.sample_id} \
            --include_all_ctgs \
            --platform {params.platform} \
            --qual {params.qual} \
            --snp_min_af {params.snp_min_af} \
            --indel_min_af {params.indel_min_af} \
            --min_coverage {params.min_coverage} \
            --min_mq {params.min_mq} \
            --conda_prefix /opt/conda/envs/clair3_rna \
            &>> {log}
        """


rule call_variants_filter_variants:
    input:
        variants=f"results/call_variants/run_clair3_rna/{config['sample_id']}/{config['sample_id']}.vcf.gz",
        transcriptome=f"results/call_orfs/filter_orfanage/orfanage_cds_{config['sample_id']}.gtf",
    output:
        output_path=f"results/call_variants/filter_variants/{config['sample_id']}.vcf",
    params:
        keep_indels=config["call_variants_keep_indels"],
        min_gq=config["call_variants_min_gq"],
    log:
        f"logs/call_variants/filter_variants/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/filter_variants.R"


rule call_variants_compress_filtered_variants:
    input:
        f"results/call_variants/filter_variants/{config['sample_id']}.vcf",
    output:
        f"results/call_variants/compress_filtered_variants/{config['sample_id']}.vcf.gz",
    params:
        keep_indels=config["call_variants_keep_indels"],
    log:
        f"logs/call_variants/compress_filtered_variants/{config['sample_id']}.log",
    conda:
        "../envs/standalone/bcftools.yaml"
    shell:
        """
        conda list &> {log};
        if [ {params.keep_indels} = "true" ]; then
            bcftools view -m2 -M2 -v snps,indels {input} | 
                bcftools sort - -o {output} -O z &>> {log};
        else
            echo "not keeping indels";
            bcftools view -m2 -M2 -v snps {input} | 
                bcftools sort - -o {output} -O z &>> {log};
        fi
        tabix {output} &>> {log}
        """
