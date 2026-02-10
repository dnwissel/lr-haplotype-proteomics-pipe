

rule phase_variants_run_whatshap:
    input:
        variants=f"results/call_variants/compress_filtered_variants/{config['sample_id']}.vcf.gz",
        aligned_reads=f"results/align/run_minimap2_genome/{config['sample_id']}.sorted.bam",
        genome="results/prepare/unify_references/genome.fna",
    output:
        f"results/phase_variants/run_whatshap/{config['sample_id']}.vcf.gz",
    log:
        f"logs/phase_variants/run_whatshap/{config['sample_id']}.log",
    conda:
        "../envs/standalone/whatshap.yaml"
    shell:
        """
        conda list &> {log};
        whatshap phase --ignore-read-groups -o {output} \
            --reference={input.genome} {input.variants} {input.aligned_reads} \
            --distrust-genotypes \
            &>> {log}
        """


rule phase_variants_postprocess_phasing:
    input:
        variants=f"results/phase_variants/run_whatshap/{config['sample_id']}.vcf.gz",
        query_transcriptome=f"results/call_orfs/filter_orfanage/orfanage_cds_{config['sample_id']}.gtf",
    output:
        resolved_vcf=f"results/phase_variants/postprocess_phasing/{config['sample_id']}.vcf",
    log:
        f"logs/phase_variants/postprocess_phasing/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/postprocess_phasing.R"


rule phase_variants_compress_postprocessed_phasing:
    input:
        f"results/phase_variants/postprocess_phasing/{config['sample_id']}.vcf",
    output:
        f"results/phase_variants/compress_postprocessed_phasing/{config['sample_id']}.vcf.gz",
    log:
        f"logs/phase_variants/compress_postprocessed_phasing/{config['sample_id']}.log",
    conda:
        "../envs/standalone/bcftools.yaml"
    shell:
        """
        conda list &> {log};
        bcftools sort {input} -o {output} -O z &>> {log};
        tabix {output} &>> {log}
        """
