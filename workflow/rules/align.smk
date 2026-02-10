

rule align_convert_gtf_to_bed:
    input:
        "results/prepare/unify_references/transcriptome.gtf",
    output:
        "results/align/convert_gtf_to_bed/transcriptome.bed",
    log:
        "logs/align/convert_gtf_to_bed.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    shell:
        """
        conda list &> {log};
        paftools.js gff2bed {input} > {output} 2>> {log}
        """


rule align_run_minimap2_genome:
    input:
        reads=config["rna_data_path"],
        transcriptome="results/align/convert_gtf_to_bed/transcriptome.bed",
        genome="results/prepare/unify_references/genome.fna",
    output:
        f"results/align/run_minimap2_genome/{config['sample_id']}.sorted.bam",
    params:
        align_map_bam_threads=config["parallel_threads"],
        align_sort_bam_threads=config["align_sort_bam_threads"],
        align_sort_bam_memory_gb=config["align_sort_bam_memory_gb"],
        sequencing_tech=config["technology"],
    threads: config["align_sort_bam_threads"] + config["parallel_threads"]
    log:
        f"logs/align/run_minimap2_genome/{config['sample_id']}.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    shell:
        """
        conda list &> {log};
        if [ {params.sequencing_tech} = "pb_iso_seq" ] ||  [ {params.sequencing_tech} = "pb_mas_seq" ]; then
            minimap2 -ax splice:hq -uf --junc-bed {input.transcriptome} \
                -t {params.align_map_bam_threads}  \
                {input.genome} {input.reads} 2>> {log} | samtools sort \
                -@ {params.align_map_bam_threads} \
                -m{params.align_sort_bam_memory_gb}g \
                -o {output} --write-index - &>> {log}
        elif [ {params.sequencing_tech} = "ont_r9_cdna" ]; then
            minimap2 -ax splice --junc-bed {input.transcriptome} \
                -t {params.align_map_bam_threads}  \
                {input.genome} {input.reads} 2>> {log} | samtools sort \
                -@ {params.align_map_bam_threads} \
                -m{params.align_sort_bam_memory_gb}g \
                -o {output} --write-index - &>> {log}
        elif [ {params.sequencing_tech} = "ont_r9_drna" ]; then
            minimap2 -ax splice -uf -k14 --junc-bed {input.transcriptome} \
                -t {params.align_map_bam_threads}  \
                {input.genome} {input.reads} 2>> {log} | samtools sort \
                -@ {params.align_map_bam_threads} \
                -m{params.align_sort_bam_memory_gb}g \
                -o {output} --write-index - &>> {log}
        fi
        """
