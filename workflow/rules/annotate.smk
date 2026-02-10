

rule annotate_calculate_per_isoform_statistics:
    input:
        sqanti_protein=f"results/annotate/postprocess_sqanti_protein/{config['sample_id']}/{config['sample_id']}.protein_classification.tsv",
        haplosaurus_json=f"results/create_mass_spec_dbs/run_haplosaurus/{config['sample_id']}.txt",
    output:
        f"results/annotate/calculate_per_isoform_statistics/{config['sample_id']}.tsv",
    log:
        f"logs/annotate/calculate_per_isoform_statistics/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/calculate_per_isoform_statistics.R"


rule annotate_differentiate_per_isoform_statistics:
    input:
        input_path=f"results/annotate/calculate_per_isoform_statistics/{config['sample_id']}.tsv",
        sqanti_protein=f"results/annotate/postprocess_sqanti_protein/{config['sample_id']}/{config['sample_id']}.protein_classification.tsv",
    output:
        f"results/annotate/differentiate_per_isoform_statistics/{config['sample_id']}.tsv",
    log:
        f"logs/annotate/differentiate_per_isoform_statistics/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/differentiate_per_isoform_statistics.R"


rule annotate_extract_splice_cds:
    input:
        genome="results/prepare/unify_references/genome.fna",
        transcriptome=f"results/call_orfs/filter_orfanage/orfanage_cds_{config['sample_id']}.gtf",
    output:
        "results/annotate/extract_splice_cds/splice.fasta",
    log:
        f"logs/annotate/extract_splice_cds/{config['sample_id']}.log",
    conda:
        "../envs/standalone/gffread.yaml"
    shell:
        """
        conda list &> {log};
        gffread -x {output} -g {input.genome} \
            {input.transcriptome} &>> {log}
        """


rule annotate_annotate_peptides:
    input:
        input_peptides=f"results/annotate/extract_peptide_sequences/{config['sample_id']}_{{hap}}.fasta",
        input_cds_homo=f"results/annotate/extract_cds_mapping_sequences/{config['sample_id']}_homo.fasta",
        input_cds_first=f"results/annotate/extract_cds_mapping_sequences/{config['sample_id']}_first.fasta",
        input_cds_second=f"results/annotate/extract_cds_mapping_sequences/{config['sample_id']}_second.fasta",
        reference_transcriptome=f"results/call_orfs/filter_orfanage/orfanage_cds_{config['sample_id']}.gtf",
        reference_extracted_cds=f"results/create_mass_spec_dbs/extract_splice_protein_db/splice_cds.fasta",
        variants=f"results/phase_variants/compress_postprocessed_phasing/{config['sample_id']}.vcf.gz",        
    output:
        output_path=f"results/annotate/annotate_peptides/{config['sample_id']}_{{hap}}.tsv",
    log:
        f"logs/annotate/annotate_peptides/{config['sample_id']}_{{hap}}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/annotate_peptides.R"


rule annotate_run_protein_inference:
    input:
        "results/prepare/install_rcppgreedysetcover/done.txt",
        contaminants_path=config["contaminants_path"],
        sage_tsv=f"results/search_mass_spec/run_sage/{config['sample_id']}/results.sage.tsv",
    output:
        output_path=f"results/annotate/run_protein_inference/{config['sample_id']}/observed_protein_groups.tsv",
    params:
        fdr_cutoff=config["search_mass_spec_fdr"],
        seed=config["seed"],
    log:
        f"logs/annotate/run_protein_inference/{config['sample_id']}.log",
    conda:
        "../envs/r/protein_inference.yaml"
    script:
        "../scripts/r/run_protein_inference.R"


rule annotate_extract_cds_mapping_sequences:
    input:
        "results/prepare/install_transmogr/done.txt",
        sqanti_protein=f"results/annotate/postprocess_sqanti_protein/{config['sample_id']}/{config['sample_id']}.protein_classification.tsv",
        haplosaurus_json=f"results/create_mass_spec_dbs/run_haplosaurus/{config['sample_id']}.txt",
        splice_fasta="results/annotate/extract_splice_cds/splice.fasta",
        genome="results/prepare/unify_references/genome.fna",
        transcriptome=f"results/call_orfs/filter_orfanage/orfanage_cds_{config['sample_id']}.gtf",
        variants=f"results/phase_variants/compress_postprocessed_phasing/{config['sample_id']}.vcf.gz",
    output:
        output_path_homozygous=f"results/annotate/extract_cds_mapping_sequences/{config['sample_id']}_homo.fasta",
        output_path_first_haplotype=f"results/annotate/extract_cds_mapping_sequences/{config['sample_id']}_first.fasta",
        output_path_second_haplotype=f"results/annotate/extract_cds_mapping_sequences/{config['sample_id']}_second.fasta",
    log:
        f"logs/annotate/extract_cds_mapping_sequences/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/extract_cds_mapping_sequences.R"


rule annotate_extract_peptide_sequences:
    input:
        search_output=f"results/search_mass_spec/run_sage/{config['sample_id']}/results.sage.tsv",
        proteome_fasta=f"results/create_mass_spec_dbs/disambiguate_proteome/{config['sample_id']}.fasta",
        homozygous_isoforms_fasta=f"results/annotate/extract_cds_mapping_sequences/{config['sample_id']}_homo.fasta",
        first_haplotype_isoforms_fasta=f"results/annotate/extract_cds_mapping_sequences/{config['sample_id']}_first.fasta",
        second_haplotype_isoforms_fasta=f"results/annotate/extract_cds_mapping_sequences/{config['sample_id']}_second.fasta",
        contaminants_path=config["contaminants_path"],
    output:
        output_path_homozygous=f"results/annotate/extract_peptide_sequences/{config['sample_id']}_homo.fasta",
        output_path_first_haplotype=f"results/annotate/extract_peptide_sequences/{config['sample_id']}_first.fasta",
        output_path_second_haplotype=f"results/annotate/extract_peptide_sequences/{config['sample_id']}_second.fasta",
    params:
        fdr=config["search_mass_spec_fdr"],
    log:
        f"logs/annotate/extract_peptide_sequences/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/extract_peptide_sequences.R"


rule annotate_prepare_peptide_mapping:
    input:
        genome="results/prepare/unify_references/genome.fna",
        transcriptome=f"results/call_orfs/filter_orfanage/orfanage_cds_{config['sample_id']}.gtf",
    output:
        output_path=directory(
            f"results/annotate/prepare_peptide_mapping/{config['sample_id']}"
        ),
    params:
        sjdbOverhang=config["annotate_map_peptides_sjdbOverhang"],
    threads: config["parallel_threads"]
    log:
        f"logs/annotate/prepare_peptide_mapping/{config['sample_id']}.log",
    conda:
        "../envs/standalone/star.yaml"
    shell:
        """
        conda list &> {log};
        STAR --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.transcriptome} \
            --sjdbOverhang {params.sjdbOverhang} &>> {log}
        """


rule annotate_map_isoforms:
    input:
        reads=f"results/annotate/extract_cds_mapping_sequences/{config['sample_id']}_{{haplotype}}.fasta",
        genome="results/prepare/unify_references/genome.fna",
    output:
        f"results/annotate/map_isoforms/{config['sample_id']}/{{haplotype}}_isoforms.sorted.bam",
    params:
        align_map_bam_threads=config["parallel_threads"],
        align_sort_bam_threads=config["align_sort_bam_threads"],
        align_sort_bam_memory_gb=config["align_sort_bam_memory_gb"],
    threads: config["align_sort_bam_threads"] + config["parallel_threads"]
    log:
        f"logs/annotate/map_isoforms/{{haplotype}}/{config['sample_id']}.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    shell:
        """
        conda list &> {log};
        minimap2 -ax splice:hq -uf \
            -t {params.align_map_bam_threads}  \
            {input.genome} {input.reads} 2>> {log} | samtools sort \
            -@ {params.align_map_bam_threads} \
            -o {output} \
            -m{params.align_sort_bam_memory_gb}g \
            --write-index - &>> {log}
        """


rule annotate_map_peptides:
    input:
        reads=f"results/annotate/extract_peptide_sequences/{config['sample_id']}_{{haplotype}}.fasta",
        index_dir=f"results/annotate/prepare_peptide_mapping/{config['sample_id']}",
    output:
        output_path=f"results/annotate/map_peptides/{config['sample_id']}/{{haplotype}}_Aligned.sortedByCoord.out.bam",
    params:
        output_path=f"results/annotate/map_peptides/{config['sample_id']}/{{haplotype}}_",
    threads: config["parallel_threads"]
    log:
        f"logs/annotate/map_peptides/{config['sample_id']}/{{haplotype}}.log",
    conda:
        "../envs/standalone/star.yaml"
    shell:
        """
        conda list &> {log};
        STAR --runMode alignReads \
            --genomeDir {input.index_dir} \
            --runThreadN {threads} \
            --readFilesIn {input.reads} \
            --outFileNamePrefix {params.output_path} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outSAMattributes Standard &>> {log}
        """


rule annotate_index_bam_files:
    input:
        peptide=f"results/annotate/map_peptides/{config['sample_id']}/{{haplotype}}_Aligned.sortedByCoord.out.bam",
        isoform=f"results/annotate/map_isoforms/{config['sample_id']}/{{haplotype}}_isoforms.sorted.bam",
    output:
        peptide=f"results/annotate/map_peptides/{config['sample_id']}/{{haplotype}}_Aligned.sortedByCoord.out.bam.bai",
        isoform=f"results/annotate/map_isoforms/{config['sample_id']}/{{haplotype}}_isoforms.sorted.bam.bai",
    log:
        f"logs/annotate/index_bam_files/{config['sample_id']}/{{haplotype}}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        conda list &> {log};
        samtools index {input.peptide} &>> {log};
        samtools index {input.isoform} &>> {log}
        """


rule annotate_prepare_transcriptome_sqanti_protein:
    input:
        tsv=f"results/create_mass_spec_dbs/disambiguate_proteome/{config['sample_id']}.tsv",
        transcriptome=f"results/call_orfs/filter_orfanage/orfanage_cds_{config['sample_id']}.gtf",
    output:
        f"results/annotate/prepare_transcriptome_sqanti_protein/{config['sample_id']}.gtf",
    log:
        f"logs/annotate/prepare_transcriptome_sqanti_protein/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/adjust_transcriptome_ids.R"


rule annotate_prepare_sqanti_protein_gffread:
    input:
        f"results/annotate/prepare_transcriptome_sqanti_protein/{config['sample_id']}.gtf",
    output:
        f"results/annotate/prepare_sqanti_protein_gffread/{config['sample_id']}.gtf",
    log:
        f"logs/annotate/prepare_sqanti_protein_gffread/{config['sample_id']}.log",
    conda:
        "../envs/standalone/gffread.yaml"
    shell:
        """
        conda list &> {log};
        gffread -E -T {input} -o {output} &>> {log}
        """


rule annotate_prepare_sqanti_protein:
    input:
        query=f"results/annotate/prepare_sqanti_protein_gffread/{config['sample_id']}.gtf",
        reference=f"results/discover_isoforms/discover_isoforms_postprocess_transcriptome_length/transcriptome.gtf",
    output:
        query=f"results/annotate/prepare_sqanti_protein/{config['sample_id']}.gtf",
        reference=f"results/annotate/prepare_sqanti_protein/transcriptome.gtf",
    log:
        f"logs/annotate/prepare_sqanti_protein/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/prepare_sqanti_protein.R"


rule annotate_prepare_sqanti_protein_construct_gene_features:
    input:
        query=f"results/annotate/prepare_sqanti_protein/{config['sample_id']}.gtf",
        reference=f"results/annotate/prepare_sqanti_protein/transcriptome.gtf",
    output:
        query=f"results/annotate/prepare_sqanti_protein_construct_gene_features/{config['sample_id']}.gtf",
        reference=f"results/annotate/prepare_sqanti_protein_construct_gene_features/transcriptome.gtf",
    log:
        f"logs/annotate/prepare_sqanti_protein_construct_gene_features/{config['sample_id']}.log",
    conda:
        "../envs/standalone/gffread.yaml"
    shell:
        """
        conda list &> {log};
        gffread -E --keep-genes {input.query} -o- 2>> {log} | \
            gffread -E - -T -o {output.query} &>> {log};
        gffread -E --keep-genes {input.reference} -o- 2>> {log} | \
            gffread -E - -T -o {output.reference} &>> {log}
        """


# Adapted from: https://github.com/sheynkman-lab/Long-Read-Proteogenomics/blob/460d53ae716419438b31a3a8052f3fce2e6175ce/main.nf
rule annotate_run_sqanti_protein:
    input:
        query=f"results/annotate/prepare_transcriptome_sqanti_protein/{config['sample_id']}.gtf",
        reference=f"results/discover_isoforms/discover_isoforms_postprocess_transcriptome_length/transcriptome.gtf",
        query_cds=f"results/annotate/prepare_sqanti_protein_construct_gene_features/{config['sample_id']}.gtf",
        reference_cds="results/annotate/prepare_sqanti_protein_construct_gene_features/transcriptome.gtf",
    output:
        f"results/annotate/run_sqanti_protein/{config['sample_id']}/{config['sample_id']}.sqanti_protein_classification.tsv",
        outdir=directory(f"results/annotate/run_sqanti_protein/{config['sample_id']}"),
    params:
        output_folder=f"results/annotate/run_sqanti_protein/{config['sample_id']}",
        sample_id=config["sample_id"],
    log:
        f"logs/annotate/run_sqanti_protein/{config['sample_id']}.log",
    conda:
        "../envs/python/sqanti.yaml"
    shell:
        """
        conda list &> {log};
        python workflow/scripts/py/sqanti3_protein.py \
                {input.query} \
                {input.query_cds} \
                {input.reference} \
                {input.reference_cds} \
                -d {output.outdir} \
                -p {params.sample_id} &>> {log}
        """


rule annotate_postprocess_sqanti_protein:
    input:
        protein_classification=f"results/annotate/run_sqanti_protein/{config['sample_id']}/{config['sample_id']}.sqanti_protein_classification.tsv",
        query=f"results/annotate/prepare_sqanti_protein_gffread/{config['sample_id']}.gtf",
        reference=f"results/discover_isoforms/discover_isoforms_postprocess_transcriptome_length/transcriptome.gtf",
    output:
        f"results/annotate/postprocess_sqanti_protein/{config['sample_id']}/{config['sample_id']}.protein_classification.tsv",
        outdir=directory(
            f"results/annotate/postprocess_sqanti_protein/{config['sample_id']}"
        ),
    params:
        sample_id=config["sample_id"],
    log:
        f"logs/annotate/postprocess_sqanti_protein/{config['sample_id']}.log",
    conda:
        "../envs/python/sqanti_postprocess.yaml"
    shell:
        """
        conda list &> {log};
        python workflow/scripts/py/get_gc_exon_and_5utr_info.py \
                --gencode_gtf {input.reference} \
                --odir {output.outdir} &>> {log};
        python workflow/scripts/py/classify_5utr_status.py \
                --gencode_exons_bed {output.outdir}/gencode_exons_for_cds_containing_ensts.bed \
                --gencode_exons_chain {output.outdir}/gc_exon_chain_strings_for_cds_containing_transcripts.tsv \
                --sample_cds_gtf {input.query} \
                --odir {output.outdir} &>> {log};
        python workflow/scripts/py/merge_5utr_info_to_pclass_table.py \
            --name {params.sample_id} \
            --utr_info {output.outdir}/pb_5utr_categories.tsv \
            --sqanti_protein_classification {input.protein_classification} \
            --odir {output.outdir} &>> {log};
        python workflow/scripts/py/protein_classification.py \
            --sqanti_protein {output.outdir}/{params.sample_id}.sqanti_protein_classification_w_5utr_info.tsv \
            --name {params.sample_id} \
            --dest_dir {output.outdir} &>> {log}
        """


rule annotate_run_vafator:
    input:
        vcf=f"results/phase_variants/compress_postprocessed_phasing/{config['sample_id']}.vcf.gz",
        homo_mapped_peptides=f"results/annotate/map_peptides/{config['sample_id']}/homo_Aligned.sortedByCoord.out.bam",
        first_mapped_peptides=f"results/annotate/map_peptides/{config['sample_id']}/first_Aligned.sortedByCoord.out.bam",
        second_mapped_peptides=f"results/annotate/map_peptides/{config['sample_id']}/second_Aligned.sortedByCoord.out.bam",
    output:
        f"results/annotate/run_vafator/{config['sample_id']}/variants_annoated.vcf",
    log:
        f"logs/annotate/run_vafator/{config['sample_id']}.log",
    conda:
        "../envs/standalone/vafator.yaml"
    shell:
        """
        conda list &> {log};
        vafator --input-vcf {input.vcf} \
            --output-vcf {output} \
            --bam homo {input.homo_mapped_peptides} \
            --bam first {input.first_mapped_peptides} \
            --bam second {input.second_mapped_peptides} \
            &>> {log}
        """


rule annotate_get_directly_observed_and_linked_variants:
    input:
        variants=f"results/annotate/run_vafator/{config['sample_id']}/variants_annoated.vcf",
        inferred_proteins=f"results/annotate/run_protein_inference/{config['sample_id']}/observed_protein_groups.tsv",
        variant_annotation=f"results/annotate/differentiate_per_isoform_statistics/{config['sample_id']}.tsv",
    output:
        linked_variants_path=f"results/annotate/get_directly_observed_and_linked_variants/{config['sample_id']}/linked_variants.vcf",
        observed_variants_path=f"results/annotate/get_directly_observed_and_linked_variants/{config['sample_id']}/observed_variants.vcf",
    log:
        f"logs/annotate/get_directly_observed_and_linked_variants/{config['sample_id']}.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/get_directly_observed_and_linked_variants.R"
