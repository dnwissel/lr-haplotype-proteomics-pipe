from snakemake.utils import min_version


configfile: "config/pipeline_config.yaml"
configfile: "config/required_config.yaml"
configfile: "config/optional_config.yaml"


min_version(config["snakemake_min_version"])


container: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"


include: "workflow/rules/prepare.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/discover_isoforms.smk"
include: "workflow/rules/call_orfs.smk"
include: "workflow/rules/call_variants.smk"
include: "workflow/rules/phase_variants.smk"
include: "workflow/rules/create_mass_spec_db.smk"
include: "workflow/rules/search_mass_spec.smk"
include: "workflow/rules/annotate.smk"


rule all:
    input:
        f"results/annotate/calculate_per_isoform_statistics/{config['sample_id']}.tsv",
        f"results/annotate/run_protein_inference/{config['sample_id']}/observed_protein_groups.tsv",
        f"results/annotate/differentiate_per_isoform_statistics/{config['sample_id']}.tsv",
        f"results/annotate/annotate_peptides/{config['sample_id']}_homo.tsv",
        f"results/annotate/annotate_peptides/{config['sample_id']}_first.tsv",
        f"results/annotate/annotate_peptides/{config['sample_id']}_second.tsv",
        peptide=expand(
            f"results/annotate/map_peptides/{config['sample_id']}/{{haplotype}}_Aligned.sortedByCoord.out.bam.bai",
            haplotype=["first", "second", "homo"],
        ),
        isoform=expand(
            f"results/annotate/map_isoforms/{config['sample_id']}/{{haplotype}}_isoforms.sorted.bam.bai",
            haplotype=["first", "second", "homo"],
        ),
        linked_variants_path=f"results/annotate/get_directly_observed_and_linked_variants/{config['sample_id']}/linked_variants.vcf",
        observed_variants_path=f"results/annotate/get_directly_observed_and_linked_variants/{config['sample_id']}/observed_variants.vcf",
