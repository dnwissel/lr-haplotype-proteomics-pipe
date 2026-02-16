# Long-read RNA-seq haplotype-resolved proteogenomics Snakemake pipeline

## Overview

This Snakemake workflow performs sample-specific, haplotype-resolved protein isoform characterization by integrating lrRNA-seq with MS data. It was originally developed as part of our paper [1].

## What the workflow does

The pipeline, broadly, performs the following steps:

1. **Alignment**: Aligns long-read RNA-seq data to the reference genome using minimap2 [2].
2. **Transcript discovery**: Discovers novel transcript isoforms using Bambu [3].
3. **ORF calling**: Predicts open reading frames (ORFs) from discovered transcripts using orfanage [4].
4. **Variant calling**: Identifies genetic variants from long-read data using Clair3-RNA [5].
5. **Phasing**: Phases variants using WhatsHap [6].
6. **Proteome construction**: Builds haplotype-resolved, sample-specific protein databases using haplosaurus [7], and DecoyPYrat [8].
7. **MS search**: Searches mass spectrometry data against the custom proteome using sage [9].
8. **Annotation**: Annotates peptides, performs protein inference, and identifies variant-linked peptides using vafator.

## Requirements

### Dependencies

- **snakemake** >= 8.2.0
- **conda**
- **singularity** or **apptainer**

### Input data

- **Long-read RNA-seq data**: lrRNA-seq FASTQ files (PB Iso-Seq or Kinnex, ONT cDNA or dRNA) (please note that lrRNA-seq data should be ready for alignment, i.e., for PB data you may want to pass FLNC reads and for ONT cDNA data, you may run a tool such as pychopper before using our pipeline).
- **Mass spectrometry data**: mzML files.
- **Reference genome**: FASTA file.
- **Reference transcriptome**: GTF file.
- **Contaminants database**: FASTA file (e.g., cRAP database [10]).

## Quick start

### 1. Clone the repo

```bash
git clone https://github.com/dnwissel/lr-haplotype-proteomics-pipe
cd lr-haplotype-proteomics-pipe
```

### 2. Configure the workflow

The pipeline uses two configuration files:

- **`config/required_config.yaml`**: Required parameters that must be set.
- **`config/optional_config.yaml`**: Optional parameters with sensible defaults (can typically be left unchanged).
- **`config/pipeline_config.yaml`**: Internal pipeline settings (can typically be let unchanged).
- **`config/sage_config.json`**: Sage search configuration that must be set.

### 3. Run the workflow

```bash
snakemake --use-conda --use-apptainer --cores <number_of_cores>
```

## Configuration guide

### Required configuration

These parameters **must** be configured before running the pipeline:

| Parameter                                  | Description                                                                                                                                                                                                            | Example                                          |
| :----------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :----------------------------------------------- |
| `sample_id`                                | Unique identifier for your sample                                                                                                                                                                                      | `"SRR187"`                                       |
| `technology`                               | Sequencing technology: `"pb_iso_seq"`, `"pb_kinnex"`, `"ont_r9_cdna"`, `"ont_002_drna"`, `"ont_r10_cdna"`, or `"ont_004_drna"`.                                                                                        | `"pb_iso_seq"`                                   |
| `rna_data_path`                            | Path to long-read RNA-seq data (BAM or FASTQ file)                                                                                                                                                                     | `"/data/rna/sample.fastq.gz"`                    |
| `ms_data_path`                             | List of paths to mzML files (one or more)                                                                                                                                                                              | `["/data/ms/file1.mzML", "/data/ms/file2.mzML"]` |
| `genome_path`                              | Full path to reference genome FASTA file                                                                                                                                                                               | `"/data/refs/genome.fna"`                        |
| `transcriptome_path`                       | Full path to reference transcriptome GTF file                                                                                                                                                                          | `"/data/refs/transcriptome.gtf"`                 |
| `contaminants_path`                        | Full path to contaminants FASTA file (e.g., cRAP database)                                                                                                                                                             | `"/data/refs/crap.fasta"`                        |
| `search_mass_spec_sage_config_path`        | Path to sage configuration JSON file. This file **must** be configured before running. Please refer to the [sage documentation](https://sage-docs.vercel.app/docs/configuration#file) for guidance on the sage config. | `"config/sage_config.json"`                      |
| `search_mass_spec_ms_enzyme`               | Protease enzyme name (e.g., `"Trypsin"`). **Must match** the enzyme setting in your SAGE config file.                                                                                                                  | `"Trypsin"`                                      |
| `search_mass_spec_ms_cleavage_position`    | Cleavage position (`"c"` for C-terminal). **Must match** the cleavage position in your SAGE config file.                                                                                                               | `"c"`                                            |
| `search_mass_spec_ms_max_missed_cleavages` | Maximum allowed missed cleavages. **Must match** the max missed cleavages setting in your SAGE config file.                                                                                                            | `2`                                              |

### Optional configuration

These parameters have sensible defaults and can typically be left unchanged:

### Resource configuration

| Parameter                  | Description                               | Default |
| -------------------------- | ----------------------------------------- | ------- |
| `parallel_threads`         | Number of threads for parallel processing | `12`    |
| `align_sort_bam_memory_gb` | Memory (GB) for BAM sorting               | `4`     |
| `align_sort_bam_threads`   | Threads for BAM sorting                   | `4`     |

### Software version configuration

These parameters are defined in `config/pipeline_config.yaml` and control which tool versions the pipeline expects.

**Note:** The Bambu version is pinned in the `workflow/envs/r/bambu.yaml` conda environment file, so it must be changed directly there.

**Note:** `clair3_rna_version` is used directly as the tag for the Clair3-RNA Docker image (`docker://hkubal/clair3-rna`), so it must match the Docker container version you want to use.

| Parameter                    | Description                            | Default    |
| ---------------------------- | -------------------------------------- | ---------- |
| `transmogr_version`          | Transmogr version                      | `"1.6.0"`  |
| `rcppgreedysetcover_version` | RcppGreedySetCover version             | `"0.1.0"`  |
| `clair3_rna_version`         | Clair3-RNA version for variant calling | `"v0.2.2"` |

### ORF calling configuration

| Parameter                         | Description                               | Default  |
| --------------------------------- | ----------------------------------------- | -------- |
| `call_orfs_orfanage_mode`         | ORFanage prediction mode (e.g., `"BEST"`) | `"BEST"` |
| `call_orfs_minimum_orf_length_nt` | Minimum ORF length in nucleotides         | `100`    |

### Variant calling configuration

| Parameter                    | Description                                    | Default  |
| ---------------------------- | ---------------------------------------------- | -------- |
| `call_variants_min_qual`     | Minimum variant quality score                  | `0`      |
| `call_variants_snp_min_af`   | Minimum allele frequency for SNPs              | `0.08`   |
| `call_variants_indel_min_af` | Minimum allele frequency for indels            | `0.15`   |
| `call_variants_min_coverage` | Minimum read coverage                          | `3`      |
| `call_variants_min_mq`       | Minimum mapping quality                        | `5`      |
| `call_variants_keep_indels`  | Whether to keep indels (`"true"` or `"false"`) | `"true"` |
| `call_variants_min_gq`       | Minimum genotype quality                       | `10`     |

### MS search configuration

| Parameter                                      | Description                             | Default |
| ---------------------------------------------- | --------------------------------------- | ------- |
| `search_mass_spec_ms_decoy_min_peptide_length` | Minimum decoy peptide length            | `5`     |
| `search_mass_spec_ms_decoy_max_peptide_length` | Maximum decoy peptide length            | `100`   |
| `search_mass_spec_ms_decoy_max_iterations`     | Maximum iterations for decoy generation | `100`   |
| `search_mass_spec_fdr`                         | False discovery rate threshold          | `0.01`  |

#### Annotation configuration

| Parameter                            | Description                                   | Default |
| ------------------------------------ | --------------------------------------------- | ------- |
| `annotate_map_peptides_sjdbOverhang` | STAR splice junction database overhang length | `100`   |

## Demo

We first download the repo and get appropriate reference files:

```bash
git clone https://github.com/dnwissel/lr-haplotype-proteomics-pipe
cd lr-haplotype-proteomics-pipe

cd references
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz
wget ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip gencode.v49.primary_assembly.annotation.gtf.gz

cd ../data
wget --content-disposition "https://zenodo.org/records/18588655/files/SRR187_subsampled.fastq.gz?download=1"
wget --content-disposition "https://zenodo.org/records/18588655/files/210728_WTC11Blue1_Tryp4,6_HCDorbi_4hr_B4_20210801032622.mzML?download=1"
```

Then, we configure `config/required_config.yaml` and `config/sage_config.yaml` appropriately for our data:

```yaml
sample_id: "SRR187"
technology: "pb_iso_seq"

rna_data_path: "data/SRR187_subsampled.fastq.gz"
ms_data_path:
  ["data/210728_WTC11Blue1_Tryp4,6_HCDorbi_4hr_B4_20210801032622.mzML"]

genome_path: "references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
transcriptome_path: "references/gencode.v49.primary_assembly.annotation.gtf"
contaminants_path: "references/crap.fasta"

search_mass_spec_sage_config_path: "config/sage_config.json"
search_mass_spec_ms_enzyme: "Trypsin"
search_mass_spec_ms_cleavage_position: "c"
search_mass_spec_ms_max_missed_cleavages: 2
```

```yaml
{
  "database":
    {
      "bucket_size": 32768,
      "enzyme":
        {
          "missed_cleavages": 2,
          "min_len": 5,
          "max_len": 100,
          "cleave_at": "KR",
          "restrict": "P",
          "c_terminal": true,
        },
      "fragment_min_mz": 150.0,
      "fragment_max_mz": 2000.0,
      "peptide_min_mass": 500.0,
      "peptide_max_mass": 5000.0,
      "static_mods": { "C": 57.0215 },
      "variable_mods": { "M": [15.9949] },
      "max_variable_mods": 2,
      "decoy_tag": "decoy_",
      "generate_decoys": false,
    },
  "precursor_tol": { "ppm": [-10, 10] },
  "fragment_tol": { "ppm": [-10, 10] },
  "precursor_charge": [2, 4],
  "predict_rt": true,
  "min_peaks": 15,
  "max_peaks": 500,
  "min_matched_peaks": 4,
  "report_psms": 1,
}
```

Lastly, since we are only searching a small amount of MS data, let's set the FDR in `config/optional_config.yaml` to a larger number to ensure we get some outputs:

```yaml
parallel_threads: 12
align_sort_bam_memory_gb: 4
align_sort_bam_threads: 4

call_orfs_orfanage_mode: "BEST"
call_orfs_minimum_orf_length_nt: 100

call_variants_min_qual: 0
call_variants_snp_min_af: 0.08
call_variants_indel_min_af: 0.15
call_variants_min_coverage: 3
call_variants_min_mq: 5
call_variants_keep_indels: "true"
call_variants_min_gq: 10

search_mass_spec_ms_decoy_min_peptide_length: 5
search_mass_spec_ms_decoy_max_peptide_length: 100
search_mass_spec_ms_decoy_max_iterations: 100
# MAKE CHANGE HERE
search_mass_spec_fdr: 0.2

annotate_map_peptides_sjdbOverhang: 100
```

Then we run:

```bash
snakemake --use-conda --use-apptainer --cores 12
```

## Output files

The pipeline generates results in the `results/` directory. Main output files are described below:

### Annotation results (primary outputs)

#### `results/annotate/run_protein_inference/{sample_id}/observed_protein_groups.tsv`

Protein groups identified from MS data using a greedy set cover algorithm based on peptide evidence from the search.

```bash
protein_group
ENST00000244217.6
ENST00000377276.5
ENST00000308683.3
ENST00000346932.9A-ENST00000346932.9B-ENST00000350527.7A-ENST00000350527.7B-ENST00000537485.5A-ENST00000537485.5B-ENST00000617533.5A-ENST00000617533.5B
ENST00000308724.9
ENST00000372652.5
ENST00000308736.7B
ENST00000367669.8
ENST00000544844.6
```

#### `results/annotate/differentiate_per_isoform_statistics/{sample_id}.tsv`

Main output file, providing information on both splice and genetic variants affecting each protein isoform (if any).

```bash
isoform_id	protein_splice_category	n_variants	n_snvs	n_indels	n_heterozygous_variants	n_heterozygous_snvs	n_heterozygous_indels	n_homozygous_variants	n_homozygous_indels	n_homozygous_snvs	contributing_variants	diff	has_indel	has_changed_stop	has_frameshift	has_resolved_frameshift
ENST00000327044.7	pFSM	1	1	0	0	0	0	1	0	1	chr1:953279_T/C	300I>V	0	0	0	0
ENST00000338591.8	pFSM	0	0	0	0	0	0	0	0	0	NA	NA	0	0	0	0
ENST00000379410.8	pFSM	1	1	0	0	0	0	1	0	1	chr1:973858_G/C	487R>P	0	0	0	0
ENST00000379407.7	pFSM	1	1	0	0	0	0	1	0	1	chr1:973858_G/C	452R>P	0	0	0	0
ENST00000379409.6	pFSM	1	1	0	0	0	0	1	0	1	chr1:973858_G/C	539R>P	0	0	0	0
ENST00000379370.7	pFSM	0	0	0	0	0	0	0	0	0	NA	NA	0	0	0	0
ENST00000620552.4	pFSM	0	0	0	0	0	0	0	0	0	NA	NA	0	0	0	0
ENST00000652369.1	pFSM	0	0	0	0	0	0	0	0	0	NA	NA	0	0	0	0
ENST00000651234.1	pFSM	0	0	0	0	0	0	0	0	0	NA	NA	0	0	0	0
```

#### `results/annotate/get_directly_observed_and_linked_variants/{sample_id}/linked_variants.vcf`

VCF file containing variants that are on the same haplotype as a variant for which direct peptide evidence was found in the MS serach.

#### `results/annotate/get_directly_observed_and_linked_variants/{sample_id}/observed_variants.vcf`

VCF file containing variants for which direct peptide evidence was found in the MS serach.

### MS search results

#### `results/search_mass_spec/run_sage/{sample_id}/results.sage.tsv`

sage search results containing peptide identifications.

#### `results/annotate/map_peptides/{sample_id}/{haplotype}_Aligned.sortedByCoord.out.bam`

BAM files containing peptide sequence alignments to the reference genome for each haplotype (`homo` for peptides only containing homozygous variants, `first`, `second` for heterozygous, respectively). We note that for alignment, peptides were back-translated to genomic sequences, with silent variants ignored.

#### `results/annotate/map_isoforms/{sample_id}/{haplotype}_isoforms.sorted.bam`

BAM files containing isoform sequence alignments to the reference genome for each haplotype (`homo` for isoforms only containing homozygous variants, `first`, `second` for heterozygous, respectively). We note that for alignment, isoform sequences were back-translated to genomic sequences, with silent variants ignored.

## Limitations

Currently, we do not sure ensure that two A haplotypes (e.g., `ENST123A` and `ENST124A` necessarily inherit the same variants, i.e., it is possible that while they share a variant, the `A` notation does not ensure this in all cases).

## Citation

If you use this pipeline, please cite:

```bibtex
@article{wissel2025samplespecific,
  title={Sample-specific haplotype-resolved protein isoform characterization via long-read RNA-seq-based proteogenomics},
  author={Wissel, David and Sheynkman, Gloria M and Robinson, Mark D},
  journal={bioRxiv},
  pages={2025--11},
  year={2025},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Bibliography

[1] Wissel, David, Gloria M. Sheynkman, and Mark D. Robinson. "Sample-specific haplotype-resolved protein isoform characterization via long-read RNA-seq-based proteogenomics." bioRxiv (2025): 2025-11.

[2] Li, Heng. "Minimap2: pairwise alignment for nucleotide sequences." Bioinformatics 34.18 (2018): 3094-3100.

[3] Chen, Ying, et al. "Context-aware transcript quantification from long-read RNA-seq data with Bambu." Nature methods 20.8 (2023): 1187-1195.

[4] Varabyou, Ales, et al. "Investigating open reading frames in known and novel transcripts using ORFanage." Nature computational science 3.8 (2023): 700-708.

[5] Zheng, Zhenxian, et al. "Clair3-RNA: A deep learning-based small variant caller for long-read RNA sequencing data." bioRxiv (2025): 2024-11.

[6] Martin, Marcel, et al. "WhatsHap: fast and accurate read-based phasing." BioRxiv (2016): 085050.

[7] Spooner, William, et al. "Haplosaurus computes protein haplotypes for use in precision drug design." Nature Communications 9.1 (2018): 4128.

[8] Wright, James C., and Jyoti S. Choudhary. "DecoyPyrat: fast non-redundant hybrid decoy sequence generation for large scale proteomics." Journal of proteomics & bioinformatics 9.6 (2016): 176.

[9] Lazear, Michael R. "Sage: an open-source tool for fast proteomics searching and quantification at scale." Journal of Proteome Research 22.11 (2023): 3652-3659.

[10] Craig, Robertson, John P. Cortens, and Ronald C. Beavis. "Open source system for analyzing, validating, and storing protein identification data." Journal of proteome research 3.6 (2004): 1234-1242.

## Contact

For questions or issues, please contact [David Wissel](mailto:dwissel@ethz.ch) or open an issue in this repository.

## License

MIT
