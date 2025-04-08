# ensoma-trace-isa-nextflow

A Nextflow pipeline for analyzing trace insertion site analysis (ISA) data for EnSoma's genetic payload delivery system.

## Overview

This pipeline processes sequencing data from insertion site analysis experiments to identify where genetic payloads have integrated into a host genome. It handles paired-end FASTQ files, performs quality control, trimming, alignment, and ultimately identifies integration sites.

## Usage

The pipeline can be run with either a local installation of the required software or using Docker (preferred method).

### Running with Docker (Recommended)

```bash
nextflow run main.nf -profile docker --sample_sheet path/to/sample_sheet.csv --fastq_dir path/to/fastq/files --genome_dir path/to/genome/files --outdir results
```

### Running with Local Installation

```bash
nextflow run main.nf --sample_sheet path/to/sample_sheet.csv --fastq_dir path/to/fastq/files --genome_dir path/to/genome/files --outdir results
```

## Parameters

| Parameter | Description | Required | Input Type |
|-----------|-------------|----------|------------|
| **Input/Output Options** | | | |
| experiment_name | The name of the TRACE ISA experiment | Yes | String (default: TRACE_ISA) |
| outdir | The output directory where the results will be saved | Yes | Directory path (default: results) |
| fastq_dir | The directory containing the FASTQ files | Yes | Directory path |
| genome_dir | The directory containing the payload and genome FASTA files | Yes | Directory path |
| sample_sheet | The CSV sample sheet | Yes | File path (CSV) |
| **FASTQ Trimming Options** | | | |
| R1_adapter | The sequence of the R1 adapter in the reads | No | String (default: NNNNCGAGTTTTAATGACTCCAACT) |
| R2_adapter | The sequence of the R2 adapter in the reads | No | String (default: AGTGGCACAGCAGTTAGGNNNNNNNNAGATGTGTATAAGAGACAG) |
| IR_seq | The sequence of the IR in the reads | No | String (default: TAAGTGTATGTAAACTTCCGACTTCAACTG) |
| min_length | Minimum FASTQ read length allowed | No | Integer (default: 20, range: 1-500) |
| allowed_errors_r1_adapter | Allowed errors for expected R1 adapter sequence | No | Integer (default: 4, range: 0-10) |
| allowed_errors_ir | Allowed errors for expected IR sequence | No | Integer (default: 4, range: 0-10) |
| allowed_errors_adapter_readthrough | Allowed errors to detect adapter readthrough | No | Integer (default: 4, range: 0-10) |
| allowed_errors_ir_readthrough | Allowed errors to detect IR readthrough | No | Integer (default: 3, range: 0-10) |
| **Alignment Filtering Options** | | | |
| exclude_flags | Flags for samtools to exclude from alignments | No | Integer (default: 3852, range: 0-10000) |
| required_flags | Flags for samtools to require for alignments | No | Integer (default: 3, range: 0-10000) |
| min_mapq | The minimum mapping quality score to retain a read | No | Integer (default: 30, range: 0-60) |
| softclip_threshold | The maximum number of allowed 5' softclipped bases | No | Integer (default: 5, range: 0-100) |
| **Insertion Site Options** | | | |
| grouping_distance | Group integration sites within this distance | No | Integer (default: 5, range: 1-100) |
| grouping_mode | Whether to take the mean or median position for integration sites with the same score | No | String (default: median, options: mean/median) |

## Sample Sheet Format

The sample sheet is a CSV file that must contain the following columns:

1. **fastq_id**: Identifier matching the prefix of FASTQ files (e.g., 'A2063-1' for 'A2063-1_S21_L001_R1_001.fastq.gz')
2. **sample_name**: The name of the sample
3. **technical_group**: Group identifier for technical replicates
4. **payload_vector_fasta**: Filename of the payload vector FASTA (located in genome_dir)
5. **integration_vector_fasta**: Filename of the integration vector FASTA (located in genome_dir)
6. **genome_fasta**: Filename of the reference genome FASTA (located in genome_dir)

Example sample sheet:

```csv
fastq_id,sample_name,technical_group,payload_vector_fasta,integration_vector_fasta,genome_fasta
sample-1,donor1_treatment1_rep1,donor1_treatment1,payload.fa,vector.fa,hg38.fa
sample-2,donor1_treatment1_rep2,donor1_treatment1,payload.fa,vector.fa,hg38.fa
sample-3,donor2_treatment1,donor2_treatment1,payload2.fa,vector.fa,hg38.fa
```

### How the Sample Sheet Works

- **fastq_id**: This identifier is used to match sample information to the corresponding FASTQ files. The pipeline expects paired-end reads with filenames following the pattern `{fastq_id}_*_R{1,2}_*.fastq.gz`. For example, if fastq_id is 'A2063-1', it will match files like 'A2063-1_S21_L001_R1_001.fastq.gz' and 'A2063-1_S21_L001_R2_001.fastq.gz'.

- **References**: The three FASTA file columns (payload_vector_fasta, integration_vector_fasta, genome_fasta) point to files within the specified genome_dir. These files are concatenated to create a composite reference for alignment, allowing the pipeline to identify reads that map to payload vectors, integration vectors, or the host genome.

## Pipeline Steps

1. **Extract UMIs**: Extracts unique molecular identifiers from reads
2. **Trim FASTQs**: Removes adapters and filters reads based on quality and length
3. **Reference Preparation**: Concatenates reference FASTA files
4. **Alignment**: Aligns trimmed reads to the concatenated reference using BWA
5. **Filtering**: Filters alignments based on mapping quality and other criteria
6. **Integration Site Calling**: Identifies potential integration sites
7. **Aggregation**: Aggregates and summarizes integration sites across samples

## Outputs

The pipeline produces the following main outputs in the specified output directory:

### Trimming Results
- `trimmed_fastqs/`: Contains trimmed FASTQ files after adapter removal
- `umi_extracted/`: FASTQ files with extracted UMI information

### Alignment Results
- `alignments/`: Contains BAM files of reads aligned to the reference
  - Raw unfiltered alignments
  - Filtered alignments based on quality criteria

### Integration Site Results
- `integration_sites/`: Contains identified integration sites
  - Individual BED files for each sample
  - Aggregated BED file combining all samples
  - Summary statistics and visualization files

### Reference Files
- `genome/`: Contains processed reference files
  - Concatenated reference sequences
  - BWA indices

## Requirements

- Nextflow 22.10.0 or later
- Docker (if using the docker profile)
- Alternatively, the following tools installed locally:
  - BWA
  - SAMtools
  - BEDTools
  - FastQC
  - Cutadapt
  - Python 3.7+

## License

[License information]
