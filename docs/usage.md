# phac-nml/viralassembly: Usage

## Introduction

This nextflow pipeline is intended for the reference-based assembly and analysis of viral sequencing data generated using either the Oxford Nanopore or Illumina sequencing platforms. The pipeline supports both amplicon and basic whole genome sequencing data for both platforms. It performs read processing, amplicon primer trimming, variant calling, consensus generation, variant annotation, and final quality-control reporting.

For Oxford Nanopore Technology (ONT) data, consensus variants are called using Clair3. Medaka and Nanopolish were deprecated as variant callers in [`v2.0.0`](https://github.com/phac-nml/measeq/releases/tag/2.0.0). Illumina data calls variants using FreeBayes by default with iVar as an option using `--use_ivar`.

For Amplicon Sequencing data it is at minimum required to:

1. Specify a path to the reads/input file
2. Specify a path to the reference genome
3. Specify an output directory
4. Specify the sequencing platform as `illumina` or `nanopore`
5. Specify the primer bed file that matches the amplicon sequencing scheme and reference genome

For Basic NGS Sequencing data it is at minimum required to:

1. Specify a path to the reads/input file
2. Specify a path to the reference genome
3. Specify an output directory
4. Specify the sequencing platform as `illumina` or `nanopore`

> [!NOTE]
> For nanopore data, the pipeline uses `r1041_e82_400bps_sup_v420` as the clair3 model by default. This is a required option when running nanopore data and can be adjusted using `--model <MODEL_NAME>` or by specifying a directory for a local model using `--local_model <PATH/TO/MODEL>`

## Index

- [phac-nml/viralassembly: Usage](#phac-nmlviralassembly-usage)
  - [Introduction](#introduction)
  - [Index](#index)
  - [Profiles](#profiles)
  - [Data Inputs](#data-inputs)
    - [Fastq Pass Directory (--fastq_pass)](#fastq-pass-directory---fastq_pass)
      - [Nanopore Input](#nanopore-input)
      - [Illumina Input](#illumina-input)
    - [Input CSV (--input)](#input-csv---input)
      - [Nanopore Example](#nanopore-example)
      - [Illumina Paired-End Example](#illumina-paired-end-example)
  - [Running the pipeline](#running-the-pipeline)
    - [Nanopore Data](#nanopore-data)
      - [Amplicon](#amplicon)
      - [Non-Amplicon](#non-amplicon)
    - [Illumina Data](#illumina-data)
      - [Amplicon](#amplicon-1)
      - [Non-Amplicon](#non-amplicon-1)
    - [Variant Callers](#variant-callers)
      - [Clair3 (Nanopore)](#clair3-nanopore)
      - [FreeBayes (Illumina)](#freebayes-illumina)
      - [iVar (Illumina)](#ivar-illumina)
    - [Other Run Notes](#other-run-notes)
    - [Updating the pipeline](#updating-the-pipeline)
    - [Reproducibility](#reproducibility)
  - [Input Parameters](#input-parameters)
    - [All Parameters](#all-parameters)
    - [Schemes and Reference](#schemes-and-reference)
    - [Metadata](#metadata)
    - [SnpEff](#snpeff)
    - [Virus Specification](#virus-specification)
      - [Currently supported viruses](#currently-supported-viruses)
    - [Virus-specific Processes](#virus-specific-processes)
    - [Nextclade](#nextclade)
      - [Explicitly Specifying Nextclade Dataset](#explicitly-specifying-nextclade-dataset)
      - [Dataset Selection Precedence](#dataset-selection-precedence)
  - [Core Nextflow Arguments](#core-nextflow-arguments)
    - [`-resume`](#-resume)
    - [`-c`](#-c)

## Profiles

Profiles are used to specify dependency installation, resources, and how to handle pipeline jobs. You can specify more than one profile but _avoid_ passing in more than one dependency managment profiles. They can be passed with `-profile <PROFILE>`

Available:

- `conda`: Utilize conda to install dependencies and environment management
- `mamba`: Utilize mamba to install dependencies and environment management
- `singularity`: Utilize singularity for dependencies and environment management
- `docker`: Utilize docker for dependencies and environment management

## Data Inputs

Two options for fastq data input: `--fastq_pass <FASTQ_PASS/>` or `--input <INPUT.csv>`

### Fastq Pass Directory (--fastq_pass)

Specify fastq data to input based on a given directory. The expected directory structure depends on the selected sequencing platform.

#### Nanopore Input

The directory can either contain barcoded directories (barcodexx), as would be seen after demultiplexing, or it could contain sample fastq files (one fastq per sample). The barcoded fastq data will be output with the barcode number but can be renamed with a [metadata tsv](#metadata) file input. The flat fastq files will keep their basename (separated out at the first `.`). Example:

Barcoded:

```
<fastq_pass>
├── barcode01
|   └── FAR41212_pass_barcode01_7d0222ac_0.fastq
├── barcode02
|   ├── FAR41212_pass_barcode02_7d0222ac_0.fastq
|   ├── FAR41212_pass_barcode02_7d0222ac_1.fastq
|   └── FAR41212_pass_barcode02_7d0222ac_2.fastq
└── barcode03
    └── FAR41212_pass_barcode03_7d0222ac_0.fastq
```

Flat:

```
<fastq_pass>
├── sample1.fastq
├── sample2.fastq
├── sample3.fastq
├── ntc.fastq
└── pos.fastq
```

#### Illumina Input

The directory should contain paired-end fastq files. Files are paired using `_R1` and `_R2` in their filenames.

Example Illumina directory:

```
<fastq_pass>
├── sample1_R1.fastq
├── sample1_R2.fastq
├── sample2_R1.fastq
├── sample2_R2.fastq
├── ntc_R1.fastq
├── ntc_R2.fastq
├── pos_R1.fastq
└── pos_R2.fastq
```

### Input CSV (--input)

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to pass in an input CSV file containing columns as follows depending on the sequencing platform used

| Column  | Required                     | Description                                                                               |
| ------- | ---------------------------- | ----------------------------------------------------------------------------------------- |
| sample  | Yes                          | Unique sample identifier                                                                  |
| fastq_1 | Yes                          | Path to the first (or only) FastQ file (`.fastq`, `.fq`, `.fastq.gz`, `fq.gz`)            |
| fastq_2 | Only for Illumina paired-end | Path to the second FastQ file for paired-end data (`.fastq`, `.fq`, `.fastq.gz`, `fq.gz`) |

#### Nanopore Example

| sample  | fastq_1                  |
| ------- | ------------------------ |
| sample1 | /path/to/sample.fastq    |
| sample2 | /path/to/sample2-1.fastq |
| sample3 | /path/to/sample-2.fastq  |
| ntc     | /path/to/control.fastq   |
| pos     | /path/to/pos.fastq       |

#### Illumina Paired-End Example

| sample  | fastq_1                    | fastq_2                    |
| ------- | -------------------------- | -------------------------- |
| sample1 | /path/to/sample_R2.fastq   | /path/to/sample_R2.fastq   |
| sample2 | /path/to/sample2_R1.fastq  | /path/to/sample2_R2.fastq  |
| sample3 | /path/to/sample-3_R1.fastq | /path/to/sample-3_R2.fastq |
| ntc     | /path/to/control_R1.fastq  | /path/to/control_R2.fastq  |
| pos     | /path/to/pos_R1.fastq      | /path/to/pos_R2.fastq      |

> [!NOTE]
> For nanopore data, the `fastq_2` column is not required. However, it shouldn't cause any issues if the column is included but the values are empty.

## Running the pipeline

### Nanopore Data

#### Amplicon

The typical command for running the pipeline with an [amplicon scheme](#schemes-and-reference) and a non-default Clair3 model with Docker for nanopore sequenced data is as follows:

```bash
nextflow run phac-nml/viralassembly \
  -profile docker \
  --platform nanopore \
  --fastq_pass FASTQ_PASS/ \
  --model 'r1041_e82_400bps_sup_v4.3.0' \
  --reference REF.fasta \
  --primer_bed PRIMER.bed \
  --outdir results
```

This will launch the pipeline with the `docker` configuration profile and use the reference and primer files supplied. Profile information [can be found above](#profiles)

#### Non-Amplicon

The typical command for running the pipeline without an amplicon scheme with the default clair3 model for nanopore sequenced data is as follows:

```bash
nextflow run phac-nml/viralassembly \
  -profile singularity \
  --platform nanopore \
  --fastq_pass FASTQ_PASS/ \
  --reference REF.fa \
  --outdir ./results
```

This will launch the pipeline with the `singularity` configuration profile and the specified reference. Profile information [can be found above](#profiles)

### Illumina Data

#### Amplicon

The typical command for running the pipeline with an [amplicon scheme](#schemes-and-reference) with Docker for illumina sequenced data is as follows:

```bash
nextflow run phac-nml/viralassembly \
  -profile docker \
  --platform illumina \
  --fastq_pass FASTQ_PASS/ \
  --reference REF.fasta \
  --primer_bed PRIMER.bed \
  --outdir results
```

This will launch the pipeline with the `docker` configuration profile and use the reference and primer files supplied. Profile information [can be found above](#profiles)

#### Non-Amplicon

The typical command for running the pipeline without an amplicon scheme with Docker for illumina sequenced data is as follows:

```bash
nextflow run phac-nml/viralassembly \
  -profile docker \
  --platform illumina \
  --fastq_pass FASTQ_PASS/ \
  --reference REF.fasta \
  --outdir results
```

This will launch the pipeline with the `docker` configuration profile and use the reference supplied. Profile information [can be found above](#profiles)

> [!TIP]
> More example commands are available in the [example commands document](./example_commands.md).

### Variant Callers

The pipeline currently supports three different variant callers, one Nanopore specific (Clair3) and two Illumina specific (FreeBayes & iVar). As of [`v2.0.0`](https://github.com/phac-nml/measeq/releases/tag/2.0.0), Medaka and Nanopolish were deprecated as variant callers and Clair3 has been set as the default variant caller for Nanopore data.

#### [Clair3](https://github.com/HKU-BAL/Clair3) (Nanopore)

Clair3 is a germline small variant caller for long-reads.

Running the pipeline with Nanopore data supports the following optional parameters related to Clair3:

- `--model <MODEL>`: Specify the base clair3 model
- `--local_model </PATH/TO/downloaded_clair3_model>`: Specify the path to a local downloaded model directory
- `--no_pool_split`: Do not split reads and variant calls by amplicon primer pool, instead call all variants at once

Clair3 comes with some models available and is defaulted to `r1041_e82_400bps_sup_v420`. Additional models available using the `--model` command will be automatically downloaded from source. You can also manually download models and then specify them with `--local_model </PATH/TO/downloaded_clair3_model>` to save having to download it each time or to use other models. Remember to pick a model that best represents the data!

#### [FreeBayes](https://github.com/freebayes/freebayes) (Illumina)

FreeBayes is a Bayesian genetic variant detector designed to find small polymorphisms. FreeBayes is the default variant caller for Illumina data and requires no additional parameters.

#### [iVar](https://github.com/andersen-lab/ivar) (Illumina)

iVar variants is part of the iVAR computational package for viral sequencing functions. iVar is used as an alternate variant caller and can be invoked with the `--use_ivar` parameter.

### Other Run Notes

Note that both analysis methods of the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

For example, an amplicon based nanopore sequencing pipeline run specified with a params file in yaml format:

```bash
nextflow run phac-nml/viralassembly -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
platform: "nanopore"
fastq_pass: "./fastq_pass"
reference: "reference.fa"
outdir: "./results/"
model: "r1041_e82_400bps_sup_v4.3.0"
primer_bed: "PRIMER.bed"
```

### Updating the pipeline

When you install the pipeline, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull phac-nml/viralassembly
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [phac-nml/viralassembly releases page](https://github.com/phac-nml/viralassembly/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Input Parameters

Use `--help` to see all options formatted on the command line

Use `--version` to see version information

### All Parameters

| Parameter                                                 | Description                                                                        | Type    | Default                     | Notes                                                                                                                                            |
| --------------------------------------------------------- | ---------------------------------------------------------------------------------- | ------- | --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Input and Output Parameters**                           |                                                                                    |         |                             |                                                                                                                                                  |
| --fastq_pass                                              | Path to directory containing `barcode##` subdirectories OR `*.fastq*` files        | Path    | null                        | See [Option for input params](#input-parameters)                                                                                                 |
| --input                                                   | Path to samplesheet with information about the samples you would like to analyse   | Path    | null                        | See [Option for input params](#input-parameters)                                                                                                 |
| --outdir                                                  | Directory name to output results to                                                | String  | null                        | Required                                                                                                                                         |
| --virus_name                                              | Virus name that sets virus specific processes and nextclade runs                   | String  | null                        | See [Virus Specification](#virus-specification)                                                                                                  |
| **Required Parameters**                                   |                                                                                    |         |                             |                                                                                                                                                  |
| --platform                                                | Sequencing platform                                                                | String  | null                        | Required: `nanopore` or `illumina`                                                                                                               |
| --reference                                               | Path to a reference FASTA file to run pipeline                                     | Path    | null                        | Required: See [Schemes and Reference](#schemes-and-reference)                                                                                    |
| **Amplicon Parameter**                                    |                                                                                    |         |                             |                                                                                                                                                  |
| --primer_bed                                              | Path to an amplicon primer bed file associated with the reference                  | Path    | null                        | See [Schemes and Reference](#schemes-and-reference)                                                                                              |
| **Nanopore Variant-Calling Parameters (Clair3)**          |                                                                                    |         |                             |                                                                                                                                                  |
| --model                                                   | Clair3 model name                                                                  | String  | 'r1041_e82_400bps_sup_v420' | Default model will not work the best for all inputs. [See clair3 docs](https://github.com/HKU-BAL/Clair3#pre-trained-models) for additional info |
| --local_model                                             | Path to clair3 local model directory to use instead of `--model`                   | Path    | ''                          | Default model will not work the best for all inputs. [See clair3 docs](https://github.com/HKU-BAL/Clair3#pre-trained-models) for additional info |
| --no_pool_split                                           | Do not split reads into separate pools                                             | Boolean | False                       | Nanopore amplicon sequencing only                                                                                                                |
| --min_qual_clair3                                         | Minumum Clair3 variant quality to keep a variant                                   | Integer | 7                           |                                                                                                                                                  |
| --min_frameshift_qual                                     | Minumum Clair3 variant quality to keep a frameshift variant                        | Integer | 15                          | Frameshift defined as not divisible by 3                                                                                                         |
| --min_allele_freq                                         | Minimum allele frequency to call a variant                                         | Number  | 0.60                        |                                                                                                                                                  |
| --min_mask_freq                                           | Minimum allele frequency at which a position may be masked with N                  | Number  | 0.25                        |                                                                                                                                                  |
| **Nanopore Minor Variant-Calling Parameters (ClairS-To)** |                                                                                    |         |                             |                                                                                                                                                  |
| --minor_variants                                          | Enable Nanopore minor variant calling with ClairS-TO                               | Boolean | False                       |                                                                                                                                                  |
| --clairsto_model                                          | ClairS-TO model name                                                               | String  | 'r1041_e82_400bps_sup_v420' | Default model will not work the best for all inputs.                                                                                             |
| --min_snv_af_ClairS                                       | Minimum allele frequency for minor single nuclotide variants                       | Number  | 0.05                        |                                                                                                                                                  |
| --min_indel_af_ClairS                                     | Minimum allele frequency for minor insertions and deletions                        | Number  | 0.15                        |                                                                                                                                                  |
| --min_cov_ClairS                                          | Minimum coverage required for minor variants                                       | Integer | 20                          |                                                                                                                                                  |
| --min_qual_ClairS                                         | Minimum quality required for minor variants                                        | Integer | 5                           |                                                                                                                                                  |
| **Illumina Variant-Calling Parameters**                   |                                                                                    |         |                             |                                                                                                                                                  |
| --use_ivar                                                | Use iVar as a variant caller instead of Freebayes                                  | Boolean | False                       |                                                                                                                                                  |
| --min_ambiguity_threshold                                 | Minimum allele frequency for IUPAC ambiguity reporting instead of reference allele | Number  | 0.25                        |                                                                                                                                                  |
| --max_ambiguity_threshold                                 | Maximum allele frequency for IUPAC ambiguity reporting instead of majority allele  | Number  | 0.75                        |                                                                                                                                                  |
| --min_indel_threshold                                     | Minimum allele frequency threshold for retaining an indel                          | Number  | 0.60                        |                                                                                                                                                  |
| --min_alt_threshold_illumina                              | Minimum fraction of observations supporting an alt allele to evaluate position     | Number  | 0.10                        |                                                                                                                                                  |
| --min_variant_qual_freebayes                              | Minumum FreeBayes variant quality to keep a variant                                | Integer | 20                          |                                                                                                                                                  |
| **Read Filtering Parameters**                             |                                                                                    |         |                             |                                                                                                                                                  |
| --min_length                                              | Maximum read length to be kept                                                     | Integer | 8000                        |                                                                                                                                                  |
| --max_length                                              | Maximum read length to be kept                                                     | Integer | Illumina: 50, Nanopore: 200 |                                                                                                                                                  |
| --min_reads                                               | Minimum number of reads required for a sample after filtering                      | Integer | 20                          |                                                                                                                                                  |
| --fastp_adapter                                           | Path to FASTA file containing adapter sequences for Fastp                          | Path    | null                        |                                                                                                                                                  |
| **General Analysis Parameters**                           |                                                                                    |         |                             |                                                                                                                                                  |
| --min_depth                                               | Minimum depth required to call a consensus position                                | Integer | 20                          | Positions with lower depth are masked with an N                                                                                                  |
| --metadata                                                | Path to metadata TSV file with columns 'sample' and 'barcode'                      | Path    | null                        | See [Metadata](#metadata) for more info                                                                                                          |
| --use_artic_tool                                          | Run the artic tool itself instead of nextflow implementation for nanopore data     | Boolean | False                       |                                                                                                                                                  |
| --normalise                                               | Target amplicon coverage used for normalization                                    | Integer | 2000                        | Entering `0` turns off normalisation. Only for amplicon sequencing                                                                               |
| --no_frameshift                                           | Filter INDEL variants that are not divisible by 3                                  | Boolean | False                       | Simple `%3 == 0` check for variants                                                                                                              |
| **SnpEff Parameters**                                     |                                                                                    |         |                             |                                                                                                                                                  |
| --skip_snpeff                                             | Skip running SnpEff                                                                | Boolean | False                       |                                                                                                                                                  |
| --gff                                                     | Path to gff3 formatted file to use in SnpEff database build                        | Path    | ''                          | Not required to run [SnpEff](#snpeff). See below for details                                                                                     |
| **Nextclade Parameters**                                  |                                                                                    |         |                             |                                                                                                                                                  |
| --skip_nextclade                                          | Skip running Nextclade                                                             | Boolean | False                       |                                                                                                                                                  |
| --nextclade_dataset_dir                                   | Path to local Nextclade dataset directory                                          | Path    | null                        | Not required to run [Nextclade](#nextclade). See below for details                                                                               |
| --nextclade_dataset_name                                  | Name of the Nextclade dataset to use                                               | Sting   | null                        | Not required to run [Nextclade](#nextclade). See below for details                                                                               |
| --nextclade_dataset_tag                                   | Dataset tag or version                                                             | String  | null                        | Not required to run [Nextclade](#nextclade). See below for details                                                                               |
| **Virus-Specific Parameters**                             |                                                                                    |         |                             |                                                                                                                                                  |
| --skip_pangolin                                           | Skip Pangolin analysis for SARS-Cov-2                                              | Boolean | False                       |                                                                                                                                                  |
| --pango_database                                          | Path to local Pangolin data directory                                              | Path    | null                        | Not required to run [Pangolin](#virus-specific-processes). See below for details                                                                 |
| **Quality-Control Parameters**                            |                                                                                    |         |                             |                                                                                                                                                  |
| --skip_qc                                                 | Skip running all QC and reporting steps                                            | Boolean | False                       |                                                                                                                                                  |
| --pcr_primer_bed                                          | Path to PCR primer bed file to check for mutations against                         | Path    | ''                          | For output QC checks                                                                                                                             |
| --neg_control_threshold                                   | Coverage threshold at which to fail negative control samples                       | Number  | 0.10                        |                                                                                                                                                  |
| --neg_ctrl_substrings                                     | Negative control sample substrings separated by a `,`                              | String  | 'ntc,neg,blank,water'       |                                                                                                                                                  |
| **Reporting Parameters**                                  |                                                                                    |         |                             |                                                                                                                                                  |
| --multiqc_report                                          | Run MultiQC report over custom report                                              | Boolean | False                       |                                                                                                                                                  |

### Schemes and Reference

Amplicon schemes are a highly targeted approach to sequencing focusing on a specific target genome. If using an amplicon scheme with this pipeline, a 7 column primer bed file is required along with the reference fasta file. This primer file will be used to trim the BAM file so that variants in primer regions are not masked out by the primers themselves.

A primer bed file titled should be organized to fit the following specifications based on using the [ARTIC/Primalbedtols v3.0.0 specs](https://github.com/artic-network/primerscheme-specs/blob/20816ff7cd53bdfaab7a605dee06de1c80be759f/pdf/primerscheme.pdf) which is minimally:

- Minimum of 7 columns
- Chrom
- Start
- End
- Primer pairs with names containing `_LEFT` and `_RIGHT`
- Primer pool numbers (1, 2, 3, etc.)
- Primer direction (+ / -)
- Primer Sequence

Example primer file format:

| MN908947.3 | 30    | 54   | nCoV-2019_1_LEFT  | 1           | +         | ATCCCGATTT |
| ---------- | ----- | ---- | ----------------- | ----------- | --------- | ---------- |
| MN908947.3 | 1183  | 1205 | nCoV-2019_1_RIGHT | 1           | -         | TTAAGCGCGC |
| MN908947.3 | 1100  | 1128 | nCoV-2019_2_LEFT  | 2           | +         | AGGGTCAGCA |
| MN908947.3 | 2244  | 2266 | nCoV-2019_2_RIGHT | 2           | -         | CCTAAGCTAG |
| ...        | ...   | ...  | ...               | ...         | ...       | ...        |
| REF ID     | Start | Stop | Primer Name       | Primer Pool | Direction | Primer Seq |

### Metadata

Input metadata is used to rename barcoded fastq files along with adding additional lines to the final overall QC csv file. Note that the metadata input is expected to be of a `TSV` format

Structure for example `metadata.tsv` file:

| sample | barcode | \<Anything else you want to add > |
| ------ | ------- | --------------------------------- |
| SR-1   | 1       | X                                 |
| SR-2   | 02      | Y                                 |
| NTC-12 | 12      | Z                                 |

### SnpEff

SnpEff is run by default on all non-segmented viruses (due to current implementation) by using the reference sequence ID to either:

1. Check if there is a SnpEff database available to download
2. Build a SnpEff database by downloading the sequence genbank file from NCBI

Instead of relying on the reference ID to build/download a database, you can instead specify a gff3 file with `--gff <PATH/TO/file.gff>` to be used with the reference sequence to create the SnpEff database

If building/downloading a database fails, the pipeline will skip over running SnpEff instead of failing out completely.

SnpEff can also be skipped entirely by passing the `--skip_snpeff` parameter

### Virus Specification

While the pipeline will run on any viral data, it also currently supports specifying a virus name to run specific analyses. Presently, it is used mostly for Nextclade dataset configuration. However, future versions of this pipeline aim to use the virus specfication as a way to set pipeline defaults and run virus specific processes.

#### Currently supported viruses

The following is a list of viruses supported by the pipeline for automatic nextclade dataset configuration. The list includes the full virus name and the abbreviation to be used with the `--virus_name` argument.

| Virus Name                    | Abbreviation for `--virus_name` |
| ----------------------------- | ------------------------------- |
| SARS-CoV-2                    | covid                           |
| Respiratory syncytial virus A | rsv_a                           |
| Respiratory syncytial virus B | rsv_b                           |
| Mpox                          | mpox                            |
| Ebola                         | ebola                           |
| Bundibugyo ebolavirus         | bsbv                            |
| Sudan                         | sudan                           |
| Measles                       | measles                         |
| Dengue                        | dengue                          |
| Yellow Fever                  | yfv                             |
| Human metapneumovirus         | hmpv                            |
| Varicella-Zoster              | vzv                             |
| Rubella                       | rubella                         |
| Mumps                         | mumps                           |
| West Nile                     | wnv                             |

### Virus-specific Processes

As indicated above, the pipeline aims to use the `--virus_name` parameter to run virus-specific processes. This is currently in the development phase and will be added as more virus-specific processes are identified based on needs at the National Microbiology Laboratory. We currently support [Pangolin](https://github.com/cov-lineages/pangolin) as a virus-specific process when the pipeline is invoked with `--virus_name covid` as a parameter. More details and processes will be added in later versions.

### Nextclade

Nextclade provides clade assignment, mutation calling, and consensus quality reporting. By default, the pipeline uses [`nextclade sort`](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/reference.html#nextclade-sort) to detect the most appropriate Nextclade dataset for each sample. For segmented viruses, dataset detection is performed seperately for each segment.

#### Explicitly Specifying Nextclade Dataset

The Nextclade dataset can also be specified explicitly in multiple ways. This is useful in cases where you want to ensure a particular dataset is used. For example, SARS-CoV-2 has multiple datasets available in the Nextclade datasets repository and `nextclade sort` may select an alternative SARS-CoV-2 dataset when processing a COVID-19 sample. As such, a dataset can be explicitly specified through one of the following methods:

1. A dataset can be downloaded directly by specifying:

   ```bash
   --nextclade_dataset_name  <DATASET_NAME>
   --nextclade_dataset_tag   <TAG>           ##Optional
   ```

   > [!NOTE]
   > The optional dataset tag is used to download a specific version of the dataset from the Nextclade datasets repository. The dataset tag can only be used together with the `--nextclade_dataset_name` parameter. If there is no tag specified, then the pipeline will download the latest version of the dataset specified.

2. A locally downloaded or custom Nextclade dataset can be supplied by specifying the directory containg the dataset:

   ```bash
   --nextclade_dataset_dir <PATH/TO/DATASET>
   ```

3. A virus name value with the `--virus_name` parameter

   ```bash
   --virus_name <VIRUS_NAME>
   ```

   Specifying the `--virus_name` parameter selects the pipeline's configured default nextclade dataset for that virus. [See above for more information](#currently-supported-viruses).

#### Dataset Selection Precedence

> [!WARNING]
> The `--nextclade_dataset_dir` and `--nextclade_dataset_name` parameters are mutually exclusive and can't be specified in the same run. They have the same level of precedence when determining which Nextclade dataset is used.

The dataset selection precedence is:

1. `--nextclade_dataset_dir` or `--nextclade_dataset_name`
2. Dataset configured through `--virus_name`
3. Automatic dataset detection using `nextclade sort`

`--virus_name` can be used together with either `--nextclade_dataset_dir` or `--nextclade_dataset_name`. In this case, the explicitly supplied nextclade dataset takes precedence only for Nextclade dataset selection, while the other virus-specific processes associated with `--virus_name` will continue to run normally.

For example:

```bash
--virus_name covid
--nextclade_dataset_name <DATASET_NAME>
```

will use `<DATASET_NAME>` for nextclade instead of the dataset configured by `--virus_name covid`, while any other SARS-CoV-2 virus-specific processes enabled by `--virus_name` will still be performed.

> [!TIP]
> Nextclade can be skipped entirely by passing the `--skip_nextclade` parameter.

## Core Nextflow Arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.
