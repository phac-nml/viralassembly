# phac-nml/viralassembly: Usage

## Introduction

This pipeline is intended to be run on either Nanopore Amplicon Sequencing data or Basic Nanopore NGS Sequencing data that can utilize a reference genome for read mapping, variant calling, and other downstream analyses. It generates variant calls, consensus sequences, and quality control information based on the reference. To do this, there are three different variant callers that can be utilized which includes: `clair3`, `medaka`, and `nanopolish` (For R9.4.1 flowcells and below only). By default, `clair3` is used and highly recommended over the other two (as of December 2025). Their inclusion being for completeness with eventual removal from the workflow.

For Amplicon Sequencing data it is at minimum required to:

1. Specify a path to the reads/input file
2. Specify the reference genome
3. Specify the primer bed file that matches the amplicon sequencing scheme and reference genome
4. Pick a variant caller and caller model

For Basic NGS Sequencing data it is at minimum required to:

1. Specify a path to the reads/input file
2. Specify a path to the reference genome
3. Pick a variant caller and caller model

## Index

- [Profiles](#profiles)
- [Data Inputs](#data-inputs)
  - [Fastq Pass Directory](#fastq-pass-directory---fastq_pass)
  - [Input CSV](#input-csv---input)
- [Variant Callers](#variant-callers)
  - [Clair3](#clair3)
  - [Medaka](#medaka)
  - [Nanopolish](#nanopolish)
- [Running the Pipeline](#running-the-pipeline)
  - [Amplicon](#amplicon)
  - [Non-Amplicon](#non-amplicon)
  - [Other Run Note](#other-run-notes)
  - [Updating the Pipeline](#updating-the-pipeline)
  - [Reproducibility](#reproducibility)
- [Input Parameters](#input-parameters)
  - [All Parameters](#all-parameters)
  - [Schemes and Reference](#schemes-and-reference)
  - [Metadata](#metadata)
  - [SnpEff](#snpeff)
- [Core Nextflow Arguments](#core-nextflow-arguments)

## Profiles

Profiles are used to specify dependency installation, resources, and how to handle pipeline jobs. You can specify more than one profile but _avoid_ passing in more than one dependency managment profiles. They can be passed with `-profile <PROFILE>`

Available:

- `conda`: Utilize conda to install dependencies and environment management
- `mamba`: Utilize mamba to install dependencies and environment management
- `singularity`: Utilize singularity for dependencies and environment management
- `docker`: Utilize docker to for dependencies and environment management

## Data Inputs

Two options for fastq data input: `--fastq_pass <FASTQ_PASS/>` or `--input <INPUT.csv>`

### Fastq Pass Directory (--fastq_pass)

Specify fastq data to input based on a given directory. The directory can either contain barcoded directories (barcodexx), as would be seen after demultiplexing, or it could contain sample fastq files (one fastq per sample). The barcoded fastq data will be output with the barcode number but can be renamed with a [metadata tsv](#metadata) file input. The flat fastq files will keep their basename (separated out at the first `.`). Example:

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

### Input CSV (--input)

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to pass in an input CSV file containing 2 columns, `sample`, and `fastq_1` where:

- `sample` is the sample name to use
- `fastq_1` is the path to one fastq file per sample in `.fastq*` format

Ex.
| sample | fastq_1 |
| - | - |
| sample1 | /path/to/sample.fastq |
| sample2 | /path/to/sample2-1.fastq |
| sample2 | /path/to/sample-2.fastq |
| ntc | /path/to/control.fastq |
| pos | /path/to/pos.fastq |

## Variant Callers

Three different variant callers are available with slightly different options regarding running with them. For the most accurate results when running with `Clair3` or `medaka` pick a model that best matches the input data!! There is a default set for both but know your data and pick what suits it best!

### [Clair3](https://github.com/HKU-BAL/Clair3)

Clair3 is a germline small variant caller for long-reads.

Running with `clair3` requires the following parameters:

- `--variant_caller clair3`: Sets clair3 as the variant caller

And has the optional parameters of:

- `--clair3_model <MODEL>`: Specify the base clair3 model
- `--clair3_local_model </PATH/TO/downloaded_clair3_model>`: Specify the path to a local downloaded model directory
- `--clair3_no_pool_split`: Do not split inputs into pools

Clair3 comes with some models available and is defaulted to `r1041_e82_400bps_sup_v420`. Additional models can be downloaded from [ONT Rerio](https://github.com/nanoporetech/rerio/tree/master) and then specified in the `--clair3_local_model </PATH/TO/downloaded_clair3_model>` parameter shown above. Remember to pick a model that best represents the data!

### [Medaka](https://github.com/nanoporetech/medaka)

Medaka is a tool to create consensus sequences and variant calls from nanopore sequencing data using neural networks and provied by ONT.

Running with `medaka` requires the following parameters:

- `--variant_caller medaka`: Sets medaka as the variant caller

And has the optional parameters of:
`--medaka_model <MODEL>`: Specify the wanted medaka model

Medaka models come built in with the tool itself with the default set to `r941_min_hac_g507` which can be changed with `--medaka_model <MODEL>` parameter shown above. More information on models [can be found here](https://github.com/nanoporetech/medaka#models). Remember to pick a model that best represents the data!

### [Nanopolish](https://github.com/jts/nanopolish)

Nanopolish is a software package for signal-level analysis of Oxford Nanopore sequencing data. It _does not presently support the R10.4 flowcells_ so as a variant caller it should only be used with R9.4 flowcells.

Running with `nanopolish` requires the following parameters:

- `--variant_caller nanopolish`
- `--fast5_pass <FAST5_PASS/>`
- `--sequencing_summary <SEQ_SUM.txt>`

Nanopolish requires the fast5 directory along with the sequencing summary file to be used as input instead of a model. As such, nanopolish requires that the read ids in the fastq files are linked by the sequencing summary file to their signal-level data in the fast5 files. This makes it **a lot** easier to run using barcoded directories but it can be done with individual read files

## Running the pipeline

### Amplicon

The typical command for running the pipeline with an [amplicon scheme](#schemes-and-reference) using Clair3 (default variant caller) with Docker is as follows:

```bash
nextflow run phac-nml/viralassembly \
  -profile docker \
  --fastq_pass FASTQ_PASS/ \
  --clair3_model 'r1041_e82_400bps_sup_v4.3.0' \
  --reference REF.fasta \
  --primer_bed PRIMER.bed \
  --outdir results
```

This will launch the pipeline with the `docker` configuration profile, the `Clair3` variant caller, and use the reference and primer files supplied. Profile information [can be found above](#profiles)

### Non-Amplicon

The typical command for running the pipeline without an amplicon scheme using `Clair3` and a different clair3 model is as follows:

```bash
nextflow run phac-nml/viralassembly \
  -profile singularity \
  --fastq_pass FASTQ_PASS/ \
  --medaka_model 'r1041_e82_400bps_sup_v4.3.0' \
  --reference REF.fa \
  --outdir ./results
```

This will launch the pipeline with the `singularity` configuration profile, the `Clair3` variant caller, and the specified reference. Profile information [can be found above](#profiles)

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

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run phac-nml/viralassembly -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
fastq_pass: "./fastq_pass"
variant_caller: "clair3"
clair3_model: "r1041_e82_400bps_sup_v4.3.0"
reference: "reference.fa"
outdir: "./results/"
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull phac-nml/viralassembly
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [phac-nml/viralassembly releases page](https://github.com/phac-nml/viralassembly/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Input Parameters

Use `--help` to see all options formatted on the command line

Use `--version` to see version information

### All Parameters

| Parameter               | Description                                                                                              | Type    | Default                     | Notes                                                                                                                                            |
| ----------------------- | -------------------------------------------------------------------------------------------------------- | ------- | --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| --fastq_pass            | Path to directory containing `barcode##` subdirectories OR Path to directory containing `*.fastq*` files | Path    | null                        | [Option for input params](#input-parameters)                                                                                                     |
| --input                 | Path to samplesheet with information about the samples you would like to analyse                         | Path    | null                        | [Option for input params](#input-parameters)                                                                                                     |
| --variant_caller        | Pick from the 3 variant callers: 'clair3', 'medaka', 'nanopolish'                                        | Choice  | 'clair3'                    | [Details above](#introduction)                                                                                                                   |
| --reference             | Required: Specify the path to a reference fasta file to run pipeline                                     | Path    | ''                          | See [schemes and reference](#schemes-and-reference)                                                                                              |
| --primer_bed            | Specify the path to an amplicon primer bed file associated with the reference                            | Path    | ''                          | See [schemes and reference](#schemes-and-reference)                                                                                              |
| --clair3_model          | Clair3 base model to be used in the pipeline                                                             | String  | 'r1041_e82_400bps_sup_v420' | Default model will not work the best for all inputs. [See clair3 docs](https://github.com/HKU-BAL/Clair3#pre-trained-models) for additional info |
| --clair3_local_model    | Path to clair3 additional model directory to use instead of a base model                                 | Path    | ''                          | Default model will not work the best for all inputs. [See clair3 docs](https://github.com/HKU-BAL/Clair3#pre-trained-models) for additional info |
| --clair3_no_pool_split  | Do not split reads into separate pools                                                                   | Boolean | False                       | Clair3 amplicon sequencing only                                                                                                                  |
| --min_qual_clair3       | Minumum Clair3 variant quality to keep a variant                                                         | Integer | 8                           |                                                                                                                                                  |
| --min_frameshift_qual   | Minumum Clair3 variant quality to keep a frameshift variant                                              | Integer | 15                          | Frameshift defined as not divisible by 3                                                                                                         |
| --min_allele_freq       | Minimum allele frequency to call a variant                                                               | Number  | 0.65                        |                                                                                                                                                  |
| --medaka_model          | Medaka model to be used in the pipeline                                                                  | String  | 'r941_min_hac_g507'         | Default model will not work the best for all inputs. [See medaka docs](https://github.com/nanoporetech/medaka#models) for additional info        |
| --skip_longshot         | Skip running longshot in medaka workflow                                                                 | Boolean | false                       | Medaka only                                                                                                                                      |
| --fast5_pass            | Path to directory containing `barcode##` fast5 subdirectories                                            | Path    | null                        | Only for nanopolish                                                                                                                              |
| --sequencing_summary    | Path to run `sequencing_summary*.txt` file                                                               | Path    | null                        | Only for nanopolish                                                                                                                              |
| --min_length            | Minimum read length to be kept                                                                           | Integer | 200                         | For artic guppyplex                                                                                                                              |
| --max_length            | Maximum read length to be kept                                                                           | Integer | 3000                        | For artic guppyplex                                                                                                                              |
| --min_reads             | Minimum size selected reads to be used in pipeline                                                       | Integer | 20                          |                                                                                                                                                  |
| --metadata              | Path to metadata TSV file with columns 'sample' and 'barcode'                                            | Path    | null                        | See [metadata](#metadata) for more info                                                                                                          |
| --use_artic_tool        | Run the artic tool itself instead of nextflow implementation                                             | Bool    | False                       | Not available with clair3                                                                                                                        |
| --normalise             | Artic minion normalise coverage option                                                                   | Integer | 1000                        | Entering `0` turns off normalisation. Only for amplicon sequencing                                                                               |
| --no_frameshift         | Filter INDEL variants that are not divisible by 3                                                        | Boolean | False                       | Simple `%3 == 0` check for variants                                                                                                              |
| --skip_snpeff           | Skip running SnpEff                                                                                      | Boolean | False                       |                                                                                                                                                  |
| --gff                   | Path to gff3 formatted file to use in SnpEff database build                                              | Path    | False                       | Not required to run [SnpEff](#snpeff). See below for details                                                                                     |
| --skip_qc               | Skip running all QC and reporting steps                                                                  | Boolean | false                       |                                                                                                                                                  |
| --multiqc_report        | Run MultiQC report over custom report report                                                             | Boolean | false                       |                                                                                                                                                  |
| --pcr_primer_bed        | Path to PCR primer bed file to check for mutations against                                               | Path    | null                        | For output QC checks                                                                                                                             |
| --neg_control_threshold | Coverage threshold at which to fail negative control samples                                             | Number  | 0.10                        |                                                                                                                                                  |
| --neg_ctrl_substrings   | Negative control sample substrings separated by a `,`                                                    | String  | 'ntc,neg,blank,water'       |                                                                                                                                                  |
| --outdir                | Directory name to output results to                                                                      | String  | null                        | Required                                                                                                                                         |

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
| ...        | ...   | ...  | ...               | ...         | ...       | CCCTAGAAA  |
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

## Core Nextflow Arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.
