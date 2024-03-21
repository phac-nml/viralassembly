# phac-nml/viralassembly: Usage

## Introduction
This pipeline is intended to be run on either Nanopore Amplicon Sequencing data or Basic Nanopore NGS Sequencing data that can utilize a reference genome for mapping variant calling, and other downstream analyses. It generates variant calls, consensus sequences, and quality control information based on the reference. To do this, there are three different variant callers that can be utilized which includes: `clair3`, `medaka`, and `nanopolish` (For R9.4.1 flowcells and below only).

For Amplicon Sequencing data it is at minimum required to:
1. Specify a path to the reads/input file
2. Specify the scheme name
3. Specify the scheme version
4. Pick a variant caller and caller model

For Basic NGS sequencing data it is required to:
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
Profiles are used to specify dependency installation, resources, and how to handle pipeline jobs. You can specify more than one profile but avoid passing in more than one dependency managment profiles. They can be passed with `-profile <PROFILE>`

Available:
- `conda`: Utilize conda to install dependencies and environment management
- `mamba`: Utilize mamba to install dependencies and environment management
- `singularity`: Utilize singularity for dependencies and environment management
- `docker`: Utilize docker to for dependencies and environment management

## Data Inputs
Two options for fastq data input: `--fastq_pass <FASTQ_PASS/>` or `--input <INPUT.csv>`

### Fastq Pass Directory (--fastq_pass)
Specify fastq data to input based on a given directory. The directory can either be barcoded, as would be seen after demultiplexing, or it could be a flat input of fastq files. The barcoded fastq data will be output with the barcode number but can be renamed with a [metadata csv]() input. The flat fastq files will keep their basename (separated out at the first `.`).

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
You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to pass in an input CSV file containing 2 columns, `sample`, and `reads` where:
- `sample` is the sample name to use
- `reads` is the path to the barcode directory containing reads

Ex.
| sample | reads |
| - | - |
| sample1 | /path/to/barcode01 |
| ntc | /path/to/barcode02 |
| pos | /path/to/barcode03 |

This will be expanded upon in future releases to allow more varied inputs for the input sheet.

## Variant Callers
Three different variant callers are available with slightly different options regarding running with them. For the most accurate results when running with `clair3` or `medaka` pick a model that best matches the input data!!

### [Clair3](https://github.com/HKU-BAL/Clair3)
Clair3 is a germline small variant caller for long-reads.

Running with `clair3` requires the following parameters:
- `--variant_caller clair3`: Sets clair3 as the variant caller

And has the optional parameters of:
- `--clair3_model <MODEL>`: Specify the base clair3 model
- `--clair3_user_variant_model </PATH/TO/downloaded_clair3_model>`:  Specify the path to an additionally downloaded model directory
- `clair3_no_pool_split`: Do not split inputs into pools

Clair3 comes with some models available and is defaulted to `r941_prom_sup_g5014`. Additional models can be downloaded from [ONT Rerio](https://github.com/nanoporetech/rerio/tree/master) and then specified in the `--clair3_user_variant_model </PATH/TO/downloaded_clair3_model>` parameter shown above. Remember to pick a model that best represents the data!

### [Medaka](https://github.com/nanoporetech/medaka)
Medaka is a tool to create consensus sequences and variant calls from nanopore sequencing data using neural networks and provied by ONT.

Running with `medaka` requires the following parameters:
- `--variant_caller medaka`: Sets medaka as the variant caller

And has the optional parameters of:
`--medaka_model <MODEL>`: Specify the wanted medaka model

Medaka models come built in with the tool itself with the default set to `r941_min_hac_g507` which can be changed with `--medaka_model <MODEL>` parameter shown above. More information on models [can be found here](https://github.com/nanoporetech/medaka#models). Remember to pick a model that best represents the data!

### [Nanopolish](https://github.com/jts/nanopolish)
Nanopolish is a software package for signal-level analysis of Oxford Nanopore sequencing data. It does not presently support the R10.4 flowcells so as a variant caller it should only be used with R9.4 flowcells. It also requires that the fastq data is in barcoded directories to work correctly.

Running with `nanopolish` requires the following parameters:
- `--variant_caller nanopolish`
- `--fast5_pass <FAST5_PASS/>`
- `--sequencing_summary <SEQ_SUM.txt>`

Nanopolish requires the fast5 directory along with the sequencing summary file to be used as input instead of a model.

## Running the pipeline

### Amplicon
The typical command for running the pipeline with an amplicon scheme using medaka and a different medaka model is as follows:

```bash
nextflow run phac-nml/viralassembly \
  -profile docker \
  --fastq_pass FASTQ_PASS/ \
  --variant_caller medaka \
  --medaka_model 'r1041_e82_400bps_sup_v4.3.0' \
  --scheme 'nCoV-2019' \
  --scheme_version 'V5.3.2' \
  --outdir ./results
```

This will launch the pipeline with the `docker` configuration profile, the `medaka` variant caller, and the `nCoV-2019` version `V5.3.2` primer scheme from https://github.com/artic-network/primer-schemes/tree/master/nCoV-2019 (default scheme repo to pull). Profile information [can be found above](#profiles)

### Non-Amplicon
The typical command for running the pipeline without an amplicon scheme using medaka and a different medaka model is as follows:

```bash
nextflow run phac-nml/viralassembly \
  -profile singularity \
  --fastq_pass FASTQ_PASS/ \
  --variant_caller medaka \
  --medaka_model 'r1041_e82_400bps_sup_v4.3.0' \
  --reference REF.fa \
  --outdir ./results
```

This will launch the pipeline with the `singularity` configuration profile, the `medaka` variant caller, and the specified reference. Profile information [can be found above](#profiles)

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
fastq_pass: './fastq_pass'
variant_caller: 'medaka'
medaka_model: 'r1041_e82_400bps_sup_v4.3.0'
reference: 'reference.fa'
outdir: './results/'
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
| Parameter | Description | Type | Default | Notes |
| - | - | - | - | - |
| --variant_caller | Pick from the 3 variant callers: 'clair3', 'medaka', 'nanopolish' | Choice | '' | Details above |
| --clair3_model | Clair3 base model to be used in the pipeline | Str | 'r941_prom_sup_g5014' | Default model will not work the best for all inputs. [See clair3 docs](https://github.com/HKU-BAL/Clair3#pre-trained-models) for additional info |
| --clair3_user_variant_model | Path to clair3 additional model directory to use instead of a base model | Path | '' | Default model will not work the best for all inputs. [See clair3 docs](https://github.com/HKU-BAL/Clair3#pre-trained-models) for additional info |
| --clair3_no_pool_split | Do not split reads into separate pools | Bool | False | Clair3 amplicon sequencing only |
| --medaka_model | Medaka model to be used in the pipeline | Str | 'r941_min_hac_g507' | Default model will not work the best for all inputs. [See medaka docs](https://github.com/nanoporetech/medaka#models) for additional info |
| --fastq_pass | Path to directory containing `barcode##` subdirectories | Path | null |  |
| --fast5_pass | Path to directory containing `barcode##` fast5 subdirectories | Path | null | Only for nanopolish |
| --sequencing_summary | Path to run `sequencing_summary*.txt` file | Path | null | Only for nanopolish |
| --min_length | Minimum read length to be kept | Int | 200 | For artic guppyplex |
| --max_length | Maximum read length to be kept | Int | 3000 | For artic guppyplex |
| --min_reads | Minimum size selected reads to be used in pipeline | Int | 20 |  |
| --reference | Specify the path to a reference fasta file to run pipeline without a primer scheme | Path | '' | Ignores all scheme inputs. See [schemes and reference](#schemes-and-reference) |
| --scheme | Name of the primer scheme to use | Str | '' | See [schemes and reference](#schemes-and-reference) |
| --scheme_version | Version name of primer scheme to use | Str | ''  | See [schemes and reference](#schemes-and-reference) |
| --scheme_repo | Github repository URL to download scheme from | Str | 'https://github.com/artic-network/primer-schemes.git' | See [schemes and reference](#schemes-and-reference) |
| --local_scheme | Path to directory containing local scheme files | Path | null | See [schemes and reference](#schemes-and-reference) |
| --metadata | Path to metadata TSV file with columns 'sample' and 'barcode' | Path | null | See [metadata](#metadata) for more info |
| --use_artic_tool | Run the artic tool itself instead of nextflow implementation | Bool | False | Not available with clair3 |
| --normalise | Artic minion normalise coverage option | Int | 1000 | Entering `0` turns off normalisation. Only for amplicon sequencing |
| --no_frameshift | Use the Artic minion no frameshift vcf filter | Bool | False | Simple `%3 == 0` check for variants |
| --use_bwa | Use BWA instead of minimap2 for read mapping | Bool | False |  |
| --skip_longshot | When running with `medaka`, skip running longshot | Bool | False | Medaka only!! |
| --skip_snpeff | Skip running SnpEff | Bool | False |  |
| --gff | Path to gff3 formatted file to use in SnpEff database build | Path | False | Not required to run [SnpEff](#snpeff). See below for details |
| --skip_qc | Skip running all QC and reporting steps | Bool | false |  |
| --custom_report | Run the custom HTML report | Bool | false | Currently requires the use of conda as there is not a singularity container yet |
| --pcr_primer_bed | Path to PCR primer bed file to check for mutations against | Path | null | For output QC checks |
| --neg_control_threshold | Coverage threshold at which to fail negative control samples | Float | 0.10 |  |
| --neg_ctrl_substrings | Negative control sample substrings separated by a `,` | Str | 'ntc,neg,blank' |  |
| --outdir | Directory name to output results to | Str | 'results' |  |
| --cache | Specify a location to store conda/singularity envs/containers for reuse | Path | null |  |

### Schemes and Reference
Amplicon schemes are a highly targeted approach to sequencing focusing on a specific target genome. If using an amplicon scheme with this pipeline, either a local directory or a URL that contains the wanted primer scheme formatted according to the below information must be provided.

If not running with an amplicon scheme, pass the `--reference <PATH/TO/reference.fasta>` argument with a reference fasta file and the pipeline will run without amplicon specific checks/outputs.

The primer scheme must contain:
- A reference genome fasta sequence titled `*reference.fasta`
- A primer bed file titled `*primer.bed`
    - Minimum of 6 columns
    - Primer pairs with names containing `_LEFT` and `_RIGHT`
    - Primer pools

Example Primer file:

| MN908947.3 | 30 | 54 | nCoV-2019_1_LEFT | 1 | + |
| - | - | - | - | - | - |
| MN908947.3 | 1183 | 1205 | nCoV-2019_1_RIGHT | 1 | - |
| MN908947.3 | 1100 | 1128 | nCoV-2019_2_LEFT | 2 | + |
| MN908947.3 | 2244 | 2266 | nCoV-2019_2_RIGHT | 2 | - |
| ... | ... | ... | ... | ... | ... |
| REF ID | Start | Stop | Primer Name | Primer Pool | Direction


The directory structure must follow the basic structure as follows:
```
primer-schemes
└── <SCHEME>
    └── <SCHEME VERSION>
        ├── reference.fasta
        └── scheme.bed
```

Example for Sars-CoV2:
```
primer-schemes
└── nCoV-2019
    ├── midnight
    |   ├── nCoV-2019.reference.fasta
    |   └── nCoV-2019.scheme.bed
    └── V1
        ├── reference.fasta
        └── scheme.bed
```

### Metadata
Input metadata is used to rename barcoded fastq files along with adding additional lines to the final overall QC csv file. Note that the metadata input is expected to be of a `TSV` format

Structure for example `metadata.tsv` file:

| sample | barcode | \<Anything else you want to add > |
| - | - | - |
| SR-1 | 1 | X |
| SR-2 | 02 | Y |
| NTC-12 | 12 | Z |

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

