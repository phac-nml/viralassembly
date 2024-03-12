# viralassembly
A generic viral assembly and QC pipeline which utilises a re-implementation of the [artic pipeline](https://github.com/artic-network/fieldbioinformatics/tree/master/artic) to separate out the individual steps allowing greater control on tool versions along with how data is run through the processes. This pipeline can be used as a starting point for analyses on viruses without dedicated workflows already available.

Some of the goals of this pipeline are:
1. Rework the artic nanopore pipeline steps as nextflow modules to deal with specific bugs and version incompatibilities
    - Example: BCFtools consensus error seen in artic pipeline sometimes
    - Allows adding in clair3 as a new variant calling tool
    - Potentially eventually work to remove `artic` as a dependency
2. Allow the pipeline to be used on other viruses with or without amplicon schemes
    - Due to the QC steps there is unfortunately a current limitation at working with segmented viruses
        - The pipeline will automatically exit after assembly and not generate QC and Reports for these at this time
        - This will be fully implemented at some point
3. Provide `Run` level and `Sample` level final reports

## Index
- [Installation](#installation)
- [Profiles](#profiles)
- [Base Run Commands](#running-commands)
    - [Nanopore - Nanopolish](#nanopore---nanopolish)
    - [Nanopore - Medaka](#nanopore---medaka)
    - [Nanopore - Clair3](#nanopore---clair3)
- [Inputs](#inputs)
    - [Input Data Options](#input-data-options)
    - [Nanopolish Pipeline Required Parameters](#nanopolish-pipeline-required-parameters)
    - [Medaka Pipeline Required Parameters](#medaka-pipeline-required-parameters)
    - [Clair3 Pipeline Required Parameters](#clair3-pipeline-required-parameters)
    - [Optional Inputs](#optional-inputs)
    - [Metadata TSV](#metadata)
    - [Primer Schemes](#schemes)
    - [SnpEff](#snpeff-1)
- [Outputs](#outputs)
- [Limitations](#limitations)
- [Citations](#citations)
- [Contributing](#contributing)
- [Legal](#legal)

## Installation
1. Download and install nextflow
    1. Download and install with [conda](https://docs.conda.io/en/latest/miniconda.html)
        - Conda command: `conda create on nextflow -c conda-forge -c bioconda nextflow`
    2. Install with the instructions at https://www.nextflow.io/

2. Run the pipeline with a profile to handle dependencies:
    - `conda`
    - `mamba`
    - `singularity`
    - `docker`

## Profiles
Profiles are used to specify dependency installation, resources, and how to handle pipeline jobs. You can specify more than one profile but avoid passing in more than one dependency managment profiles.

Available:
- `conda`: Utilize conda to install dependencies and environment management
- `mamba`: Utilize mamba to install dependencies and environment management
- `singularity`: Utilize singularity for dependencies and environment management
- `docker`: Utilize docker to for dependencies and environment management

Custom configs can be created and specified to be used as an additional config file with `-c <CONFIG>`. More info [available here](https://www.nextflow.io/docs/latest/config.html) or [here](https://training.nextflow.io/basic_training/config/#config-syntax)

## Running Commands
Simple commands to run input data. Input data can be done in three different ways:
1. Passing `--fastq_pass </PATH/TO/fastq_pass>` where `fastq_pass` is a directory containing `barcode##` subdirectories with fastq files
2. Passing `--fastq_pass </PATH/TO/fastqs>` where `fastqs` is a directory containing `.fastq*` files
3. Passing `--input <samplesheet.csv>` where `samplesheet.csv` is a CSV file with two columns
    1. `sample` - The name of the sample
    2. `reads`  - Path to a directory containing reads for the sample in `.fastq*` format

The basic examples will show how to run the pipeline using the `--fastq_pass` input but it could be subbed in for the `--input` CSV file if wanted

### Nanopore - Nanopolish
Running the pipeline with [nanopolish](https://github.com/jts/nanopolish) for variant calls requires fastq files, fast5 files, and the sequencing summary file. When running, the pipeline will look for subdirectories off of the input directory called `barcode##` to be used in the pipeline.

Basic command:
```bash
nextflow run /PATH/TO/artic-generic-nf/main.nf \
    -profile <PROFILE(s)> \
    --variant_caller 'nanopolish' \
    --fastq_pass </PATH/TO/fastq_pass> \
    --fast5_pass </PATH/TO/fast5_pass> \
    --sequencing_summary </PATH/TO/sequencing_summary.txt>
    <OPTIONAL INPUTS>
```

[Optional inputs](#inputs) could include:
- Different schemes
- Metadata
- Filtering options
- Using base `artic minion` instead of nextflow implementation
- Running SnpEff for variant consequence prediction
- Output reporting options

### Nanopore - Medaka
Running the pipeline with [medaka](https://github.com/nanoporetech/medaka) for variant calls requires fastq files and a medaka model. When running, the pipeline will either:
- Look for subdirectories off of the input "--fastq_pass" directory called `barcode##` to be used in the pipeline
- Look for fastq files in the input "--fastq_pass" directory called `*.fastq*` to be used in the pipeline

Basic command:
```bash
nextflow run /PATH/TO/artic-generic-nf/main.nf \
    -profile <PROFILE(s)> \
    --variant_caller 'medaka' \
    --fastq_pass </PATH/TO/fastq_pass> \
    --medaka_model <Medaka Model>
    <OPTIONAL INPUTS>
```

[Optional inputs](#optional-inputs) could include:
- Different schemes
- Metadata
- Filtering options
- Using base `artic minion` instead of nextflow implementation
- Running SnpEff for variant consequence prediction
- Output reporting options

Medaka model information [can be found here](https://github.com/nanoporetech/medaka#models)

### Nanopore - Clair3
Running the pipeline with [Clair3](https://github.com/HKU-BAL/Clair3) for variant calls requires fastq files and a clair3 model. When running, the pipeline will either:
- Look for subdirectories off of the input "--fastq_pass" directory called `barcode##` to be used in the pipeline
- Look for fastq files in the input "--fastq_pass" directory called `*.fastq*` to be used in the pipeline

This pipeline utilizes the same steps as the artic fieldbioinformatics minion pipeline but with each step run using nextflow to allow clair3 to be easily slotted in.

Basic command:
```bash
nextflow run /PATH/TO/artic-generic-nf/main.nf \
    -profile <PROFILE(s)> \
    --variant_caller 'clair3' \
    --fastq_pass </PATH/TO/fastq_pass> \
    --clair3_model <Clair3 Model or /PATH/TO/clair3_model>
    <OPTIONAL INPUTS>
```

[Optional inputs](#optional-inputs) could include:
- Different schemes
- Metadata
- Filtering options
- Running SnpEff for variant consequence prediction
- Output reporting options

## Inputs
The required inputs are based off of what variant caller is being used. The caller is specified with `--variant_caller [nanopolish, medaka, clair3]` and then the remaining required parameters are based on that although the only real difference is that nanopolish requires additional files.

Most of the optional parameters are available for all three pipelines

### Input Data Options
There are two ways to specify how fastq files are input into the pipeline:

1. `--fastq_dir <PATH/TO/FASTQs>`
    - Using this you can either provide the path to the barcode directories found when demultiplexing nanopore data or to a flat directory of named fastq files
    - If using a barcoded directory, you can also rename the fastq barcode by providing a `--metadata metadata.tsv` file

2. `--input <input.csv>`
    - Pass in an input CSV file containing 2 columns, `sample`, and `reads` where:
        - `sample` is the sample name to use
        - `reads` is the path to the barcode directory containing reads (for now it is just directory, will update later to work on flat fastqs too)
    - This will find the reads in the given directory and attach the sample name to them

The below examples all use method 1 but swapping out `--fastq_dir <FASTQs>` for `--input <input.csv>` would also work for all of them

### Nanopolish Pipeline Required Parameters

| Parameter | Description | Type | Default | Notes |
| - | - | - | - | - |
| --variant_caller 'nanopolish' | Set nanopolish for variant calls | Choice |  |  |
| --fastq_pass | Path to directory containing `barcode##` subdirectories | Path | null |  |
| --fast5_pass | Path to directory containing `barcode##` fast5 subdirectories | Path | null |  |
| --sequencing_summary | Path to run `sequencing_summary*.txt` file | Path | null |  |

### Medaka Pipeline Required Parameters

| Parameter | Description | Type | Default | Notes |
| - | - | - | - | - |
| --variant_caller 'medaka' | Set medaka for variant calls | Choice |  |  |
| --fastq_pass | Path to directory containing `barcode##` subdirectories or a directory containing `*.fastq*` files | Path | null |  |
| --medaka_model | Medaka model to be used in the pipeline | Str | 'r941_min_hac_g507' | Default model will not work the best for all inputs. [See medaka docs](https://github.com/nanoporetech/medaka#models) for additional info |

### Clair3 Pipeline Required Parameters

| Parameter | Description | Type | Default | Notes |
| - | - | - | - | - |
| --variant_caller 'clair3' | Set clair3 for variant calls | Choice |  |  |
| --fastq_pass | Path to directory containing `barcode##` subdirectories or a directory containing `*.fastq*` files | Path | null |  |
| --clair3_model | Clair3 model to be used in the pipeline | Str/Path | 'r941_prom_sup_g5014' | Default model will not work the best for all inputs. [See clair3 docs](https://github.com/HKU-BAL/Clair3#pre-trained-models) for additional info |

### Optional Inputs
Use `--help` to see all options formatted on the command line

Use `--version` to see version information

#### Read Filtering
| Parameter | Description | Type |  Default | Notes |
| - | - | - | - | - |
| --min_length | Minimum read length to be kept | Int | 200 | For artic guppyplex |
| --max_length | Maximum read length to be kept | Int | 3000 | For artic guppyplex |
| --min_reads | Minimum size selected reads to be used in pipeline | Int | 20 |  |

#### Scheme
| Parameter | Description | Type |  Default | Notes |
| - | - | - | - | - |
| --reference_no_scheme | Specify the path to a reference fasta file to run pipeline without a primer scheme | Path | '' | Ignores all scheme inputs |
| --scheme | Name of the primer scheme to use | Str | 'nCoV-2019' | See [schemes](#schemes) for more info |
| --scheme_version | Version name of primer scheme to use | Str | 'freed_V2_nml'  | See [schemes](#schemes) for more info |
| --scheme_repo | Github repository URL to download scheme from | Str | 'https://github.com/DarianHole/primer-schemes.git' | See [schemes](#schemes) for more info |
| --local_scheme | Path to directory containing local scheme files | Path | null | See [schemes](#schemes) for more info |

#### Pipeline Options
| Parameter | Description | Type |  Default | Notes |
| - | - | - | - | - |
| --metadata | Path to metadata TSV file with columns 'sample' and 'barcode' | Path | null | See [metadata](#metadata) for more info |
| --use_artic_tool | Run the artic tool itself instead of nextflow implementation | Bool | False | Not available with clair3 |
| --normalise | Artic minion normalise coverage option | Int | 1000 | Entering `0` turns off normalisation |
| --no_frameshift | Use the Artic minion no frameshift vcf filter | Bool | False | Simple `%3 == 0` check for variants |
| --use_bwa | Use BWA instead of minimap2 for read mapping | Bool | False |  |
| --skip_longshot | When running with `medaka`, skip running longshot | Bool | False | Medaka only!! |

#### SnpEff
| Parameter | Description | Type |  Default | Notes |
| - | - | - | - | - |
| --skip_snpeff | Skip running SnpEff | Bool | False |  |
| --gff | Path to gff3 formatted file to use in SnpEff database build | Path | False | If not given, the Reference ID will be used to attempt to pull the genbank file |

#### QC and Reporting
| Parameter | Description | Type |  Default | Notes |
| - | - | - | - | - |
| --skip_qc | Skip running all QC and reporting steps | Bool | false |  |
| --custom_report | Run the custom HTML report | Bool | false | Currently requires the use of conda as there is not a singularity container yet |
| --pcr_primer_bed | Path to PCR primer bed file to check for mutations against | Path | null | For output QC checks |
| --neg_control_threshold | Coverage threshold at which to fail negative control samples | Float | 0.10 |  |
| --neg_ctrl_substrings | Negative control sample substrings separated by a `,` | Str | 'ntc,neg,blank' |  |

#### Other Generic Options
| Parameter | Description | Type |  Default | Notes |
| - | - | - | - | - |
| --outdir | Directory name to output results to | Str | 'results' |  |
| --cache | Specify a location to store conda/singularity envs/containers for reuse | Path | null |  |

### Metadata
Input metadata is used to rename barcoded fastq files along with adding additional lines to the final overall QC csv file. Note that the metadata input is expected to be of a `TSV` format

Structure for example `metadata.tsv` file:

| sample | barcode | \<Anything else you want to add > |
| - | - | - |
| SR-1 | 1 | X |
| SR-2 | 02 | Y |
| NTC-12 | 12 | Z |

### Schemes
Amplicon schemes are a highly targeted approach to sequencing focusing on a specific target genome. If using an amplicon scheme with this pipeline, either a local directory or a URL that contains the wanted primer scheme formatted according to the below information must be provided.

If not running with an amplicon scheme, pass the `--reference_no_scheme <PATH/TO/reference.fasta>` argument with a reference fasta file and the pipeline will run without amplicon specific checks/outputs.

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

Example for covid:
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

### SnpEff
SnpEff is run by default on all non-segmented viruses by using the reference sequence ID to either:
1. Check if there is a SnpEff database available to download
2. Build a SnpEff database by downloading the sequence genbank file from NCBI

Instead of relying on the reference ID to build/download a database, you can instead specify a gff3 file with `--gff <PATH/TO/file.gff>` to be used with the reference sequence to create a database

If building/downloading a database fails, the pipeline will skip over running SnpEff instead of failing out completely.

SnpEff can also be skipped entirely by passing the `--skip_snpeff` parameter

## Outputs
Outputs are separated based off of their tool or file format in the `results/` directory by default.

### Consensus Sequences
Final consensus sequence output by the pipeline found in the `consensus/` directory as: `SAMPLE.consensus.fasta`

### Bam Files
BAM files in the `bam/` directory include:
- `SAMPLE.primertrimmed.rg.sorted.bam`
    - If using a scheme, these are the bam files used for variant calling and downstream QC steps
- `SAMPLE.sorted.bam`
    - If not using a scheme, these are the bam files used for all further steps
        - Still output if using a scheme though

### VCF Files
VCF files found in the `vcf/` directory include:
- `SAMPLE.pass.vcf.gz`
    - Variants passing variant filtering
- `SAMPLE.pass.norm.vcf.gz`
    - BCFtools normalized variants passing variant filtering

### SnpEff Outputs
Output files from SnpEff found in the `snpeff/` directory including:
- `SAMPLE.ann.vcf`
    - SnpEff annotated VCF file
- `SAMPLE.csv`
    - SnpEff CSV summary

### Variation Files
Positional variation information based on the pileups and formatted as a CSV file is available in the `variation_csvs/` directory

### Sample Summary Files
- `sample_csvs/`
    - Directory containing the individual csv file for each sample
- `sample_mqc/`
    - Directory containing MultiQC HTML reports for each sample

### Overall outputs
- `overall.qc.csv`
    - Contains all of the calculated per sample metadata from the run along with giving a qc and negative control status
- `Overall-Run-MultiQC-Report_multiqc_report.html`
    - MultiQC HTML report containing all sample if running multiQC output
- `reportDashboard.html`
    - Custom HTML report containing all samples if running the custom output
    - Note that it may take a bit to load if running a lot of samples
    ![custom-report](./pictures/custom_report.png)
- `pipeline_info/`
    - Directory containing nextflow run information

## Limitations
Current limitations include:

1. Nanopore data only at this time
2. Currently runs for viruses using a reference genome
    - Segmented viruses will exit before the QC section for now
3. Custom report can only work when running with `conda`
4. SnpEff issues in running and database building/downloading
    - Database building/downloading requires one of three things:
        - The reference ID is in the SnpEff database
            - This allows the database to be downloaded
        - A gff3 file
            - This is used with the reference sequence to build a database
        - A well annotated NCBI genome matching the reference ID 
            - This will pull the genbank file and use that to build a database
    - Running SnpEff with singularity sometimes leads to a lock issue which is hopefully fixed

## Citations
This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
 
> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
> In addition, references of tools and data used in this pipeline are as follows:

Detailed citations for utilized tools are found in [citations.md](./citations.md)

## Contributing
Contributions are welcome through creating PRs or Issues

## Legal
Copyright 2023 Government of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
