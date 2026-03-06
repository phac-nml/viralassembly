# viralassembly

A generic viral assembly and QC pipeline which utilises a re-implementation of the [ARTIC pipeline](https://github.com/artic-network/fieldbioinformatics/tree/master/artic) to separate out the individual steps allowing greater control on tool versions along with how data is run through the processes. This pipeline can be used as a starting point for analyses on viruses without dedicated workflows already available.

This pipeline is intended to be run on either Nanopore Amplicon Sequencing data or Basic Nanopore NGS Sequencing data that can utilize a reference genome for read mapping, variant calling, and other downstream analyses. It generates variant calls, consensus sequences, and quality control information based on the reference. To do this, there are three different variant callers that can be utilized which includes: `clair3`, `medaka`, and `nanopolish` (For R9.4.1 flowcells and below only). By default, `clair3` is used and highly recommended over the other two (as of December 2025). Their inclusion being for completeness with eventual removal from the workflow.

The goals of this pipeline are:

1. Provide a generic viral pipeline for the NML Surviellance Platform IRIDA-Next
2. Provide detailed and useful `Run` and `Sample` level final reports
3. Allow the pipeline to be used on other viruses with or without amplicon schemes

## Index

- [Installation](#installation)
- [Base Run Commands](#running-commands)
  - [Nanopore - Clair3](#nanopore---clair3---default)
  - [Nanopore - Medaka](#nanopore---medaka)
  - [Nanopore - Nanopolish](#nanopore---nanopolish)
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

2. Run the pipeline with one of the following profiles to handle dependencies (or use your own ``-profile` if you have one!):
   - `conda`
   - `mamba`
   - `singularity`
   - `docker`

## Running Commands

Simple commands to run input data. Data can be ingested by the pipeline in three different ways:

1. Passing `--fastq_pass </PATH/TO/fastq_pass>` where `fastq_pass` is a directory containing `barcode##` subdirectories with fastq files
   - This is the format that matches data which has been basecalled with `dorado`
   - Samples can be renamed by passing in the [`--metadata` parameter](./docs/usage.md#metadata) with a TSV file mapping the barcode to the sample name
2. Passing `--fastq_pass </PATH/TO/fastqs>` where `fastqs` is a directory containing named `*.fastq*` files
   - Sample names are based off of the file names
3. Passing `--input <samplesheet.csv>` where `samplesheet.csv` is a CSV file with two columns
   1. `sample` - The name of the sample
   2. `fastq_1` - Path to one fastq file per sample in `.fastq*` format

The basic examples will show how to run the pipeline using the `--fastq_pass` input but it could be subbed in for the `--input` CSV file if wanted.

_All detailed running information is available in the [usage docs](./docs/usage.md)_

### Nanopore - Clair3 - Default

Running the pipeline with [Clair3](https://github.com/HKU-BAL/Clair3) (default and recommended variant caller) requires fastq files and a Clair3 model. When running with `--fastq_pass`, the pipeline will either:

- Look for subdirectories off of the input "--fastq_pass" directory called `barcode##` to be used in the pipeline
- Look for fastq files in the input "--fastq_pass" directory called `*.fastq*` to be used in the pipeline

This pipeline utilizes the same steps as the artic fieldbioinformatics minion pipeline but with each step run using nextflow to allow Clair3 to be easily slotted in. See the [Clair3 Section](./docs/usage.md#clair3) of the usage docs for more information

Basic command:

```bash
nextflow run /PATH/TO/artic-generic-nf/main.nf \
    -profile <PROFILE(s)> \
    --fastq_pass </PATH/TO/fastq_pass> \
    --reference <REF.fa> \
    --clair3_model <MODEL> \
    <OPTIONAL INPUTS>
```

[Optional inputs](./docs/usage.md#all-parameters) could include:

- [Amplicon scheme](./docs/usage.md#schemes-and-reference) instead of just a reference fasta file
- Metadata
- Filtering options
- Running SnpEff for variant consequence prediction
- Output reporting options

The pipeline could also be run with `--input CSV` to pass in an input CSV file with the sample names and fastq pathes

### Nanopore - Medaka

Running the pipeline with [medaka](https://github.com/nanoporetech/medaka) for variant calling requires fastq files and a medaka model. When running the pipeline will either:

- Look for subdirectories off of the input "--fastq_pass" directory called `barcode##` to be used in the pipeline
- Look for fastq files in the input "--fastq_pass" directory called `*.fastq*` to be used in the pipeline

See the [medaka section](./docs/usage.md#medaka) of the usage docs for more information

Basic command:

```bash
nextflow run /PATH/TO/artic-generic-nf/main.nf \
    -profile <PROFILE(s)> \
    --variant_caller 'medaka' \
    --fastq_pass </PATH/TO/fastq_pass> \
    --medaka_model <Medaka Model> \
    --reference <REF.fa> \
    <OPTIONAL INPUTS>
```

[Optional inputs](./docs/usage.md#all-parameters) could include:

- [Amplicon scheme](./docs/usage.md#schemes-and-reference) instead of just a reference fasta file
- Metadata
- Filtering options
- Running SnpEff for variant consequence prediction
- Output reporting options

Medaka model information [can be found here](https://github.com/nanoporetech/medaka#models)

### Nanopore - Nanopolish

Running the pipeline with [nanopolish](https://github.com/jts/nanopolish) for variant calls requires fastq files, fast5 files, and the sequencing summary file instead of providing a model. As such, nanopolish requires that the read ids in the fastq files are linked by the sequencing summary file to their signal-level data in the fast5 files. This makes it **a lot** easier to run using barcoded directories but it can be run with individual read files

See the [nanopolish section](./docs/usage.md#nanopolish) of the usage docs for more information

Basic command:

```bash
nextflow run /PATH/TO/artic-generic-nf/main.nf \
    -profile <PROFILE(s)> \
    --variant_caller 'nanopolish' \
    --fastq_pass </PATH/TO/fastq_pass> \
    --fast5_pass </PATH/TO/fast5_pass> \
    --sequencing_summary </PATH/TO/sequencing_summary.txt> \
    --reference <REF.fa>
    <OPTIONAL INPUTS>
```

[Optional inputs](./docs/usage.md#all-parameters) could include:

- [Amplicon scheme](./docs/usage.md#schemes-and-reference) instead of just a reference fasta file
- Metadata
- Filtering options
- Running SnpEff for variant consequence prediction
- Output reporting options

## Outputs

Outputs are separated based off of their tool or file format and found in the `results/` directory by default.

Outputs include:

- Consensus fasta files
- VCF files
- Bam files
- HTML summary files (either custom or MultiQC)

_More output information on pipeline steps and output files can be found in the [output docs](./docs/output.md)_

## Limitations

Current limitations include:

1. Nanopore data only at this time
2. Currently runs for viruses using a reference genome
   - Segmented viruses will exit before the QC section for now while looking into how to best report them
3. SnpEff and database building/downloading can be finicky
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

Licensed under the MIT License (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

https://opensource.org/license/mit/

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
