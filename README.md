# viralassembly
A generic viral assembly and QC pipeline which utilises a re-implementation of the [artic pipeline](https://github.com/artic-network/fieldbioinformatics/tree/master/artic) to separate out the individual steps allowing greater control on tool versions along with how data is run through the processes. This pipeline can be used as a starting point for analyses on viruses without dedicated workflows already available.

This pipeline is currently intended to be run on either Nanopore Amplicon Sequencing data or Basic Nanopore NGS Sequencing data that can utilize a reference genome for mapping. It generates variant calls, consensus sequences, and quality control information based on the reference. To do this, there are three different variant callers that can be utilized which includes: `clair3`, `medaka`, and `nanopolish` (For R9.4.1 flowcells and below only),

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
- [Base Run Commands](#running-commands)
    - [Nanopore - Clair3](#nanopore---clair3)
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

2. Run the pipeline with a profile to handle dependencies:
    - `conda`
    - `mamba`
    - `singularity`
    - `docker`

## Running Commands
Simple commands to run input data. Input data can be done in three different ways:
1. Passing `--fastq_pass </PATH/TO/fastq_pass>` where `fastq_pass` is a directory containing `barcode##` subdirectories with fastq files
2. Passing `--fastq_pass </PATH/TO/fastqs>` where `fastqs` is a directory containing `.fastq*` files
3. Passing `--input <samplesheet.csv>` where `samplesheet.csv` is a CSV file with two columns
    1. `sample` - The name of the sample
    2. `reads`  - Path to a directory containing reads for the sample in `.fastq*` format

The basic examples will show how to run the pipeline using the `--fastq_pass` input but it could be subbed in for the `--input` CSV file if wanted.

*All detailed running information is available in the [usage docs](./docs/usage.md)*

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
    --clair3_model <Clair3 Model or /PATH/TO/clair3_model> \
    --reference <REF.fa>
    <OPTIONAL INPUTS>
```

[Optional inputs](./docs/usage.md#all-parameters) could include:
- Different schemes
- Metadata
- Filtering options
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

[Optional inputs](./docs/usage.md#all-parameters) could include:
- Different schemes
- Metadata
- Filtering options
- Using base `artic minion` instead of nextflow implementation
- Running SnpEff for variant consequence prediction
- Output reporting options

Medaka model information [can be found here](https://github.com/nanoporetech/medaka#models)

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

[Optional inputs](./docs/usage.md#all-parameters) could include:
- Different schemes
- Metadata
- Filtering options
- Using base `artic minion` instead of nextflow implementation
- Running SnpEff for variant consequence prediction
- Output reporting options

## Outputs
Outputs are separated based off of their tool or file format and found in the `results/` directory by default.

Outputs include:
- Consensus fasta files
- VCF files
- Bam files
- HTML summary files (either custom or multiqc)

*More output information can be found in the [output docs](./docs/output.md)

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
