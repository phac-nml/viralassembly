# phac-nml/viralassembly: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from either the MultiQC report or the custom report, which both summarise the results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps (where a `*` indicates a final output kept in the top level results directory):

- [Preprocessing](#preprocessing)
  - [Reference Stats](#reference-stats)* - Get reference genome information needed for variant calling and QC
  - [Artic Guppyplex](#artic-guppyplex) - Read length filtering
  - [Chopper](#chopper) - Additional Read QC
  - [Nanostat](#nanostat) - Read statistics

- [Variant Calling](#variant-calling)
  - [Minimap2](#minimap2)* - Read mapping
  - [Artic Align Trim](#artic-align_trim)* - Primer trimming and normalisation
  - [Clair3](#clair3) - Determine initial variants with clair3
  - [Medaka](#medaka) - Determine initial variants with medaka
  - [Nanopolish](#nanopolish) - Determine initial variants with nanopolish
  - [Longshot](#longshot)* - Genotype and phase called medaka variants
  - [Variant Filter](#variant-filter) - Filter variants not matching required criteria

- [Consensus Generation](#consensus-generation)
  - [Artic Mask](#artic-mask) - Mask failing variants and low depth sites in preparation for consensus generation
  - [BCFtools Norm](#bcftools-norm)* - Left-align and normalize indels along with make sure the reference alleles match
  - [BCFtools Consensus](#bcftools-consensus)* - Create consensus sequence from VCF variants and Masked sites

- [QC and Reporting]()
  - [SnpEff]() - Variant annotation and functional prediction
  - Mosdepth - Genome and amplicon coverage QC plots
  - Qualimap BAMQC - Alignment quality and metrics
  - Samtools - Flagstat QC information
  - BCFtools - Variant quality and statistics
  - Variation CSV - Custom reporting script for finding and calculating variation in the BAM pileups
  - Amplicon Completeness - Custom reporting script for calculating amplicon completeness based on bedtools output
  - [QC compilation]() - Custom reporting scripts for each sample and the overall run
  - [MultiQC]() - Sample and Run HTML visual report
  - [Custom Report]() - Custom HTML Run visual report

Additionally [Pipeline information](#pipeline-information) which includes report metrics generated during the workflow execution can also be found

### Preprocessing
Initial processing steps and statistic gathering. The reference statistics are output to their own final folder while the other statistics are passed to the final multiqc report.

#### Reference Stats
<details markdown="1">
<summary>Output files</summary>

- `reference/`
  - `genome.bed`: Genomic information in bed format that has the coordiantes of the reference genome needed for nanopolish
  - `refstats.txt`: Genomic information in a format needed for clair3
  - `*.fai`: Samtools faidx fai file for reference genome
</details>

The reference files are generated with both `awk` and `samtools` and are needed as different inputs for downstream tools.

#### Artic Guppyplex
Select reads by size and generate size selected fastq files.

#### Chopper
Filter and trim fastq reads by quality and length (Non-amplicon only).

#### Nanostat
[Nanostat](https://github.com/wdecoster/nanostat) generates plots and statistics on trimmed fastq files for the final multiqc reports (Non-amplicon only).

----------

### Variant Calling
Read mapping and variant calling. Note that only one of `clair3`, `medaka`, and `nanopolish` is used. In the end, final normalized passing and failing variants are output along with the BAM files to their respective folders.

#### Minimap2
<details markdown="1">
<summary>Output files</summary>

- `bam/`
  - `*.sorted.bam`: Sorted bam file from minimap2 and samtools
</details>

The sorted BAM file from minimap2 and samtools.

#### Artic Align_Trim
*Amplicon only*
<details markdown="1">
<summary>Output files</summary>

- `bam/`
  - `*.trimmed.rg.sorted.bam`: Artic align_trim output which normalises coverage and assigns reads to amplicons
  - `*.primertrimmed.rg.sorted.bam`: Artic align_trim output which normalises coverage and assigns reads to amplicons along with softmasking the primer sequences
    - The primertrimmed file is used for subsequent variant calling
</details>

See [the artic core pipeline](https://artic.readthedocs.io/en/latest/minion/#core-pipeline) for more info on how `align_trim` trims the BAM files.

#### Clair3
Run clair3 variant caller on BAM files to create initial variant calls in VCF format.

#### Medaka
Run medaka variant caller on BAM files to create initial variant calls in VCF format.

#### Nanopolish
Run nanopolish variant caller on BAM files, fast5 files, and the sequencing summary file to create initial variant calls in VCF format.

#### Longshot
<details markdown="1">
<summary>Output files</summary>

- `vcf/`
  - `*.longshot.merged.vcf`: Longshot phased VCF file
</details>

Genotype and phase the variants from the initial medaka VCF variant file. [Longshot](https://github.com/pjedge/longshot)

#### Variant Filter
<details markdown="1">
<summary>Output files</summary>

- `vcf/`
  - `*.pass.vcf.gz`: VCF file containing variants passing quality filters
  - `*.pass.vcf.gz.tbi`: VCF index file containing variants passing quality filters
  - `*.fail.vcf`: VCF file containing variants failing quality filters
</details>

Pass/Fail variants based on quality for the final consensus sequence generation.

----------

### Consensus Generation
Final consensus sequence generation based on passing/failing variants and sequencing depth.

#### Artic Mask
Mask low depth and failing variants to create a preconsensus sequence for BCFtools consensus.

#### BCFtools Norm
<details markdown="1">
<summary>Output files</summary>

- `vcf/`
  - `*.pass.norm.vcf.gz`: VCF file containing variants passing quality filters that have their indels normalized and reference positions fixed
    - Reference positions may need to be fixed if there are overlapping variants
</details>

BCFtools norm is utilized to fix locations in which one two variants overlap which during BCFtools consensus would crash the pipeline previously. [BCFtools](https://samtools.github.io/bcftools/bcftools.html#norm)

#### BCFtools Consensus
<details markdown="1">
<summary>Output files</summary>

- `consensus/`
  - `*..consensus.fasta`: Fasta file containing the final output consensus sequence with applied variants and masked sites
</details>

Final output consensus sequence for the sample with variants applied and low coverage/failing variants masked with N's. [BCFtools](https://samtools.github.io/bcftools/bcftools.html#norm)

----------

### 










Pass/Fail variants based on quality for the final consensus sequence generation

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

:::note
The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.
:::

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
