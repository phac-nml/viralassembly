# phac-nml/viralassembly: Example Commands

A variety of example commands using different parameter options to display how to use each

## Index

- [phac-nml/viralassembly: Example Commands](#phac-nmlviralassembly-example-commands)
  - [Index](#index)
  - [Nanopore](#nanopore)
    - [Non-Amplicon](#non-amplicon)
    - [Amplicon](#amplicon)
  - [Illumina](#illumina)
    - [Non-Amplicon](#non-amplicon-1)
    - [Amplicon](#amplicon-1)
  - [Parameter File](#parameter-file)

## Nanopore

### Non-Amplicon

Minimal input Nanopore run - Containerization with Singularity, input reads with samplesheet, calls variants using clair3 with a default model of `r1041_e82_400bps_sup_v420`, attempts to run SnpEff and Nextclade using output consensus sequences, and creates a custom final report

```bash
nextflow run phac-nml/viralassembly \
  -profile singularity \
  --input SAMPLESHEET.csv \
  --platform nanopore \
  --reference REF.fa \
  --outdir ./results
```

Input samplesheet run, containerization with Singularity, calls variants using a non-default Clair3 model, runs SnpEff with gff file and Nextclade using a local nextclade dataset, and creates a custom final report

```bash
nextflow run phac-nml/viralassembly \
  -profile singularity \
  --input SAMPLESHEET.csv \
  --platform nanopore \
  --reference REF.fa \
  --model r1041_e82_400bps_sup_v410 \
  --gff REF.gff \
  --nextclade_dataset_dir PATH/TO/LOCAL_DIR \
  --outdir ./results
```

### Amplicon

Fastq directory input run, Conda as a dependency manager, calls variants using a local Clair3 model, attempts to run SnpEff, runs Nextclade using the latest version of a specific Nextclade dataset downloaded from Nextstrain, and creates a final multiqc report output

```bash
nextflow run phac-nml/viralassembly \
  -profile conda \
  --fastq_pass FASTQ_PASS/ \
  --platform nanopore \
  --local_model PATH/TO/r1041_e82_400bps_sup_v420 \
  --reference REF.fa \
  --primer_bed PRIMER.bed \
  --nextclade_dataset_name DATASET_NAME \
  --multiqc_report \
  --outdir ./results
```

Fastq directory input run for barcoded directories with a metadata file for sample names, containerization with Docker, calls variants using clair3 with a default model of `r1041_e82_400bps_sup_v420` with Clair3 primer-pool splitting disabled, attempts to run SnpEff, minor variant calling enabled, Nextclade skipped, and creates a custom final report

```bash
nextflow run phac-nml/viralassembly \
  -profile docker \
  --fastq_pass FASTQ_PASS/ \
  --platform nanopore \
  --reference REF.fa \
  --primer_bed PRIMER.bed \
  --no_pool_split \
  --minor_varaints \
  --skip_nextclade \
  --metadata FILE.tsv \
  --outdir ./results
```

---

## Illumina

### Non-Amplicon

Minimal input Illumina run - Containerization with Singularity, input reads with samplesheet, calls variants using FreeBayes as the default variant caller for Illumina data, attempts to run SnpEff and Nextclade using output consensus sequences, and creates a custom final report

```bash
nextflow run phac-nml/viralassembly \
  -profile singularity \
  --input SAMPLESHEET.csv \
  --platform illumina \
  --reference REF.fa \
  --outdir ./results
```

Input samplesheet run, containerization with Singularity, calls variants using FreeBayes (default), runs Nextclade using output consensus sequences, skips SnpEff annotation, skips QC, and creates a final multiqc report output

```bash
nextflow run phac-nml/viralassembly \
  -profile singularity \
  --input SAMPLESHEET.csv \
  --platform illumina \
  --reference REF.fa \
  --skip_qc \
  --skip_snpeff \
  --multiqc_report \
  --outdir ./results
```

### Amplicon

Input samplesheet run, containerization with Docker, virus name for nextclade dataset specification and virus specific processes, calls variants using FreeBayes (default), attempts to run SnpEff, and creates a final multiqc report output

```bash
nextflow run phac-nml/viralassembly \
  -profile docker \
  --input SAMPLESHEET.csv \
  --platform illumina \
  --reference REF.fa \
  --primer_bed PRIMER.bed \
  --virus_name VIRUS_NAME \
  --multiqc_report
  --outdir ./results
```

Fastq directory input run, Conda as a dependency manager, calls variants using iVar variant caller, runs SnpEff with gff file, runs Nextclade using a specific version of a specific nextclade dataset downloaded from Nextstrain and creates custom final report

```bash
nextflow run phac-nml/viralassembly \
  -profile conda \
  --fastq_pass FASTQ_PASS/ \
  --platform illumina \
  --reference REF.fa \
  --primer_bed PRIMER.bed \
  --gff REF.gff \
  --use_ivar \
  --nextclade_dataset_name DATASET_NAME \
  --nextclade_dataset_tag DATASET_TAG \
  --outdir ./results
```

---

## Parameter File

Instead of supplying all parameters on the command line every time, they can be stored in a YAML file and invoked as follows:

```bash
nextflow run phac-nml/viralassembly \
  -profile docker \
  -params-file params.yaml
```

Example `params.yaml` for Nanopore amplicon run:

```yaml
fastq_pass: "FASTQ_PASS/"
platfrom: "nanopore"
reference: "REF.fa"
primer_bed: "PRIMER.bed"
virus_name: "VIRUS_NAME"
gff: "REF.gff"
metadata: "FILE.tsv"
local_model: "PATH/TO/LOCAL_MODEL"
outdir: "./results"
```

Example `params.yaml` for Illumina non-amplicon run:

```yaml
input: "SAMPLESHEET.csv"
platfrom: "illumina"
reference: "REF.fa"
virus_name: "VIRUS_NAME"
gff: "REF.gff"
outdir: "./results"
```
