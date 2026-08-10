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

Minimal Nanopore Run

```bash
nextflow run phac-nml/viralassembly \
  -profile singularity \
  --input SAMPLESHEET.csv \
  --platform nanopore \
  --reference REF.fa \
  --outdir ./results
```

Run with input samplesheet, different Clair3 model, gff for SnpEff, local nextclade dataset, singularity and custom final report

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

Run with a local Clair3 model, fastq directory, a specific nextclade dataset, conda and the multiqc report output

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

Run with Clair3 primer-pool splitting disabled, minor variant calling enabled, nextclade skipped, metadata file, docker, and custom final report

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

Minimal Illumina Run

```bash
nextflow run phac-nml/viralassembly \
  -profile singularity \
  --input SAMPLESHEET.csv \
  --platform illumina \
  --reference REF.fa \
  --outdir ./results
```

Run with input samplesheet, automatic nextclade dataset detection, freebayes variant caller (default), skipping SnpEff annotation, skipping QC, singularity, and multiqc report

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

Run with input samplesheet, virus name for nextclade dataset specification and virus specific processes, freebayes variant caller (default), docker, and multiqc report

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

Run with fastq_directory, gff for SnpEff, iVar variant caller, nextclade dataset with name and tag, conda, and custom report

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
