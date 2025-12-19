# phac-nml/viralassembly: Example Commands

A variety of example commands using different parameter options to display how to use each

## Amplicon

### Clair3

Clair3 with a local model, primer file, fastq directory, conda, and the multiqc report output

```bash
nextflow run phac-nml/viralassembly \
  -profile conda \
  --fastq_pass FASTQ_PASS/ \
  --clair3_model PATH/TO/r1041_e82_400bps_sup_v420 \
  --reference REF.fa \
  --primer_bed PRIMER.bed \
  --multiqc_report \
  --outdir ./results
```

### Medaka

Minimal input medaka with conda, an input csv file for data, and just a reference file

```bash
nextflow run phac-nml/viralassembly \
  -profile conda \
  --input INPUT.csv \
  --variant_caller 'medaka' \
  --reference REF.fa \
  --primer_bed PRIMER.bed \
  --outdir ./results
```

### Nanopolish

Nanopolish run using singularity

```bash
nextflow run phac-nml/viralassembly \
  -profile singularity \
  --input INPUT.csv \
  --fast5_pass FAST5_PASS/ \
  --sequencing_summart SEQ_SUM.txt \
  --variant_caller 'nanopolish' \
  --reference REF.fa \
  --primer_bed PRIMER.bed \
  --outdir ./results
```

---

## Non-Amplicon

### Clair3

Minimal clair3 with docker using a fastq input directory along wth a gff3 reference file for SnpEff

```bash
nextflow run phac-nml/viralassembly \
  -profile docker \
  --fastq_pass FASTQ_PASS/ \
  --reference ./REFERENCE.fa \
  --gff ./REFERENCE.gff
```

### Medaka

Medaka with conda skipping QC and SnpEff

```bash
nextflow run phac-nml/viralassembly \
  -profile conda \
  --input INPUT.csv \
  --variant_caller 'medaka' \
  --reference ./REFERENCE.fa \
  --skip_qc \
  --skip_snpeff
```

### Nanopolish

Nanopolish running with conda, filtering the read lengths to be shorter, and creating a custom report

```bash
nextflow run phac-nml/viralassembly \
  -profile conda \
  --input INPUT.csv \
  --fast5_pass FAST5_PASS/ \
  --sequencing_summart SEQ_SUM.txt \
  --variant_caller 'nanopolish' \
  --reference ./REFERENCE.fa \
  --min_length 100 \
  --max_length 600 \
  --outdir ./results
```
