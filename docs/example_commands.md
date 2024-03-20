# Example Commands
A variety of example commands using different parameter options to display how to use each

## Amplicon

### Clair3
Clair3 with a local model, local scheme, fastq directory, conda, and the custom report output

```bash
nextflow run phac-nml/viralassembly \
  -profile conda \
  --fastq_pass FASTQ_PASS/ \
  --variant_caller 'clair3' \
  --clair3_model ./r1041_e82_400bps_sup_v420 \
  --local_scheme ./primer_schemes \
  --scheme 'hCMV' \
  --scheme_version 'V1' \
  --custom_report \
  --outdir ./results 
```

### Medaka
Minimal input medaka with conda, an input csv file for data, and the nCoV-2019 scheme

```bash
nextflow run phac-nml/viralassembly \
  -profile conda \
  --input INPUT.csv \
  --variant_caller 'medaka' \
  --scheme 'nCoV-2019' \
  --scheme_version 'V5.3.2' \
  --outdir ./results 
```

### Nanopolish
Nanopolish run using singularity and the base artic command line tool (instead of the default nextflow implementation)

```bash
nextflow run phac-nml/viralassembly \
  -profile singularity \
  --input INPUT.csv \
  --fast5_pass FAST5_PASS/ \
  --sequencing_summart SEQ_SUM.txt \
  --variant_caller 'nanopolish' \
  --scheme 'nCoV-2019' \
  --scheme_version 'V5.3.2' \
  --use_artic_tool \
  --outdir ./results 
```

--------------------------

## Non-Amplicon

### Clair3
Minimal clair3 with docker using a fastq input directory along wth a gff3 reference file for SnpEff

```bash
nextflow run phac-nml/viralassembly \
  -profile docker \
  --fastq_pass FASTQ_PASS/ \
  --variant_caller 'clair3' \
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
