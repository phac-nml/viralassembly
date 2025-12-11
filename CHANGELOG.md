# phac-nml/viralassembly: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0 - [Unreleased]

### `Added`

- Input schema JSON and validation
- FORMAT_INPUT workflow
  - Handles the input data now
- `nf-schema@2.0.0` plugin

### `Changed`

- `--input SAMPLESHEET_CSV` header
  - Went from `reads` with path to barcode directories to `fastq_1` with path to fastq files
- Fixed bug so that SNPEff will now work with given gff files
  - Issue was typo related in the build module
- Fixed bug with `calc_bam_variation` caused by genome case
- Log and error statements
- Fixed the cache directory statements

## v1.0.0 - [2024-03-22]

Initial release of `phac-nml/viralassembly`, created from combining the [nf-core](https://nf-co.re/) template with the artic steps.

### `Added`

- All initial pipeline features and logic
- All initial docs and images
