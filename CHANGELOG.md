# phac-nml/viralassembly: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.2.0-dev] - Unreleased

Large update get ready to go into IRIDA-Next surveillance platform. Mostly focusing on best practices with not too many logic changes overall. `Clair3` has been made the default and recommended variant caller with most of the changes focusing on it. Primer schemes were also changed to just require a primer bed file and a reference file to make them easier to run

Parameters have been added and adjusted so that is something to be aware of. This is a dev release as we move to release 2.0.0 with all the changes coming later.

### `Added`

- Clair3 specific variant calling parameters
  - `--min_qual_clair3`
  - `--min_frameshift_qual`
  - `--min_allele_freq`
- Generic variant calling parameters added
  - `--min_depth`
- Primer bed parameter added to run primer schemes
  - `--primer_bed`
- Basic nf-tests added
- `nf-iridanext` plugin and config added

### `Removed`

- Specifics relating to amplicon primer-scheme repos and formatting to just require a bed and reference file.
  - This includes the following parameters:
    - `--scheme`
    - `--scheme_version`
    - `--scheme_repo`
    - `--local_scheme`
  - Along with the workflows and processes associated with downloading and checking:
    - `DOWNLOAD_SCHEME`
    - `SIMPLE_SCHEME_VALIDATE`
- Removed the `CAT_FASTQ` module and process when creating the sample channel with the input list
  - Change made just as how it is being done needs to be reevaluated

### `Changed`

- Minimum nextflow version bumped to `24.10.0`
- `nf-validation` replaces `nf-schema`
- The `lib/` folder has been replaced by `nf-validation` along with the nf-core utils subworkflows along with the local initialization subworkflow
- `--variant_caller` defaults to 'clair3' and isn't required to be set
- `--clair3_model` defaults to 'r1041_e82_400bps_sup_v420' now
- `--neg_ctrl_substrings` defaults to 'neg,ntc,blank,water' now (added in ',water')

### `ToDo`

- Fix SnpEff workflow
- Reevaluate more of the best practices for modules, parameters, workflows, and the modules.config file

## v1.1.0 - Unreleased

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

## [v1.0.0] - 2024-03-22

Initial release of `phac-nml/viralassembly`, created from combining the [nf-core](https://nf-co.re/) template with the artic steps.

### `Added`

- All initial pipeline features and logic
- All initial docs and images

[v1.0.0]: https://github.com/phac-nml/measeq/releases/tag/1.0.0
