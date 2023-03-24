# Changelog

## 0.4.0 — 03-24-2023

### Added

* Adds GRCh38FullAnalysisSetWithDecoyHLA (used by 1000 Genomes).
* `ngs convert`: adds conversion between many common next-generation sequencing
  file formats.
* `ngs index`: adds CRAM indexing.
* `ngs list`: adds `ngs list plots`.
* `ngs plot`: adds `VariantAlleleDistributionPlot`.
* `ngs qc`: adds VAF file output to Edits facet.
* `ngs view`: adds `ngs view`.

### Revised

* Makes `--verbose` and `--quiet` global arguments.
* Updates `noodles` to 0.34.0.
* `ngs list`: changes `ngs list reference-genomes` -> `ngs list genomes`.
* `ngs plot`: adds `--only` flag to `ngs plot sample` and `ngs plot cohort`.
* `ngs qc`: adds `--only` argument to qc command.
* `ngs qc`: now preserves struct order of JSON outputs.
* `ngs qc`: improves logging for qc subcommand.
* `ngs qc`: adds coverage bins to coverage qc facet.
* `ngs qc`: updates second-pass qc facets to respect `-n` records

## 0.3.0 — 10-10-2022

### Added

* Adds reference genome support for `GRCh38_no_alt_AnalysisSet`, `hs37d5`,
  `hg38m1x`, and `T2T-CHM13`.
* `ngs index`: adds `ngs index command to index common bioinformatics formats.
* `ngs list`: adds `ngs list` command to list out particular subjects that are
  supported by the `ngs` command line tool.
* `ngs plot`: adds `ngs plot` command to visualize the output of data from `ngs
  qc`. In the initial sample-level implementation, we support graphs for GC
  Content Distribution & Quality Score Distribution. In the initial cohort-level
  implementation, we just support GC Content Distribution.
* `ngs qc`: adds coverage and edits quality check in a second pass.
* `ngs qc`: adds mate mismatched sequence id and CIGAR accumulation to the
  General quality control facet.

### Revised

* Unifies command line arguments for number of records (`-n`).
* `ngs generate`: now supports better read names (the location where read one
  originated from is now in the read name).
* `ngs qc`: all results are aggregated into a single file now.

### Fixed

* `ngs generate`: fixed off by one error when generating records (one too many
  records was being generated).


### Important Chores

* Minimum supported Rust version is now 1.64.0.
* Updates license to be either MIT or Apache 2.0 licensed (at the user's
  discretion).
* Updates dependencies as of 09/29/2022.
* Adds lint groups for documentation, Rust 2021 compatability, and Rust 2018
  idioms. This caused a few changes in the code, as well as a massive
  improvement in documentation.
* The code was reorganized by subcommand (in terms of file system structure).


## 0.2.0 — 09-18-2022

### Added

* `ngs generate`: adds support to generate files resulting from one or more reference FASTAs.
* `ngs qc`: adds `ngs qc` command with the initial four modules:
  * _General metrics_ reports general statistics about the records contained within the file.
  * _GC content_ reports statistics regarding the GC content for records within the file.
  * _Genomic features_ reports statistics regarding genomic features contained within a GFF file.
  * _Template length_ reports statistics regarding the template lengths contained within the file.

### Removed

* `ngs flagstat`: `ngs flagstat` has been removed in favor of the general metrics module in `ngs qc`.

## 0.1.0 — 05-10-2022

### Added

* `ngs derive`: added experimental support for `ngs derive instrument`.
* `ngs flagstat`: added experimental support for `ngs flagstat`.
