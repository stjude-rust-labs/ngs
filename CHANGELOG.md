# Changelog

## 0.2.0 — 09-18-2022

### Added

* `ngs qc`: adds `ngs qc` command with the initial four modules:
  * _General metrics_ reports general statistics about the records contained within the file.
  * _GC content_ reports statistics regarding the GC content for records within the file.
  * _Genomic features_ reports statistics regarding genomic features contained within a GFF file.
  * _Template length_ reports statistics regarding the template lengths contained within the file.
* `ngs generate`: adds support to generate files resulting from one or more reference FASTAs.

### Removed

* `ngs flagstat`: `ngs flagstat` has been removed in favor of the general metrics module in `ngs qc`.

## 0.1.0 — 05-10-2022

### Added

* `ngs derive`: added experimental support for `ngs derive instrument`.
* `ngs flagstat`: added experimental support for `ngs flagstat`.
