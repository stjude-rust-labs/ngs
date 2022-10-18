//! Extensions to and utilities concerning [`PathBuf`]s.
//!
//! # Overview
//!
//! [`PathBuf`][std::path::PathBuf] leaves a decent amount to be desired by way
//! of ergonomics. This module implements some commonly used methods within the
//! `ngs` command line tool. We recommend you read through the provided
//! functionality carefully and use these methods in place of rolling your own.
//!
//! # Examples
//!
//! For example, one very common task when working with genomic file formats is
//! to create an index for that file. By convention, the indices are generally
//! of the form `<filename>.<extra index>`. Thus, it's desirable to have a
//! method to easily append to an existing [`PathBuf`][std::path::PathBuf].
//!
//! In this module, we provide such a function ([`AppendExtension`] and the
//! respective [`append_extension`][AppendExtension::append_extension] function).
//!
//! ```
//! use std::path::PathBuf;
//! // Trait must be in scope to use it.
//! use ngs::utils::pathbuf::AppendExtension;
//!
//! assert_eq!(
//!     PathBuf::from("hello.txt")
//!         .append_extension("world")
//!         .unwrap(),
//!     PathBuf::from("hello.txt.world"))
//! ```

use std::ffi::OsStr;
use std::path::PathBuf;

use anyhow::bail;

/// A trait that is intended to add a
/// [`append_extension`][AppendExtension::append_extension] method to
/// [`PathBuf`]. This makes it significantly more ergonomic to work with things
/// like index files where the filename is simply the target file name with some
/// extra extension.
pub trait AppendExtension {
    /// Appends an extension with the specified futher extension.
    ///
    /// ```
    /// use std::path::PathBuf;
    /// use ngs::utils::pathbuf::AppendExtension;
    ///
    /// let bam = PathBuf::from("~/test.bam");
    /// let bai = bam.append_extension("bai").unwrap();
    /// assert_eq!(bai.file_name().unwrap(), "test.bam.bai");
    /// ```
    fn append_extension<P>(self, ext: P) -> anyhow::Result<Self>
    where
        Self: Sized,
        P: AsRef<OsStr>;
}

impl AppendExtension for PathBuf {
    fn append_extension<P>(mut self, ext: P) -> anyhow::Result<Self>
    where
        P: AsRef<OsStr>,
    {
        let mut new_ext = match self.extension() {
            Some(ext) => ext.to_os_string(),
            None => bail!("path did not have an extension: {}", self.display()),
        };

        new_ext.push(".");
        new_ext.push(ext);

        self.set_extension(new_ext);
        Ok(self)
    }
}
