//! Utilities related to path manipulation.

use std::{ffi::OsStr, path::PathBuf};

use anyhow::bail;

/// A trait that is intended to add a `append_extension` method to [`PathBuf`].
/// This makes it significantly more ergonomic to work with things like index
/// files where the filename is simply the target file name with some extra
/// extension.
pub trait AppendExtension {
    /// Appends an extension with the specified futher extension.
    ///
    /// ```
    /// use std::path::PathBuf;
    /// use ngs::utils::path::AppendExtension;
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
