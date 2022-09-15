use noodles_fasta as fasta;
use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

pub fn open<P>(src: P) -> io::Result<fasta::Reader<BufReader<File>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = File::open(path);

    match path.extension().and_then(|x| x.to_str()) {
        Some("gz") => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "This command does not yet support gzipped FASTA files. \
                       Please unzip your FASTA file and try again.",
        )),
        Some("fa") | Some("fna") | Some("fasta") => {
            let reader = file.map(BufReader::new)?;
            Ok(fasta::Reader::new(reader))
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Unknown extension for FASTA file.",
        )),
    }
}
