use anyhow::bail;
use noodles_fasta as fasta;
use std::{fs::File, io::BufReader, path::Path};

pub fn open<P>(src: P) -> anyhow::Result<fasta::Reader<BufReader<File>>>
where
    P: AsRef<Path>,
{
    let path = src.as_ref();
    let file = File::open(path);

    match path.extension().and_then(|x| x.to_str()) {
        Some("gz") => bail!(
            "This command does not yet support gzipped FASTA files. \
                       Please unzip your FASTA file and try again."
        ),
        Some("fa") | Some("fna") | Some("fasta") => {
            let reader = file.map(BufReader::new)?;
            Ok(fasta::Reader::new(reader))
        }
        _ => bail!("Unknown extension for FASTA file."),
    }
}
