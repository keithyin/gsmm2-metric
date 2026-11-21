//! gskits `fastx_reader` 的 noodles-fasta/fastq 移植。
//!
//! 行为与 gskits(bio crate 版)保持一致:
//! - fasta:name 取 header 首个空白分隔 token,无 qual;
//! - fastq:name 取首 token,qual 为 phred(即 ASCII - 33)。

use std::{fs::File, io::BufReader};

use noodles::fasta;
use noodles::fastq;

use crate::kits::read_info::ReadInfo;

pub struct FastaFileReader {
    fname: String,
    reader: fasta::io::Reader<BufReader<File>>,
}

impl FastaFileReader {
    pub fn new(fname: String) -> Self {
        let file = File::open(&fname).unwrap_or_else(|e| panic!("open file error: {fname}: {e}"));
        let reader = fasta::io::Reader::new(BufReader::new(file));
        Self { fname, reader }
    }

    pub fn get_fname(&self) -> &str {
        &self.fname
    }
}

impl Iterator for FastaFileReader {
    type Item = ReadInfo;

    fn next(&mut self) -> Option<Self::Item> {
        let record = match self.reader.records().next() {
            Some(Ok(record)) => record,
            Some(Err(e)) => panic!("read fasta error {}: {e}", self.fname),
            None => return None,
        };

        Some(ReadInfo::new_fa_record(
            String::from_utf8_lossy(record.name()).into_owned(),
            String::from_utf8_lossy(record.sequence().as_ref()).into_owned(),
        ))
    }
}

pub struct FastqReader {
    fname: String,
    reader: fastq::io::Reader<BufReader<File>>,
}

impl FastqReader {
    pub fn new(fname: String) -> Self {
        let file = File::open(&fname).unwrap_or_else(|e| panic!("open {fname} error: {e}"));
        let reader = fastq::io::Reader::new(BufReader::new(file));
        Self { fname, reader }
    }

    pub fn get_fname(&self) -> &str {
        &self.fname
    }
}

impl Iterator for FastqReader {
    type Item = ReadInfo;

    fn next(&mut self) -> Option<Self::Item> {
        let record = match self.reader.records().next() {
            Some(Ok(record)) => record,
            Some(Err(e)) => panic!("read fastq record error. {}: {e}", self.fname),
            None => return None,
        };

        Some(ReadInfo::new_fq_record(
            String::from_utf8_lossy(record.definition().name()).into_owned(),
            String::from_utf8_lossy(record.sequence()).into_owned(),
            record.quality_scores().iter().map(|v| v - 33).collect(),
        ))
    }
}

pub fn read_fastx<T: Iterator<Item = ReadInfo>>(reader: T) -> Vec<ReadInfo> {
    reader.into_iter().collect()
}
