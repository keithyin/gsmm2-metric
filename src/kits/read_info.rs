//! gskits `ds::ReadInfo` 的本地移植。
//!
//! `from_bam_record` 由 rust-htslib 版本改为接收 noodles `RecordBuf`,
//! 字段提取语义(见 `bam_ext::BamRecordExt`)与 gskits 保持一致。

use noodles::sam::alignment::RecordBuf;

use crate::kits::bam_ext::BamRecordExt;

#[derive(Debug, Default)]
pub struct ReadInfo {
    pub name: String,
    pub seq: String,
    /// BAM np tag（pass 数），输入过滤使用。
    pub np: Option<u32>,
    /// read 质量（fq 路径由 phreq 计算，BAM 路径取 rq tag），输入过滤使用。
    pub rq: Option<f32>,
    pub qual: Option<Vec<u8>>, // phreq, no offset
}

impl ReadInfo {
    pub fn new_fa_record(name: String, seq: String) -> Self {
        Self {
            name,
            seq,
            ..Default::default()
        }
    }

    pub fn new_fq_record(name: String, seq: String, qual: Vec<u8>) -> Self {
        let mut res = ReadInfo::new_fa_record(name, seq);
        res.qual = Some(qual);
        res
    }

    /// gskits `ReadInfo::from_bam_record` 的 noodles 移植。
    pub fn from_bam_record(record: &RecordBuf, qname_suffix: Option<&str>) -> Self {
        let record_ext = BamRecordExt::new(record);
        let mut qname = record_ext.get_qname();
        if let Some(suffix) = qname_suffix {
            qname.push_str(suffix);
        }

        Self {
            name: qname,
            seq: record_ext.get_seq(),
            np: record_ext.get_np(),
            rq: record_ext.get_rq(),
            qual: Some(record_ext.get_qual().to_vec()),
        }
    }
}
