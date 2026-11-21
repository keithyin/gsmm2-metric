//! gskits `ds::ReadInfo` 的本地移植。
//!
//! `from_bam_record` 由 rust-htslib 版本改为接收 noodles `RecordBuf`,
//! 字段提取语义(见 `bam_ext::BamRecordExt`)与 gskits 保持一致。

use std::collections::HashSet;

use noodles::sam::alignment::RecordBuf;

use crate::kits::bam_ext::BamRecordExt;

#[derive(Debug, Default)]
pub struct ReadInfo {
    pub name: String,
    pub seq: String,
    pub cx: Option<u8>,
    pub ch: Option<u32>,
    pub np: Option<u32>,
    pub rq: Option<f32>,
    pub qual: Option<Vec<u8>>, // phreq, no offset
    pub dw: Option<Vec<u8>>,
    pub ar: Option<Vec<u8>>,
    pub cr: Option<Vec<u8>>,
    pub be: Option<Vec<u32>>,
    pub nn: Option<Vec<u8>>,
    pub wd: Option<Vec<u8>>, // linker width
    pub sd: Option<Vec<u8>>, // standard devition
    pub sp: Option<Vec<u8>>, // slope
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
    ///
    /// 修复:gskits 0.15.15 原实现里 `sp` tag 实际读的是 `b"sd"`(笔误),
    /// 这里按字面语义读 `b"sp"`。
    pub fn from_bam_record(
        record: &RecordBuf,
        qname_suffix: Option<&str>,
        tags: &HashSet<String>,
    ) -> Self {
        let record_ext = BamRecordExt::new(record);
        let mut qname = record_ext.get_qname();
        if let Some(suffix) = qname_suffix {
            qname.push_str(suffix);
        }

        let to_u8_list = |v: Option<Vec<u32>>| -> Option<Vec<u8>> {
            v.map(|v| v.into_iter().map(|v| v as u8).collect())
        };

        let dw = if tags.contains("dw") {
            to_u8_list(record_ext.get_dw())
        } else {
            None
        };

        let ar = if tags.contains("ar") {
            to_u8_list(record_ext.get_ar())
        } else {
            None
        };

        let cr = if tags.contains("cr") {
            to_u8_list(record_ext.get_cr())
        } else {
            None
        };

        let nn = if tags.contains("nn") {
            to_u8_list(record_ext.get_nn())
        } else {
            None
        };

        let wd = if tags.contains("wd") {
            to_u8_list(record_ext.get_uint_list(b"wd"))
        } else {
            None
        };

        let sd = if tags.contains("sd") {
            to_u8_list(record_ext.get_uint_list(b"sd"))
        } else {
            None
        };

        let sp = if tags.contains("sp") {
            to_u8_list(record_ext.get_uint_list(b"sp")) // 修复:原 gskits 读的是 b"sd"
        } else {
            None
        };

        Self {
            name: qname,
            seq: record_ext.get_seq(),
            cx: record_ext.get_cx(),
            ch: record_ext.get_ch(),
            np: record_ext.get_np(),
            rq: record_ext.get_rq(),
            qual: Some(record_ext.get_qual().to_vec()),
            dw,
            ar,
            cr,
            be: record_ext.get_be(),
            nn,
            wd,
            sd,
            sp,
        }
    }
}
