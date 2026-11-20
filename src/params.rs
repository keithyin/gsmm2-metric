//! 本地参数结构，字段与方法名与被移除的 mm2 0.25.6 `params.rs` 保持一致，
//! cli.rs / 测试代码只改 use 路径即可。bandmap 侧的映射见 `aligner.rs`。

use std::collections::HashSet;

use gskits::{
    gsbam::bam_record_ext::{BamRecord, BamRecordExt},
    utils::Range,
};

#[derive(Default)]
pub struct InputFilterParams {
    pub np_range: Option<Range<i32>>,
    pub rq_range: Option<Range<f32>>,
    pub qname_suffix: Option<String>,
}

impl InputFilterParams {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_np_range(mut self, np_range_str: &str) -> Self {
        self.np_range = Some(Range::<i32>::new(np_range_str));
        self
    }

    pub fn set_rq_range(mut self, rq_range_str: &str) -> Self {
        self.rq_range = Some(Range::<f32>::new(rq_range_str));
        self
    }

    pub fn valid(&self, record: &BamRecord) -> bool {
        let record_ext = BamRecordExt::new(record);

        if let Some(np_range_) = &self.np_range {
            if let Some(np) = record_ext.get_np() {
                if !np_range_.within_range(np as i32) {
                    return false;
                }
            }
        }

        if let Some(rq_range_) = &self.rq_range {
            if let Some(rq) = record_ext.get_rq() {
                if !rq_range_.within_range(rq) {
                    return false;
                }
            }
        }

        if let Some(qname_suffix) = &self.qname_suffix {
            if !record_ext.get_qname().ends_with(qname_suffix) {
                return false;
            }
        }

        true
    }
}

#[derive(Debug, Clone, Default)]
pub struct AlignParams {
    pub matching_score: Option<i32>,
    pub mismatch_penalty: Option<i32>,
    pub gap_open_penalty: Option<String>,
    pub gap_extension_penalty: Option<String>,

    pub min_cnt: Option<i32>,
    pub min_dp_max: Option<i32>,
    pub min_chain_score: Option<i32>,
    pub min_ksw_len: Option<i32>,
}

impl AlignParams {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_m_score(mut self, m_score: i32) -> Self {
        self.matching_score = Some(m_score);
        self
    }

    pub fn set_mm_score(mut self, mm_score: i32) -> Self {
        self.mismatch_penalty = Some(mm_score);
        self
    }

    pub fn set_gap_open_penalty(mut self, go: String) -> Self {
        self.gap_open_penalty = Some(go);
        self
    }

    pub fn set_gap_extension_penalty(mut self, ge: String) -> Self {
        self.gap_extension_penalty = Some(ge);
        self
    }

    /// 与 mm2 `TOverrideAlignerParam::modify_aligner` 等价的 bandmap 映射：
    /// - a -> scoring.match_score, b -> scoring.mismatch_penalty
    /// - "q,q2" -> scoring.gap_open / gap_open2
    /// - "e,e2" -> scoring.gap_extend / gap_extend2
    /// - min_cnt -> chaining.min_cnt, min_chain_score -> chaining.min_chain_score
    /// - min_dp_max -> alignment.min_dp_max, min_ksw_len -> alignment.min_dp_len
    pub fn apply_to(&self, opt: &mut bandmap::align::map::MapOptions) {
        if let Some(m) = self.matching_score {
            opt.scoring.match_score = m;
        }

        if let Some(mm) = self.mismatch_penalty {
            opt.scoring.mismatch_penalty = mm;
        }

        if let Some(gap_o) = &self.gap_open_penalty {
            if gap_o.contains(",") {
                let (o1, o2) = gap_o.rsplit_once(",").unwrap();
                opt.scoring.gap_open = o1.parse::<i32>().unwrap();
                opt.scoring.gap_open2 = o2.parse::<i32>().unwrap();
            } else {
                opt.scoring.gap_open = gap_o.parse::<i32>().unwrap();
            }
        }

        if let Some(gap_e) = &self.gap_extension_penalty {
            if gap_e.contains(",") {
                let (e1, e2) = gap_e.rsplit_once(",").unwrap();
                opt.scoring.gap_extend = e1.parse::<i32>().unwrap();
                opt.scoring.gap_extend2 = e2.parse::<i32>().unwrap();
            } else {
                opt.scoring.gap_extend = gap_e.parse::<i32>().unwrap();
            }
        }

        if let Some(min_cnt) = self.min_cnt {
            opt.chaining.min_cnt = min_cnt;
        }

        if let Some(min_dp_max) = self.min_dp_max {
            opt.alignment.min_dp_max = min_dp_max;
        }

        if let Some(min_chain_score) = self.min_chain_score {
            opt.chaining.min_chain_score = min_chain_score;
        }
        if let Some(min_ksw_len) = self.min_ksw_len {
            opt.alignment.min_dp_len = min_ksw_len;
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct MapParams {}

#[derive(Debug, Clone, Copy, Default)]
pub struct IndexParams {
    pub kmer: Option<usize>,
    pub wins: Option<usize>,
}

impl IndexParams {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_kmer(mut self, kmer: usize) -> Self {
        self.kmer = Some(kmer);
        self
    }

    pub fn set_wins(mut self, wins: usize) -> Self {
        self.wins = Some(wins);
        self
    }
}

#[derive(Debug, Clone)]
pub struct OupParams {
    pub discard_secondary: bool,
    pub discard_supplementary: bool,
    pub oup_identity_threshold: f32,
    pub oup_coverage_threshold: f32,
    pub discard_multi_align_reads: bool,
    pub pass_through_tags: HashSet<String>,
}

impl Default for OupParams {
    fn default() -> Self {
        Self {
            discard_secondary: false,
            discard_supplementary: false,
            oup_identity_threshold: -1.0,
            oup_coverage_threshold: -1.0,
            discard_multi_align_reads: false,
            pass_through_tags: HashSet::new(),
        }
    }
}

impl OupParams {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_discard_secondary(mut self, discard_secondary: bool) -> Self {
        if discard_secondary {
            assert!(self.discard_multi_align_reads == false);
        }
        self.discard_secondary = discard_secondary;
        self
    }

    pub fn set_discard_supplementary(mut self, discard_supplementary: bool) -> Self {
        if discard_supplementary {
            assert!(self.discard_multi_align_reads == false);
        }
        self.discard_supplementary = discard_supplementary;
        self
    }
    pub fn set_oup_identity_threshold(mut self, oup_identity_threshold: f32) -> Self {
        self.oup_identity_threshold = oup_identity_threshold;
        self
    }
    pub fn set_oup_coverage_threshold(mut self, oup_coverage_threshold: f32) -> Self {
        self.oup_coverage_threshold = oup_coverage_threshold;
        self
    }
    pub fn set_discard_multi_align_reads(mut self, discard_multi_align_reads: bool) -> Self {
        if discard_multi_align_reads {
            assert!(self.discard_secondary == false);
            assert!(self.discard_supplementary == false);
        }

        self.discard_multi_align_reads = discard_multi_align_reads;
        self
    }
}
