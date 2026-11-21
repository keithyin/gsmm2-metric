//! gskits `gsbam::bam_record_ext::BamRecordExt` 的 noodles 移植,以及
//! rust-htslib `bam::ext::BamRecordExtensions::aligned_pairs_full` 的语义复刻。

use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::cigar::op::{Kind, Op};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;

pub struct BamRecordExt<'a> {
    record: &'a RecordBuf,
}

impl<'a> BamRecordExt<'a> {
    pub fn new(record: &'a RecordBuf) -> Self {
        BamRecordExt { record }
    }

    pub fn get_qname(&self) -> String {
        self.record
            .name()
            .map(|v| String::from_utf8_lossy(v).into_owned())
            .unwrap_or_default()
    }

    pub fn get_seq(&self) -> String {
        String::from_utf8_lossy(self.record.sequence().as_ref()).into_owned()
    }

    pub fn get_qual(&self) -> &[u8] {
        self.record.quality_scores().as_ref()
    }

    /// 0-based 左端参考位置(unmapped/no-pos 时 panic,对应 gskits 的 assert)。
    pub fn reference_start(&self) -> usize {
        let pos = self
            .record
            .alignment_start()
            .expect("record has no alignment_start (unmapped?)");
        usize::from(pos) - 1
    }

    /// 0-based 右端(开区间)参考位置 = reference_start + sum(M/D/N/=/X)。
    pub fn reference_end(&self) -> usize {
        let aligned_rlen: usize = self
            .record
            .cigar()
            .as_ref()
            .iter()
            .filter(|op| op.kind().consumes_reference())
            .map(|op| op.len())
            .sum();

        self.reference_start() + aligned_rlen
    }

    /// 首个 op 为 soft clip 时取其长度,否则 0。
    pub fn query_alignment_start(&self) -> usize {
        match self.record.cigar().as_ref().first() {
            Some(op) if op.kind() == Kind::SoftClip => op.len(),
            _ => 0,
        }
    }

    /// query_alignment_start + sum(M/I/=/X)。
    pub fn query_alignment_end(&self) -> usize {
        let aligned_qlen: usize = self
            .record
            .cigar()
            .as_ref()
            .iter()
            .filter(|op| {
                matches!(
                    op.kind(),
                    Kind::Match | Kind::Insertion | Kind::SequenceMatch | Kind::SequenceMismatch
                )
            })
            .map(|op| op.len())
            .sum();

        self.query_alignment_start() + aligned_qlen
    }

    pub fn get_cx(&self) -> Option<u8> {
        self.get_int(b"cx").map(|v| v as u8)
    }

    pub fn get_ch(&self) -> Option<u32> {
        self.get_int(b"ch")
    }

    pub fn get_np(&self) -> Option<u32> {
        self.get_int(b"np")
    }

    pub fn get_rq(&self) -> Option<f32> {
        self.get_float(b"rq")
    }

    pub fn get_dw(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"dw")
    }

    pub fn get_ar(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"ar")
    }

    pub fn get_cr(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"cr")
    }

    pub fn get_nn(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"nn")
    }

    pub fn get_be(&self) -> Option<Vec<u32>> {
        self.get_uint_list(b"be")
    }

    /// int 标量(I8/U8/I16/U16/I32/U32)→ u32,语义与 gskits `get_int` 一致。
    pub fn get_int(&self, tag: &[u8; 2]) -> Option<u32> {
        let tag = Tag::new(tag[0], tag[1]);
        match self.record.data().get(&tag)? {
            Value::Int8(v) => Some(*v as u32),
            Value::UInt8(v) => Some(*v as u32),
            Value::Int16(v) => Some(*v as u32),
            Value::UInt16(v) => Some(*v as u32),
            Value::Int32(v) => Some(*v as u32),
            Value::UInt32(v) => Some(*v as u32),
            _ => None,
        }
    }

    /// float 标量 → f32,语义与 gskits `get_float` 一致。
    pub fn get_float(&self, tag: &[u8; 2]) -> Option<f32> {
        let tag = Tag::new(tag[0], tag[1]);
        match self.record.data().get(&tag)? {
            Value::Float(v) => Some(*v),
            _ => None,
        }
    }

    /// uint 数组(B:C/B:S/B:I)→ Vec<u32>,语义与 gskits `get_uint_list` 一致。
    pub fn get_uint_list(&self, tag: &[u8; 2]) -> Option<Vec<u32>> {
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let tag = Tag::new(tag[0], tag[1]);
        match self.record.data().get(&tag)? {
            Value::Array(array) => match array {
                Array::UInt8(v) => Some(v.iter().map(|&v| v as u32).collect()),
                Array::UInt16(v) => Some(v.iter().map(|&v| v as u32).collect()),
                Array::UInt32(v) => Some(v.iter().map(|&v| v as u32).collect()),
                _ => None,
            },
            _ => None,
        }
    }
}

/// rust-htslib `aligned_pairs_full` 的语义复刻:逐碱基产出 `[qpos, rpos]`。
///
/// - qpos:0-based,覆盖整条 read(含 soft clip 前缀);rpos:绝对 0-based 参考位置
/// - `M`/`=`/`X` → `[Some(qpos), Some(rpos)]`,双指针同时推进
/// - `I`/`S` → `[Some(qpos), None]`,只推进 read 指针
/// - `D`/`N` → `[None, Some(rpos)]`,只推进参考指针
/// - `H` 跳过;`P` 不支持(panic,与 rust-htslib 行为一致)
pub fn aligned_pairs_full(record: &RecordBuf) -> AlignedPairsFull<'_> {
    AlignedPairsFull {
        ops: record.cigar().as_ref().iter().copied(),
        genome_pos: match record.alignment_start() {
            Some(pos) => usize::from(pos) as i64 - 1,
            None => -1,
        },
        read_pos: 0,
        remaining_match_bp: 0,
        remaining_ins_bp: 0,
        remaining_del_bp: 0,
    }
}

pub struct AlignedPairsFull<'a> {
    ops: std::iter::Copied<std::slice::Iter<'a, Op>>,
    genome_pos: i64,
    read_pos: i64,
    remaining_match_bp: u32,
    remaining_ins_bp: u32,
    remaining_del_bp: u32,
}

impl Iterator for AlignedPairsFull<'_> {
    type Item = [Option<i64>; 2];

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.remaining_match_bp > 0 {
                self.remaining_match_bp -= 1;
                self.genome_pos += 1;
                self.read_pos += 1;
                return Some([Some(self.read_pos - 1), Some(self.genome_pos - 1)]);
            }
            if self.remaining_ins_bp > 0 {
                self.remaining_ins_bp -= 1;
                self.read_pos += 1;
                return Some([Some(self.read_pos - 1), None]);
            }
            if self.remaining_del_bp > 0 {
                self.remaining_del_bp -= 1;
                self.genome_pos += 1;
                return Some([None, Some(self.genome_pos - 1)]);
            }

            let op = self.ops.next()?;
            let len = op.len() as u32;
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    self.genome_pos += 1;
                    self.read_pos += 1;
                    self.remaining_match_bp = len - 1;
                    return Some([Some(self.read_pos - 1), Some(self.genome_pos - 1)]);
                }
                Kind::Insertion | Kind::SoftClip => {
                    self.read_pos += 1;
                    self.remaining_ins_bp = len - 1;
                    return Some([Some(self.read_pos - 1), None]);
                }
                Kind::Deletion | Kind::Skip => {
                    self.genome_pos += 1;
                    self.remaining_del_bp = len - 1;
                    return Some([None, Some(self.genome_pos - 1)]);
                }
                Kind::HardClip => {
                    // no advance
                }
                Kind::Pad => panic!("Padding (Cigar::Pad) is not supported."),
            }
        }
    }
}

#[cfg(test)]
mod test {
    use noodles::sam::alignment::record::cigar::op::Op;
    use noodles::sam::alignment::record_buf::Cigar as RecordBufCigar;

    use super::*;

    // 复刻 rust-htslib aligned_pairs_full 语义的对照用例。
    // cigar: 5S 3= 2X 2I 1D 2= 1N 2= (seq len = 5+3+2+2+2+2 = 16)
    fn fixture() -> RecordBuf {
        let mut record = RecordBuf::default();
        *record.alignment_start_mut() = noodles::core::Position::try_from(101).ok(); // 0-based 100
        *record.cigar_mut() = RecordBufCigar::from(vec![
            Op::new(Kind::SoftClip, 5),
            Op::new(Kind::SequenceMatch, 3),
            Op::new(Kind::SequenceMismatch, 2),
            Op::new(Kind::Insertion, 2),
            Op::new(Kind::Deletion, 1),
            Op::new(Kind::SequenceMatch, 2),
            Op::new(Kind::Skip, 1),
            Op::new(Kind::SequenceMatch, 2),
        ]);
        record
    }

    #[test]
    fn test_aligned_pairs_full() {
        let pairs: Vec<_> = aligned_pairs_full(&fixture()).collect();

        // 5S → qpos 0..5, rpos None
        assert_eq!(pairs[0], [Some(0), None]);
        assert_eq!(pairs[4], [Some(4), None]);
        // 3= → q 5..8, r 100..103
        assert_eq!(pairs[5], [Some(5), Some(100)]);
        // 2X → q 8..10, r 103..105
        assert_eq!(pairs[8], [Some(8), Some(103)]);
        // 2I → q 10..12, rpos None
        assert_eq!(pairs[10], [Some(10), None]);
        // 1D → r 105
        assert_eq!(pairs[12], [None, Some(105)]);
        // 2= → q 12..14, r 106..108
        assert_eq!(pairs[13], [Some(12), Some(106)]);
        // 1N → r 108
        assert_eq!(pairs[15], [None, Some(108)]);
        // 2= → q 14..16, r 109..111
        assert_eq!(pairs[16], [Some(14), Some(109)]);
        assert_eq!(pairs.len(), 18);
    }

    #[test]
    fn test_hard_clip_skipped() {
        let mut record = RecordBuf::default();
        *record.alignment_start_mut() = noodles::core::Position::try_from(1).ok(); // 0-based 0
        *record.cigar_mut() = RecordBufCigar::from(vec![
            Op::new(Kind::HardClip, 5),
            Op::new(Kind::SequenceMatch, 3),
            Op::new(Kind::HardClip, 5),
        ]);
        let pairs: Vec<_> = aligned_pairs_full(&record).collect();
        assert_eq!(
            pairs,
            vec![
                [Some(0), Some(0)],
                [Some(1), Some(1)],
                [Some(2), Some(2)]
            ]
        );
    }
}
