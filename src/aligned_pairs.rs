use crate::adapter::Mapping;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AlignOp {
    Match(u32),
    Ins(u32),
    Del(u32),
    RefSkip(u32),
    SoftClip(u32),
    HardClip(u32),
    Pad(u32),
    Equal(u32),
    Diff(u32),
}

impl From<u8> for AlignOp {
    fn from(value: u8) -> Self {
        match value {
            0 => AlignOp::Match(1),
            1 => AlignOp::Ins(1),
            2 => AlignOp::Del(1),
            3 => AlignOp::RefSkip(1),
            4 => AlignOp::SoftClip(1),
            5 => AlignOp::HardClip(1),
            6 => AlignOp::Pad(1),
            7 => AlignOp::Equal(1),
            8 => AlignOp::Diff(1),
            _ => panic!("invalid code. {}", value),
        }
    }
}

pub trait TAlignedPairs {
    fn target_start(&self) -> i64;
    fn target_end(&self) -> i64;
    fn query_start(&self) -> i64; // 原始序列的开始
    fn query_end(&self) -> i64; // 原始序列的结束
    fn is_reverse(&self) -> bool;
    fn cigars(&self) -> &[(u32, u8)]; // 如果是 reverse，那么是 reverse 后的 cigar。

    /// 如果是 reverse的，返回结果中，query的位置是倒序的
    fn aligned_pairs(&self) -> impl Iterator<Item = (Option<i64>, Option<i64>, AlignOp)> {
        let mut query_cursor = if self.is_reverse() {
            self.query_end() - 1
        } else {
            self.query_start()
        };

        let query_corsor_move_direction = if self.is_reverse() { -1 } else { 1 };

        let mut target_cursor = self.target_start();

        self.cigars().iter().copied().flat_map(move |(cnt, op)| {
            let cnt = cnt as i64;
            match op {
                0 => panic!("eqx required"),
                1 => {
                    // ins
                    let old_query_cursor = query_cursor;
                    query_cursor += query_corsor_move_direction * cnt;
                    Box::new((0..cnt).into_iter().map(move |shift| {
                        (
                            Some(old_query_cursor + query_corsor_move_direction * shift),
                            None,
                            op.into(),
                        )
                    }))
                        as Box<dyn Iterator<Item = (Option<i64>, Option<i64>, AlignOp)>>
                }
                2 => {
                    // del
                    let old_target_cursor = target_cursor;
                    target_cursor += cnt;
                    Box::new(
                        (0..cnt)
                            .into_iter()
                            .map(move |shift| (None, Some(old_target_cursor + shift), op.into())),
                    )
                        as Box<dyn Iterator<Item = (Option<i64>, Option<i64>, AlignOp)>>
                }
                7 | 8 => {
                    let old_query_cursor = query_cursor;
                    query_cursor += query_corsor_move_direction * cnt;
                    let old_target_cursor = target_cursor;
                    target_cursor += cnt;

                    Box::new((0..cnt).into_iter().map(move |shift| {
                        (
                            Some(old_query_cursor + query_corsor_move_direction * shift),
                            Some(old_target_cursor + shift),
                            op.into(),
                        )
                    }))
                        as Box<dyn Iterator<Item = (Option<i64>, Option<i64>, AlignOp)>>
                }
                op => panic!("not a valid op:{}", op),
            }
        })
    }
}

impl TAlignedPairs for Mapping {
    fn cigars(&self) -> &[(u32, u8)] {
        self.cigar.as_ref().unwrap()
    }

    fn is_reverse(&self) -> bool {
        self.is_rev()
    }

    fn query_end(&self) -> i64 {
        self.query_end as i64
    }
    fn query_start(&self) -> i64 {
        self.query_start as i64
    }
    fn target_end(&self) -> i64 {
        self.target_end as i64
    }
    fn target_start(&self) -> i64 {
        self.target_start as i64
    }
}

#[cfg(test)]
mod test {
    use crate::adapter::{Mapping, Strand};
    use crate::aligned_pairs::{AlignOp, TAlignedPairs};

    #[test]
    fn test_align_info() {
        let align_info = Mapping {
            query_name: None,
            query_start: 1,
            query_end: 12,
            strand: Strand::Reverse,
            target_name: None,
            target_len: 22,
            target_start: 0,
            target_end: 10,
            mapq: 0,
            is_primary: true,
            is_supplementary: false,
            cigar: Some(vec![(4, 7), (2, 1), (4, 7), (1, 2), (1, 7)]),
            alignment_score: None,
        };
        let pairs: Vec<_> = align_info.aligned_pairs().collect();
        // reverse：query 从 11 倒走；cigar 4= 2I 4= 1D 1=
        assert_eq!(pairs[0], (Some(11), Some(0), AlignOp::Equal(1)));
        assert_eq!(pairs[4], (Some(7), None, AlignOp::Ins(1)));
        assert_eq!(pairs[6], (Some(5), Some(4), AlignOp::Equal(1)));
        assert_eq!(pairs[10], (None, Some(8), AlignOp::Del(1)));
        assert_eq!(pairs[11], (Some(1), Some(9), AlignOp::Equal(1)));
    }
}
