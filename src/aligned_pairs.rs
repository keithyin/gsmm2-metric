/*

0 => Cigar::Match(cnt),
                1 => Cigar::Ins(cnt),
                2 => Cigar::Del(cnt),
                3 => Cigar::RefSkip(cnt),
                4 => Cigar::SoftClip(cnt),
                5 => Cigar::HardClip(cnt),
                6 => Cigar::Pad(cnt),
                7 => Cigar::Equal(cnt),
                8 => Cigar::Diff(cnt),

*/

use mm2::minimap2::Mapping;

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
        &self.alignment.as_ref().unwrap().cigar.as_ref().unwrap()
    }

    fn is_reverse(&self) -> bool {
        match self.strand {
            mm2::minimap2::Strand::Forward => false,
            mm2::minimap2::Strand::Reverse => true,
        }
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
    use mm2::minimap2::{Alignment, Mapping};

    use crate::aligned_pairs::TAlignedPairs;

    #[test]
    fn test_align_info() {
        let align_info = Mapping {
            query_name: None,
            query_len: None,
            query_start: 1,
            query_end: 12,
            strand: mm2::minimap2::Strand::Reverse,
            target_name: None,
            target_len: 22,
            target_start: 0,
            target_end: 10,
            match_len: 0,
            block_len: 0,
            mapq: 0,
            is_primary: true,
            is_supplementary: false,
            alignment: Some(Alignment {
                nm: 0,
                cigar: Some(vec![(4, 7), (2, 1), (4, 7), (1, 2), (1, 7)]),
                cigar_str: None,
                md: None,
                cs: None,
                alignment_score: None,
            }),
        };
        for align_info in align_info.aligned_pairs() {
            println!("align_info:{:?}", align_info);
        }
    }
}
