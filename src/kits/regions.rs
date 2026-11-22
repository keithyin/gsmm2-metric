//! gskits `ds::region` 的本地移植（aligned metric 需要 `ovlp_ratio` / `gaps`）。
//!
//! 语义与 gskits 原版逐行一致，勿随意“优化”：
//! - `Regions` 构造不排序、不合并，完全保留输入顺序；
//! - `ovlp_ratio` 返回 f32（分子 = 未合并区间上的扫描线重叠长度，
//!   端点相接不算重叠；分母 = 合并后总跨度）；
//! - `gaps` 返回 `Vec<i64>`，可为负（重叠/乱序时不 clamp）。

use std::ops::{Deref, DerefMut};

#[derive(Debug, Clone)]
pub struct Region {
    start: usize,
    end: usize,
}

impl Region {
    pub fn new(s: usize, e: usize) -> Self {
        Self { start: s, end: e }
    }
    pub fn length(&self) -> usize {
        self.end - self.start
    }
}

pub struct Regions(Vec<Region>);

impl Deref for Regions {
    type Target = Vec<Region>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Regions {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl From<&Vec<(usize, usize)>> for Regions {
    fn from(value: &Vec<(usize, usize)>) -> Self {
        let regions = value
            .iter()
            .map(|v| Region {
                start: v.0,
                end: v.1,
            })
            .collect::<Vec<_>>();
        Self(regions)
    }
}

impl Regions {
    pub fn total_length(&self) -> usize {
        self.iter().map(|v| v.length()).sum()
    }

    pub fn merge_regions(&self) -> Self {
        let mut regions = self.0.clone();
        regions.sort_by_key(|v| v.start);
        let mut merged_regions: Vec<Region> = Vec::new();

        for region in regions {
            if let Some(last) = merged_regions.last_mut() {
                if region.start <= last.end {
                    last.end = last.end.max(region.end);
                } else {
                    merged_regions.push(region);
                }
            } else {
                merged_regions.push(region);
            }
        }
        Self(merged_regions)
    }

    pub fn ovlp_length(&self) -> usize {
        let mut total_overlap = 0;
        let mut events = Vec::new();

        for region in &self.0 {
            events.push((region.start, 1));
            events.push((region.end, -1));
        }

        events.sort();

        let mut current_depth = 0;
        let mut prev_pos = events[0].0;

        for (pos, change) in events {
            if current_depth > 1 {
                total_overlap += pos - prev_pos;
            }
            current_depth += change;
            prev_pos = pos;
        }

        total_overlap
    }

    /// 覆盖 >=2 次的坐标长度 / 合并后总跨度。
    ///
    /// 前置条件：`self` 非空（gskits 原版在空集合上会因 `ovlp_length` 的
    /// `events[0]` 越界而 panic；调用方沿用 mm2 的「非空才调」模式）。
    pub fn ovlp_ratio(&self) -> f32 {
        let merged_regions = self.merge_regions();
        let len_after_merge = merged_regions.total_length();
        let ovlp_len = self.ovlp_length();
        if len_after_merge == 0 {
            0.0
        } else {
            ovlp_len as f32 / len_after_merge as f32
        }
    }

    /// 相邻区间之间的 gap 值（当前 start - 前一个 end），不 clamp、可为负。
    /// 前置条件：`self` 已按 start 排序。
    /// `init_pos` 计入首个 gap、`end_pos` 计入末尾 gap；空集合返回 `vec![]`。
    pub fn gaps(&self, init_pos: Option<usize>, end_pos: Option<usize>) -> Vec<i64> {
        if self.is_empty() {
            return vec![];
        }
        let mut gaps = vec![];

        let mut pre_pos = init_pos;
        self.iter().for_each(|region| {
            if let Some(pre_pos) = pre_pos {
                let gap = region.start as i64 - pre_pos as i64;
                gaps.push(gap);
            }
            pre_pos = Some(region.end);
        });

        if let Some(end_pos) = end_pos {
            let gap = end_pos as i64 - pre_pos.unwrap() as i64;
            gaps.push(gap);
        }

        gaps
    }
}

#[cfg(test)]
mod test {
    use super::Regions;

    /// gskits 自带用例：完全重叠的两个区间，ratio = 1。
    #[test]
    fn test_ovlp_ratio_full() {
        let regions: Regions = (&vec![(1, 5), (1, 5)]).into();
        assert_eq!(regions.ovlp_length(), 4);
        assert_eq!(regions.merge_regions().total_length(), 4);
        assert_eq!(regions.ovlp_ratio(), 1.0);
    }

    /// gskits 自带用例：端点相接不算重叠，ratio = 0。
    #[test]
    fn test_ovlp_ratio_adjacent() {
        let regions: Regions = (&vec![(1, 5), (5, 10)]).into();
        assert_eq!(regions.ovlp_length(), 0);
        assert_eq!(regions.merge_regions().total_length(), 9);
        assert_eq!(regions.ovlp_ratio(), 0.0);
    }

    /// 部分重叠：[(1,5),(3,8)] 重叠 2bp，合并后 [1,8) 长 7 → 2/7。
    #[test]
    fn test_ovlp_ratio_partial() {
        let regions: Regions = (&vec![(1, 5), (3, 8)]).into();
        assert_eq!(regions.ovlp_length(), 2);
        assert_eq!(regions.merge_regions().total_length(), 7);
        assert_eq!(regions.ovlp_ratio(), 2.0 / 7.0);
    }

    /// gaps：qlen=100，regions=[(1,5),(5,10)] → [1, 0, 90]。
    #[test]
    fn test_gaps() {
        let regions: Regions = (&vec![(1, 5), (5, 10)]).into();
        assert_eq!(regions.gaps(Some(0), Some(100)), vec![1, 0, 90]);
    }

    /// gaps：重叠区间产生负 gap（不 clamp）：[(1,5),(1,5)] → [1, -4, 95]。
    #[test]
    fn test_gaps_negative() {
        let regions: Regions = (&vec![(1, 5), (1, 5)]).into();
        assert_eq!(regions.gaps(Some(0), Some(100)), vec![1, -4, 95]);
    }
}
