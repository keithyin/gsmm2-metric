//! gskits `utils::Range` 的本地移植(仅保留本项目用到的部分)。

use std::{fmt::Debug, str::FromStr};

/// parse the range str to Vector of tuple
/// range fmt: b1:e1,b2:e2,be3,b4:e4. corresponds to [[b1, e1], [b2, e2], [be3, be3], [b4, e4]]. left inclusive, right inclusive
fn range_parser<T>(range_str: &str) -> Vec<(T, T)>
where
    T: FromStr,
    <T as FromStr>::Err: Debug,
{
    range_str
        .trim()
        .split(",")
        .map(|single_range| {
            let items = single_range
                .trim()
                .split(":")
                .filter(|item| !item.trim().is_empty())
                .collect::<Vec<_>>();
            if items.len() == 2 {
                (
                    items[0].parse::<T>().expect(items[0]),
                    items[1].parse::<T>().expect(items[1]),
                )
            } else {
                (
                    items[0].parse::<T>().expect(items[0]),
                    items[0].parse::<T>().expect(items[0]),
                )
            }
        })
        .collect()
}

#[derive(Debug, Clone)]
pub struct Range<T> {
    range: Vec<(T, T)>,
}

impl<T> Range<T>
where
    T: PartialOrd,
    T: FromStr,
    <T as FromStr>::Err: Debug,
{
    pub fn new(range_str: &str) -> Self {
        Self {
            range: range_parser::<T>(range_str),
        }
    }

    pub fn within_range(&self, v: T) -> bool {
        let mut result = false;
        for (b, e) in &self.range {
            if *b <= v && v <= *e {
                result = true;
                break;
            }
        }
        result
    }
}

#[cfg(test)]
mod test {
    use super::Range;

    #[test]
    fn test_range() {
        let range = Range::<i32>::new("1:3,5,7:9");
        assert!(range.within_range(1));
        assert!(range.within_range(3));
        assert!(range.within_range(5));
        assert!(range.within_range(9));
        assert!(!range.within_range(4));
        assert!(!range.within_range(6));
    }
}
