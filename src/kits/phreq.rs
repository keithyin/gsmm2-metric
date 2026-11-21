//! gskits `phreq` 模块的本地移植(仅保留本项目用到的 phreq_list_2_quality)。

pub fn phreq_list_2_quality(baseq: &[u8]) -> Option<f32> {
    if baseq.len() == 0 {
        return None;
    }
    let length = baseq.len() as f64;
    let err_rate = baseq
        .iter()
        .map(|v| *v as f64)
        .map(|v| 10.0_f64.powf(v / -10.0_f64))
        .reduce(|acc, v| acc + v)
        .unwrap()
        / length;
    Some((1.0_f64 - err_rate) as f32)
}

#[cfg(test)]
mod test {
    use super::phreq_list_2_quality;

    #[test]
    fn test_phreq_list_2_quality() {
        let quality = phreq_list_2_quality(&[20, 20, 20, 20]);
        assert!((quality.unwrap() - 0.99) < 1e-3);

        let quality = phreq_list_2_quality(&[30, 30, 30, 30]);
        assert!((quality.unwrap() - 0.999) < 1e-4);
    }
}
