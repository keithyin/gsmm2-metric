use std::{collections::HashMap, sync::Arc};

use bandmap::Preset;
use gskits::fastx_reader::{fasta_reader::FastaFileReader, read_fastx};

use crate::{
    cli::MetricArgs,
    params::{AlignParams, IndexParams},
};

/// 把 CLI 的 preset 字符串映射为 bandmap 的 Preset。
/// mm2 原来只支持 map-ont / map-pb / map-hifi；这里沿用 bandmap 的完整集合。
pub fn parse_preset(preset: &str) -> Option<Preset> {
    match preset {
        "lr" | "map-ont" => Some(Preset::MapOnt),
        "map10k" | "map-pb" => Some(Preset::MapPb),
        "map-hifi" | "lr:hq" => Some(Preset::MapHifi),
        "lr:hqae" => Some(Preset::LrHqae),
        "map-iclr" => Some(Preset::MapIclr),
        "sr" | "short" => Some(Preset::Sr),
        "splice" => Some(Preset::Splice),
        "splice:hq" => Some(Preset::SpliceHq),
        "splice:sr" => Some(Preset::SpliceSr),
        "cdna" => Some(Preset::Cdna),
        "asm5" => Some(Preset::Asm5),
        "asm10" => Some(Preset::Asm10),
        "asm20" => Some(Preset::Asm20),
        "ava-ont" => Some(Preset::AvaOnt),
        "ava-pb" => Some(Preset::AvaPb),
        _ => None,
    }
}

/// 用参数结构构建 bandmap Aligner（单索引覆盖全部 target）。
///
/// 与 mm2 per-target aligner 数组不同，bandmap 一次建全索引；`Send + Sync`，
/// 跨线程共享一个实例即可。
pub fn build_aligner_from_params(
    preset: &str,
    index_params: &IndexParams,
    align_params: &AlignParams,
    targets: &[gskits::ds::ReadInfo],
) -> bandmap::Aligner {
    let preset = parse_preset(preset).unwrap_or_else(|| panic!("unknown preset: {}", preset));

    let seqs: Vec<(String, Vec<u8>)> = targets
        .iter()
        .map(|v| (v.name.clone(), v.seq.as_bytes().to_vec()))
        .collect();

    let mut builder = bandmap::Aligner::builder(preset);
    if let Some(kmer) = index_params.kmer {
        builder = builder.k(kmer);
    }
    if let Some(wins) = index_params.wins {
        builder = builder.w(wins);
    }
    // mid_occ：保持 bandmap 默认的按索引校准（mid_occ_frac=2e-4，与 CLI minimap2 的
    // mm_mapopt_update 语义一致）。注意 mm2 的内存索引路径（minimap2-rs
    // with_seq_and_id）写死 mapopt.mid_occ = 1000，不走校准，属于该 crate 的实现
    // 怪癖；实测在 S. aureus 合成数据上钉死 1000 反而破坏一致性（poly-A 高频 seed
    // 进入链合导致重对齐窗口漂移），校准值（MG1655 实测 12）一致性更好，故不模仿。

    let mut aligner = builder.from_seqs(seqs).unwrap();

    align_params.apply_to(aligner.options_mut());

    // 与 mm2 版本的 hack 保持一致：保留尽量多的 secondary/alternative hit，
    // 过滤交给 main.rs 的 identity>0.75 / 同 target 逻辑。
    aligner.options_mut().filtering.best_n = 10000;
    aligner.options_mut().filtering.pri_ratio = 0.2;

    // identity / aligned_pairs 依赖 EQX cigar（op 7/8），必须显式打开。
    aligner.output_config_mut().eqx = true;

    aligner
}

pub fn build_aligner(
    args: &MetricArgs,
    preset: &str,
    _tot_threads: Option<usize>,
) -> (HashMap<Arc<String>, Arc<String>>, bandmap::Aligner) {
    let target_filename = &args.io_args.target;

    let fa_iter = FastaFileReader::new(target_filename.to_string());
    let targets = read_fastx(fa_iter);
    let targetname2seq = targets
        .iter()
        .map(|v| (Arc::new(v.name.clone()), Arc::new(v.seq.clone())))
        .collect::<HashMap<_, _>>();

    let index_params = args.index_args.to_index_params();
    let align_params = args.align_args.to_align_params();

    let aligner = build_aligner_from_params(preset, &index_params, &align_params, &targets);

    (targetname2seq, aligner)
}
