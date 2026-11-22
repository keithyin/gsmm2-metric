use std::{
    fs::File,
    io::{BufWriter, Write},
    sync::Arc,
    thread,
};

use clap::Parser;
use crossbeam::channel::{Receiver, Sender};
use noodles::{bam, bgzf};
use once_cell::sync::Lazy;

use crate::{
    adapter::Mapping,
    aligner::build_aligner,
    cli::{Cli, MetricArgs, Mode},
    global_data::GlobalData,
    kits::{
        fastx::{FastaFileReader, FastqReader},
        pbar::{self, DEFAULT_INTERVAL},
        phreq::phreq_list_2_quality,
        read_info::ReadInfo,
    },
    metrics::{
        Metric, MetricContext, MetricOutcome, hp_metric_v2::HpMetricV2,
        hp_tr_metric_v2::HpTrMetricV2, hp_tr_metric_v3::HpTrMetricV3,
    },
    msa::msa_entrance,
    params::{InputFilterParams, OupParams},
};

pub mod adapter;
pub mod aligned_pairs;
pub mod aligner;
pub mod cli;
pub mod global_data;
pub mod kits;
pub mod metrics;
pub mod msa;
pub mod params;

/// metric 管线的最低总线程数（1 读入 + N-2 计算/预留 + 1 落盘）。
/// README 的「建议 ≥10 线程」与此对应。
const MIN_THREADS: usize = 10;
/// 读入/落盘两端的预留线程数；计算 worker 数 = 总线程 - 预留。
const OUTPUT_THREAD_RESERVE: usize = 4;
/// 读入 → 计算 → 落盘 各级 channel 的容量（背压上限）。
const CHANNEL_CAP: usize = 1000;
/// BAM 输入的 BGZF 解压 worker 数（等价 rust-htslib set_threads(4)）。
const BAM_BGZF_WORKERS: usize = 4;
/// 参与 metric 统计的 hit 最低 identity（same-target 过滤后应用）。
const MIN_PRIMARY_IDENTITY: f32 = 0.75;

fn main() {
    let mut cli = Cli::parse();
    cli.post_process_param();
    let log_filepath = cli.metric_args.get_log_file(cli.mode);
    let log_file = std::fs::File::create(&log_filepath)
        .unwrap_or_else(|_| panic!("create run log file error: {}", log_filepath));
    let (non_blocking, _guard) = tracing_appender::non_blocking(log_file);

    let time_fmt = time::format_description::parse(
        "[year]-[month padding:zero]-[day padding:zero]#[hour]:[minute]:[second]",
    )
    .unwrap();

    let time_offset = time::UtcOffset::current_local_offset().unwrap_or(time::UtcOffset::UTC);
    let timer = tracing_subscriber::fmt::time::OffsetTime::new(time_offset, time_fmt.clone());
    tracing_subscriber::fmt::fmt()
        .with_timer(timer)
        .with_ansi(false)
        .with_writer(non_blocking)
        .init();

    // --noSeco 为 mm2 版 CLI 兼容旗标：本管线中重比对产生的 secondary hit 在
    // compute_metric 的保留过滤中被无条件丢弃（恒等价于 --noSeco 开启），
    // 旗标本身无可观测效果，仅接受并留痕。
    if cli.metric_args.oup_args.discard_secondary {
        tracing::info!(
            "noSeco: secondary hits are always excluded in this pipeline; flag accepted for mm2 CLI compatibility"
        );
    }

    match cli.mode {
        Mode::HpTrV2 => {
            metric_entrance::<HpTrMetricV2>(&cli.preset, cli.threads, cli.mode, &cli.metric_args);
        }

        Mode::HpTrV3 => {
            metric_entrance::<HpTrMetricV3>(&cli.preset, cli.threads, cli.mode, &cli.metric_args);
        }

        Mode::HpV2 => {
            metric_entrance::<HpMetricV2>(&cli.preset, cli.threads, cli.mode, &cli.metric_args);
        }

        Mode::Bam2Msa => {
            msa_entrance(&cli);
        }
    }
}

fn metric_entrance<M>(preset: &str, tot_threads: Option<usize>, mode: Mode, args: &MetricArgs)
where
    M: Metric,
{
    let tot_threads = tot_threads.unwrap_or(num_cpus::get());
    assert!(
        tot_threads >= MIN_THREADS,
        "at least {MIN_THREADS} threads are needed"
    );

    let (targetname2seq, aligner) = build_aligner(args, preset, Some(tot_threads));

    let global_data = Arc::new(GlobalData::new(targetname2seq));

    let oup_params = args.oup_args.to_oup_params();

    let ref_anchored = args.oup_args.reference_anchored;

    let inp_filter_params = args.io_args.to_input_filter_params();
    thread::scope(|thread_scope| {
        let aligner = &aligner;
        let input_filter_params = &inp_filter_params;
        let oup_params = &oup_params;
        let (qs_sender, qs_recv) = crossbeam::channel::bounded(CHANNEL_CAP);

        thread_scope
            .spawn(move || query_seq_sender(&args.io_args.query, qs_sender, input_filter_params));
        let metric_threads = tot_threads - OUTPUT_THREAD_RESERVE;
        let (metric_sender, metric_recv) = crossbeam::channel::bounded(CHANNEL_CAP);

        for _ in 0..metric_threads {
            thread_scope.spawn({
                let qs_recv_ = qs_recv.clone();
                let metric_sender_ = metric_sender.clone();
                let global_data_ = global_data.clone();
                move || {
                    compute_metric_worker::<M>(
                        qs_recv_,
                        metric_sender_,
                        global_data_,
                        aligner,
                        oup_params,
                        ref_anchored,
                    );
                }
            });
        }
        drop(qs_recv);
        drop(metric_sender);
        let oup_filename = args.get_oup_file(mode);
        let header = M::COLUMNS.join("\t");
        dump_metric_worker(&header, metric_recv, &oup_filename, true);
    });
}

/// mm2 0.25.x / 0.26.1 的 `query_seq_sender` 在 fq/fastq 分支把 qname 后缀错误地
/// 追加到了 `record.seq` 而不是 `record.name`，会污染 read 序列且无法区分多文件
/// 来源。这里在本仓库内重新实现一份并修复该问题；fa / bam 分支行为与 mm2 保持一致。
pub fn query_seq_sender(
    filenames: &[String],
    sender: Sender<ReadInfo>,
    input_filter_params: &InputFilterParams,
) {
    for (file_idx, filename) in filenames.iter().enumerate() {
        let qname_suffix = if filenames.len() == 1 {
            None
        } else {
            Some(format!("__{}", file_idx))
        };
        if filename.ends_with("fa") || filename.ends_with("fasta") || filename.ends_with("fna") {
            let fasta_reader = FastaFileReader::new(filename.clone());

            for mut record in fasta_reader {
                if let Some(suffix) = &qname_suffix {
                    record.name.push_str(suffix);
                }
                sender.send(record).unwrap();
            }
        } else if filename.ends_with("bam") {
            let file = File::open(filename).unwrap();
            let worker_count = std::num::NonZero::new(BAM_BGZF_WORKERS).unwrap();
            let bgzf_reader = bgzf::io::MultithreadedReader::with_worker_count(worker_count, file);
            let mut bam_h = bam::io::Reader::from(bgzf_reader);
            let header = bam_h.read_header().unwrap();
            for record in bam_h.record_bufs(&header) {
                let record = record.unwrap();
                let mut read_info = ReadInfo::from_bam_record(&record, None);
                if input_filter_params.valid(&read_info) {
                    if let Some(suffix) = &qname_suffix {
                        read_info.name.push_str(suffix);
                    }
                    sender.send(read_info).unwrap();
                }
            }
        } else if filename.ends_with("fq") || filename.ends_with("fastq") {
            let fq_iter = FastqReader::new(filename.to_string());
            for mut record in fq_iter {
                if let Some(suffix) = &qname_suffix {
                    record.name.push_str(suffix);
                }
                let rq = phreq_list_2_quality(record.qual.as_ref().unwrap()).unwrap_or(0.);
                record.rq = Some(rq);
                sender.send(record).unwrap();
            }
        } else {
            panic!(
                "invalid file format {}. bam/fa/fasta/fna supported",
                filename
            );
        }
    }
}

pub fn compute_metric_worker<M>(
    receiver: Receiver<ReadInfo>,
    sender: Sender<MetricOutcome>,
    global_data: Arc<GlobalData>,
    aligner: &bandmap::Aligner,
    oup_params: &OupParams,
    ref_anchored: bool,
) where
    M: Metric,
{
    for query_record in receiver {
        let hits = adapter::align_single_query_to_targets(aligner, &query_record);
        sender
            .send(compute_metric::<M>(
                hits,
                query_record,
                oup_params,
                global_data.clone(),
                ref_anchored,
            ))
            .unwrap();
    }
}

/// --noMar 判定：存在 supplementary/secondary chain 即视为 multi-mapping query。
/// 语义对应 mm2 align_worker 的 `records.len() > 1 -> continue`：bandmap 恰有
/// 一个 `is_primary && !is_supplementary` 的 sam_pri mapping，任何其余 mapping
/// （supplementary / secondary，含 `is_primary && is_supplementary` 的原生旗标
/// 形态）都是额外 chain。
fn has_multiple_chains(hits: &[Mapping]) -> bool {
    hits.iter().any(|v| !v.is_primary || v.is_supplementary)
}

pub fn compute_metric<M>(
    hits: Vec<Mapping>,
    read_info: ReadInfo,
    oup_params: &OupParams,
    global_data: Arc<GlobalData>,
    ref_anchored: bool,
) -> MetricOutcome
where
    M: Metric,
{
    dump_hits_if_enabled(&hits);

    // --noMar：query 比对产生不止一个 chain（存在 supplementary/secondary）时，
    // 整条 query 过滤掉，不产出 metric 行。判定必须在所有下游过滤之前——
    // secondary 会在下面的 primary/supp 保留过滤中被丢弃，same-target 过滤会丢
    // 跨 target 的 alt chain，但它们都必须参与 chain 判定。语义对应 mm2
    // align_worker 的 `records.len() > 1 -> continue`（存在 supp/secondary 即扔）。
    if oup_params.discard_multi_align_reads && has_multiple_chains(&hits) {
        tracing::debug!(
            "noMar: query {} has supplementary/secondary chains, filtered",
            read_info.name
        );
        return MetricOutcome::Skipped;
    }

    let hits = hits
        .into_iter()
        .filter(|v| v.is_primary || (v.is_supplementary && !oup_params.discard_supplementary))
        .collect::<Vec<_>>();

    // 取第一个 primary 的 target 名；hits 非空但无 primary 时回退为空串
    // （不再越界 panic，见 bug-report.md Bug 7）。此时下游过滤后必为空，
    // 不会产出任何行，回退值不影响输出。
    let target_name = hits
        .iter()
        .find(|v| v.is_primary)
        .and_then(|v| v.target_name.clone())
        .unwrap_or_default();

    let hits = hits
        .into_iter()
        .filter(|v| v.target_name.as_ref().unwrap() == &target_name)
        .filter(|v| v.identity() > MIN_PRIMARY_IDENTITY)
        .collect::<Vec<_>>();

    let mut metric = M::default();
    metric.compute(MetricContext {
        read: read_info,
        target_name,
        mappings: hits,
        global_data,
        reference_anchored: ref_anchored,
    })
}

/// 调试开关（环境变量解析一次，避免每条 read 探测）。
static DUMP_HITS_ENABLED: Lazy<bool> = Lazy::new(|| std::env::var_os("GSM2_DUMP_HITS").is_some());

/// 调试用：设置环境变量 GSM2_DUMP_HITS 后，把每条 read 的全部 hits
/// （过滤前）打到 stderr，便于端到端对比时定位差异来自比对层还是 metric 层。
fn dump_hits_if_enabled(hits: &[Mapping]) {
    if !*DUMP_HITS_ENABLED {
        return;
    }
    for hit in hits {
        eprintln!(
            "HIT\t{}\t{}\t{:?}\t{}-{}\t{}-{}\tprimary={}\tsupp={}\tidentity={:.4}\tscore={:?}\tcigar={}",
            hit.query_name.as_ref().map(|v| v.as_str()).unwrap_or(""),
            hit.target_name.as_ref().map(|v| v.as_str()).unwrap_or(""),
            hit.strand,
            hit.query_start,
            hit.query_end,
            hit.target_start,
            hit.target_end,
            hit.is_primary,
            hit.is_supplementary,
            hit.identity(),
            hit.alignment_score,
            hit.cigar
                .as_ref()
                .map(|c| c
                    .iter()
                    .map(|(n, op)| format!("{}{}", n, op))
                    .collect::<String>())
                .unwrap_or_default(),
        );
    }
}

pub fn dump_metric_worker(
    header: &str,
    metric_recv: Receiver<MetricOutcome>,
    fname: &str,
    enable_pb: bool,
) {
    let mut writer = BufWriter::new(
        File::create(fname).unwrap_or_else(|_| panic!("create file error: {}", fname)),
    );

    writeln!(&mut writer, "{}", header).unwrap();

    let pb = if enable_pb {
        Some(pbar::get_spin_pb(
            format!("gsmm2-metric: writing metric to {}", fname),
            DEFAULT_INTERVAL,
        ))
    } else {
        None
    };

    for outcome in metric_recv {
        if let Some(v) = pb.as_ref() {
            v.inc(1)
        }
        if let MetricOutcome::Rows(rows) = outcome {
            for row in rows {
                writeln!(&mut writer, "{}", row).unwrap();
            }
        }
    }

    // 显式 flush 并检查错误：drop-flush 会静默吞掉磁盘满等错误，导致
    // 「成功」退出却丢尾部数据
    writer
        .flush()
        .unwrap_or_else(|e| panic!("flush metric output {fname} failed: {e}"));

    if let Some(v) = pb.as_ref() {
        v.finish()
    }
}

#[cfg(test)]
mod test {
    use super::has_multiple_chains;
    use crate::adapter::{Mapping, Strand};

    /// has_multiple_chains 只看旗标，不依赖 cigar/identity，dummy 字段即可。
    fn mapping(is_primary: bool, is_supplementary: bool) -> Mapping {
        Mapping {
            query_name: None,
            query_start: 0,
            query_end: 10,
            strand: Strand::Forward,
            target_name: None,
            target_len: 100,
            target_start: 0,
            target_end: 10,
            mapq: 60,
            is_primary,
            is_supplementary,
            cigar: None,
            alignment_score: None,
        }
    }

    #[test]
    fn test_has_multiple_chains() {
        // 空 hits（unmapped query）：不算 multi，走原有空 metric 逻辑
        assert!(!has_multiple_chains(&[]));

        // 只有 primary：唯一 chain，保留
        assert!(!has_multiple_chains(&[mapping(true, false)]));

        // primary + supplementary（分裂比对）：过滤
        assert!(has_multiple_chains(&[
            mapping(true, false),
            mapping(false, true)
        ]));

        // primary + secondary（多位点 alt chain）：过滤
        assert!(has_multiple_chains(&[
            mapping(true, false),
            mapping(false, false)
        ]));

        // reorder_like_mm2 归一后的「假 primary」（bandmap is_primary&&is_supp
        // 已被归一为非 primary）：同样视为多余 chain
        assert!(has_multiple_chains(&[
            mapping(true, false),
            mapping(true, true)
        ]));
    }
}
