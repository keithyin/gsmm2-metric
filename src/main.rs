use std::{
    fs::File,
    io::{BufWriter, Write},
    sync::Arc,
    thread,
};

use clap::Parser;
use crossbeam::channel::{Receiver, Sender};
use gskits::{
    ds::ReadInfo,
    fastx_reader::{fasta_reader::FastaFileReader, fastq_reader::FastqReader},
    pbar::{self, DEFAULT_INTERVAL},
    phreq::phreq_list_2_quality,
};
use rust_htslib::bam::Read;

use crate::{
    adapter::Mapping,
    aligner::build_aligner,
    cli::{Cli, MetricArgs, Mode},
    global_data::GlobalData,
    metrics::{
        TMetric, hp_metric_v2::HpMetricV2, hp_tr_metric_v2::HpTrMetricV2,
        hp_tr_metric_v3::HpTrMetricV3,
    },
    msa::msa_entrence,
    params::{InputFilterParams, OupParams},
};

pub mod adapter;
pub mod aligned_pairs;
pub mod aligner;
pub mod cli;
pub mod global_data;
pub mod metrics;
pub mod msa;
pub mod params;

fn main() {
    let cli = Cli::parse();
    let log_filepath = cli.metric_args.get_log_file(cli.mode);
    let log_file = std::fs::File::create(&log_filepath)
        .expect(&format!("create run log file error: {}", &log_filepath));
    let (non_blocking, _guard) = tracing_appender::non_blocking(log_file);

    let time_fmt = time::format_description::parse(
        "[year]-[month padding:zero]-[day padding:zero]#[hour]:[minute]:[second]",
    )
    .unwrap();

    let time_offset =
        time::UtcOffset::current_local_offset().unwrap_or_else(|_| time::UtcOffset::UTC);
    let timer = tracing_subscriber::fmt::time::OffsetTime::new(time_offset, time_fmt.clone());
    tracing_subscriber::fmt::fmt()
        .with_timer(timer)
        .with_ansi(false)
        .with_writer(non_blocking)
        .init();

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
            msa_entrence(&cli);
        }
    }
}

fn metric_entrance<M>(preset: &str, tot_threads: Option<usize>, mode: Mode, args: &MetricArgs)
where
    M: TMetric,
{
    let tot_threads = tot_threads.unwrap_or(num_cpus::get());
    assert!(tot_threads >= 10, "at least 10 threads are needed");

    let (targetname2seq, aligner) = build_aligner(args, preset, Some(tot_threads));

    let global_data = Arc::new(GlobalData::new(targetname2seq));

    let oup_params = args.oup_args.to_oup_params();

    let ref_anchored = args.oup_args.reference_anchored;

    let inp_filter_params = args.io_args.to_input_filter_params();
    thread::scope(|thread_scope| {
        let aligner = &aligner;
        let input_filter_params = &inp_filter_params;
        let oup_params = &oup_params;
        let (qs_sender, qs_recv) = crossbeam::channel::bounded(1000);

        thread_scope.spawn(move || {
            query_seq_sender(
                &args.io_args.query,
                qs_sender,
                input_filter_params,
                oup_params,
            )
        });
        let metric_threads = tot_threads - 4;
        let (metric_sender, metric_recv) = crossbeam::channel::bounded(1000);

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
        dump_metric_worker::<M>(metric_recv, &oup_filename, true);
    });
}

/// mm2 0.25.x / 0.26.1 的 `query_seq_sender` 在 fq/fastq 分支把 qname 后缀错误地
/// 追加到了 `record.seq` 而不是 `record.name`，会污染 read 序列且无法区分多文件
/// 来源。这里在本仓库内重新实现一份并修复该问题；fa / bam 分支行为与 mm2 保持一致。
pub fn query_seq_sender(
    filenames: &Vec<String>,
    sender: Sender<ReadInfo>,
    input_filter_params: &InputFilterParams,
    oup_params: &OupParams,
) {
    let mut file_idx = 0;
    for filename in filenames {
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
            let mut bam_h = rust_htslib::bam::Reader::from_path(filename).unwrap();
            bam_h.set_threads(4).unwrap();
            for record in bam_h.records() {
                let record = record.unwrap();
                if input_filter_params.valid(&record) {
                    sender
                        .send(ReadInfo::from_bam_record(
                            &record,
                            qname_suffix.as_ref().map(|v| v.as_str()),
                            &oup_params.pass_through_tags,
                        ))
                        .unwrap();
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

        file_idx += 1;
    }
}

pub fn compute_metric_worker<M>(
    receiver: Receiver<gskits::ds::ReadInfo>,
    sender: Sender<M>,
    global_data: Arc<GlobalData>,
    aligner: &bandmap::Aligner,
    oup_params: &OupParams,
    ref_anchored: bool,
) where
    M: TMetric,
{
    for query_record in receiver {
        let hits = adapter::align_single_query_to_targets(aligner, &query_record);
        sender
            .send(compute_metric(
                hits,
                &query_record,
                oup_params,
                global_data.clone(),
                ref_anchored,
            ))
            .unwrap();
    }
}

pub fn compute_metric<M>(
    hits: Vec<Mapping>,
    read_info: &gskits::ds::ReadInfo,
    oup_params: &OupParams,
    global_data: Arc<GlobalData>,
    ref_anchored: bool,
) -> M
where
    M: TMetric,
{
    let mut metric = M::default();
    metric.set_qname(read_info.name.clone());

    dump_hits_if_enabled(&hits);

    let hits = hits
        .into_iter()
        .filter(|v| v.is_primary || (v.is_supplementary && !oup_params.discard_supplementary))
        .collect::<Vec<_>>();

    let target_name = if hits.len() > 0 {
        hits.iter()
            .filter(|v| v.is_primary)
            .take(1)
            .map(|v| v.target_name.clone())
            .collect::<Vec<_>>()[0]
            .take()
            .unwrap()
    } else {
        Arc::new("".to_string())
    };

    metric.set_target_name(target_name.clone());

    let hits = hits
        .into_iter()
        .filter(|v| v.target_name.as_ref().unwrap() == &target_name)
        .filter(|v| v.identity() > 0.75)
        .collect::<Vec<_>>();

    metric.set_global_data(global_data);
    metric.set_mappings(hits);
    metric.compute_metric(&read_info, ref_anchored);
    metric.set_metric_str();
    metric
}

/// 调试用：设置环境变量 GSM2_DUMP_HITS 后，把每条 read 的全部 hits
/// （过滤前）打到 stderr，便于端到端对比时定位差异来自比对层还是 metric 层。
fn dump_hits_if_enabled(hits: &[Mapping]) {
    if std::env::var_os("GSM2_DUMP_HITS").is_none() {
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

pub fn dump_metric_worker<M>(metric_recv: Receiver<M>, fname: &str, enable_pb: bool)
where
    M: TMetric,
{
    let mut writer =
        BufWriter::new(File::create(fname).expect(&format!("create file error: {}", fname)));

    writeln!(&mut writer, "{}", M::csv_header()).unwrap();

    let pb = if enable_pb {
        Some(pbar::get_spin_pb(
            format!("gsmm2-metric: writing metric to {}", fname),
            DEFAULT_INTERVAL,
        ))
    } else {
        None
    };

    for metric in metric_recv {
        pb.as_ref().map(|v| v.inc(1));
        let metric_str = metric.get_metric_str().as_ref().unwrap();
        if !metric_str.is_empty() {
            writeln!(&mut writer, "{}", metric_str).unwrap();
        }
    }

    pb.as_ref().map(|v| v.finish());
}
