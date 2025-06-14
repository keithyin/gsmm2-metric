use std::{
    fs::File,
    io::{BufWriter, Write},
    sync::Arc,
    thread,
};

use clap::Parser;
use crossbeam::channel::{Receiver, Sender};
use gskits::pbar::{self, DEFAULT_INTERVAL};
use mm2::query_seq_sender;

use crate::{
    aligner::build_aligner,
    cli::{Cli, MetricArgs, Mode},
    global_data::GlobalData,
    metrics::{TMetric, hp_metric::HpMetric, hp_tr_metric::HpTrMetric},
};

pub mod aligned_pairs;
pub mod aligner;
pub mod cli;
pub mod global_data;
pub mod metrics;

fn main() {
    let cli = Cli::parse();
    match cli.mode {
        Mode::HpTr => {
            metric_entrance::<HpTrMetric>(&cli.preset, cli.threads, cli.mode, &cli.metric_args);
        }
        Mode::Hp => {
            metric_entrance::<HpMetric>(&cli.preset, cli.threads, cli.mode, &cli.metric_args);
        }
    }
}

fn metric_entrance<M>(preset: &str, tot_threads: Option<usize>, mode: Mode, args: &MetricArgs)
where
    M: TMetric,
{
    let tot_threads = tot_threads.unwrap_or(num_cpus::get());
    assert!(tot_threads >= 10, "at least 10 threads are needed");

    let (targetname2seq, aligners) = build_aligner(args, preset, Some(tot_threads));

    let global_data = Arc::new(GlobalData::new(targetname2seq));

    let oup_params = args.oup_args.to_oup_params();
    let inp_filter_params = args.io_args.to_input_filter_params();
    thread::scope(|thread_scope| {
        let aligners = &aligners;
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
                        aligners,
                        oup_params,
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

pub fn compute_metric_worker<M>(
    receiver: Receiver<mm2::gskits::ds::ReadInfo>,
    sender: Sender<M>,
    global_data: Arc<GlobalData>,
    aligners: &Vec<mm2::NoMemLeakAligner>,
    oup_params: &mm2::params::OupParams,
) where
    M: TMetric,
{
    for query_record in receiver {
        let hits = mm2::align_single_query_to_targets(&query_record, aligners);
        sender
            .send(compute_metric(
                hits,
                &query_record,
                oup_params,
                global_data.clone(),
            ))
            .unwrap();
    }
}

pub fn compute_metric<M>(
    hits: Vec<mm2::minimap2::Mapping>,
    read_info: &mm2::gskits::ds::ReadInfo,
    oup_params: &mm2::params::OupParams,
    global_data: Arc<GlobalData>,
) -> M
where
    M: TMetric,
{
    let mut metric = M::default();
    metric.set_qname(read_info.name.clone());

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
        .collect::<Vec<_>>();

    metric.set_global_data(global_data);
    metric.set_mappings(hits);
    metric.compute_metric(&read_info);
    metric.set_metric_str();
    metric
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
            format!("gsmm2-time-err: writing metric to {}", fname),
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
