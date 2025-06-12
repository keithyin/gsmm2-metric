use std::{collections::HashMap, sync::Arc};

use mm2::{
    NoMemLeakAligner,
    gskits::fastx_reader::{fasta_reader::FastaFileReader, read_fastx},
};

use crate::cli::MetricArgs;

pub fn build_aligner(
    args: &MetricArgs,
    preset: &str,
    tot_threads: Option<usize>,
) -> (HashMap<Arc<String>, Arc<String>>, Vec<NoMemLeakAligner>) {
    let tot_threads = tot_threads.unwrap_or(num_cpus::get());
    let target_filename = &args.io_args.target;

    let fa_iter = FastaFileReader::new(target_filename.to_string());
    let targets = read_fastx(fa_iter);
    let targetname2seq = targets
        .iter()
        .map(|v| (Arc::new(v.name.clone()), Arc::new(v.seq.clone())))
        .collect::<HashMap<_, _>>();

    let index_params = args.index_args.to_index_params();
    let map_params = args.map_args.to_map_params();
    let align_params = args.align_args.to_align_params();
    let oup_params = args.oup_args.to_oup_params();

    let mut aligners = mm2::build_aligner(
        preset,
        &index_params,
        &map_params,
        &align_params,
        &oup_params,
        &targets,
        tot_threads,
    );
    aligners.iter_mut().for_each(|aligner| {
        aligner.mapopt.best_n = 10000;
        aligner.mapopt.pri_ratio = 0.2;
    });

    (targetname2seq, aligners)
}
