use std::{fs, path, str::FromStr};

use clap::{self, Args, Parser};
use mm2::params::{AlignParams, IndexParams, InputFilterParams, MapParams, OupParams};

#[derive(Debug, Parser, Clone)]
#[command(version, about, long_about=None)]
pub struct Cli {
    #[arg(long = "threads")]
    pub threads: Option<usize>,

    #[arg(long="preset", default_value_t=String::from_str("map-ont").unwrap(), 
        help="read https://lh3.github.io/minimap2/minimap2.html for more details")]
    pub preset: String,

    #[arg(value_enum, long="mode")]
    pub mode: Mode,

    #[arg(
        long = "short-aln",
        help = "if the 30 < query_len OR target <200, use short aln mode"
    )]
    pub short_aln: bool,

    #[command(flatten)]
    pub metric_args: MetricArgs,
}
impl Cli {
    pub fn post_process_param(&mut self) {
        if self.short_aln {
            if self.metric_args.index_args.kmer.is_none() {
                self.metric_args.index_args.kmer = Some(5);
            }
            if self.metric_args.index_args.wins.is_none() {
                self.metric_args.index_args.wins = Some(3);
            }

            if self.metric_args.align_args.min_cnt.is_none() {
                self.metric_args.align_args.min_cnt = Some(3);
            }
            if self.metric_args.align_args.min_dp_max.is_none() {
                self.metric_args.align_args.min_dp_max = Some(20);
            }

            if self.metric_args.align_args.min_chain_score.is_none() {
                self.metric_args.align_args.min_chain_score = Some(5);
            }

            if self.metric_args.align_args.min_ksw_len.is_none() {
                self.metric_args.align_args.min_ksw_len = Some(0);
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, clap::ValueEnum)]
pub enum Mode {
    /// homopolymer and tandem repeats
    HpTr,

    HpTrV2,

    /// homopolymer only
    Hp,

    HpV2,
}

#[derive(Debug, Args, Clone, Copy, Default)]
pub struct IndexArgs {
    #[arg(long, help = "minimizer kmer")]
    pub kmer: Option<usize>,

    #[arg(long, help = "minimizer window size")]
    pub wins: Option<usize>,
}

impl IndexArgs {
    pub fn to_index_params(&self) -> IndexParams {
        let mut param = IndexParams::new();
        param = if let Some(kmer) = self.kmer {
            param.set_kmer(kmer)
        } else {
            param
        };

        param = if let Some(wins) = self.wins {
            param.set_wins(wins)
        } else {
            param
        };

        param
    }
}

#[derive(Debug, Args, Clone)]
pub struct MetricArgs {
    #[command(flatten)]
    pub io_args: IoArgs,

    #[command(flatten)]
    pub index_args: IndexArgs,

    #[command(flatten)]
    pub map_args: MapArgs,

    #[command(flatten)]
    pub align_args: AlignArgs,

    #[command(flatten)]
    pub oup_args: OupArgs,
}

impl MetricArgs {
    pub fn get_oup_file(&self, mode: Mode) -> String {
        self.io_args.get_oup_filename(mode)
    }
    pub fn get_log_file(&self, mode: Mode) -> String {
        format!("{}.log", self.get_oup_file(mode))
    }
}

#[derive(Debug, Args, Clone)]
pub struct IoArgs {
    #[arg(
        short = 'q',
        help = "query file paths, 
    if multiple query_filepath are provided, 
    their query_names will be rewriten to ___0, ___1, and so on, 
    based on the order of the filenames. valid file format bam/fa/fq",
        required = true
    )]
    pub query: Vec<String>,

    /// target.fasta
    #[arg(long = "target", short = 't')]
    pub target: String,

    #[arg(long = "out", help = "output filename")]
    pub outfn: Option<String>,

    #[arg(
        long = "np-range",
        help = "1:3,5,7:9 means [[1, 3], [5, 5], [7, 9]]. only valid for bam input that contains np field"
    )]
    pub np_range: Option<String>,

    #[arg(
        long = "rq-range",
        help = "0.9:1.1 means 0.9<=rq<=1.1. only valid for bam input that contains rq field"
    )]
    pub rq_range: Option<String>,

    #[arg(long = "qname-suffix", help = "suffix for query name")]
    pub qname_suffix: Option<String>,
}

impl IoArgs {
    pub fn get_oup_filename(&self, mode: Mode) -> String {
        let out_fn = if let Some(oup) = &self.outfn {
            oup.to_string()
        } else {
            format!(
                "{}.{:?}-fact.csv",
                self.query[0].rsplit_once(".").unwrap().0,
                mode
            )
        };
        if let Some(dir) = path::Path::new(&out_fn).to_path_buf().parent() {
            if !dir.exists() {
                fs::create_dir_all(dir).unwrap();
            }
        }

        out_fn
    }

    pub fn to_input_filter_params(&self) -> InputFilterParams {
        let mut param = InputFilterParams::new();
        param = if let Some(ref np_range_str) = self.np_range {
            param.set_np_range(np_range_str)
        } else {
            param
        };

        param = if let Some(ref rq_range_str) = self.rq_range {
            param.set_rq_range(rq_range_str)
        } else {
            param
        };

        param.qname_suffix = self.qname_suffix.clone();

        param
    }
}

#[derive(Debug, Args, Clone, Default)]
pub struct MapArgs {}

impl MapArgs {
    pub fn to_map_params(&self) -> MapParams {
        MapParams::default()
    }
}

#[derive(Debug, Args, Clone, Default)]
pub struct AlignArgs {
    #[arg(short = 'm', help = "matching_score>=0, recommend 2")]
    matching_score: Option<i32>,

    #[arg(short = 'M', help = "mismatch_penalty >=0, recommend 4")]
    mismatch_penalty: Option<i32>,

    #[arg(short = 'o', help = "gap_open_penalty >=0, recommend 4,24")]
    gap_open_penalty: Option<String>,

    #[arg(short = 'e', help = "gap_extension_penalty >=0, recommend 2,1")]
    gap_extension_penalty: Option<String>,

    #[arg(long = "min-cnt", help = "min_cnt")]
    pub min_cnt: Option<i32>,
    #[arg(long = "min-dp-max", help = "min dp max")]
    pub min_dp_max: Option<i32>,
    #[arg(long = "min-chain-score", help = "min chain score")]
    pub min_chain_score: Option<i32>,
    #[arg(long = "min-ksw-len", help = "min ksw len")]
    pub min_ksw_len: Option<i32>,
}

impl AlignArgs {
    pub fn to_align_params(&self) -> AlignParams {
        let mut param = AlignParams::new();
        param.matching_score = self.matching_score;
        param.mismatch_penalty = self.mismatch_penalty;

        param = if let Some(ref go) = self.gap_open_penalty {
            param.set_gap_open_penalty(go.to_string())
        } else {
            param
        };

        param = if let Some(ref ge) = self.gap_extension_penalty {
            param.set_gap_extension_penalty(ge.to_string())
        } else {
            param
        };

        param.min_cnt = self.min_cnt;
        param.min_dp_max = self.min_dp_max;
        param.min_chain_score = self.min_chain_score;
        param.min_ksw_len = self.min_ksw_len;

        param
    }
}

#[derive(Debug, Args, Clone, Default)]
pub struct OupArgs {
    #[arg(
        long = "noSupp",
        help = "discard supplementary alignment, only keep the primary alignemnt for metric"
    )]
    pub discard_supplementary: bool,

    #[arg(
        long = "noMar",
        help = "discard multi-mapping reads, only keep the unique mapping reads for metric"
    )]
    pub discard_multi_mapping_reads: bool,
}

impl OupArgs {
    pub fn to_oup_params(&self) -> OupParams {
        let mut param = OupParams::new();
        param = param
            .set_discard_secondary(true)
            .set_discard_supplementary(self.discard_supplementary)
            .set_oup_identity_threshold(-1.0)
            .set_oup_coverage_threshold(-1.0)
            .set_discard_multi_align_reads(self.discard_multi_mapping_reads);
        param
    }
}
